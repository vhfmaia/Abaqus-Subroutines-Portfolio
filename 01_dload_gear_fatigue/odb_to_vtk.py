# -*- coding: utf-8 -*-
"""
Abaqus ODB to VTK ASCII Exporter
Pure 7-bit ASCII format compatible with Abaqus Python (2.7 & 3.x).
"""
import os
import sys
from odbAccess import openOdb

def export_odb_to_vtk(odb_name, step_name=None, output_prefix="gear_vtk"):
    if not os.path.exists(odb_name):
        print("[ERROR] ODB file not found: %s" % odb_name)
        return

    print("[INFO] Opening ODB: %s" % odb_name)
    odb = openOdb(odb_name, readOnly=True)

    # 1. Select analysis step
    step_keys = list(odb.steps.keys())
    if not step_keys:
        print("[ERROR] No steps found in ODB.")
        odb.close()
        return

    if step_name is None:
        step_name = step_keys[0]
    step = odb.steps[step_name]
    print("[INFO] Processing Step: %s" % step_name)

    # 2. Select instance / assembly
    assembly = odb.rootAssembly
    if len(assembly.instances) > 0:
        instance = list(assembly.instances.values())[0]
    else:
        instance = assembly

    nodes = instance.nodes
    elements = instance.elements

    # Map node labels to 0-based VTK point indices
    node_map = {}
    node_coords = []
    for idx, node in enumerate(nodes):
        node_map[node.label] = idx
        node_coords.append(node.coordinates)

    # Abaqus Element Type to VTK Cell Type mapping
    vtk_cell_types = {
        'C3D4': 10,   # VTK_TETRA (Linear tetrahedron)
        'C3D10': 24,  # VTK_QUADRATIC_TETRA (Quadratic tetrahedron)
        'C3D8': 12,   # VTK_HEXAHEDRON (Linear brick)
        'C3D8R': 12,  # VTK_HEXAHEDRON (Reduced integration brick)
        'C3D6': 13    # VTK_WEDGE (Linear prism)
    }

    # Extract valid mesh cells and connectivity
    cell_data = []
    cell_types = []
    total_cell_int_count = 0
    active_element_labels = []

    for elem in elements:
        elem_type = elem.type.upper()
        if elem_type in vtk_cell_types:
            c_type = vtk_cell_types[elem_type]
            cell_types.append(c_type)
            conn = [node_map[nl] for nl in elem.connectivity]
            cell_data.append([len(conn)] + conn)
            total_cell_int_count += len(conn) + 1
            active_element_labels.append(elem.label)

    num_nodes = len(nodes)
    num_cells = len(cell_types)
    print("[INFO] Mesh extracted: %d nodes, %d elements." % (num_nodes, num_cells))

    # 3. Process field outputs frame by frame
    for f_idx, frame in enumerate(step.frames):
        vtk_filename = "%s_%04d.vtk" % (output_prefix, f_idx)
        print("[INFO] Exporting Frame %d (Time = %.4f s) -> %s" % (f_idx, frame.frameValue, vtk_filename))

        # Extract nodal displacements (U)
        u_dict = {}
        if 'U' in frame.fieldOutputs:
            u_field = frame.fieldOutputs['U']
            for val in u_field.values:
                u_dict[val.nodeLabel] = val.data

        # Extract element von Mises stress (averaged per element)
        mises_sum = {}
        mises_cnt = {}
        if 'S' in frame.fieldOutputs:
            s_field = frame.fieldOutputs['S']
            for val in s_field.values:
                el_lbl = val.elementLabel
                if el_lbl not in mises_sum:
                    mises_sum[el_lbl] = val.mises
                    mises_cnt[el_lbl] = 1
                else:
                    mises_sum[el_lbl] += val.mises
                    mises_cnt[el_lbl] += 1

        # Write Legacy VTK ASCII file
        with open(vtk_filename, 'w') as f:
            f.write("# vtk DataFile Version 3.0\n")
            f.write("Abaqus ODB Export Frame %d Time %f\n" % (f_idx, frame.frameValue))
            f.write("ASCII\n")
            f.write("DATASET UNSTRUCTURED_GRID\n")

            # Point coordinates
            f.write("POINTS %d float\n" % num_nodes)
            for pt in node_coords:
                f.write("%.6f %.6f %.6f\n" % (pt[0], pt[1], pt[2]))

            # Cell connectivity
            f.write("\nCELLS %d %d\n" % (num_cells, total_cell_int_count))
            for c in cell_data:
                f.write(" ".join(str(val) for val in c) + "\n")

            # Cell types
            f.write("\nCELL_TYPES %d\n" % num_cells)
            for ct in cell_types:
                f.write("%d\n" % ct)

            # Point data: Displacement vectors
            f.write("\nPOINT_DATA %d\n" % num_nodes)
            f.write("VECTORS Displacement float\n")
            for node in nodes:
                u = u_dict.get(node.label, (0.0, 0.0, 0.0))
                f.write("%.6e %.6e %.6e\n" % (u[0], u[1], u[2]))

            # Cell data: Von Mises scalar field
            f.write("\nCELL_DATA %d\n" % num_cells)
            f.write("SCALARS Von_Mises float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for el_lbl in active_element_labels:
                if el_lbl in mises_sum and mises_cnt.get(el_lbl, 0) > 0:
                    avg_mises = mises_sum[el_lbl] / float(mises_cnt[el_lbl])
                else:
                    avg_mises = 0.0
                f.write("%.6e\n" % avg_mises)

    odb.close()
    print("[INFO] Export completed successfully.")

if __name__ == "__main__":
    import glob
    target_odb = "Gear_torque.odb"
    if len(sys.argv) > 1 and sys.argv[1].endswith(".odb"):
        target_odb = sys.argv[1]
    elif not os.path.exists(target_odb):
        found_odbs = glob.glob("*.odb")
        if found_odbs:
            target_odb = found_odbs[0]
            print("[INFO] Target ODB not specified, auto-detected: %s" % target_odb)

    export_odb_to_vtk(target_odb, output_prefix="gear_helical")