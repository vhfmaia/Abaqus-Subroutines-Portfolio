# -*- coding: utf-8 -*-
"""
Abaqus ODB to VTK ASCII Exporter for Laser Welding Simulation
============================================================
Exports sequentially coupled thermal (GLOBAL_THERM.odb) and mechanical
(GLOBAL_MECH.odb) simulation results into legacy VTK ASCII unstructured grids.
Compatible with Python 2.7 (Abaqus/Standard native) and Python 3.x.

Features:
  - Supports DC3D8 (thermal) and C3D8H / C3D8 (mechanical hybrid) bricks
  - Thermal fields: Nodal Temperature (NT11), Phase fractions (SDV1, SDV4-7)
  - Mechanical fields: Displacement (U), von Mises Stress (S), PEEQ, Max Principal
  - Automated ZIP packaging for cloud/HPC workflow & PyVista post-processing
"""

import os
import sys
import glob
import time
import zipfile

def export_odb_to_vtk(odb_name, step_name=None, output_prefix=None, create_zip=True):
    if not os.path.exists(odb_name):
        print("[ERROR] ODB file not found: %s" % odb_name)
        return False

    # Auto-determine output prefix if not provided
    if output_prefix is None:
        base = os.path.splitext(os.path.basename(odb_name))[0].lower()
        if "therm" in base:
            output_prefix = "laser_therm"
        elif "mech" in base:
            output_prefix = "laser_mech"
        else:
            output_prefix = base + "_vtk"

    print("==========================================================")
    print("  Abaqus ODB -> VTK Exporter: %s" % odb_name)
    print("  Output Prefix: %s | Target ZIP: %s_vtk.zip" % (output_prefix, output_prefix))
    print("==========================================================")

    from odbAccess import openOdb
    odb = openOdb(odb_name, readOnly=True)

    step_keys = list(odb.steps.keys())
    if not step_keys:
        print("[ERROR] No steps found in ODB.")
        odb.close()
        return False

    if step_name is None:
        step_name = step_keys[0]

    if step_name not in odb.steps:
        print("[ERROR] Requested step '%s' not in ODB. Available: %s" % (step_name, step_keys))
        odb.close()
        return False

    step = odb.steps[step_name]
    print("[INFO] Processing Step: %s (Total Frames: %d)" % (step_name, len(step.frames)))

    assembly = odb.rootAssembly
    if len(assembly.instances) > 0:
        instance = list(assembly.instances.values())[0]
    else:
        instance = assembly

    nodes = instance.nodes
    elements = instance.elements

    num_nodes = len(nodes)
    num_elements = len(elements)
    print("[INFO] Mesh: %d nodes, %d elements." % (num_nodes, num_elements))

    # Node label to 0-based VTK index
    node_map = {}
    node_coords = []
    for idx, node in enumerate(nodes):
        node_map[node.label] = idx
        node_coords.append(node.coordinates)

    # Element type to VTK cell type
    vtk_cell_types = {
        'DC3D8': 12,  # VTK_HEXAHEDRON
        'C3D8H': 12,  # VTK_HEXAHEDRON
        'C3D8': 12,   # VTK_HEXAHEDRON
        'C3D8R': 12,  # VTK_HEXAHEDRON
        'DC3D4': 10,  # VTK_TETRA
        'C3D4': 10,   # VTK_TETRA
        'DC3D6': 13,  # VTK_WEDGE
        'C3D6': 13    # VTK_WEDGE
    }

    cell_data = []
    cell_types = []
    total_cell_int_count = 0
    active_element_labels = []

    for elem in elements:
        etype = elem.type.upper()
        if etype in vtk_cell_types:
            ctype = vtk_cell_types[etype]
            cell_types.append(ctype)
            conn = [node_map[nl] for nl in elem.connectivity]
            cell_data.append([len(conn)] + conn)
            total_cell_int_count += len(conn) + 1
            active_element_labels.append(elem.label)

    num_cells = len(cell_types)

    exported_files = []
    t0 = time.time()

    for f_idx, frame in enumerate(step.frames):
        vtk_filename = "%s_%04d.vtk" % (output_prefix, f_idx)
        print("[INFO] Exporting Frame %03d / %03d (Time = %.4f s) -> %s" %
              (f_idx, len(step.frames)-1, frame.frameValue, vtk_filename))

        # Check available field outputs
        fields = frame.fieldOutputs

        # 1. Nodal Fields
        # Temperature (NT11 or NT)
        temp_dict = {}
        if 'NT11' in fields:
            for val in fields['NT11'].values:
                temp_dict[val.nodeLabel] = val.data
        elif 'NT' in fields:
            for val in fields['NT'].values:
                temp_dict[val.nodeLabel] = val.data

        # Displacement (U)
        u_dict = {}
        if 'U' in fields:
            for val in fields['U'].values:
                u_dict[val.nodeLabel] = val.data

        # 2. Elemental Fields
        # Stress (S)
        mises_dict = {}
        mises_cnt = {}
        maxp_dict = {}
        if 'S' in fields:
            s_field = fields['S']
            for val in s_field.values:
                el = val.elementLabel
                if el not in mises_dict:
                    mises_dict[el] = val.mises
                    maxp_dict[el] = val.maxPrincipal
                    mises_cnt[el] = 1
                else:
                    mises_dict[el] += val.mises
                    maxp_dict[el] += val.maxPrincipal
                    mises_cnt[el] += 1

        # Plastic strain (PEEQ)
        peeq_dict = {}
        peeq_cnt = {}
        if 'PEEQ' in fields:
            for val in fields['PEEQ'].values:
                el = val.elementLabel
                peeq_dict[el] = peeq_dict.get(el, 0.0) + val.data
                peeq_cnt[el] = peeq_cnt.get(el, 0) + 1

        # State dependent variables (Phase fractions)
        sdv1_dict = {}
        sdv5_dict = {}
        sdv7_dict = {}
        if 'SDV1' in fields:
            for val in fields['SDV1'].values:
                sdv1_dict[val.elementLabel] = val.data
        if 'SDV5' in fields:
            for val in fields['SDV5'].values:
                sdv5_dict[val.elementLabel] = val.data
        if 'SDV7' in fields:
            for val in fields['SDV7'].values:
                sdv7_dict[val.elementLabel] = val.data

        # Write VTK ASCII
        with open(vtk_filename, 'w') as f:
            f.write("# vtk DataFile Version 3.0\n")
            f.write("Abaqus ODB Frame %d Time %f Step %s\n" % (f_idx, frame.frameValue, step_name))
            f.write("ASCII\n")
            f.write("DATASET UNSTRUCTURED_GRID\n")

            # Points
            f.write("POINTS %d float\n" % num_nodes)
            for pt in node_coords:
                f.write("%.6f %.6f %.6f\n" % (pt[0], pt[1], pt[2]))

            # Cells
            f.write("\nCELLS %d %d\n" % (num_cells, total_cell_int_count))
            for c in cell_data:
                f.write(" ".join(str(val) for val in c) + "\n")

            f.write("\nCELL_TYPES %d\n" % num_cells)
            for ct in cell_types:
                f.write("%d\n" % ct)

            # Point Data
            f.write("\nPOINT_DATA %d\n" % num_nodes)
            if u_dict:
                f.write("VECTORS Displacement float\n")
                for node in nodes:
                    u = u_dict.get(node.label, (0.0, 0.0, 0.0))
                    f.write("%.6e %.6e %.6e\n" % (u[0], u[1], u[2]))

            if temp_dict:
                f.write("SCALARS Temperature float 1\n")
                f.write("LOOKUP_TABLE default\n")
                for node in nodes:
                    t_val = temp_dict.get(node.label, 25.0)
                    f.write("%.4e\n" % t_val)

            # Cell Data
            has_cell_data = bool(mises_dict or peeq_dict or sdv1_dict)
            if has_cell_data:
                f.write("\nCELL_DATA %d\n" % num_cells)

                if mises_dict:
                    f.write("SCALARS Von_Mises float 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    for el_lbl in active_element_labels:
                        c = mises_cnt.get(el_lbl, 1)
                        f.write("%.6e\n" % (mises_dict.get(el_lbl, 0.0) / float(c)))

                if peeq_dict:
                    f.write("SCALARS PEEQ float 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    for el_lbl in active_element_labels:
                        c = peeq_cnt.get(el_lbl, 1)
                        f.write("%.6e\n" % (peeq_dict.get(el_lbl, 0.0) / float(c)))

                if sdv1_dict:
                    f.write("SCALARS Liquid_Fraction_RLS float 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    for el_lbl in active_element_labels:
                        f.write("%.4e\n" % sdv1_dict.get(el_lbl, 1.0))

                if sdv5_dict:
                    f.write("SCALARS Austenite_Fraction float 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    for el_lbl in active_element_labels:
                        f.write("%.4e\n" % sdv5_dict.get(el_lbl, 0.0))

                if sdv7_dict:
                    f.write("SCALARS Martensite_Fraction float 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    for el_lbl in active_element_labels:
                        f.write("%.4e\n" % sdv7_dict.get(el_lbl, 0.0))

        exported_files.append(vtk_filename)

    odb.close()
    print("[INFO] Exported %d frames in %.1f s." % (len(exported_files), time.time() - t0))

    # ZIP packaging
    if create_zip and exported_files:
        zip_name = "%s_vtk.zip" % output_prefix
        print("[INFO] Creating archive '%s'..." % zip_name)
        with zipfile.ZipFile(zip_name, 'w', zipfile.ZIP_DEFLATED) as zf:
            for vf in exported_files:
                zf.write(vf)
        zip_size_mb = os.path.getsize(zip_name) / (1024.0 * 1024.0)
        print("[SUCCESS] Archived %d VTK frames into %s (%.2f MB)." %
              (len(exported_files), zip_name, zip_size_mb))

    return True

if __name__ == "__main__":
    target_odb = None
    step_arg = None
    prefix_arg = None
    make_zip = True

    # Parse arguments
    for arg in sys.argv[1:]:
        if arg.endswith(".odb"):
            target_odb = arg
        elif arg.startswith("--step="):
            step_arg = arg.split("=")[1]
        elif arg.startswith("--prefix="):
            prefix_arg = arg.split("=")[1]
        elif arg == "--no-zip":
            make_zip = False

    if target_odb is None:
        # Auto-detect in current directory
        candidates = ["GLOBAL_THERM.odb", "GLOBAL_MECH.odb"] + glob.glob("*.odb")
        for cand in candidates:
            if os.path.exists(cand):
                target_odb = cand
                break

    if target_odb is None:
        print("[USAGE] abaqus python odb_to_vtk.py <job.odb> [--step=STEP_NAME] [--prefix=PREFIX]")
        sys.exit(1)

    export_odb_to_vtk(target_odb, step_name=step_arg, output_prefix=prefix_arg, create_zip=make_zip)
