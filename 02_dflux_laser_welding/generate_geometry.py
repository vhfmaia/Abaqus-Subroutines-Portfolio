"""
generate_geometry.py
====================
Generates a generic, fully parametric two-part axisymmetric assembly for a
circumferential weld simulation in Abaqus/Standard (heat transfer, DC3D8):

    SHAFT : hollow cylinder  r in [R_BORE, R_JOINT],  z in [0, H_SHAFT]
    HUB   : ring             r in [R_JOINT, R_OUT],   z in [Z_HUB, H_SHAFT]

The two parts share the cylindrical joint surface r = R_JOINT, but are meshed
as independent bodies (duplicate nodes) and connected in the main input file
with a surface-to-surface *TIE. The top faces are flush at z = H_SHAFT, where
the beam (DFLUX.for) travels along the joint line.

The structured hexahedral mesh is graded so that the elements are smallest
(~H_FINE) in the weld zone and grow geometrically away from it.

Output: GEOMETRY_THERM.inp  (nodes, elements, element/node sets, surfaces)

Usage:  python generate_geometry.py [output.inp]
"""
import sys
import numpy as np

# ---------------------------------------------------------------------------
# Geometry (mm) - must be consistent with RADIUS and Z0 in DFLUX.for
# ---------------------------------------------------------------------------
R_BORE = 8.0     # shaft bore radius
R_JOINT = 15.0   # joint radius  (= weld path RADIUS in DFLUX.for)
R_OUT = 28.0     # hub outer radius
H_SHAFT = 20.0   # shaft height = top surface (= Z0 in DFLUX.for)
Z_HUB = 8.0      # hub bottom face

# ---------------------------------------------------------------------------
# Mesh controls
# ---------------------------------------------------------------------------
H_FINE = 0.30        # target element size in the weld zone (mm)
FINE_R = 2.7         # half-width of the fine radial band around R_JOINT (mm)
FINE_Z = 4.0         # depth of the fine band below the top surface (mm)
GROWTH = 1.25        # max. geometric growth ratio outside the fine band
N_THETA = 312        # circumferential divisions (arc ~0.30 mm at R_JOINT)

# DFLUX element set: generous envelope around the heat source
DFLUX_DR = 2.7       # |r - R_JOINT| <= DFLUX_DR
DFLUX_DZ = 3.6       # z >= H_SHAFT - DFLUX_DZ


def uniform(a, b, h):
    """Uniform spacing from a to b with element size close to h."""
    n = max(1, int(round(abs(b - a) / h)))
    return np.linspace(a, b, n + 1)


def graded(a, b, h0, g):
    """Coordinates from a to b; first element (at a) ~h0, growing by <= g."""
    L = abs(b - a)
    sizes = []
    s = h0
    while sum(sizes) + s < L:
        sizes.append(s)
        s *= g
    sizes.append(s)
    sizes = np.array(sizes) * L / np.sum(sizes)    # rescale to fit exactly
    x = a + np.sign(b - a) * np.concatenate([[0.0], np.cumsum(sizes)])
    x[-1] = b
    return x


def join(*arrs):
    """Concatenate monotonic coordinate arrays, removing duplicated ends."""
    out = [arrs[0]]
    for a in arrs[1:]:
        out.append(a[1:])
    return np.concatenate(out)


# --- 1D grids ---------------------------------------------------------------
r_shaft = join(graded(R_JOINT - FINE_R, R_BORE, H_FINE, GROWTH)[::-1],
               uniform(R_JOINT - FINE_R, R_JOINT, H_FINE))
r_hub = join(uniform(R_JOINT, R_JOINT + FINE_R, H_FINE),
             graded(R_JOINT + FINE_R, R_OUT, H_FINE, GROWTH))
z_fine = uniform(H_SHAFT - FINE_Z, H_SHAFT, H_FINE)
z_hub = join(graded(H_SHAFT - FINE_Z, Z_HUB, H_FINE, GROWTH)[::-1], z_fine)
z_shaft = join(graded(Z_HUB, 0.0, z_hub[1] - z_hub[0], GROWTH)[::-1], z_hub)
theta = np.linspace(0.0, 2.0 * np.pi, N_THETA + 1)[:-1]   # periodic


def build_part(r, z, node0, elem0):
    """Structured cylindrical hex mesh. Returns nodes, elements, index map."""
    nr, nz, nt = len(r), len(z), len(theta)
    R, T, Z = np.meshgrid(r, theta, z, indexing="ij")          # (nr, nt, nz)
    ids = node0 + np.arange(nr * nt * nz).reshape(nr, nt, nz)
    xyz = np.column_stack([(R * np.cos(T)).ravel(),
                           (R * np.sin(T)).ravel(), Z.ravel()])
    i, j, k = np.meshgrid(np.arange(nr - 1), np.arange(nt),
                          np.arange(nz - 1), indexing="ij")
    i, j, k = i.ravel(), j.ravel(), k.ravel()
    jp = (j + 1) % nt
    # DC3D8 ordering: face 1-2-3-4 at k, 5-6-7-8 at k+1 (right-handed r,t,z)
    conn = np.column_stack([ids[i, j, k], ids[i + 1, j, k],
                            ids[i + 1, jp, k], ids[i, jp, k],
                            ids[i, j, k + 1], ids[i + 1, j, k + 1],
                            ids[i + 1, jp, k + 1], ids[i, jp, k + 1]])
    eids = elem0 + np.arange(len(conn))
    return ids.ravel(), xyz, eids, conn, (i, j, k)


def centroid(xyz_all, node0, conn):
    return xyz_all[conn - node0].mean(axis=1)


def write_ids(f, ids, per_line=16):
    ids = np.asarray(ids)
    for s in range(0, len(ids), per_line):
        f.write(", ".join(str(v) for v in ids[s:s + per_line]) + "\n")


def main(out="GEOMETRY_THERM.inp"):
    s_nid, s_xyz, s_eid, s_conn, (si, sj, sk) = build_part(r_shaft, z_shaft, 1, 1)
    h_node0 = s_nid[-1] + 1
    h_elem0 = s_eid[-1] + 1
    h_nid, h_xyz, h_eid, h_conn, (hi, hj, hk) = build_part(r_hub, z_hub,
                                                           h_node0, h_elem0)
    nid = np.concatenate([s_nid, h_nid])
    xyz = np.vstack([s_xyz, h_xyz])
    eid = np.concatenate([s_eid, h_eid])
    conn = np.vstack([s_conn, h_conn])

    # --- element sets ------------------------------------------------------
    c = centroid(xyz, 1, conn)
    rc = np.hypot(c[:, 0], c[:, 1])
    dflux = eid[(np.abs(rc - R_JOINT) <= DFLUX_DR) &
                (c[:, 2] >= H_SHAFT - DFLUX_DZ)]

    # --- surfaces (element-based, DC3D8 face numbering) --------------------
    # S1: z-  S2: z+  S3: theta-  S4: r+  S5: theta+  S6: r-
    s_nr, s_nz = len(r_shaft) - 1, len(z_shaft) - 1
    h_nr, h_nz = len(r_hub) - 1, len(z_hub) - 1
    z_hub_bot = z_shaft[len(z_shaft) - len(z_hub)]
    s_zc = 0.5 * (z_shaft[sk] + z_shaft[sk + 1])
    surf = {
        "SURF_SHAFT_JOINT": [(s_eid[(si == s_nr - 1) & (s_zc > z_hub_bot)], "S4")],
        "SURF_HUB_JOINT":   [(h_eid[hi == 0], "S6")],
        "SURF_FILM": [
            (s_eid[si == 0], "S6"),                                  # bore
            (s_eid[sk == 0], "S1"),                                  # shaft bottom
            (s_eid[sk == s_nz - 1], "S2"),                           # shaft top
            (s_eid[(si == s_nr - 1) & (s_zc < z_hub_bot)], "S4"),    # shaft outer
            (h_eid[hi == h_nr - 1], "S4"),                           # hub outer
            (h_eid[hk == 0], "S1"),                                  # hub bottom
            (h_eid[hk == h_nz - 1], "S2"),                           # hub top
        ],
    }

    with open(out, "w") as f:
        f.write("** Generic two-part axisymmetric weld assembly\n")
        f.write("** Generated by generate_geometry.py - DO NOT EDIT BY HAND\n")
        f.write(f"** R_BORE={R_BORE} R_JOINT={R_JOINT} R_OUT={R_OUT} "
                f"H_SHAFT={H_SHAFT} Z_HUB={Z_HUB} H_FINE={H_FINE} "
                f"N_THETA={N_THETA}\n")
        f.write("*NODE\n")
        for n, (x, y, z) in zip(nid, xyz):
            f.write(f"{n:d}, {x:.10g}, {y:.10g}, {z:.10g}\n")
        for name, ids, cn in (("SHAFT_SECTION", s_eid, s_conn),
                              ("HUB_SECTION", h_eid, h_conn)):
            f.write(f"*ELEMENT, TYPE=DC3D8, ELSET={name}\n")
            for e, row in zip(ids, cn):
                f.write(f"{e:d}, " + ", ".join(str(v) for v in row) + "\n")
        for name, rows in surf.items():
            f.write(f"*SURFACE, NAME={name}, TYPE=ELEMENT\n")
            for ids, face in rows:
                for e in ids:
                    f.write(f"{e:d}, {face}\n")
        f.write("*ELSET, ELSET=ELSET_ALL, GENERATE\n")
        f.write(f"{eid[0]}, {eid[-1]}, 1\n")
        f.write("*NSET, NSET=NSET_ALL, GENERATE\n")
        f.write(f"{nid[0]}, {nid[-1]}, 1\n")
        f.write("*ELSET, ELSET=DFLUX_ELSET\n")
        write_ids(f, dflux)

    print(f"Written {out}")
    print(f"  nodes    : {len(nid):,}")
    print(f"  elements : {len(eid):,}  (shaft {len(s_eid):,}, hub {len(h_eid):,})")
    print(f"  DFLUX_ELSET elements : {len(dflux):,}")
    print(f"  radial divisions shaft/hub : {s_nr}/{h_nr},"
          f"  axial shaft/hub : {s_nz}/{h_nz},  theta : {N_THETA}")
    print(f"  arc length at joint : {2*np.pi*R_JOINT/N_THETA:.3f} mm")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "GEOMETRY_THERM.inp")
