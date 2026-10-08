"""
verify_model.py
===============
Pre-processing checks that can be run without Abaqus:

 1. Reads the heat-source parameters directly from DFLUX.for (PARAMETER
    statements), so the check always matches the compiled subroutine.
 2. Parses GEOMETRY_THERM.inp (nodes, DC3D8 elements, DFLUX_ELSET).
 3. Checks the Jacobian of every DFLUX element at the 2x2x2 Gauss points.
 4. Integrates the volumetric flux over DFLUX_ELSET with the same 2x2x2
    Gauss rule used by Abaqus for DC3D8 and compares it with ETA*QTOT
    at several beam positions -> mesh-discretisation energy error.
 5. Optionally (--plot) writes docs/heat_source.png.

Usage:  python verify_model.py [GEOMETRY_THERM.inp] [--plot]
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# 1. DFLUX parameters
# ---------------------------------------------------------------------------
def read_dflux_parameters(path=None):
    if path is None:
        path = HERE / "DFLUX.f" if (HERE / "DFLUX.f").exists() else HERE / "DFLUX.for"
    p = {
        "QTOT": 1.8e6, "ETA": 0.60, "R_E": 1.0, "R_I": 0.75, "ZI": 3.0,
        "RADIUS": 15.0, "VBEAM": 1200.0, "RAMP_U_DEG": 10.0, "WELD_DEG": 360.0,
        "RAMP_D_DEG": 10.0, "Z0": 20.0, "TOL": 1e-8
    }
    pat_f77 = re.compile(r"^\s+PARAMETER\s*\(\s*(\w+)\s*=\s*([-+0-9.EeDd]+)\s*\)")
    pat_f90 = re.compile(r"^\s*real\([^)]+\),\s*parameter\s*::\s*(\w+)\s*=\s*([-+0-9.EeDd]+)", re.IGNORECASE)
    
    mapping = {
        "QTOT_NOMINAL": "QTOT", "ETA_ABSORB": "ETA", "R_TOP": "R_E",
        "R_BOTTOM": "R_I", "DEPTH_ZI": "ZI", "WELD_RADIUS": "RADIUS",
        "TRAVEL_SPEED": "VBEAM", "DEG_RAMP_UP": "RAMP_U_DEG",
        "DEG_WELD": "WELD_DEG", "DEG_RAMP_DOWN": "RAMP_D_DEG",
        "Z_SURFACE": "Z0", "CUTOFF_TOL": "TOL"
    }
    
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    for line in text.splitlines():
        m = pat_f77.match(line)
        if m:
            p[m.group(1).upper()] = float(m.group(2).replace("D", "E").replace("d", "e"))
        m90 = pat_f90.match(line)
        if m90:
            var_name = m90.group(1).upper()
            val = float(m90.group(2).replace("D", "E").replace("d", "e"))
            if var_name in mapping:
                p[mapping[var_name]] = val
            else:
                p[var_name] = val

    p["OMEGA"] = (p["VBEAM"] / 60.0) / p["RADIUS"]
    e3 = math.exp(3.0)
    p["E3_RATIO"] = e3 / (e3 - 1.0)
    p["DENOM"] = p["ZI"] * (p["R_I"]**2 + p["R_E"] * p["R_I"] + p["R_E"]**2)
    p["PREFAC"] = (9.0 * p["ETA"] * p["QTOT"] * p["E3_RATIO"]) / (math.pi * p["DENOM"])
    p["TOTDEG"] = p["RAMP_U_DEG"] + p["WELD_DEG"] + p["RAMP_D_DEG"]
    return p


def flux(P, xyz, t):
    """Python port of DFLUX.for (vectorised). xyz: (..., 3) array."""
    th = P["OMEGA"] * t
    thd = math.degrees(th)
    if thd > P["TOTDEG"]:
        return np.zeros(xyz.shape[:-1])
    if thd < P["RAMP_U_DEG"]:
        amp = thd / P["RAMP_U_DEG"]
    elif thd < P["RAMP_U_DEG"] + P["WELD_DEG"]:
        amp = 1.0
    elif P["RAMP_D_DEG"] > 0:
        amp = (P["TOTDEG"] - thd) / P["RAMP_D_DEG"]
    else:
        amp = 0.0
    zr = np.clip(P["Z0"] - xyz[..., 2], 0.0, None)
    r0 = P["R_E"] + (P["R_I"] - P["R_E"]) * zr / P["ZI"]
    rl = np.hypot(xyz[..., 0] - P["RADIUS"] * math.cos(th),
                  xyz[..., 1] - P["RADIUS"] * math.sin(th))
    rcut = r0
    q = P["PREFAC"] * np.exp(-3.0 * rl**2 / r0**2) * amp
    q[(zr > P["ZI"]) | (rl > rcut)] = 0.0
    return q


# ---------------------------------------------------------------------------
# 2. Mesh parser (only what is needed)
# ---------------------------------------------------------------------------
def read_mesh(path):
    nodes, elems, elsets = {}, {}, {}
    mode, name, gen = None, None, False
    with open(path) as f:
        for line in f:
            if line.startswith("**"):
                continue
            if line.startswith("*"):
                kw = line.upper().replace(" ", "")
                mode = None
                if kw.startswith("*NODE") and not kw.startswith("*NODEOUTPUT"):
                    mode = "n"
                elif kw.startswith("*ELEMENT,"):
                    mode = "e"
                    name = re.search(r"ELSET=([^,\n]+)", kw).group(1)
                    elsets.setdefault(name, [])
                elif kw.startswith("*ELSET"):
                    mode = "s"
                    name = re.search(r"ELSET=([^,\n]+)", kw).group(1)
                    gen = "GENERATE" in kw
                    elsets.setdefault(name, [])
                continue
            v = line.replace(",", " ").split()
            if not v:
                continue
            if mode == "n":
                nodes[int(v[0])] = tuple(map(float, v[1:4]))
            elif mode == "e":
                elems[int(v[0])] = list(map(int, v[1:9]))
                elsets[name].append(int(v[0]))
            elif mode == "s":
                if gen:
                    a, b = int(v[0]), int(v[1])
                    s = int(v[2]) if len(v) > 2 else 1
                    elsets[name].extend(range(a, b + 1, s))
                else:
                    elsets[name].extend(map(int, v))
    return nodes, elems, elsets


# ---------------------------------------------------------------------------
# 3-4. Gauss integration on DC3D8
# ---------------------------------------------------------------------------
SG = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
               [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)


def shape(xi):
    N = 0.125 * np.prod(1 + SG * xi, axis=1)
    dN = np.empty((8, 3))
    for a in range(3):
        o = [b for b in range(3) if b != a]
        dN[:, a] = 0.125 * SG[:, a] * (1 + SG[:, o[0]] * xi[o[0]]) * \
            (1 + SG[:, o[1]] * xi[o[1]])
    return N, dN


def gauss_points(X):
    """X: (ne, 8, 3) -> integration point coords (ne, 8, 3), weights (ne, 8)."""
    g = 1.0 / math.sqrt(3.0)
    ip = np.empty((len(X), 8, 3))
    w = np.empty((len(X), 8))
    for q, xi in enumerate(SG * g):
        N, dN = shape(xi)
        ip[:, q] = np.einsum("i,eij->ej", N, X)
        w[:, q] = np.linalg.det(np.einsum("ia,eib->eab", dN, X))
    return ip, w


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    mesh = args[0] if args else str(HERE / "GEOMETRY_THERM.inp")
    P = read_dflux_parameters()
    print("DFLUX parameters read from DFLUX.for:")
    for k in ("QTOT", "ETA", "R_E", "R_I", "ZI", "RADIUS", "VBEAM",
              "RAMP_U_DEG", "WELD_DEG", "RAMP_D_DEG", "Z0", "TOL"):
        print(f"  {k:<11}= {P[k]:g}")
    t_on = math.radians(P["TOTDEG"]) / P["OMEGA"]
    print(f"  -> absorbed power ETA*QTOT = {P['ETA']*P['QTOT']/1e3:.1f} W")
    print(f"  -> peak flux Q0            = {P['PREFAC']:.4e} mW/mm3")
    print(f"  -> angular speed           = {P['OMEGA']:.4f} rad/s")
    print(f"  -> beam-on time            = {t_on:.4f} s ({P['TOTDEG']:g} deg)")

    nodes, elems, elsets = read_mesh(mesh)
    ids = elsets["DFLUX_ELSET"]
    X = np.array([[nodes[n] for n in elems[e]] for e in ids])
    ip, w = gauss_points(X)
    print(f"\nMesh: {len(nodes):,} nodes, {len(elems):,} elements, "
          f"DFLUX_ELSET: {len(ids):,} elements")
    print(f"  min detJ in DFLUX_ELSET : {w.min():.3e} "
          f"({'OK' if w.min() > 0 else 'NEGATIVE JACOBIAN!'})")
    rn = np.hypot(X[..., 0], X[..., 1])
    print(f"  DFLUX_ELSET extent r=[{rn.min():.2f}, {rn.max():.2f}] "
          f"z=[{X[..., 2].min():.2f}, {X[..., 2].max():.2f}] mm")
    rcut = math.sqrt(-math.log(P["TOL"]) / 3.0) * max(P["R_E"], P["R_I"])
    need = (P["RADIUS"] - rcut, P["RADIUS"] + rcut, P["Z0"] - P["ZI"], P["Z0"])
    ok = (rn.min() <= need[0] and rn.max() >= need[1] and
          X[..., 2].min() <= need[2] and X[..., 2].max() >= need[3] - 1e-9)
    print(f"  required envelope       r=[{need[0]:.2f}, {need[1]:.2f}] "
          f"z=[{need[2]:.2f}, {need[3]:.2f}] mm -> {'OK' if ok else 'TOO SMALL'}")
    print(f"  integration points above Z0 : {(ip[..., 2] > P['Z0'] + 1e-9).sum()}")

    print("\nEnergy check (2x2x2 Gauss, as Abaqus DC3D8):")
    print(f"  {'t [s]':>7} {'theta [deg]':>11} {'amp':>5} "
          f"{'integrated [W]':>15} {'target [W]':>11} {'ratio':>7}")
    for frac in (0.0, 0.13, 0.37, 0.5, 0.71, 0.9, 0.97):
        t = frac * t_on
        thd = math.degrees(P["OMEGA"] * t)
        E = (flux(P, ip, t) * w).sum() / 1e3
        q_ref = flux(P, np.array([[P["RADIUS"] * math.cos(P["OMEGA"] * t),
                                   P["RADIUS"] * math.sin(P["OMEGA"] * t),
                                   P["Z0"]]]), t)[0]
        amp = q_ref / P["PREFAC"]
        tgt = P["ETA"] * P["QTOT"] / 1e3 * amp
        ratio = E / tgt if tgt > 0 else float("nan")
        print(f"  {t:7.3f} {thd:11.2f} {amp:5.2f} {E:15.2f} {tgt:11.2f} {ratio:7.4f}")

    if "--plot" in sys.argv:
        plot(P)


def plot(P):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = HERE / "docs"
    out.mkdir(exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    # (a) cross-section through the beam axis at t=0: x = RADIUS + s, z
    s = np.linspace(-2.0, 2.0, 401)
    z = np.linspace(P["Z0"] - P["ZI"] - 0.5, P["Z0"], 301)
    S, Z = np.meshgrid(s, z)
    xyz = np.stack([P["RADIUS"] + S, np.zeros_like(S), Z], axis=-1)
    q = flux(P, xyz, 0.0)
    cs = ax[0].contourf(S, P["Z0"] - Z, q / 1e3, levels=30, cmap="inferno")
    ax[0].invert_yaxis()
    ax[0].set_xlabel("radial offset from beam axis [mm]")
    ax[0].set_ylabel("depth below top surface [mm]")
    ax[0].set_title("Conical Gaussian source, q [W/mm$^3$]")
    fig.colorbar(cs, ax=ax[0])
    # (b) power amplitude vs travelled angle
    t_on = math.radians(P["TOTDEG"]) / P["OMEGA"]
    tt = np.linspace(0, 1.15 * t_on, 600)
    amp = []
    for t in tt:
        p0 = np.array([[P["RADIUS"] * math.cos(P["OMEGA"] * t),
                        P["RADIUS"] * math.sin(P["OMEGA"] * t), P["Z0"]]])
        amp.append(flux(P, p0, t)[0] / P["PREFAC"])
    ax[1].plot(tt, np.array(amp) * P["ETA"] * P["QTOT"] / 1e3, lw=2)
    ax[1].set_xlabel("step time [s]")
    ax[1].set_ylabel("absorbed power [W]")
    ax[1].set_title("Power schedule (ramp-up / weld / ramp-down)")
    ax[1].grid(alpha=0.3)
    sec = ax[1].secondary_xaxis(
        "top", functions=(lambda t: np.degrees(P["OMEGA"] * t),
                          lambda d: np.radians(d) / P["OMEGA"]))
    sec.set_xlabel("beam angle [deg]")
    fig.tight_layout()
    fig.savefig(out / "heat_source.png", dpi=150)
    print(f"\nWritten {out / 'heat_source.png'}")


if __name__ == "__main__":
    main()
