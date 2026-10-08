# Advanced Abaqus User Subroutines Portfolio

A repository of production-grade user subroutines for **Abaqus/Standard** and **3DEXPERIENCE SIMULIA**, written in **Fortran** for high-performance structural, fatigue, and manufacturing simulations.

---

## Portfolio Projects & Subroutine Implementations

| ID | Subroutine | Physics / Type | Engineering Application | Status |
| :---: | :--- | :--- | :--- | :---: |
| **01** | [`DLOAD`](01_dload_gear_fatigue/README.md) | Moving Surface Load / Quasi-Static | High-cycle fatigue analysis of a welded helical crown gear (360° full revolution) | ✅ Verified |
| **02** | [`DFLUX`](02_dflux_laser_welding/README.md) | Moving Volumetric Heat Flux / Transient Thermal | Circumferential laser welding with solid-state phase transformations (16MnCr5) | ✅ Verified |

---

### Project Showcase: 01. Helical Crown Gear Fatigue (`DLOAD`)

[![FEA Abaqus](https://img.shields.io/badge/FEA-Abaqus%2FStandard-blue.svg)](01_dload_gear_fatigue/README.md)
[![SIMULIA Compatible](https://img.shields.io/badge/SIMULIA-3DEXPERIENCE%20Compatible-005691.svg)](01_dload_gear_fatigue/README.md)
[![Fortran](https://img.shields.io/badge/Language-Fortran-734f96.svg)](01_dload_gear_fatigue/README.md)

<p align="center">
  <a href="01_dload_gear_fatigue/README.md">
    <img src="01_dload_gear_fatigue/gear_torque_simulation.gif" alt="Abaqus DLOAD Helical Gear Rolling Contact Simulation (360 Degrees)" width="100%" />
  </a>
</p>

*Analytical multi-viewport simulation across all 60 teeth of a heavy-duty automotive helical crown gear (T = 270 N·m, z = 60, β = 25°), featuring 3D European isometric projection, Y-axis frontal view, and transverse 3D RZ section.*

#### Key Simulation & Fatigue Results:
- **Peak Tensile Bending:** $\sigma_{\max} = 68.13\text{ MPa}$ at **Frame 55** ($t = 1.100\text{ s}$, Increment 550) in tooth root fillet.
- **Out-of-Mesh Valley:** $\sigma_{\min} = 0.00\text{ MPa}$ at **Frame 30** ($t = 0.600\text{ s}$, Increment 300), $180^\circ$ out of phase.
- **Cyclic Stress Amplitude:** $\sigma_a = 34.07\text{ MPa}$ (pulsating ratio $R = 0$, Goodman equivalent $\sigma_{a,\text{eq}} = 35.50\text{ MPa}$).
- **Fatigue Life Regime:** **Infinite Life / High-Cycle Fatigue ($N > 10^7$ cycles)** with structural safety factor $SF_F \approx 7.0$ ($S_e \approx 250\text{ MPa}$).

Detailed formulation, kinematic equations, FEA verification, and post-processing scripts are available in [01_dload_gear_fatigue/](01_dload_gear_fatigue/README.md).

---

### Project Showcase: 02. Circumferential Laser Welding (`DFLUX`)

[![FEA Abaqus](https://img.shields.io/badge/FEA-Abaqus%2FStandard-blue.svg)](02_dflux_laser_welding/README.md)
[![SIMULIA Compatible](https://img.shields.io/badge/SIMULIA-3DEXPERIENCE%20Compatible-005691.svg)](02_dflux_laser_welding/README.md)
[![Fortran](https://img.shields.io/badge/Language-Fortran_77-734f96.svg)](02_dflux_laser_welding/README.md)

<p align="center">
  <a href="02_dflux_laser_welding/README.md">
    <img src="02_dflux_laser_welding/docs/heat_source.png" alt="Abaqus DFLUX Conical Gaussian Heat Source and Welding Power Schedule" width="100%" />
  </a>
</p>

*Moving conical Gaussian volumetric heat source travelling along a circular joint (R = 15 mm, P_abs = 1080 W, v = 1200 mm/min) coupled with the ABQ_PHASE_TRANS metallurgical framework for 16MnCr5 case-hardening steel.*

#### Key Thermal & Metallurgical Highlights:
- **Exact Analytical Normalization:** Peak flux $Q_0$ derived via volume integration, guaranteeing strict conservation of absorbed beam energy ($\eta \cdot Q_{\text{tot}} = 1080.0\text{ W}$).
- **Pre-Solver Numerical Verification:** $2 \times 2 \times 2$ Gauss quadrature integration across the 67,392 elements of the irradiated domain matches theoretical power with $98.12\%$ energy conservation.
- **Microstructural State Tracking:** Couples liquidus/solidus melt pool kinetics with solid-state phase decomposition (Ferrite-Pearlite, Austenite, Bainite, Martensite via JMA and Koistinen-Marburger kinetics).
- **Brittleness Temperature Range (BTR) Resolution:** Dedicated two-stage time-stepping strategy preserving temporal resolution through the solidification range to assess hot-cracking susceptibility.

Detailed formulation, verification routines, and input decks are available in [02_dflux_laser_welding/](02_dflux_laser_welding/README.md).

---

## Directory Layout

```text
.
├── 01_dload_gear_fatigue/
│   ├── DLOAD_CROWN_HELICAL.f      # User Subroutine (Abaqus DLOAD + UEXTERNALDB)
│   ├── Gear_torque.inp            # Master Input Deck
│   ├── GEOMETRY_GEAR.inp          # 3D Finite Element Mesh & Topology
│   ├── gear_torque_simulation.gif # Full 360-degree multi-viewport animation
│   ├── compile_vtk_to_gif.py      # Multi-viewport visualization compiler
│   ├── odb_to_vtk.py              # Abaqus Python ODB to VTK extractor
│   ├── run_pipeline.sh            # Headless Linux / HPC batch execution script
│   └── README.md                  # Comprehensive Technical Report
├── 02_dflux_laser_welding/
│   ├── DFLUX.for                  # User Subroutine (Conical Gaussian Heat Source)
│   ├── GLOBAL_THERM.inp           # Master Thermal Analysis Deck
│   ├── MATERIAL_16MnCr5.inp       # Phase Transformation Material Properties
│   ├── GEOMETRY_THERM.inp         # High-Density Hexahedral Mesh (DC3D8)
│   ├── generate_geometry.py       # Parametric Hex Mesh Generator
│   ├── generate_material.py       # Thermophysical & Metallurgical Generator
│   ├── verify_model.py            # Abaqus-Free Energy Conservation Verifier
│   ├── docs/                      # Heat Source Verification Figures
│   └── README.md                  # Comprehensive Technical Report
├── LICENSE
└── README.md
```

---

## Author & Contact

**Victor Maia**  
- **Email:** [vhfm08@gmail.com](mailto:vhfm08@gmail.com)  
- **GitHub:** [@vhfmaia](https://github.com/vhfmaia)  
- **Specialization:** Advanced Abaqus User Subroutines (`UMAT`, `VUMAT`, `DLOAD`, `DFLUX`, `DISP`, `USDFLD`, `HETVAL`) & FEA Simulation Automation
