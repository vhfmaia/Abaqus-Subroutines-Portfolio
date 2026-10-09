# Advanced Abaqus User Subroutines Portfolio

<div align="justify">
A repository of production-grade user subroutines for <b>Abaqus/Standard</b> and <b>3DEXPERIENCE SIMULIA</b>, written in <b>Fortran</b> for high-performance structural, fatigue, and manufacturing simulations.
</div>

---

## Portfolio Projects & Subroutine Implementations

| ID | Subroutine | Physics / Type | Engineering Application | Status |
| :---: | :--- | :--- | :--- | :---: |
| **01** | [`DLOAD`](01_dload_gear_fatigue/README.md) | Moving Surface Load / Quasi-Static | High-cycle fatigue analysis of a welded helical crown gear (360° full revolution) | ✅ Verified |
| **02** | [`DFLUX`](02_dflux_laser_welding/README.md) | Moving Conical Flux / Thermo-Mechanical | Sequentially coupled laser welding with phase transformations and residual stress | ✅ Verified |
| **03** | [`UINTER`](03_uinter_welding_interface/README.md) | Surface Interaction / BTR Hot Cracking | 3-State contact interface (Raw $\rightarrow$ Welded $\rightarrow$ Cracked) with BTR criteria | ⚠️ Exploratory (R&D Backlog) |
| **04** | [`UMAT`](04_umat_cohesive_welding/README.md) | Cohesive Zone Modeling / Fracture Energy $G_c$ | Mesh-independent 3-state thermo-mechanical weld interface (`COH3D8`) | ⚠️ In Work (Testing Phase) |

---

### Project Showcase: 01. Helical Crown Gear Fatigue (`DLOAD`)

[![FEA Abaqus](https://img.shields.io/badge/FEA-Abaqus%2FStandard-blue.svg)](01_dload_gear_fatigue/README.md)
[![SIMULIA Compatible](https://img.shields.io/badge/SIMULIA-3DEXPERIENCE%20Compatible-005691.svg)](01_dload_gear_fatigue/README.md)
[![Fortran](https://img.shields.io/badge/Language-Fortran_2008%2F90-734f96.svg)](01_dload_gear_fatigue/README.md)

<p align="center">
  <a href="01_dload_gear_fatigue/README.md">
    <img src="01_dload_gear_fatigue/gear_torque_simulation.gif" alt="Abaqus DLOAD Helical Gear Rolling Contact Simulation (360 Degrees)" width="100%" />
  </a>
</p>

<div align="justify">
Analytical multi-viewport simulation across all 60 teeth of a heavy-duty automotive helical crown gear ($T = 270\text{ N}\cdot\text{m}$, $z = 60$, $\beta = 25^\circ$), featuring 3D European isometric projection, Y-axis frontal view, and transverse 3D RZ section.
</div>

#### Key Simulation & Fatigue Results:
- **Peak Tensile Bending:** $\sigma_{\max} = 68.13\text{ MPa}$ at **Frame 55** ($t = 1.100\text{ s}$, Increment 550) in tooth root fillet.
- **Out-of-Mesh Valley:** $\sigma_{\min} = 0.00\text{ MPa}$ at **Frame 30** ($t = 0.600\text{ s}$, Increment 300), $180^\circ$ out of phase.
- **Cyclic Stress Amplitude:** $\sigma_a = 34.07\text{ MPa}$ (pulsating ratio $R = 0$, Goodman equivalent $\sigma_{a,\text{eq}} = 35.50\text{ MPa}$).
- **Fatigue Life Regime:** **Infinite Life / High-Cycle Fatigue** ($N > 10^7\text{ cycles}$) with structural safety factor $SF_F \approx 7.0$ ($S_e \approx 250\text{ MPa}$).

<div align="justify">
Detailed formulation, kinematic equations, FEA verification, and post-processing scripts are available in <a href="01_dload_gear_fatigue/README.md">01_dload_gear_fatigue/</a>.
</div>

---

### Project Showcase: 02. Circumferential Laser Welding (`DFLUX`)

[![FEA Abaqus](https://img.shields.io/badge/FEA-Abaqus%2FStandard-blue.svg)](02_dflux_laser_welding/README.md)
[![SIMULIA Compatible](https://img.shields.io/badge/SIMULIA-3DEXPERIENCE%20Compatible-005691.svg)](02_dflux_laser_welding/README.md)
[![Fortran](https://img.shields.io/badge/Language-Fortran_2008%2F90-734f96.svg)](02_dflux_laser_welding/README.md)

<p align="center">
  <a href="02_dflux_laser_welding/README.md">
    <img src="02_dflux_laser_welding/conical_heat_source_improved_plot.svg" alt="Abaqus DFLUX 3D Conical Gaussian Heat Source and Welding Power Schedule" width="100%" />
  </a>
</p>

<div align="justify">
Sequentially coupled thermo-mechanical analysis of circumferential laser welding ($R = 15\text{ mm}$, $P_{\text{abs}} = 1080\text{ W}$, $v = 1200\text{ mm/min}$) with moving conical Gaussian heat source, metallurgical phase transformations (<code>ABQ_PHASE_TRANS</code>), hybrid elements (<code>C3D8H</code>), continuous annealing at $1530^\circ\text{C}$, and residual stress prediction for 16MnCr5 case-hardening steel.
</div>

#### Key Thermal, Metallurgical & Mechanical Highlights:
- **3D Conical Gaussian Heat Source:** Analytically normalized peak flux $Q_0 = 4.6935 \times 10^5\text{ mW/mm}^3$ with $10^\circ$ linear power ramp-up, $360^\circ$ steady weld, and $10^\circ$ crater-filling ramp-down overlap.
- **Unified Material Architecture:** Single `Material_16MnCr5.inp` shared by both thermal and mechanical simulations, ensuring seamless data consistency.
- **Metallurgical Kinetics & Solidification (BTR):** Tracks liquidus/solidus melt pool kinetics and solid-state phase decomposition (Ferrite-Pearlite, Austenite, Bainite, Martensite via JMA and Koistinen-Marburger kinetics).
- **Sequentially Coupled Structural Analysis:** Mitigates high-temperature volumetric locking via hybrid formulation elements (`C3D8H`), hydrostatic gravity stabilization, and line search non-linear controls.

<div align="justify">
Detailed formulation, verification routines, simulation input decks, and automated post-processing pipeline are available in <a href="02_dflux_laser_welding/README.md">02_dflux_laser_welding/</a>.
</div>

---

### Project Showcase: 03 & 04. Advanced Interface Modeling Evolution

<div align="justify">
The portfolio demonstrates a complete engineering R&D cycle comparing two distinct methodologies to capture the 3 physical states of a welded interface (<b>Raw $\rightarrow$ Welded $\rightarrow$ Cracked</b>):
</div>

```text
┌──────────────────────────────────────────────┐       Identified Bottlenecks:       ┌──────────────────────────────────────────────┐
│  03_uinter_welding_interface (UINTER)        │   • Pathological Mesh Dependency    │  04_umat_cohesive_welding (UMAT + CZM)       │
│  • Surface Interaction Routine               │ ─────────────────────────────────>  │  • Cohesive Zone Elements (COH3D8)           │
│  • Contact Master/Slave formulation          │   • Severe Discontinuity (SDI)      │  • Fracture Energy Gc Regularization         │
│  • Status: ⚠️ Exploratory (R&D Backlog)      │   • MPI Domain Split Conflicts      │  • Status: ⚠️ In Work (Testing Phase)        │
└──────────────────────────────────────────────┘                                     └──────────────────────────────────────────────┘
```

<div align="justify">
<ul>
  <li><b><a href="03_uinter_welding_interface/README.md">03. Exploratory Interface (UINTER)</a>:</b> Formulates the 3-state state machine on contact surfaces. Benchmarking identified critical mesh-dependency due to displacement jumps vs strain rates, prompting its classification as an exploratory R&D backlog module.</li>
  <li><b><a href="04_umat_cohesive_welding/README.md">04. Cohesive Zone Interface (UMAT)</a>:</b> Under active development and testing to evaluate mesh-independent fracture energy $G_c$ regularization. Uses an exact $3 \times 3$ analytical $\mathbf{DDSDDE}$ Jacobian to investigate convergence without contact chattering or SDI cutbacks (strictly a research prototype under evaluation, not for production use).</li>
</ul>
</div>

---

## Directory Layout

```text
.
├── 01_dload_gear_fatigue/
│   ├── DLOAD_CROWN_HELICAL.f                  # User Subroutine (Modern Fortran DLOAD + UEXTERNALDB)
│   ├── Gear_torque.inp                        # Master Input Deck
│   ├── GEOMETRY_GEAR.inp                      # 3D Finite Element Mesh & Topology
│   ├── gear_torque_simulation.gif             # Full 360-degree multi-viewport animation (11 FPS)
│   ├── compile_vtk_to_gif.py                  # Multi-viewport visualization compiler
│   ├── odb_to_vtk.py                          # Abaqus Python ODB to VTK extractor
│   ├── run_pipeline.sh                        # Headless Linux / HPC batch execution script
│   └── README.md                              # Comprehensive Technical Report
├── 02_dflux_laser_welding/
│   ├── Disk_heatsource_TH.inp                 # Master Thermal Analysis Deck (DC3D8)
│   ├── Disk_heatsource_ME.inp                 # Master Mechanical Analysis Deck (C3D8H, Sequentially Coupled)
│   ├── Geometry_TH.inp                        # Thermal Mesh Deck (DC3D8 Solid Bricks)
│   ├── Geometry_ME.inp                        # Mechanical Hybrid Mesh Deck (C3D8H, NSET_BASE)
│   ├── Material_16MnCr5.inp                   # Unified Thermo-Elasto-Plastic & Phase Transformation Deck
│   ├── dflux_disk_conical_gaussian.f          # User Subroutine (Modern Fortran TDC Conical Model)
│   ├── conical_heat_source_improved_plot.svg  # 4-Panel 3D Conical Heat Source & Trajectory Plot (SVG)
│   ├── sequentially_coupled_architecture.svg  # Multi-Physics Data Flow & Pipeline Architecture (SVG)
│   ├── odb_to_vtk.py                          # Abaqus Python ODB to VTK Extractor & ZIP Packager
│   ├── compile_vtk_to_gif.py                  # Multi-Viewport Animated Report Compiler
│   ├── run_pipeline.sh                        # Headless Linux / HPC Batch Execution Pipeline
│   └── README.md                              # Comprehensive Technical Report
├── 03_uinter_welding_interface/
│   ├── uinter_btr_hot_cracking.f              # User Subroutine (UINTER State Machine in Fortran)
│   ├── dummy_btr_contact_patch.inp            # 2-Element Contact Patch Benchmark Deck
│   └── README.md                              # Technical Report & R&D Backlog Documentation
├── 04_umat_cohesive_welding/
│   ├── umat_cohesive_welding_btr.f            # Prototype User Material Subroutine (Modern Fortran)
│   ├── cohesive_weld_btr_verification.inp     # 3-Element Cohesive Zone Benchmark Deck (COH3D8)
│   └── README.md                              # Comprehensive Technical Report & Formulation
├── LICENSE
└── README.md
```

---

## Disclaimer

<div align="justify">
This repository is published for professional portfolio and academic demonstration purposes. All CAD dimensions, process parameters, and metallurgical kinetics presented herein are synthetic, open-literature benchmarks designed to illustrate advanced FEA and Fortran programming methodologies, and do not represent any confidential data.
</div>

---

## Author & Contact

**Victor Maia**  
- **Email:** [vhfm08@gmail.com](mailto:vhfm08@gmail.com)  
- **GitHub:** [@vhfmaia](https://github.com/vhfmaia)  
- **Specialization:** Advanced Abaqus User Subroutines (`UMAT`, `VUMAT`, `DLOAD`, `DFLUX`, `UINTER`, `DISP`, `USDFLD`, `HETVAL`) & FEA Simulation Automation
