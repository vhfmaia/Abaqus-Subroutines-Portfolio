# Circumferential Laser Welding: Sequentially Coupled Thermo-Mechanical FEA & Phase Transformations

[![FEA: Abaqus/Standard](https://img.shields.io/badge/Solver-Abaqus%2FStandard-005a9c.svg)](https://www.3ds.com/products-services/simulia/products/abaqus/)
[![SIMULIA Compatible](https://img.shields.io/badge/SIMULIA-3DEXPERIENCE%20Compatible-005691.svg)](https://www.3ds.com/)
[![Language: Modern Fortran](https://img.shields.io/badge/Language-Fortran_2008%2F90-734f96.svg)](https://en.wikipedia.org/wiki/Fortran)
[![Automation: Pipeline](https://img.shields.io/badge/Automation-Shell%20%26%20Python-3776ab.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](../LICENSE)

An industrial-grade, benchmark finite element analysis (FEA) framework for simulating **circumferential laser welding** of axisymmetric powertrain assemblies (shaft-to-hub cylindrical joint) using Abaqus/Standard and SIMULIA 3DEXPERIENCE.

The project implements a **sequentially coupled thermo-mechanical analysis**:
1. **Thermal Stage:** Solves transient heat conduction using user subroutine `DFLUX.f` (moving conical Gaussian volumetric heat source with angular power scheduling) and metallurgical phase transformations (`ABQ_PHASE_TRANS`) in low-alloy steel (16MnCr5).
2. **Mechanical Stage:** Maps the transient thermal history onto a structural hybrid mesh (`C3D8H`), computing thermal distortions, plastic strain accumulation, continuous high-temperature annealing at $1500^\circ\text{C}$, and post-weld residual stresses.
3. **Automated Reporting:** Includes an automated pipeline (`run_pipeline.sh`, `odb_to_vtk.py`, `compile_vtk_to_gif.py`) extracting `.odb` results into VTK unstructured grids and compiling multi-viewport animated GIF visual reports.

---

## 1. Physics & Mathematical Formulation

### 1.1 Conical Gaussian Heat Source ($q(r, z)$)

The laser volumetric heat flux is formulated in a local cylindrical frame tracking the instantaneous position of the beam axis:

$$q(r, z) = Q_0 \cdot \exp\left( -3 \frac{r^2}{R_0(z)^2} \right)$$

where:
- $r = \sqrt{(X - X_c)^2 + (Y - Y_c)^2}$ is the radial distance from the beam focal center in the horizontal plane.
- $z$ is the depth measured from the top irradiated joint surface ($z = Z_0 - Z_{\text{coord}}$, with $0 \le z \le z_i$).
- $R_0(z)$ is the local conical beam radius at depth $z$, interpolating linearly between top surface radius $r_e$ and root radius $r_i$:

$$R_0(z) = r_e + (r_i - r_e) \frac{z}{z_i}$$

```
                Top Irradiated Surface (z = 0)
         |---------------- 2 * r_e ----------------|
         \                                         /
          \                                       /  |
           \                                     /   | Penetration Depth z_i
            \                                   /    |
             |---------------- 2 * r_i ----------|   v
                   Bottom Root (z = z_i)
```

### 1.2 Analytical Volume Energy Conservation

To guarantee that the integrated thermal power equals the net absorbed beam power ($P_{\text{abs}} = \eta \cdot Q_{\text{tot}}$), the peak flux $Q_0$ is derived analytically by integrating over the conical envelope:

$$\int_V q(r, z) \, dV = \int_0^{z_i} \left[ \int_0^\infty Q_0 \exp\left( -3 \frac{r^2}{R_0(z)^2} \right) 2\pi r \, dr \right] dz$$

Evaluating the radial Gaussian integral:

$$\int_0^\infty \exp\left( -3 \frac{r^2}{R_0(z)^2} \right) 2\pi r \, dr = \frac{\pi R_0(z)^2}{3}$$

Integrating along the penetration depth $z$:

$$\int_0^{z_i} \frac{\pi R_0(z)^2}{3} \, dz = \frac{\pi}{3} \int_0^{z_i} \left[ r_e + (r_i - r_e) \frac{z}{z_i} \right]^2 dz = \frac{\pi z_i}{9} \left( r_e^2 + r_e r_i + r_i^2 \right)$$

Equating to $\eta \cdot Q_{\text{tot}}$ yields the exact normalisation prefactor implemented in `DFLUX.f`:

$$Q_0 = \frac{9 \, \eta Q_{\text{tot}}}{\pi z_i \left( r_e^2 + r_e r_i + r_i^2 \right)}$$

![3D Conical Heat Source & Schedule](docs/heat_source_3d.png)

### 1.3 Kinematics & Angular Power Schedule

The beam revolves along the circular joint ($R = 15.0\text{ mm}$) at linear travel speed $v = 20.0\text{ mm/s}$ ($1200\text{ mm/min}$):

$$\omega = \frac{v}{R} = 1.3333\text{ rad/s}, \qquad \theta(t) = \omega \cdot t$$
$$X_c(t) = R \cos(\theta), \qquad Y_c(t) = R \sin(\theta)$$

To prevent hot-cracking, sudden thermal shock at beam initiation, and keyhole collapse defects at weld termination, the subroutine applies a continuous 3-stage angular power schedule:

| Segment | Angular Domain | Physical Time | Power Amplitude $\text{AMP}(\theta)$ | Purpose |
| :--- | :---: | :---: | :---: | :--- |
| **1. Ramp-Up** | $0^\circ \le \theta < 10^\circ$ | $0.000 \to 0.131\text{ s}$ | $\theta / 10^\circ$ (Linear) | Smooth keyhole initiation |
| **2. Steady Weld** | $10^\circ \le \theta \le 370^\circ$ | $0.131 \to 4.843\text{ s}$ | $1.000$ ($1080\text{ W}$ absorbed) | Full $360^\circ$ circumferential joint weld |
| **3. Ramp-Down Overlap** | $370^\circ < \theta \le 380^\circ$ | $4.843 \to 4.974\text{ s}$ | $1.0 - (\theta - 370^\circ)/10^\circ$ | Crater filling & hot-crack mitigation |
| **4. Beam Off** | $\theta > 380^\circ$ | $4.974 \to 6.000\text{ s}$ | $0.000$ | Solidification & initial cooling |

Total beam active time is $t_{\text{beam}} = 4.974\text{ s}$, transferring $E_{\text{net}} = 5,236\text{ J}$ ($261.8\text{ J/mm}$ heat input along the joint circumference).

---

## 2. Sequentially Coupled Thermo-Mechanical Architecture

```
+-----------------------------------------------------------------------------------------+
|                                1. THERMAL ANALYSIS DECK                                 |
|                                                                                         |
|  GLOBAL_THERM.inp <---+--- GEOMETRY_THERM.inp (DC3D8 Linear Hex Elements)               |
|                       +--- MATERIAL_16MnCr5.inp (k(T), cp(T), rho, ABQ_PHASE_TRANS)     |
|                       +--- DFLUX.f (Conical Gaussian Subroutine: 10° Up / 360° / 10° Dn)|
|                                                                                         |
|  Step 1: Welding & Solidification (dt = 0.012 s, CFL <= 0.80)                           |
|  Step 2: Cooling to Ambient (DELTMX = 25 °C adaptive kinetics)                          |
+--------------------------------------------+--------------------------------------------+
                                             |
                                 GLOBAL_THERM.odb / .fil
                                 (Nodal Temperatures NT11)
                                             |
                                             v
+--------------------------------------------+--------------------------------------------+
|                               2. MECHANICAL ANALYSIS DECK                               |
|                                                                                         |
|  GLOBAL_MECH.inp  <---+--- GEOMETRY_MECH.inp (C3D8H Hybrid Brick Elements, NSET_BASE)   |
|                       +--- MATERIAL_16MnCr5_MECH.inp (E(T), nu(T), alpha(T), sigma_y(T))|
|                       +--- *TEMPERATURE, FILE=GLOBAL_THERM (Mapped Thermal History)     |
|                       +--- *ANNEAL TEMPERATURE = 1500.0 °C (Resets Molten PEEQ)         |
|                       +--- *DLOAD, GRAV (Hydrostatic Fluid Stabilization)               |
|                                                                                         |
|  Step 1: Transient Thermal Stress & Distortions (NLGEOM=YES, Line Search Controls)      |
|  Step 2: Post-Weld Cooling & Locked-In Residual Stress Field                            |
+--------------------------------------------+--------------------------------------------+
                                             |
                                  GLOBAL_MECH.odb
                                             |
                                             v
+--------------------------------------------+--------------------------------------------+
|                          3. AUTOMATED POST-PROCESSING & REPORT                          |
|                                                                                         |
|  odb_to_vtk.py       -->  Extracts NT11, U, S (Von Mises), PEEQ to VTK ASCII            |
|  laser_*_vtk.zip     -->  Compacts frames for cloud HPC & archive portability           |
|  compile_vtk_to_gif.py--> Multi-viewport 3DEXPERIENCE simulation report GIF             |
+-----------------------------------------------------------------------------------------+
```

### 2.1 Numerical Stability: CFL Criterion

In moving-source thermal FEA, capturing the steep temperature gradients across the laser fusion line requires synchronizing circumferential element size with time incrementation:

$$C = \frac{v \cdot \Delta t}{\Delta h} \le 1.0$$

- Welding velocity: $v = 20.0\text{ mm/s}$.
- Circumferential element pitch: $\Delta h = 2\pi R / N_\theta = 2\pi(15.0)/312 \approx 0.302\text{ mm}$.
- CFL limit: $\Delta t_{\text{max}} = 0.302 / 20.0 = 0.0151\text{ s}$.

In `GLOBAL_THERM.inp`, Step 1 enforces a fixed time increment $\Delta t = 0.012\text{ s}$, corresponding to a Courant number of $C = 0.80$. This ensures the Gaussian heat source smoothly traverses each element without artificial peak temperature oscillations or skipped nodes.

### 2.2 Mitigation of Volumetric Locking (`C3D8H` Hybrid Elements)

At temperatures near melting ($T > 1200^\circ\text{C}$) and in the plastic regime, steel behaves near-incompressibly ($\nu \to 0.5$). In standard first-order brick elements (`C3D8`), the kinematic incompressibility constraint causes severe **volumetric locking** across the steep thermal gradient of the fusion line, leading to negative Jacobians ($\det(J) \le 0$) and premature Newton-Raphson divergence.

**Engineering Solution:**
- `GEOMETRY_MECH.inp` converts all solid elements to **`C3D8H` (8-node linear brick, hybrid with constant pressure)**. Hydrostatic pressure is treated as an independently interpolated variable, eliminating volumetric locking.
- Body-force gravity loading (`*DLOAD, GRAV`) stabilizes the hydrostatic pressure field in molten regions.

### 2.3 Continuous Annealing & Plastic Strain Reset

Molten material cannot sustain shear stress or accumulate dislocation work hardening. Abaqus resets the equivalent plastic strain (`PEEQ = 0`) when an element exceeds the annealing threshold:

```inp
*ANNEAL TEMPERATURE
 1500.0
```

To preserve robust numerical convergence across the solid-liquid transition, the mechanical material card (`MATERIAL_16MnCr5_MECH.inp`) provides a smooth, monotonic decay of yield stress down to $1.5\text{ MPa}$ at $1500^\circ\text{C}$ (avoiding numerical zero-stiffness singularities).

---

## 3. Project File Structure

```text
02_dflux_laser_welding/
├── DFLUX.f                   # User subroutine (Modern Fortran 2008/90, moving conical source)
├── DFLUX.for                 # Subroutine fallback (Fixed-form Fortran 77)
├── GLOBAL_THERM.inp          # Master thermal deck (DC3D8, 2-step welding & cooling)
├── GEOMETRY_THERM.inp        # Thermal mesh deck (DC3D8 heat transfer solid bricks)
├── MATERIAL_16MnCr5.inp      # Thermal & phase transformation material deck (ABQ_PHASE_TRANS)
├── GLOBAL_MECH.inp           # Master mechanical deck (C3D8H, sequentially coupled)
├── GEOMETRY_MECH.inp         # Mechanical mesh deck (C3D8H hybrid elements, NSET_BASE)
├── MATERIAL_16MnCr5_MECH.inp # Temperature-dependent elasto-plastic & annealing material deck
├── odb_to_vtk.py             # Abaqus Python ODB to VTK ASCII exporter & ZIP bundler
├── compile_vtk_to_gif.py     # 3DEXPERIENCE multi-viewport animated GIF report compiler
├── run_pipeline.sh           # HPC / 3DEXPERIENCE automated batch runner & pipeline
├── verify_model.py           # Verification script & Gauss energy conservation check
├── docs/
│   ├── heat_source_3d.png    # 3D conical heat flux, trajectory & power schedule plot
│   └── LINKEDIN_POST.md      # Summary technical overview
└── README.md
```

---

## 4. Verification & Energy Integration

The model integrity and energy conservation can be verified standalone without an Abaqus solver licence:

```bash
python verify_model.py
```

### Numerical Energy Conservation on 278,304-Element Mesh:

Integrating $q(r, z)$ across the 67,392 elements of `DFLUX_ELSET` using $2 \times 2 \times 2$ Gauss quadrature (identical to Abaqus `DC3D8` formulation):

| Beam Angle $\theta$ | Process Status | Theoretical Power | Integrated on FE Mesh | Ratio ($E_{\text{num}} / E_{\text{exact}}$) |
| :---: | :---: | :---: | :---: | :---: |
| $0.0^\circ$ | Start ($t=0\text{ s}$) | $0.0\text{ W}$ | $0.0\text{ W}$ | — |
| $49.4^\circ$ | Steady Weld | $1080.0\text{ W}$ | $1059.7\text{ W}$ | **$98.12\%$** |
| $140.6^\circ$ | Steady Weld | $1080.0\text{ W}$ | $1059.7\text{ W}$ | **$98.12\%$** |
| $190.0^\circ$ | Steady Weld | $1080.0\text{ W}$ | $1059.7\text{ W}$ | **$98.12\%$** |
| $269.8^\circ$ | Steady Weld | $1080.0\text{ W}$ | $1059.7\text{ W}$ | **$98.12\%$** |
| $368.6^\circ$ | Overlap Ramp-Down | $1080.0\text{ W}$ | $1059.7\text{ W}$ | **$98.12\%$** |

*Result:* Minimum element Jacobian determinant is $\det(J) = 2.873 \times 10^{-3} > 0$ (no inverted elements). Discretization error is under **$1.88\%$**, validating mesh adequacy.

---

## 5. Execution Guide (Abaqus / 3DEXPERIENCE)

### Automated Pipeline Execution

Execute the entire decoupled workflow or individual stages via `run_pipeline.sh`:

```bash
# Run full decoupled pipeline (Thermal FEA -> Thermal VTK -> Mechanical FEA -> Mechanical VTK)
./run_pipeline.sh all

# Or run stage-by-stage:
./run_pipeline.sh therm      # Runs GLOBAL_THERM.inp + DFLUX.f
./run_pipeline.sh vtk_therm  # Exports laser_therm_vtk.zip
./run_pipeline.sh mech       # Runs GLOBAL_MECH.inp
./run_pipeline.sh vtk_mech   # Exports laser_mech_vtk.zip
```

### Manual Command-Line Execution

#### 1. Thermal Analysis
```bash
abaqus job=GLOBAL_THERM input=GLOBAL_THERM.inp user=DFLUX.f cpus=4 interactive
```

#### 2. Thermal Result Extraction (VTK ZIP)
```bash
abaqus python odb_to_vtk.py GLOBAL_THERM.odb --prefix=laser_therm
```

#### 3. Mechanical Analysis (Importing Temperatures)
```bash
abaqus job=GLOBAL_MECH input=GLOBAL_MECH.inp cpus=4 interactive
```

#### 4. Mechanical Result Extraction (VTK ZIP)
```bash
abaqus python odb_to_vtk.py GLOBAL_MECH.odb --prefix=laser_mech
```

#### 5. Generate Multi-Viewport Animated Visual Report
```bash
python compile_vtk_to_gif.py laser_mech_vtk.zip --fps 15 -o laser_welding_report.gif
```

---

## 6. Key Results to Evaluate

- **Melt Pool Geometry:** Isotherm $T \ge T_{\text{liquidus}} = 1515^\circ\text{C}$ identifies the molten weld bead width and penetration depth ($z_i \approx 3.0\text{ mm}$).
- **Heat Affected Zone (HAZ):** Region bounded between $Ac_1 = 734^\circ\text{C}$ and $T_{\text{solidus}} = 1490^\circ\text{C}$.
- **Phase Distributions (`SDV4`–`SDV7`):** Martensite formation in the rapidly cooled HAZ and bainite/ferrite in adjacent parent material.
- **Residual Stress State:** Peak hoop ($\sigma_{\theta\theta}$) and axial ($\sigma_{zz}$) tensile stresses locked along the weld fusion line, balanced by compressive stress in the surrounding shaft and hub body.

---

## 7. Standards & References

1. **EN 1993-1-2:2005:** *Eurocode 3: Design of steel structures — Part 1-2: General rules — Structural fire design* (temperature-dependent thermal properties).
2. **EN 10084:** *Case hardening steels — Technical delivery conditions* (chemical composition of 16MnCr5).
3. **Andrews, K. W. (1965):** *Empirical formulae for the calculation of critical temperatures in steels*, Journal of the Iron and Steel Institute (JISI), 203, 721–727.
4. **Goldak, J., Chakravarti, A., & Bibby, M. (1984):** *A new finite element model for welding heat sources*, Metallurgical Transactions B, 15(2), 299–305.
5. **Koistinen, D. P., & Marburger, R. E. (1959):** *A general equation prescribing the extent of the austenite-martensite transformation in pure iron-carbon alloys and plain carbon steels*, Acta Metallurgica, 7(1), 59–60.

---

## 8. Disclaimer

This repository is published for professional portfolio and academic demonstration purposes. All CAD dimensions, process parameters, and metallurgical kinetics presented herein are synthetic, open-literature benchmarks designed to illustrate advanced FEA and Fortran programming methodologies, and do not represent any confidential automotive powertrain production data.

---

## 9. Author & Contact

**Victor Maia**  
- **Email:** [vhfm08@gmail.com](mailto:vhfm08@gmail.com)  
- **GitHub:** [@vhfmaia](https://github.com/vhfmaia)  
- **Specialization:** Advanced Abaqus Subroutines (`UMAT`, `VUMAT`, `DLOAD`, `DFLUX`, `DISP`, `USDFLD`, `HETVAL`) & FEA Simulation Automation
