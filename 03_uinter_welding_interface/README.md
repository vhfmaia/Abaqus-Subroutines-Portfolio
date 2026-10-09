# 03. 3-State Contact Interface & BTR Hot Cracking (`UINTER`)

[![Status](https://img.shields.io/badge/Status-Exploratory%20%2F%20In%20Work%20(R%26D%20Backlog)-orange.svg)](README.md)
[![Abaqus](https://img.shields.io/badge/Abaqus%2FStandard-User%20Interaction%20(UINTER)-blue.svg)](README.md)
[![Language](https://img.shields.io/badge/Language-Fortran_90%2F2008-734f96.svg)](README.md)
[![Next-Gen Succeeded](https://img.shields.io/badge/Succeeded%20By-04__umat__cohesive__welding-success.svg)](../04_umat_cohesive_welding/README.md)

---

## Executive Summary & Engineering Status

> [!WARNING]
> **DEVELOPMENT STATUS: EXPLORATORY RESEARCH / IN WORK (R&D BACKLOG)**  
> This directory documents an exploratory implementation of an Abaqus/Standard **User Interaction Subroutine (`UINTER`)** designed to model the progressive multi-stage constitutive behavior of an industrial welded joint (**Raw $\rightarrow$ Welded $\rightarrow$ Cracked**).  
>
> While conceptually elegant, **rigorous finite-element benchmarking revealed fundamental numerical bottlenecks**—primarily **pathological mesh dependency** and **severe discontinuity iterations (SDI cutbacks)**. Consequently, this interface formulation is preserved here for research reference and is **officially succeeded by the regularized Cohesive Zone UMAT architecture** documented in [**`04_umat_cohesive_welding`**](../04_umat_cohesive_welding/README.md).

---

## 1. Physical Motivation: The 3-State Contact Interface

In conventional finite element welding analyses, joint interfaces are typically simplified via static `*TIE` constraints (as in baseline sequential thermal-mechanical pipelines). However, a `*TIE` constraint glues the entire 360° circumference from time $t = 0$, completely ignoring:
1. **Pre-weld joint gaps and press-fit mechanics** ahead of the advancing beam;
2. **Semi-solid mushy zone tearing** during cooling through the **Brittleness Temperature Range (BTR)**;
3. **Permanent weld seam rupture or hot cracking**.

To bridge this gap, an interaction model was formulated to dynamically transition between three discrete constitutive states:

```text
                     ┌───────────────────────────────────────┐
                     │         STATE 0: RAW (Virgin)         │
                     │  - Pure unilateral penalty contact    │
                     │  - Free tensile opening (Pn = 0)      │
                     │  - Coulomb friction (µ ≈ 0.15)        │
                     └───────────────────────────────────────┘
                                         │
                                         │ Beam arrival: T >= T_liquidus (1530 °C)
                                         ▼
                     ┌───────────────────────────────────────┐
                     │        BTR MUSH-ZONE SOLIDIFICATION   │
                     │         [T_solidus, T_liquidus]       │
                     │   Tracks normal opening rate Δu_n     │
                     └───────────────────────────────────────┘
                                    /         \
          Won criterion NOT violated           Won criterion violated
          (Δu_n / L0 < ε_crit)                (Δu_n / L0 >= ε_crit)
                        /                           \
                       ▼                             ▼
       ┌───────────────────────────────┐   ┌───────────────────────────────┐
       |    STATE 2: WELDED (Sound)    |   |    STATE 3: CRACKED (Tear)    |
       | - Bilateral stiffness Kn, Ks  |   | - Permanent loss of cohesion  |
       | - Transmits tension & shear   |   |   in tension (Pn = 0, un > 0) |
       | - Full continuum continuity   |   | - Compressive contact only    |
       └───────────────────────────────┘   └───────────────────────────────┘
```

---

## 2. Constitutive Formulation (`UINTER`)

The subroutine is called at every contact slave integration point. It receives the relative separation vector $\mathbf{u} = \{u_n, u_{s1}, u_{s2}\}^T$, temperature $T$, and maintains 5 internal state variables (`STATEV`):

### 2.1 State Machine Governing Equations

1. **State 0 — RAW (Press-fit Unwelded):**
   $$\begin{cases} P_n = K_{\text{comp}} \cdot u_n, \quad |\boldsymbol{\tau}| \le \mu |P_n| & \text{if } u_n < 0 \text{ (compression)} \\ P_n = 0, \quad \boldsymbol{\tau} = \mathbf{0} & \text{if } u_n \ge 0 \text{ (clearance/separation)} \end{cases}$$

2. **State 1 — BTR Mushy Zone Solidification ($T_{\text{sol}} \le T \le T_{\text{liq}}$):**
   * Solidus: $T_{\text{solidus}} = 1485.0\,^\circ\text{C}$
   * Liquidus: $T_{\text{liquidus}} = 1530.0\,^\circ\text{C}$
   * Critical hot tearing index based on the Won / Prokhorov criterion:
     $$\varepsilon_{\text{BTR}} = \frac{\Delta u_{n,\text{accum}}}{L_0}, \quad I_{\text{Won}} = \frac{\varepsilon_{\text{BTR}}}{\varepsilon_{\text{crit}}}$$
   * If $I_{\text{Won}} \ge 1.0 \implies$ Transitions to **State 3 (CRACKED)**.
   * If cooled below $T_{\text{solidus}}$ with $I_{\text{Won}} < 1.0 \implies$ Transitions to **State 2 (WELDED)**.

3. **State 2 — WELDED (Consolidated Joint):**
   $$P_n = K_{\text{bond}} \cdot u_n, \quad \boldsymbol{\tau} = G_{\text{bond}} \cdot \mathbf{u}_s \quad (\forall u_n \in (-\infty, +\infty))$$

4. **State 3 — CRACKED (Permanent Hot Tear / Delamination):**
   $$\begin{cases} P_n = K_{\text{comp}} \cdot u_n & \text{if } u_n < 0 \text{ (crack closure / contact)} \\ P_n = 0, \quad \boldsymbol{\tau} = \mathbf{0} & \text{if } u_n \ge 0 \text{ (open crack)} \end{cases}$$

---

## 3. Critical Numerical Bottlenecks Identified During Testing

Extensive patch tests (`dummy_btr_contact_patch.inp`) revealed three major obstacles that prevent robust industrial deployment of `UINTER` for hot cracking:

### 3.1 Pathological Mesh Dependency
* **The Root Cause:** A contact interaction formulation operates on discrete **displacement jumps** ($\Delta u_n$ in mm), whereas metallurgical hot tearing criteria (Won, Rappaz-Drezet-Gremaud RDG) govern **continuum strain rates** ($\dot{\varepsilon}$).
* To convert $\Delta u_n$ to strain, `UINTER` requires an arbitrary *virtual length scale* ($L_0 = 0.35\text{ mm}$ pad thickness).
* **Consequence:** Refining the finite element mesh from $0.35\text{ mm}$ to $0.10\text{ mm}$ concentrates the compliance into fewer nodes. The apparent strain shoots up artificially, triggering **widespread false-positive hot cracking predictions**. Coarser meshes fail to detect true tears.

### 3.2 Contact Chattering & Severe Discontinuity Iterations (SDI)
* When slave nodes cross the solidus line, the normal contact stiffness jumps abruptly from zero ($K = 0$) to bonded stiffness ($K_{\text{bond}} \approx 5 \times 10^5\text{ N/mm}^3$).
* In Abaqus/Standard, this instantaneous switch triggers **Severe Discontinuity Iterations (SDI)**. The Newton-Raphson solver experiences severe cutbacks (`1U`, `2U`, `3U`) and frequent solver stalls, multiplying CPU time by $5\times - 8\times$.

### 3.3 MPI Domain Decomposition Conflicts
* In multi-core HPC environments (`abaqus cpus=8..16`), master and slave surfaces are frequently partitioned across different CPU memory domains.
* Synchronizing non-local state transitions across distributed sub-domains introduces communication overhead and occasional domain-boundary locking.

### 3.4 Empirical Temporal Discretization Study (BTR Sampling Sensitivity)

To rigorously investigate the temporal sensitivity of the `UINTER` routine, a controlled benchmark was executed across 4 time-step discretizations. The thermal cooling rate was held constant at $\dot{T} = 2812.5\,^\circ\text{C/s}$ across the BTR interval ($[1485\,^\circ\text{C}, 1530\,^\circ\text{C}]$, total span $\Delta t_{\text{BTR}} = 0.0160\text{ s}$), while Step 5 was systematically partitioned into 2, 3, 4, and 5 temporal sampling points:

| Numerical Metric | 2 Points in BTR | 3 Points in BTR | 4 Points in BTR | 5 Points in BTR | Asymptotic Behavior |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Fixed $\Delta t$ in Step 5** | $0.01185\text{ s}$ | $0.00790\text{ s}$ | $0.00547\text{ s}$ | $0.00444\text{ s}$ | $2.7\times$ refinement |
| **Fracture Step 5 Inc.** | **Inc. 5** | **Inc. 6** | **Inc. 8** | **Inc. 10** | Proportional tracking |
| **Fracture Temp. ($T_{\text{frac}}$)** | **$1469.67\,^\circ\text{C}$** | **$1469.67\,^\circ\text{C}$** | **$1479.92\,^\circ\text{C}$** | **$1478.00\,^\circ\text{C}$** | Resolves solidus boundary |
| **Incremental Strain ($\Delta \varepsilon_{\text{inc}}$)** | $8.53 \times 10^{-3}$ | $9.36 \times 10^{-3}$ | $4.04 \times 10^{-3}$ | $2.58 \times 10^{-3}$ | Proportional to $\Delta t$ |
| **Won Index ($\varepsilon / \varepsilon_{\text{crit}}$)** | $1.2866$ | $1.2515$ | $1.1062$ | $1.1394$ | Converges to $\approx 1.12$ |
| **Severe Discontinuities** | 3 SDI iters | 3 SDI iters | 3 SDI iters | 3 SDI iters | Systematic across nodes 9–12 |
| **Final Gap `COPEN` ($25^\circ\text{C}$)** | **$2.0636\,\mu\text{m}$** | **$2.0663\,\mu\text{m}$** | **$2.0680\,\mu\text{m}$** | **$2.0684\,\mu\text{m}$** | **Monotonic ($< 0.23\%$ error)** |
| **Final Pressure `CPRESS`** | $-970.2\text{ MPa}$ | $-971.5\text{ MPa}$ | $-972.3\text{ MPa}$ | $-972.5\text{ MPa}$ | **Monotonic ($< 0.24\%$ error)** |

#### Key Insights from the Temporal Benchmark:
1. **Asymptotic Structural Convergence:** Once ruptured, the final room-temperature crack opening `COPEN` converges monotonically with an error below $0.23\%$ ($2.0636\,\mu\text{m} \rightarrow 2.0684\,\mu\text{m}$), demonstrating that the post-rupture contact release mechanics are structurally sound.
2. **Thermal Sampling Overshoot:** Coarse time stepping ($\le 3$ points in BTR) skips past the solidus temperature ($1485^\circ\text{C}$), registering fracture with a numerical lag at $1469.67^\circ\text{C}$. A minimum of **4 to 5 increments inside the BTR** ($\Delta t \le 0.005\text{ s}$) is strictly necessary to resolve the true physical onset of hot tearing.
3. **Severe Discontinuity Overhead:** In all cases, the sudden loss of cohesion across interface nodes triggers **3 Severe Discontinuity Iterations (SDI)** before standard Newton-Raphson equilibrium can be re-established.

---

## 4. Architectural Comparison: UINTER vs. Cohesive UMAT

To resolve these limitations, the project pivoted to a **Cohesive Zone Modeling (CZM)** approach via a thermo-mechanical `UMAT`:

| Performance Metric | Surface Interaction (`UINTER`) | Cohesive Zone (`UMAT` + `COH3D8`) |
| :--- | :--- | :--- |
| **Mathematical Domain** | 2D Contact Master/Slave Surface | Continuous Interface Finite Elements |
| **Mesh Dependency** | ❌ **High** (arbitrary $L_0$ scaling) | ✅ **Zero** (regularized by fracture energy $G_c$) |
| **Solver Convergence** | ⚠️ High SDI cutbacks & chattering | ✅ Smooth quadratic Newton-Raphson |
| **HPC / MPI Scaling** | ⚠️ Domain boundary communication bottlenecks | ✅ Perfect linear scaling (element-local integration) |
| **Industrial Readiness** | ⚠️ Exploratory Research (Backlog) | 🚀 **Production Architecture** ([Project 04](../04_umat_cohesive_welding/README.md)) |

---

## 5. File Inventory

```text
03_uinter_welding_interface/
├── uinter_btr_hot_cracking.f        # User Subroutine (UINTER state machine in Fortran)
├── dummy_btr_contact_patch.inp      # 2-Element unit verification patch test deck
└── README.md                        # Technical report & research documentation
```

---

## 6. How to Run the Unit Patch Test

To compile and verify the exploratory UINTER routine against the 2-element patch test:

```bash
abaqus job=dummy_btr_contact_patch user=uinter_btr_hot_cracking.f interactive
```

> **Next Step:** See [**`04_umat_cohesive_welding`**](../04_umat_cohesive_welding/README.md) for the production-grade, mesh-independent implementation.
