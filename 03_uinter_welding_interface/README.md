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
