#!/bin/bash
# ==============================================================================
# Circumferential Laser Welding Simulation Pipeline
# Sequentially Coupled Thermo-Mechanical FEA (Abaqus / SIMULIA / 3DEXPERIENCE)
# Author: Victor Maia (vhfm08@gmail.com)
#
# Usage:
#   ./run_pipeline.sh [all | therm | mech | vtk_therm | vtk_mech | auto]
# ==============================================================================
set -euo pipefail

CPUS="${CPUS:-4}"

# Intelligent Auto-Detection of execution target
if [ -z "${1:-}" ] || [ "${1:-}" = "auto" ]; then
    if [ -f "Disk_heatsource_TH.odb" ] && [ ! -f "Disk_heatsource_ME.inp" ]; then
        echo "[AUTO-DETECT] Detected completed thermal run (Disk_heatsource_TH.odb) in standalone thermal folder."
        echo "[AUTO-DETECT] Executing post-processing VTK extraction for thermal results..."
        TARGET="vtk_therm"
    elif [ -f "Disk_heatsource_ME.odb" ]; then
        echo "[AUTO-DETECT] Detected completed mechanical run (Disk_heatsource_ME.odb)."
        echo "[AUTO-DETECT] Executing post-processing VTK extraction for mechanical results..."
        TARGET="vtk_mech"
    elif [ -f "Disk_heatsource_TH.inp" ] && [ ! -f "Disk_heatsource_ME.inp" ]; then
        echo "[AUTO-DETECT] Standalone thermal simulation deck detected. Running thermal FEA + VTK..."
        TARGET="therm_pipeline"
    else
        TARGET="all"
    fi
else
    TARGET="$1"
fi

echo "=================================================================="
echo "  Circumferential Laser Welding Simulation Pipeline"
echo "  Mode: ${TARGET} | Parallel Cores: ${CPUS}"
echo "=================================================================="

# 1. Environment & Launcher Auto-Detection
if [ -n "${GRID_PRODUCT_PATH:-}" ] && [ -f "${GRID_PRODUCT_PATH}/linux_a64/code/bin/SMALauncher" ]; then
    echo "[INFO] 3DEXPERIENCE Grid environment detected via GRID_PRODUCT_PATH."
    LAUNCHER="${GRID_PRODUCT_PATH}/linux_a64/code/bin/SMALauncher"
    RUN_ABQ="${LAUNCHER} abaqus"
    RUN_PY="${LAUNCHER} python"
elif command -v SMALauncher >/dev/null 2>&1; then
    echo "[INFO] Native SMALauncher command detected."
    LAUNCHER="SMALauncher"
    RUN_ABQ="${LAUNCHER} abaqus"
    RUN_PY="${LAUNCHER} python"
elif command -v abaqus >/dev/null 2>&1; then
    echo "[INFO] Standard Abaqus command detected."
    RUN_ABQ="abaqus"
    RUN_PY="abaqus python"
else
    echo "[WARNING] Abaqus launcher not found in PATH. Defaulting to 'abaqus'."
    RUN_ABQ="abaqus"
    RUN_PY="abaqus python"
fi

# Function: Run Thermal Simulation
run_thermal_fea() {
    echo ""
    echo ">>> Running Thermal FEA Analysis (Disk_heatsource_TH.inp + dflux_disk_conical_gaussian.f)..."
    USER_SUB="dflux_disk_conical_gaussian.f"
    ${RUN_ABQ} job=Disk_heatsource_TH input=Disk_heatsource_TH.inp user="${USER_SUB}" cpus="${CPUS}" interactive
    echo "[SUCCESS] Thermal analysis completed: Disk_heatsource_TH.odb generated."
}

# Function: Extract Thermal VTK
run_thermal_vtk() {
    echo ""
    echo ">>> Extracting Thermal Results to VTK (Disk_heatsource_TH.odb)..."
    if [ ! -f "Disk_heatsource_TH.odb" ]; then
        echo "[ERROR] Disk_heatsource_TH.odb not found!"
        exit 1
    fi
    ${RUN_PY} odb_to_vtk.py Disk_heatsource_TH.odb --prefix=laser_therm
    echo "[SUCCESS] Thermal VTK archive ready: laser_therm_vtk.zip"
}

# Function: Run Mechanical Simulation
run_mechanical_fea() {
    echo ""
    echo ">>> Running Mechanical FEA Analysis (Disk_heatsource_ME.inp)..."
    if [ ! -f "Disk_heatsource_TH.odb" ] && [ ! -f "Disk_heatsource_TH.fil" ]; then
        echo "[ERROR] Cannot run mechanical step without thermal results (Disk_heatsource_TH.odb or .fil)!"
        exit 1
    fi
    ${RUN_ABQ} job=Disk_heatsource_ME input=Disk_heatsource_ME.inp cpus="${CPUS}" interactive
    echo "[SUCCESS] Mechanical analysis completed: Disk_heatsource_ME.odb generated."
}

# Function: Extract Mechanical VTK
run_mechanical_vtk() {
    echo ""
    echo ">>> Extracting Mechanical Results to VTK (Disk_heatsource_ME.odb)..."
    if [ ! -f "Disk_heatsource_ME.odb" ]; then
        echo "[ERROR] Disk_heatsource_ME.odb not found!"
        exit 1
    fi
    ${RUN_PY} odb_to_vtk.py Disk_heatsource_ME.odb --prefix=laser_mech
    echo "[SUCCESS] Mechanical VTK archive ready: laser_mech_vtk.zip"
}

# Execution Dispatch
case "${TARGET}" in
    therm)
        run_thermal_fea
        ;;
    vtk_therm)
        run_thermal_vtk
        ;;
    therm_pipeline)
        run_thermal_fea
        run_thermal_vtk
        ;;
    mech)
        run_mechanical_fea
        ;;
    vtk_mech)
        run_mechanical_vtk
        ;;
    mech_pipeline)
        run_mechanical_fea
        run_mechanical_vtk
        ;;
    all)
        run_thermal_fea
        run_thermal_vtk
        run_mechanical_fea
        run_mechanical_vtk
        ;;
    *)
        echo "[ERROR] Unknown mode '${TARGET}'. Valid modes: all, therm, vtk_therm, therm_pipeline, mech, vtk_mech, auto"
        exit 1
        ;;
esac

echo ""
echo "=================================================================="
echo "  Pipeline execution finished successfully!"
echo "=================================================================="
