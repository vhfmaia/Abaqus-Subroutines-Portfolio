#!/bin/bash
# ==============================================================================
# Circumferential Laser Welding Simulation Pipeline
# Sequentially Coupled Thermo-Mechanical FEA (Abaqus / SIMULIA)
# Author: Victor Maia (vhfm08@gmail.com)
#
# Usage:
#   ./run_pipeline.sh [all | therm | mech | vtk_therm | vtk_mech]
# ==============================================================================
set -euo pipefail

TARGET="${1:-all}"
CPUS="${CPUS:-4}"

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
    echo ">>> [STAGE 1/4] Running Thermal FEA Analysis (Disk_heatsource_TH.inp + dflux_disk_conical_gaussian.f)..."
    USER_SUB="dflux_disk_conical_gaussian.f"
    ${RUN_ABQ} job=Disk_heatsource_TH input=Disk_heatsource_TH.inp user="${USER_SUB}" cpus="${CPUS}" interactive
    echo "[SUCCESS] Thermal analysis completed: Disk_heatsource_TH.odb generated."
}

# Function: Extract Thermal VTK
run_thermal_vtk() {
    echo ""
    echo ">>> [STAGE 2/4] Extracting Thermal Results to VTK (Disk_heatsource_TH.odb)..."
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
    echo ">>> [STAGE 3/4] Running Mechanical FEA Analysis (Disk_heatsource_ME.inp)..."
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
    echo ">>> [STAGE 4/4] Extracting Mechanical Results to VTK (Disk_heatsource_ME.odb)..."
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
    mech)
        run_mechanical_fea
        ;;
    vtk_mech)
        run_mechanical_vtk
        ;;
    all)
        run_thermal_fea
        run_thermal_vtk
        run_mechanical_fea
        run_mechanical_vtk
        ;;
    *)
        echo "[ERROR] Unknown mode '${TARGET}'. Valid modes: all, therm, vtk_therm, mech, vtk_mech"
        exit 1
        ;;
esac

echo ""
echo "=================================================================="
echo "  Pipeline execution finished successfully!"
echo "=================================================================="
