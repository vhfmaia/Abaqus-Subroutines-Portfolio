#!/bin/bash
set -euo pipefail

echo "=========================================================="
echo "  Abaqus / 3DEXPERIENCE: Exportacao ODB -> VTK & ZIP      "
echo "=========================================================="

JOB_NAME="${1:-Gear_torque}"
ODB_FILE="${JOB_NAME%.odb}.odb"

# 1. Identificar o motor de execucao Python do Abaqus
if [ -n "${GRID_PRODUCT_PATH:-}" ] && [ -f "${GRID_PRODUCT_PATH}/linux_a64/code/bin/SMALauncher" ]; then
    echo "[INFO] Ambiente 3DEXPERIENCE Grid detetado via GRID_PRODUCT_PATH."
    LAUNCHER="${GRID_PRODUCT_PATH}/linux_a64/code/bin/SMALauncher"
elif command -v SMALauncher >/dev/null 2>&1; then
    echo "[INFO] Comando 'SMALauncher' detetado no sistema."
    LAUNCHER="SMALauncher"
elif command -v abaqus >/dev/null 2>&1; then
    echo "[INFO] Comando 'abaqus' standard detetado no sistema."
    LAUNCHER="abaqus"
else
    echo "[AVISO] Interpretador standard nao encontrado. A tentar SMALauncher direto."
    LAUNCHER="SMALauncher"
fi

# 2. Executar o script de conversao ODB -> VTK
echo "[INFO] A invocar interpretador Python: ${LAUNCHER} python odb_to_vtk.py ${ODB_FILE}"
if "${LAUNCHER}" python odb_to_vtk.py "${ODB_FILE}"; then
    echo "[INFO] Exportacao concluida com SMALauncher python."
else
    echo "[INFO] A tentar chamada direta de script SMALauncher..."
    "${LAUNCHER}" odb_to_vtk.py "${ODB_FILE}"
fi

# 3. Compactar todos os ficheiros VTK em gear_helical_vtk.zip
echo "[INFO] A compactar ficheiros VTK em 'gear_helical_vtk.zip'..."
if command -v zip >/dev/null 2>&1; then
    zip -q -9 gear_helical_vtk.zip gear_helical_*.vtk
    echo "[SUCESSO] Ficheiro 'gear_helical_vtk.zip' gerado com sucesso via zip!"
else
    "${LAUNCHER}" python -c "import zipfile, glob, os; vtks = sorted(glob.glob('gear_helical_*.vtk')); z = zipfile.ZipFile('gear_helical_vtk.zip', 'w', zipfile.ZIP_DEFLATED); [z.write(f) for f in vtks]; z.close(); print('[SUCESSO] gear_helical_vtk.zip criado (%d MB) com %d frames!' % (os.path.getsize('gear_helical_vtk.zip') // 1048576, len(vtks)))" || python -c "import zipfile, glob; vtks = sorted(glob.glob('gear_helical_*.vtk')); z = zipfile.ZipFile('gear_helical_vtk.zip', 'w', zipfile.ZIP_DEFLATED); [z.write(f) for f in vtks]; z.close()"
fi

echo "=========================================================="
echo "  Post-processing concluido com sucesso!                  "
echo "  Ficheiro unico 'gear_helical_vtk.zip' pronto para download! "
echo "=========================================================="
