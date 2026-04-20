#!/bin/bash

# ── Build ────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[Build] Configuring (vcpkg preset, Polyscope ON)..."
cmake --preset vcpkg -DQMAT_WITH_POLYSCOPE=ON -S "$SCRIPT_DIR" || { echo "[Build] Configure FAILED."; exit 1; }
echo "[Build] Building (vcpkg-release preset)..."
cmake --build --preset vcpkg-release || { echo "[Build] Build FAILED."; exit 1; }
echo "[Build] Done."
echo ""

# Paths
QMAT_ROOT="C:/Users/alirz/Projects/Graphics/QMAT_old working version  exe file/qmat_x64/qmat"
ABC_DATASET_DIR="D:/datasets/abc_full_10k/out_ABC_v6_knn_poission40_20_15_10"
OBJ_DATASET_DIR="D:/datasets/abc_full_10k/ABC_input(use scaled_sf.obj)"
MA_STRUCT_DIR="D:/datasets/abc_full_10k/out_ABC_v6_knn_poission40_20_15_10"
EXE="$QMAT_ROOT/build/Release/qmat_cli.exe"
LOG_DIR="$QMAT_ROOT/output/abc_dataset_logs"
TIMEOUT_SECS=900  # 15 minutes

# # Git Bash on Windows does not ship GNU timeout; use a no-op fallback.
# if ! command -v timeout &>/dev/null; then
#     timeout() { local t=$1; shift; "$@"; }
# fi

# mkdir -p "$LOG_DIR"

# total=0
# success=0
# failed=0
# skipped=0
# timedout=0

# echo "=============================="
# echo " QMAT ABC Batch Runner"
# echo " $(date)"
# echo " Logs -> $LOG_DIR"
# echo "=============================="
# echo ""

# # Use find instead of glob to safely enumerate dirs in a path with special chars
# while IFS= read -r -d '' mesh_dir; do

#     folder_name=$(basename "$mesh_dir")
#     mesh_name="${folder_name%.obj}"

#     # Find the .ma_struct file — filename contains a variable timestamp
#     mat_file=$(find "$MA_STRUCT_DIR/$mesh_name/mat" -maxdepth 1 -name "mat_${mesh_name}.obj__*.ma_struct" 2>/dev/null | head -1)
#     obj_file="$OBJ_DATASET_DIR/$mesh_name/${mesh_name}_scaled_sf.obj"
#     log_file="$LOG_DIR/${mesh_name}.log"

#     total=$((total + 1))

#     output_dir="$QMAT_ROOT/output/abc_dataset/$folder_name"
#     {
#         echo "=============================="
#         echo " Mesh    : $mesh_name"
#         echo " Started : $(date)"
#         echo " MAT     : $mat_file"
#         echo " OBJ     : $obj_file"
#         echo " Output  : $output_dir"
#         echo " Log     : $log_file"
#         echo "=============================="
#     } > "$log_file"

#     if [ ! -f "$mat_file" ]; then
#         msg="[SKIP] MAT file not found: $mat_file"
#         echo "$msg" | tee -a "$log_file"
#         skipped=$((skipped + 1))
#         continue
#     fi

#     if [ ! -f "$obj_file" ]; then
#         msg="[SKIP] OBJ file not found: $obj_file"
#         echo "$msg" | tee -a "$log_file"
#         skipped=$((skipped + 1))
#         continue
#     fi

#     mkdir -p "$output_dir"
#     echo "[RUN] $mesh_name"
#     echo "      Output -> $output_dir"
#     echo "      Log    -> $log_file"

#     timeout "$TIMEOUT_SECS" "$EXE" \
#         --output "$output_dir" \
#         --matstruct "$mat_file" \
#         --feature-angle 40.0 \
#         --k 0.0001 \
#         --simplify 1 \
#         "$obj_file" \
#         >> "$log_file" 2>&1

#     exit_code=$?

#     echo "" >> "$log_file"
#     echo "------------------------------" >> "$log_file"

#     if [ $exit_code -eq 124 ]; then
#         msg="[TIMEOUT] Killed after ${TIMEOUT_SECS}s -- $(date)"
#         timedout=$((timedout + 1))
#     elif [ $exit_code -eq 0 ]; then
#         msg="[SUCCESS] Exit 0 -- $(date)"
#         success=$((success + 1))
#     else
#         msg="[FAILED] Exit code $exit_code -- $(date)"
#         failed=$((failed + 1))
#     fi

#     echo "$msg" | tee -a "$log_file"
#     echo "------------------------------" >> "$log_file"
#     echo "" >> "$log_file"

# done < <(find "$ABC_DATASET_DIR" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

# echo ""
# echo "=============================="
# echo " DONE -- $(date)"
# echo " Total   : $total"
# echo " Success : $success"
# echo " Failed  : $failed"
# echo " Timeout : $timedout"
# echo " Skipped : $skipped"
# echo "=============================="
