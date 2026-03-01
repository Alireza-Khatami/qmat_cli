#!/usr/bin/env bash
# batch_qmat.sh
#
# For every .obj / .off file found in <input_folder>, runs qmat_cli.exe at
# each simplification target and organises all output files as:
#
#   <output_root>/
#     <mesh_name>/
#       <N>/
#         <mesh_name>_v_<N>.ma
#         <mesh_name>_v_<N>_sampledpoints.txt
#         <mesh_name>_v_<N>_vertex_samples.txt
#         ... (all other files written by qmat_cli with that prefix)
#
# Usage:
#   ./batch_qmat.sh <input_folder> [output_root]
#
# Examples:
#   ./batch_qmat.sh /c/Data/meshes
#   ./batch_qmat.sh /c/Data/meshes /c/Data/qmat_results

# ─── configuration ────────────────────────────────────────────────────────────

# Derive the exe path relative to this script's own directory so the
# absolute path never has to be hardcoded (and avoids Git Bash path-with-
# spaces translation issues).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QMAT_EXE="$SCRIPT_DIR/qmat_cli.exe"

K_VALUE=0.0001

# Simplification targets (vertex counts) — edit to taste
SIMP_TARGETS=(50 100 250 500 1000 2000 5000)
# SIMP_TARGETS=( 500 1000 2000 5000)

# ─── argument handling ────────────────────────────────────────────────────────

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input_folder> [output_root]"
    exit 1
fi

INPUT_FOLDER="${1%/}"          # strip trailing slash if any
OUTPUT_ROOT="${2:-${INPUT_FOLDER}_qmat_results}"

# ─── sanity checks ────────────────────────────────────────────────────────────

if [[ ! -d "$INPUT_FOLDER" ]]; then
    echo "ERROR: input folder not found: $INPUT_FOLDER"
    exit 1
fi

if [[ ! -f "$QMAT_EXE" ]]; then
    echo "ERROR: qmat_cli.exe not found at:"
    echo "  $QMAT_EXE"
    exit 1
fi

mkdir -p "$OUTPUT_ROOT"

# ─── helpers ─────────────────────────────────────────────────────────────────

# Convert a bash/MSYS path to a Windows path that qmat_cli.exe understands.
# Works with Git Bash (cygpath) and falls back to a simple slash-swap.
to_win_path() {
    if command -v cygpath &>/dev/null; then
        cygpath -w "$1"
    else
        # replace /c/ → C:/ style (basic MSYS convention)
        echo "$1" | sed 's|^/\([a-zA-Z]\)/|\1:/|' | tr '/' '\\'
    fi
}

# ─── main loop ────────────────────────────────────────────────────────────────

total_runs=0
failed_runs=0

# Collect mesh files (obj and off, case-insensitive via find)
mapfile -t MESH_FILES < <(
    find "$INPUT_FOLDER" -maxdepth 1 -type f \
        \( -iname "*.obj" -o -iname "*.off" \) | sort
)

if [[ ${#MESH_FILES[@]} -eq 0 ]]; then
    echo "No .obj or .off files found in: $INPUT_FOLDER"
    exit 0
fi

echo "Found ${#MESH_FILES[@]} mesh(es) in: $INPUT_FOLDER"
echo "Output root : $OUTPUT_ROOT"
echo "Targets     : ${SIMP_TARGETS[*]}"
echo "k value     : $K_VALUE"
echo ""

for mesh_file in "${MESH_FILES[@]}"; do
    mesh_basename=$(basename "$mesh_file")
    mesh_name="${mesh_basename%.*}"            # strip extension

    echo "══════════════════════════════════════"
    echo "Mesh: $mesh_name"
    echo "══════════════════════════════════════"

    for N in "${SIMP_TARGETS[@]}"; do
        dest_dir="$OUTPUT_ROOT/$mesh_name/$N"
        mkdir -p "$dest_dir"

        # qmat_cli.exe writes all output files using this prefix
        output_prefix="$dest_dir/${mesh_name}_v_${N}"

        # Convert paths for the Windows executable
        win_input=$(to_win_path "$mesh_file")
        win_prefix=$(to_win_path "$output_prefix")

        echo "  simplify → $N  |  $dest_dir"

        "$QMAT_EXE" "$win_input" \
            --output  "$win_prefix" \
            --simplify "$N" \
            --k "$K_VALUE"

        exit_code=$?
        (( total_runs++ )) || true

        if [[ $exit_code -ne 0 ]]; then
            echo "  !! qmat_cli.exe failed (exit $exit_code) for $mesh_name / target $N"
            (( failed_runs++ )) || true
        else
            # Count how many files landed in dest_dir
            n_files=$(find "$dest_dir" -maxdepth 1 -type f | wc -l)
            echo "     $n_files file(s) written"
        fi
    done

    echo ""
done

# ─── summary ─────────────────────────────────────────────────────────────────

echo "══════════════════════════════════════"
echo "Done."
echo "  Total runs : $total_runs"
echo "  Failed     : $failed_runs"
echo "  Results in : $OUTPUT_ROOT"
echo "══════════════════════════════════════"



#Usage

  # Basic — output goes to <input_folder>_qmat_results/
#   ./batch_qmat.sh /c/Data/meshes

#   # Explicit output folder
#   ./batch_qmat.sh /c/Data/meshes /c/Data/qmat_results

#   What it does

#   For each .obj / .off in the input folder, it loops over all 7 targets and calls:
#   qmat_cli.exe <mesh> --output <dest>/<meshname>_v_<N> --simplify <N> --k 0.0001

#   Because --output already points to the destination directory, all files land there immediately — no post-move step needed.

#   Output structure

#   qmat_results/
#     bunny/
#       50/   bunny_v_50.ma   bunny_v_50_sampledpoints.txt ...
#       100/  bunny_v_100.ma  ...
#       250/  ...
#     armadillo/
#       50/   ...

#   Key points

#   - Path conversion (to_win_path): handles the spaces in the QMAT path by using cygpath if available, otherwise falls back to a sed-based MSYS path rewrite,
#    so the Windows EXE sees proper backslash paths.
#   - Non-fatal failures: a bad run prints a warning and continues with the next target/mesh.
#   - Summary at the end shows total runs and how many failed.


# The issue is that Git Bash translates paths differently than the bash tool I used — [[ -f "C:/..." ]] can fail in some Git Bash versions with spaces in
#   the path. The cleanest fix is to derive the EXE path relative to the script's own location instead of hardcoding an absolute path, since qmat_cli.exe sits
#    in the same folder as the script: