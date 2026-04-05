#!/usr/bin/env bash
# Build qmat_cli in Release mode with Polyscope disabled.
# Uses the 'vcpkg-release' preset (binaryDir = build/, QMAT_WITH_POLYSCOPE=OFF).
# Run from any directory — paths are relative to this script's location.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Configuring (preset: vcpkg, Release, QMAT_WITH_POLYSCOPE=OFF) ==="
cmake --preset vcpkg -S "$SCRIPT_DIR" -DUSE_MATSTRUCT_INITIALIZATION=ON

echo "=== Building ==="
cmake --build --preset vcpkg-release --target qmat_cli --parallel

echo "=== Done ==="
echo "Binary: $SCRIPT_DIR/build/Release/qmat_cli.exe"
