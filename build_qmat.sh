#!/usr/bin/env bash
set -e

QNDF=/groups/czhang/axk230084/conda/envs/qndf
SRC_DIR=/groups/xguo/axk230084/src
REPO_DIR="$SRC_DIR/qmat_cli"

# Clone if not already cloned
if [ ! -d "$REPO_DIR" ]; then
    echo "=== Cloning qmat_cli ==="
    /home/axk230084/bin/git clone https://github.com/Alireza-Khatami/qmat_cli "$REPO_DIR"
else
    echo "=== Repo already exists, skipping clone ==="
fi

# Configure
echo "=== Configuring with CMake ==="
"$QNDF/bin/cmake" \
    -S "$REPO_DIR" \
    -B "$REPO_DIR/build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DQMAT_WITH_POLYSCOPE=OFF \
    -DCMAKE_PREFIX_PATH="$QNDF" \
    -DEigen3_DIR="$QNDF/share/eigen3/cmake" \
    -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON \
    -DCMAKE_CXX_COMPILER="$QNDF/bin/x86_64-conda-linux-gnu-g++"

# Build
echo "=== Building ==="
"$QNDF/bin/cmake" --build "$REPO_DIR/build" -j$(nproc)

echo "=== Done ==="
echo "Binary: $REPO_DIR/build/qmat_cli"
