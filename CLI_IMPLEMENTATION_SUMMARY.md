# QMAT CLI Implementation Summary

## Goal
Create a command-line interface to run QMAT's core functionality without Qt GUI:
- Load .OFF mesh files
- Compute Medial Axis (computeMA)
- Simplify slab mesh (simplifySlab)
- Export results as .ma files

---

## Files Created

### 1. `main_cli.cpp` (NEW)
A standalone CLI entry point that:
- Parses command-line arguments (`--simplify`, `--k`, `--output`, `--help`)
- Loads OFF mesh files using CGAL stream reader
- Computes the Medial Axis (Delaunay triangulation + MA extraction)
- Optionally simplifies the slab mesh to a target vertex count
- Exports results as `.ma` files

**CLI Usage:**
```
qmat_cli <input.off> [options]

Options:
  --simplify <N>     Simplify to N vertices (default: no simplification)
  --k <value>        K factor for slab initialization (default: 0.00001)
  --output <prefix>  Output file prefix (default: input filename)
  --help             Show help message
```

**Example:**
```
qmat_cli model.off --simplify 1000 --k 0.00001 --output result
```

### 2. `CMakeLists.txt` (NEW)
CMake build configuration that:
- Sets C++ standard to C++14 (compatible with CGAL 4.x)
- Finds dependencies via vcpkg or manual paths
- Builds the `qmat_cli` executable
- Includes only computation source files (excludes Qt/GUI files)
- Links GMP/MPFR libraries required by CGAL
- Adds necessary Windows defines (`_USE_MATH_DEFINES`, `NOMINMAX`, etc.)

### 3. `vcpkg.json` (NEW)
vcpkg manifest file specifying dependencies:
- **CGAL 4.14-3** (overridden - required for API compatibility)
- **Boost** (required by CGAL)
- **GMP/MPFR** (required by CGAL exact arithmetic)
- **Eigen3** (linear algebra)

---

## Files Modified

### 1. `ThreeDimensionalShape.cpp`
- Removed unused `#include <QString>` to eliminate Qt dependency

### 2. `Mesh.h`
- Initially modified to support newer CGAL API (circumcenter cell base)
- **Reverted** to original when switching to CGAL 4.14-3 (older API)

---

## Source Files Included in CLI Build

**Core computation (no Qt dependencies):**
- `main_cli.cpp`
- `Mesh.cpp`
- `ThreeDimensionalShape.cpp`
- `SlabMesh.cpp`
- `PrimMesh.cpp`
- `NonManifoldMesh/nonmanifoldmesh.cpp`
- `LinearAlgebra/Wm4Math.cpp`
- `LinearAlgebra/Wm4Matrix.cpp`
- `LinearAlgebra/Wm4Vector.cpp`
- `ColorRamp/ColorRamp.cpp`
- `GeometryObjects/GeometryObjects.cpp`

**Excluded (GUI-only):**
- `main.cpp` - Qt GUI entry point
- `medialaxissimplification3d.cpp` - Qt GUI wrapper
- `GLWidget.cpp` - OpenGL rendering
- `PsRender/PsRender.cpp` - PostScript rendering
- `GeneratedFiles/*` - Qt MOC/UIC generated files

---

## Dependencies

### Required (via vcpkg):
- **CGAL 4.14-3** - Computational Geometry Algorithms Library
- **Boost** - C++ libraries (required by CGAL)
- **GMP** - GNU Multiple Precision Arithmetic Library
- **MPFR** - Multiple Precision Floating-Point Reliable Library
- **Eigen3** - Linear algebra library

### NOT required (removed):
- Qt5Core, Qt5Gui, Qt5Widgets, Qt5OpenGL
- OpenGL (opengl32.lib, glu32.lib)
- qtmain.lib

---

## Build Instructions

### Step 1: Install Dependencies via vcpkg
```bash
cd "C:\Users\alirz\Projects\Graphics\QMAT_old working version  exe file\qmat_x64\qmat"
C:\Users\alirz\Projects\vcpkg\vcpkg.exe install --triplet x64-windows
```

### Step 2: Configure with CMake
```bash
cmake -B build -DCMAKE_TOOLCHAIN_FILE="C:/Users/alirz/Projects/vcpkg/scripts/buildsystems/vcpkg.cmake"
```

### Step 3: Build
```bash
cmake --build build --config Release
```

### Output
The executable will be at: `build/Release/qmat_cli.exe`

---

## Alternative Build (Manual Paths)

If not using vcpkg toolchain, the CMakeLists.txt falls back to:
- vcpkg packages directory: `C:/Users/alirz/Projects/vcpkg/packages`
- Eigen: `C:/Users/alirz/Projects/Graphics/QMAT/include/eigen-3.4.0`

---

## Key Code Flow (from main_cli.cpp)

```cpp
// 1. Load OFF file
std::ifstream stream(inputFile);
stream >> shape.input;
shape.input.computebb();
shape.input.GenerateList();
shape.input.GenerateRandomColor();
shape.input.compute_normals();

// 2. Create mesh domain for inside/outside queries
Polyhedron pol;
streampol >> pol;
Mesh_domain* domain = new Mesh_domain(pol);
shape.input.domain = domain;

// 3. Compute Delaunay Triangulation
shape.input.computedt();

// 4. Compute Medial Axis
shape.ComputeInputNMM();
// Exports raw MA to <output>.ma

// 5. If simplification requested:
shape.LoadInputNMM(maFile);  // Load MA into slab mesh
shape.LoadSlabMesh();         // Initialize for simplification
shape.slab_mesh.CleanIsolatedVertices();
shape.slab_mesh.Simplify(reductionCount);

// 6. Export simplified MA
shape.slab_mesh.Export(outputPrefix);
```

---

## Issues Encountered & Solutions

### Issue 1: Qt dependency in ThreeDimensionalShape.cpp
**Solution:** Removed unused `#include <QString>`

### Issue 2: CGAL not found
**Solution:** Created vcpkg.json manifest and CMakeLists.txt with proper paths

### Issue 3: Boost split into many vcpkg packages
**Solution:** CMakeLists.txt uses `file(GLOB ...)` to find all boost-* package include directories

### Issue 4: M_PI undefined on MSVC
**Solution:** Added `-D_USE_MATH_DEFINES` to CMake definitions

### Issue 5: CGAL API incompatibility (newer vcpkg CGAL 6.x requires C++17 and different cell base)
**Solution:** Created vcpkg.json with override to use CGAL 4.14-3 (compatible with original codebase API)

---

## Verification

After building:
1. Run: `build/Release/qmat_cli.exe model.off --simplify 500`
2. Check output `.ma` files are created
3. Optionally compare with GUI output for same input

---

## Notes

- The original project was built with CGAL 4.8.1, Boost 1.61
- vcpkg's oldest available CGAL is 4.11, with 4.14-3 being the latest 4.x version
- CGAL 5.x and 6.x have breaking API changes (requires C++17, different cell base classes)
- The vcpkg.json uses a baseline override to get CGAL 4.14-3 instead of the current 6.0.1
