# QMAT CLI - CMake Build Configuration
#
# This builds a command-line version of QMAT without Qt dependencies.
#
# Windows (vcpkg):
#   cmake -B build -DCMAKE_TOOLCHAIN_FILE=C:/Users/alirz/Projects/vcpkg/scripts/buildsystems/vcpkg.cmake
#   cmake --build build --config Release
#
# Linux / HPC (system or module-loaded libs):
#   module load gcc cmake cgal eigen boost   # adapt to your cluster's module names
#   cmake -B build -DCMAKE_BUILD_TYPE=Release -DQMAT_WITH_POLYSCOPE=OFF
#   cmake --build build -j$(nproc)
#
# Switch edge-collapse energy (default = QMAT 4D):
#   cmake -B build -DQMAT_USE_QEM_MESH_PLANE=ON      # MeshLab-style plane-quadric QEM on
#                                                    # the medial mesh (see md_files/
#                                                    # qem_mesh_plane_flow.md).
#
# Dependencies will be installed to: qmat/library/  (Windows/vcpkg only)
# The executable will be at: build/Release/qmat_cli.exe (Windows)
#                             build/qmat_cli           (Linux)

cmake_minimum_required(VERSION 3.16)

if(WIN32)
    # Enforce vcpkg toolchain on Windows — dependencies are not available otherwise.
    if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
        message(FATAL_ERROR
            "On Windows you must pass the vcpkg toolchain file:\n"
            "  cmake -B build "
            "-DCMAKE_TOOLCHAIN_FILE=C:/Users/alirz/Projects/vcpkg/scripts/buildsystems/vcpkg.cmake")
    endif()
    # vcpkg installs packages into the "library" folder (must be before project())
    set(VCPKG_INSTALLED_DIR "${CMAKE_CURRENT_SOURCE_DIR}/library")
endif()

project(qmat_cli VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Windows-specific settings
if(WIN32)
    add_definitions(-DWIN32 -DUNICODE -D_UNICODE)
    add_definitions(-D_USE_MATH_DEFINES)  # For M_PI
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-DNOMINMAX)  # Prevent Windows.h min/max macros
endif()

# ============================================================================
# Optional: Polyscope  (-DQMAT_WITH_POLYSCOPE=ON to enable)
# Downloaded via FetchContent into the build directory.
# Requires OpenGL/GLFW — not available on HPC compute nodes.
# ============================================================================

# Polyscope requires OpenGL/GLFW — not available on HPC compute nodes.
# Default ON for Windows (interactive use), OFF on Linux.
if(WIN32)
    option(QMAT_WITH_POLYSCOPE "Link Polyscope into qmat_cli" ON)
else()
    option(QMAT_WITH_POLYSCOPE "Link Polyscope into qmat_cli" OFF)
endif()

option(USE_MATSTRUCT_INITIALIZATION "Load initial MAT from external typed .ma file instead of computing from Voronoi/Delaunay" OFF)

# MeshLab-style QEM: per-vertex plane quadric K = sum_{f at v} p_f p_f^T over
# incident MEDIAL triangles (planes built from sphere centres; radii ignored).
# When ON, EvaluateEdgeCollapseCost uses this energy instead of the 4D QMAT
# slab quadric. See md_files/qem_medial_mesh_vs_slab.md and
# md_files/qem_mesh_plane_flow.md.
option(QMAT_USE_QEM_MESH_PLANE "Use MeshLab-style plane-quadric QEM on the medial mesh" OFF)

# QEM-only condition checks: when ON, the merge gates in MinCostEdgeCollapse /
# MinCostBoundaryEdgeCollapse / CanMerge are restricted to the genuine VCG/MeshLab
# collapse checks — WouldCreateNonManifold (link condition), FlipsNormals
# (normal-flip), and FailsHardQualityCheck (sliver quality). QMAT's own gates
# (cluster-type, struct-ids, topo_contractable) are skipped. This mode requires
# the QEM mesh-plane energy, so it forces QMAT_USE_QEM_MESH_PLANE ON.
option(ONLY_USE_QEM_CONDITION_CHECKS "Restrict merge gates to QEM-only; implies QMAT_USE_QEM_MESH_PLANE" OFF)
if(ONLY_USE_QEM_CONDITION_CHECKS)
    set(QMAT_USE_QEM_MESH_PLANE ON CACHE BOOL "Forced ON by ONLY_USE_QEM_CONDITION_CHECKS" FORCE)
    message(STATUS "ONLY_USE_QEM_CONDITION_CHECKS=ON -> forcing QMAT_USE_QEM_MESH_PLANE=ON")
endif()

if(QMAT_WITH_POLYSCOPE)
    include(FetchContent)
    FetchContent_Declare(
        polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG        v2.3.0
        GIT_PROGRESS   TRUE
        # Polyscope bundles glfw/imgui/glad/stb as git submodules.
        # Do NOT add GIT_SHALLOW TRUE or the submodules will be missing.
    )
    set(POLYSCOPE_MASTER_PROJECT OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(polyscope)
endif()

# ============================================================================
# Optional: vcglib (-DQMAT_WITH_VCGLIB=ON to enable; ON by default)
# Header-only, pulled via FetchContent from cnr-isti-vclab/vcglib.
# Used by src/VcgDirectSimplifier.cpp to invoke vcg::tri::TriEdgeCollapseQuadric
# (the real MeshLab simplifier) on a vcg::TriMesh built from the slab mesh.
#
# We reuse our vcpkg Eigen3 by setting VCG_ALLOW_SYSTEM_EIGEN + EIGEN3_INCLUDE_DIR
# BEFORE FetchContent_MakeAvailable so vcglib's CMakeLists.txt picks it up and
# does NOT pull the bundled Eigen (avoids two Eigens in the include path).
# ============================================================================
option(QMAT_WITH_VCGLIB "Pull vcglib via FetchContent and enable --vcg-simplify" ON)
if(QMAT_WITH_VCGLIB)
    include(FetchContent)
    FetchContent_Declare(
        vcglib
        GIT_REPOSITORY https://github.com/cnr-isti-vclab/vcglib.git
        GIT_TAG        2023.12
        GIT_PROGRESS   TRUE
    )
    # Tell vcglib to use our Eigen instead of its bundled copy.  Eigen3::Eigen
    # is set by the find_package(Eigen3) below — but FetchContent_MakeAvailable
    # runs *now*, so we resolve the include dir directly from vcpkg's layout.
    set(VCG_ALLOW_SYSTEM_EIGEN ON CACHE BOOL "" FORCE)
    set(VCG_ALLOW_BUNDLED_EIGEN OFF CACHE BOOL "" FORCE)
    set(VCG_HEADER_ONLY ON CACHE BOOL "" FORCE)
    set(VCG_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(VCG_USE_OPENMP OFF CACHE BOOL "" FORCE)
    if(WIN32)
        set(EIGEN3_INCLUDE_DIR "${VCPKG_INSTALLED_DIR}/x64-windows/include/eigen3" CACHE PATH "" FORCE)
    endif()
    FetchContent_MakeAvailable(vcglib)
endif()

# ============================================================================
# Find Dependencies
# Windows: vcpkg toolchain handles everything via find_package normally.
# Linux/HPC: use system packages or environment-module-loaded libraries.
# ============================================================================

if(WIN32)
    find_package(CGAL CONFIG REQUIRED)
    find_package(Eigen3 CONFIG REQUIRED)
    # find_package(Boost CONFIG REQUIRED)
else()
    # CGAL's own CMake config is preferred; falls back to FindCGAL module.
    find_package(CGAL REQUIRED)
    # Eigen3 is header-only; available system-wide or via module load eigen.
    find_package(Eigen3 REQUIRED)
    # GMP and MPFR are CGAL runtime deps — usually present system-wide on HPC.
    find_package(GMP QUIET)
    find_package(MPFR QUIET)
endif()

# ============================================================================
# Source Files
# ============================================================================

set(QMAT_CLI_SOURCES
    src/main_cli.cpp
    src/Mesh.cpp
    src/ThreeDimensionalShape.cpp
    src/SlabMesh.cpp
    src/VcgQuadricSimplifier.cpp
    src/QemRejectionViz.cpp
    src/PrimMesh.cpp
    src/ObjLoader.cpp
    src/NonManifoldMesh/nonmanifoldmesh.cpp
    src/LinearAlgebra/Wm4Math.cpp
    src/LinearAlgebra/Wm4Matrix.cpp
    src/LinearAlgebra/Wm4Vector.cpp
    src/ColorRamp/ColorRamp.cpp
    src/GeometryObjects/GeometryObjects.cpp
    src/matfp/sharp_feature_detection.cpp
)
if(QMAT_WITH_POLYSCOPE)
    list(APPEND QMAT_CLI_SOURCES src/MatVisualizerCommon.cpp)
endif()
if(QMAT_WITH_VCGLIB)
    list(APPEND QMAT_CLI_SOURCES src/VcgDirectSimplifier.cpp)
    if(QMAT_WITH_POLYSCOPE)
        list(APPEND QMAT_CLI_SOURCES src/VdeVisualizer.cpp)
    endif()
endif()

set(QMAT_CLI_HEADERS
    src/Mesh.h
    src/ThreeDimensionalShape.h
    src/SlabMesh.h
    src/VcgQuadricSimplifier.h
    src/QemRejection.h
    src/QemRejectionViz.h
    src/PrimMesh.h
    src/ObjLoader.h
    src/tiny_obj_loader.h
    src/NonManifoldMesh/nonmanifoldmesh.h
    src/LinearAlgebra/Wm4Math.h
    src/LinearAlgebra/Wm4Matrix.h
    src/LinearAlgebra/Wm4Vector.h
    src/ColorRamp/ColorRamp.h
    src/GeometryObjects/GeometryObjects.h
    src/matfp/sharp_feature_detection.h
    src/matfp/Args.h
    src/matfp/pre_types.h
)
if(QMAT_WITH_POLYSCOPE)
    list(APPEND QMAT_CLI_HEADERS src/MatVisualizerCommon.h)
endif()
if(QMAT_WITH_VCGLIB)
    list(APPEND QMAT_CLI_HEADERS src/VcgDirectSimplifier.h)
    if(QMAT_WITH_POLYSCOPE)
        list(APPEND QMAT_CLI_HEADERS src/VdeVisualizer.h)
    endif()
endif()

# ============================================================================
# Create Executable
# ============================================================================

add_executable(qmat_cli ${QMAT_CLI_SOURCES} ${QMAT_CLI_HEADERS})

# ============================================================================
# Include Directories
# ============================================================================

target_include_directories(qmat_cli PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/LinearAlgebra
    ${CMAKE_CURRENT_SOURCE_DIR}/src/NonManifoldMesh
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ColorRamp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/GeometryObjects
    ${CMAKE_CURRENT_SOURCE_DIR}/src/matfp
)

# ============================================================================
# Link Libraries (vcpkg handles include paths automatically)
# ============================================================================

target_link_libraries(qmat_cli PRIVATE
    CGAL::CGAL
    Eigen3::Eigen
)

if(NOT WIN32)
    # On Linux/HPC, explicitly link GMP and MPFR if found as standalone targets.
    if(GMP_FOUND)
        target_link_libraries(qmat_cli PRIVATE ${GMP_LIBRARIES})
    endif()
    if(MPFR_FOUND)
        target_link_libraries(qmat_cli PRIVATE ${MPFR_LIBRARIES})
    endif()
endif()

if(QMAT_WITH_POLYSCOPE)
    target_link_libraries(qmat_cli PRIVATE polyscope)
    target_compile_definitions(qmat_cli PRIVATE QMAT_WITH_POLYSCOPE)
endif()

if(QMAT_WITH_VCGLIB)
    target_link_libraries(qmat_cli PRIVATE vcglib)
    target_compile_definitions(qmat_cli PRIVATE QMAT_WITH_VCGLIB)
endif()

if(USE_MATSTRUCT_INITIALIZATION)
    target_compile_definitions(qmat_cli PRIVATE USE_MATSTRUCT_INITIALIZATION)
endif()

if(ONLY_USE_QEM_CONDITION_CHECKS)
    target_compile_definitions(qmat_cli PRIVATE ONLY_USE_QEM_CONDITION_CHECKS)
endif()

if(QMAT_USE_QEM_MESH_PLANE)
    target_compile_definitions(qmat_cli PRIVATE QMAT_USE_QEM_MESH_PLANE)
endif()

# tiny_obj_loader.h is in the main source directory (header-only library)
# No additional include path needed since CMAKE_CURRENT_SOURCE_DIR is already included
# ============================================================================
# Compiler Settings
# ============================================================================

if(MSVC)
    target_compile_options(qmat_cli PRIVATE /W3 /wd4267 /wd4244 /wd4305)
endif()

# ============================================================================
# Installation
# ============================================================================

install(TARGETS qmat_cli RUNTIME DESTINATION bin)

# ============================================================================
# Print Configuration Summary
# ============================================================================

message(STATUS "")
message(STATUS "QMAT CLI Configuration:")
message(STATUS "  Platform        : ${CMAKE_SYSTEM_NAME}")
message(STATUS "  C++ Standard    : ${CMAKE_CXX_STANDARD}")
if(WIN32)
message(STATUS "  vcpkg dir       : ${VCPKG_INSTALLED_DIR}")
endif()
message(STATUS "  CGAL            : Found")
message(STATUS "  Eigen3          : Found")
message(STATUS "  Polyscope       : ${QMAT_WITH_POLYSCOPE}")
message(STATUS "  vcglib          : ${QMAT_WITH_VCGLIB}  (header-only FetchContent; enables --simplifier vcg-direct)")
message(STATUS "  MatStruct init  : ${USE_MATSTRUCT_INITIALIZATION}")
message(STATUS "  Edge energy     : ${QMAT_USE_QEM_MESH_PLANE}  (ON=MeshLab plane-quadric, OFF=QMAT 4D)")
message(STATUS "  QEM-only gates  : ${ONLY_USE_QEM_CONDITION_CHECKS}  (ON=manifold+flip+quality only, skip QMAT cluster/struct/topo gates)")
message(STATUS "  matfp sources   : matfp/sharp_feature_detection.cpp")
message(STATUS "")
