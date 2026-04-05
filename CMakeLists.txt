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
    main_cli.cpp
    Mesh.cpp
    ThreeDimensionalShape.cpp
    SlabMesh.cpp
    MatCollapseHistory.cpp
    PrimMesh.cpp
    ObjLoader.cpp
    NonManifoldMesh/nonmanifoldmesh.cpp
    LinearAlgebra/Wm4Math.cpp
    LinearAlgebra/Wm4Matrix.cpp
    LinearAlgebra/Wm4Vector.cpp
    ColorRamp/ColorRamp.cpp
    GeometryObjects/GeometryObjects.cpp
    matfp/sharp_feature_detection.cpp
)

set(QMAT_CLI_HEADERS
    Mesh.h
    ThreeDimensionalShape.h
    SlabMesh.h
    PrimMesh.h
    ObjLoader.h
    tiny_obj_loader.h
    NonManifoldMesh/nonmanifoldmesh.h
    LinearAlgebra/Wm4Math.h
    LinearAlgebra/Wm4Matrix.h
    LinearAlgebra/Wm4Vector.h
    ColorRamp/ColorRamp.h
    GeometryObjects/GeometryObjects.h
    matfp/sharp_feature_detection.h
    matfp/Args.h
    matfp/pre_types.h
)

# ============================================================================
# Create Executable
# ============================================================================

add_executable(qmat_cli ${QMAT_CLI_SOURCES} ${QMAT_CLI_HEADERS})

# ============================================================================
# Include Directories
# ============================================================================

target_include_directories(qmat_cli PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebra
    ${CMAKE_CURRENT_SOURCE_DIR}/NonManifoldMesh
    ${CMAKE_CURRENT_SOURCE_DIR}/ColorRamp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeometryObjects
    ${CMAKE_CURRENT_SOURCE_DIR}/matfp
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

if(USE_MATSTRUCT_INITIALIZATION)
    target_compile_definitions(qmat_cli PRIVATE USE_MATSTRUCT_INITIALIZATION)
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
message(STATUS "  MatStruct init  : ${USE_MATSTRUCT_INITIALIZATION}")
message(STATUS "  matfp sources   : matfp/sharp_feature_detection.cpp")
message(STATUS "")
