# QMAT CLI - CMake Build Configuration
#
# This builds a command-line version of QMAT without Qt dependencies.
#
# Build instructions (using vcpkg manifest mode):
#   cd qmat
#   cmake -B build -DCMAKE_TOOLCHAIN_FILE=C:/Users/alirz/Projects/vcpkg/scripts/buildsystems/vcpkg.cmake
#   cmake --build build --config Release
#
# Dependencies will be installed to: qmat/library/
# The executable will be at: build/Release/qmat_cli.exe (Windows)

cmake_minimum_required(VERSION 3.14)

# Set vcpkg to install packages into "library" folder (must be before project())
set(VCPKG_INSTALLED_DIR "${CMAKE_CURRENT_SOURCE_DIR}/library")

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
# Find Dependencies via vcpkg (manifest mode)
# All packages installed to: ${CMAKE_CURRENT_SOURCE_DIR}/library
# ============================================================================

find_package(CGAL CONFIG REQUIRED)
# find_package(Boost CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

# ============================================================================
# Source Files
# ============================================================================

set(QMAT_CLI_SOURCES
    main_cli.cpp
    Mesh.cpp
    ThreeDimensionalShape.cpp
    SlabMesh.cpp
    PrimMesh.cpp
    NonManifoldMesh/nonmanifoldmesh.cpp
    LinearAlgebra/Wm4Math.cpp
    LinearAlgebra/Wm4Matrix.cpp
    LinearAlgebra/Wm4Vector.cpp
    ColorRamp/ColorRamp.cpp
    GeometryObjects/GeometryObjects.cpp
)

set(QMAT_CLI_HEADERS
    Mesh.h
    ThreeDimensionalShape.h
    SlabMesh.h
    PrimMesh.h
    NonManifoldMesh/nonmanifoldmesh.h
    LinearAlgebra/Wm4Math.h
    LinearAlgebra/Wm4Matrix.h
    LinearAlgebra/Wm4Vector.h
    ColorRamp/ColorRamp.h
    GeometryObjects/GeometryObjects.h
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
)

# ============================================================================
# Link Libraries (vcpkg handles include paths automatically)
# ============================================================================

target_link_libraries(qmat_cli PRIVATE
    CGAL::CGAL
    Eigen3::Eigen
    # Boost::boost
    )

#include the tiny_obj_loader.h library (header only)
target_include_directories(qmat_cli PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tinyobjloader
)
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
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  vcpkg installed dir: ${VCPKG_INSTALLED_DIR}")
message(STATUS "  CGAL: Found")
message(STATUS "  Eigen3: Found")
message(STATUS "")
