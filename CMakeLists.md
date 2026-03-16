 preset (VS2022, x64).

  The approach for Polyscope should be:
  1. Check if polyscope is a port in your local vcpkg install — do you want me to check that now?
  2. If it is: add it to vcpkg.json as a feature-flagged dependency, add the option() + find_package
  block to CMakeLists.txt, and add two new presets (vcpkg-polyscope-debug / vcpkg-polyscope-release)
  with -DQMAT_WITH_POLYSCOPE=ON in cacheVariables
  3. If it isn't in vcpkg: use FetchContent but still install into library/ to match the existing
  convention


Changed:
  - cmake_minimum_required bumped 3.14 → 3.16 (needed for FetchContent_MakeAvailable)
  - Added the option + FetchContent block (only runs when QMAT_WITH_POLYSCOPE=ON)
  - Added target_link_libraries / target_compile_definitions for polyscope (same guard)
  - Summary message now shows Polyscope status
  - CMakePresets.json gets two new presets: vcpkg-polyscope (configure) and vcpkg-polyscope-debug /
  vcpkg-polyscope-release (build)

  Unchanged:
  - Same single qmat_cli target, same sources, same entry point
  - vcpkg + library/ folder untouched — polyscope goes into build/_deps/ (FetchContent's standard location)
  - Default build (vcpkg preset, no flag) is byte-for-byte identical to before

  To build with Polyscope:
  cmake --preset vcpkg-polyscope
  cmake --build build --preset vcpkg-polyscope-release

  To build without (as before):
  cmake --preset vcpkg
  cmake --build build --preset vcpkg-release