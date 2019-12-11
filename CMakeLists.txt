# 3.9 required to use per-target c++ std requirement
cmake_minimum_required(VERSION 3.9)

project(ActsFatras LANGUAGES CXX)

# could be included in a larger project, e.g. acts-framework, that also
# includes acts-core as a subproject. in this case we do not need to
# explicitely add it here.
if(NOT TARGET ActsCore)
  find_package(Acts REQUIRED COMPONENTS Core)
endif()
# heterogeneous lookup in set-like containers requires 1.68
find_package(Boost 1.68 REQUIRED COMPONENTS unit_test_framework)

include(GNUInstallDirs)

add_subdirectory(Core)

# always build unit tests; must be set from the main CMakeLists.txt
enable_testing()
add_subdirectory(Tests)

# CMake package configuration files
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/ActsFatrasConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/ActsFatrasConfig.cmake")
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/ActsFatrasConfig.cmake"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake")
install(
  EXPORT ActsFatrasTargets
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras")
