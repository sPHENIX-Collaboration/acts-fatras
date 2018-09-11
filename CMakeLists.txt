cmake_minimum_required(VERSION 3.7)

project(Fatras LANGUAGES CXX)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(USE_GEANT4 "Build Geant4-based code" on)

# could be included in a larger project, e.g. acts-framework, that also
# includes acts-core as a subproject. in this case we do not need to
# explicitely add it here.

if(NOT TARGET ActsCore)
  find_package(Acts REQUIRED COMPONENTS Core Legacy)
endif()

include(GNUInstallDirs)
include(FatrasComponentsHelpers)
include(CTest)

# Geant4 specific configuration
if(USE_GEANT4)
  set(Geant4_CONFIG_DEBUG on)
  find_package(Geant4 REQUIRED gdml)
endif()

add_subdirectory(Core)
add_subdirectory(Tests)
if(USE_GEANT4)
    add_subdirectory(Plugins)
endif()

# CMake package configuration files
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/FatrasConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/FatrasConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake")
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/FatrasConfig.cmake"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake")
install(
  EXPORT FatrasTargets
  NAMESPACE Fatras::
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras")

