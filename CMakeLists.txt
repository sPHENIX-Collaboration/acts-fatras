cmake_minimum_required(VERSION 3.7)

project(Fatras LANGUAGES CXX)

option(USE_GEANT4 "Build Geant4-based code" on)

# could be included in a larger project, e.g. acts-framework, that also
# includes acts-core as a subproject. in this case we do not need to
# explicitely add it here.
if(NOT TARGET Acts::Core)
  find_package(Acts REQUIRED COMPONENTS Core Legacy MaterialPlugin)
endif()

# require C++14 globally
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

include(GNUInstallDirs)
include(CTest)

add_subdirectory(Core)
add_subdirectory(Tests)
  
# Geant4 specific configuration
if(USE_GEANT4)
  set(Geant4_CONFIG_DEBUG on)
  find_package(Geant4 REQUIRED gdml)
endif()
add_subdirectory(Plugins)
add_subdirectory(Tests)


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

