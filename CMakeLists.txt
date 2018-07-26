cmake_minimum_required(VERSION 3.7)

project(Fatras LANGUAGES CXX)

option(USE_GEANT4 "Build Geant4-based code" on)

# could be included in a larger project, e.g. acts-framework, that also
# includes acts-core as a subproject. in this case we do not need to
# explicitely add it here.
if(NOT TARGET Acts::Core)
  find_package(Acts REQUIRED COMPONENTS Core Legacy MaterialPlugin) #TODO: MaterialPlugin should be optional (only needed by G4)
endif()

# simplify handling of optional components
# all arguments after the path are evaluated in the if(...) context
function(add_subdirectory_if path)
  file(RELATIVE_PATH _name ${PROJECT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
  if(${ARGN})
    add_subdirectory(${path})
    message(STATUS "Enable component '${_name}'")
  else()
    message(STATUS "Disable component '${_name}'")
  endif()
endfunction()

# require C++14 globally
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

include(GNUInstallDirs)

include(CTest)
add_subdirectory(Core)


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


  
