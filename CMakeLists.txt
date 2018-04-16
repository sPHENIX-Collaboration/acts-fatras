cmake_minimum_required(VERSION 3.7)

project(ACTSFatras LANGUAGES CXX)

# could be included in a larger project, e.g. acts-framework, that also
# includes acts-core as a subproject. in this case we do not need to
# explicitely add it here.
if(NOT TARGET ACTSCore)
  find_package(ACTS REQUIRED COMPONENTS Core)
endif()

# require C++14 globally
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

add_subdirectory(Core)
add_subdirectory(Samplers)
add_subdirectory(Tests)

  
# CMake package configuration files
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/ACTSFatrasConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/ACTSFatrasConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake")
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/ACTSFatrasConfig.cmake"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake")
install(
  EXPORT ACTSFatrasTargets
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake")
  