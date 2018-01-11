cmake_minimum_required(VERSION 3.7)

project(ActsFatras LANGUAGES CXX)

# could be included in a larger project, e.g. acts-framework, that also
# includes acts-core as a subproject. in this case we do not need to
# explicitely add it here.
if(NOT TARGET ACTSCore)
  find_package(ACTS REQUIRED COMPONENTS Core)
endif()

# Fatras library
set(Fatras_SOURCES
  src/EnergyLoss.cpp
  src/RandomNumberDistributions.cpp)
add_library(Fatras SHARED ${Fatras_SOURCES})
target_include_directories(Fatras PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(Fatras PUBLIC ACTSCore)
set_target_properties(Fatras PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED TRUE)

# install library and include files to default directories
include(GNUInstallDirs)
install(
  TARGETS Fatras
  EXPORT FatrasTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/Fatras
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# CMake package configuration files
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/FatrasConfig.cmake.in "${PROJECT_BINARY_DIR}/FatrasConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras")
install(
  EXPORT FatrasTargets
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras")
install(
  FILES "${PROJECT_BINARY_DIR}/FatrasConfig.cmake"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras")
