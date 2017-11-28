cmake_minimum_required(VERSION 3.5)

project(Fatras VERSION 0.0.1 LANGUAGES CXX)

# dependencies
find_package(ACTS REQUIRED COMPONENTS Core)

include(GNUInstallDirs)

# Fatras library
set(Fatras_SOURCES
  src/EnergyLoss.cpp
  src/RandomNumberDistributions.cpp)
add_library(Fatras SHARED ${Fatras_SOURCES})
target_include_directories(Fatras PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(Fatras PUBLIC ACTS::ACTSCore)
set_property(TARGET Fatras PROPERTY CXX_STANDARD 14)
set_property(TARGET Fatras PROPERTY CXX_STANDARD_REQUIRED TRUE)

# install library and include files to default directories
install(
  TARGETS Fatras
  EXPORT FatrasTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/Fatras DESTINATION include)

# CMake package configuration files
include(CMakePackageConfigHelpers)

configure_package_config_file(
  cmake/FatrasConfig.cmake.in "${PROJECT_BINARY_DIR}/FatrasConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras")
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/FatrasConfigVersion.cmake"
  COMPATIBILITY SameMajorVersion)

install(
  EXPORT FatrasTargets
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras")
install(
  FILES "${PROJECT_BINARY_DIR}/FatrasConfig.cmake" "${PROJECT_BINARY_DIR}/FatrasConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras")
