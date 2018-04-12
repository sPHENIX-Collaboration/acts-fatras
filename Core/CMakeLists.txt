cmake_minimum_required(VERSION 3.7)

project(ACTSFatras LANGUAGES CXX)

# could be included in a larger project, e.g. acts-framework, that also
# includes acts-core as a subproject. in this case we do not need to
# explicitely add it here.
if(NOT TARGET ACTSCore)
  find_package(ACTS REQUIRED COMPONENTS Core)
endif()

# main library
add_library(ACTSFatras SHARED
  src/EnergyLoss.cpp
  src/RandomNumberDistributions.cpp)
target_include_directories(ACTSFatras PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTSFatras PUBLIC ACTSCore)
set_target_properties(ACTSFatras PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED TRUE)

# install library and include files to default directories
include(GNUInstallDirs)
install(
  TARGETS ACTSFatras
  EXPORT ACTSFatrasTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/Fatras
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
