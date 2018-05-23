cmake_minimum_required(VERSION 3.7)

project(ActsFatras LANGUAGES CXX)

# could be included in a larger project, e.g. acts-framework, that also
# includes acts-core as a subproject. in this case we do not need to
# explicitely add it here.
if(NOT TARGET Acts::Core)
  find_package(Acts REQUIRED COMPONENTS Core)
endif()

# main library
add_library(ActsFatras SHARED
  src/EnergyLoss.cpp
  src/RandomNumberDistributions.cpp)
target_include_directories(ActsFatras PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ActsFatras PUBLIC Acts::Core)
set_target_properties(ActsFatras PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED TRUE)

# install library and include files to default directories
include(GNUInstallDirs)
install(
  TARGETS ActsFatras
  EXPORT ActsFatrasTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/Fatras
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# CMake package configuration files
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/ActsFatrasConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/ActsFatrasConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake")
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/ActsFatrasConfig.cmake"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake")
install(
  EXPORT ActsFatrasTargets
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake")
