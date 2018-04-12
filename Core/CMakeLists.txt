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
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Fatras)
target_include_directories(ACTSFatras PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
  



