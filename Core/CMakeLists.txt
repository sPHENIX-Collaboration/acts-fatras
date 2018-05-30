# main library
add_library(FatrasCore SHARED src/RandomNumberDistributions.cpp)
target_include_directories(FatrasCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(FatrasCore PUBLIC ACTSCore)
set_target_properties(FatrasCore PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED TRUE)

# install library and include files to default directories
include(GNUInstallDirs)
install(
  TARGETS FatrasCore
  EXPORT FatrasCoreTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/Fatras
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Fatras)
target_include_directories(FatrasCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
  



