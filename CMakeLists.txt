cmake_minimum_required (VERSION 3.5)
project (Fatras)

# Use GNU-style hierarchy for installing build products
include (GNUInstallDirs)

# set the current version of Fatras
set (FATRAS_VERSION 0.0.1)

# convert the version into a plain number to be used in the version header file
string (REPLACE "." "" FATRAS_VERSION_NUMBER ${FATRAS_VERSION})

# default installation direction
set (CMAKE_INSTALL_PREFIX "/opt/AFTS/${AFTS_VERSION}" CACHE STRING "installation directory")

#---------------------------------------------------------

# get list of all source files
file (GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

# find the ACTS package
find_package(ACTS COMPONENTS Core)

# define library target
add_library(Fatras SHARED ${src_files})

# setup include directories
target_include_directories(Fatras PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(Fatras PUBLIC ACTS::ACTSCore)

# set compiler flags
target_compile_options(Fatras PUBLIC -std=c++14)

# set installation directories
install (TARGETS Fatras
         EXPORT FatrasTargets
	 LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})	
	 
# install CMake targets
install (EXPORT FatrasTargets
         DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras") 

install (DIRECTORY include/Fatras DESTINATION include)



#--- CMake Config Files -----------------------------------------------
# - Use CMake's module to help generating relocatable config files
include (CMakePackageConfigHelpers)

# - Versioning
write_basic_package_version_file ("${PROJECT_BINARY_DIR}/FatrasConfigVersion.cmake"
                                  VERSION ${FATRAS_VERSION}
                                  COMPATIBILITY SameMajorVersion)


# - Install time config and target files
configure_package_config_file (cmake/FatrasConfig.cmake.in "${PROJECT_BINARY_DIR}/FatrasConfig.cmake"
                               INSTALL_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras"
                               PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)

# - install and export
install (FILES "${PROJECT_BINARY_DIR}/FatrasConfigVersion.cmake" "${PROJECT_BINARY_DIR}/FatrasConfig.cmake"
         DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Fatras")
