cmake_minimum_required (VERSION 2.6)
project (AFTS)

# set the current version of AFTS
set (AFTS_VERSION 0.0.1)

# default installation direction
set (CMAKE_INSTALL_PREFIX "/opt/AFTS/${AFTS_VERSION}" CACHE STRING "installation directory")

add_subdirectory (Core)
