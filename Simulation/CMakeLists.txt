find_package(ACTS)

# get list of all source files
file (GLOB src_files "src/*.cxx")

# define library target
add_library(FATRAS SHARED ${src_files})

# setup include directories
target_include_directories(FATRAS PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(FATRAS PUBLIC ACTSCore)

# set compiler flags
target_compile_options(FATRAS PUBLIC -std=c++14)

# set installation directories
install (TARGETS FATRAS
	 LIBRARY DESTINATION lib
	 RUNTIME DESTINATION bin)

install (DIRECTORY include/FATRAS DESTINATION include)

