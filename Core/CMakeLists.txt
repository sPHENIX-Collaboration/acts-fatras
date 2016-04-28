find_package(ACTS)

# get list of all source files
file (GLOB src_files "src/*/*.cxx")

# define library target
add_library(FATRASCore SHARED ${src_files})

# setup include directories
target_include_directories(FATRASCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(FATRASCore PUBLIC ACTSCore)

# set compiler flags
target_compile_options(FATRASCore PUBLIC -std=c++14)

# set installation directories
install (TARGETS FATRASCore
	 LIBRARY DESTINATION lib
	 RUNTIME DESTINATION bin)

install (DIRECTORY include/FATRAS DESTINATION include)

