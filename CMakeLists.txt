cmake_minimum_required(VERSION 3.5)

project(ACTSFatras LANGUAGES CXX)

# require C++14 globally
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

add_subdirectory(Core)
add_subdirectory(Tests)
