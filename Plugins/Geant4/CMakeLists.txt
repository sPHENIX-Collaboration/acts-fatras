include(${Geant4_USE_FILE})

file(GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_library(FatrasGeant4Plugin SHARED ${src_files})
target_include_directories(FatrasGeant4Plugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(FatrasGeant4Plugin PUBLIC ${Geant4_INCLUDE_DIRS})
target_link_libraries(FatrasGeant4Plugin PUBLIC ActsCore)
target_link_libraries(FatrasGeant4Plugin PUBLIC ${Geant4_LIBRARIES})

install(TARGETS FatrasGeant4Plugin LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/Fatras DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
