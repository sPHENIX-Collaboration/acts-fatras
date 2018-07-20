include(${Geant4_USE_FILE})

file(GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_library(ACTSFatrasGeant4Plugin SHARED ${src_files})
target_include_directories(ACTSFatrasGeant4Plugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ACTSFatrasGeant4Plugin PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(ACTSFatrasGeant4Plugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ACTSFatrasGeant4Plugin PUBLIC Acts::Core Acts::MaterialPlugin)
target_link_libraries(ACTSFatrasGeant4Plugin PUBLIC ${Geant4_LIBRARIES})
target_link_libraries(ACTSFatrasGeant4Plugin PUBLIC ${ROOT_LIBRARIES})

install(TARGETS ACTSFatrasGeant4Plugin LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/Fatras DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
