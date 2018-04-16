add_library(ACTSFatrasSamplers INTERFACE)
target_include_directories(ACTSFatrasSamplers INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ACTSFatrasSamplers INTERFACE ACTS::ACTSCore ACTSFatras)

install(TARGETS ACTSFatrasSamplers LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/Fatras DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Fatras)
