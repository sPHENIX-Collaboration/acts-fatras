# add a unittest executable w/ default dependencies and register it
macro(add_unittest _name)
  # for now we use the same name also for the target
  set(_target "ActsFatras${_name}")
  # assume source file and target share the name
  add_executable(${_target} "${_name}.cpp" ${ARGN})
  target_include_directories(
    ${_target}
    PRIVATE "${PROJECT_SOURCE_DIR}/Tests/Common")
  target_link_libraries(
    ${_target}
    PRIVATE ActsCore ActsFatras Boost::unit_test_framework)
  # register as unittest executable
  add_test(NAME ${_target} COMMAND ${_target})
endmacro()

add_subdirectory(Kernel)
add_subdirectory(Physics)
add_subdirectory(Selectors)
