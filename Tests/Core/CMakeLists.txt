add_executable (PhysicsListTests PhysicsListTests.cpp)
target_include_directories (PhysicsListTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (PhysicsListTests PRIVATE ACTSCore)
add_test (NAME PhysicsListUnitTest COMMAND PhysicsListTests)