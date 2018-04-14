add_executable (PhysicsListTests PhysicsListTests.cpp)
target_include_directories (PhysicsListTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (PhysicsListTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME PhysicsListUnitTest COMMAND PhysicsListTests)

add_executable (FatrasDefinitionsTests FatrasDefinitionsTests.cpp)
target_include_directories (FatrasDefinitionsTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (FatrasDefinitionsTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME FatrasDefinitionsUnitTest COMMAND FatrasDefinitionsTests)
