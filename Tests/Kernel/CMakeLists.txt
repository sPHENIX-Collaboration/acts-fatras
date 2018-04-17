add_executable (FatrasDefinitionsTests FatrasDefinitionsTests.cpp)
target_include_directories (FatrasDefinitionsTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (FatrasDefinitionsTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME FatrasDefinitionsUnitTest COMMAND FatrasDefinitionsTests)

add_executable (PhysicsListTests PhysicsListTests.cpp)
target_include_directories (PhysicsListTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (PhysicsListTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME PhysicsListUnitTest COMMAND PhysicsListTests)

add_executable (SelectorListTests SelectorListTests.cpp)
target_include_directories (SelectorListTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (SelectorListTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME SelectorListUnitTest COMMAND SelectorListTests)

add_executable (ProcessTests ProcessTests.cpp)
target_include_directories (ProcessTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (ProcessTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME ProcessUnitTest COMMAND ProcessTests)



