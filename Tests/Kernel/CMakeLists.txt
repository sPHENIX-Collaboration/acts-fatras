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

add_executable (ScatteringTests ScatteringTests.cpp)
target_include_directories (ScatteringTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (ScatteringTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME ScatteringUnitTest COMMAND ScatteringTests)

add_executable (EnergyLossTests EnergyLossTests.cpp)
target_include_directories (EnergyLossTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (EnergyLossTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME EnergyLossUnitTest COMMAND EnergyLossTests)


