add_executable (ScatteringTests ScatteringTests.cpp)
target_include_directories (ScatteringTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (ScatteringTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME ScatteringUnitTest COMMAND ScatteringTests)

add_executable (EnergyLossTests EnergyLossTests.cpp)
target_include_directories (EnergyLossTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (EnergyLossTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME EnergyLossUnitTest COMMAND EnergyLossTests)


