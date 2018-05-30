add_executable (FatrasScatteringTests ScatteringTests.cpp)
target_include_directories (FatrasScatteringTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (FatrasScatteringTests PRIVATE Acts::Core Acts::Legacy Fatras::Core)
add_test (NAME FatrasScatteringUnitTest COMMAND FatrasScatteringTests)

add_executable (FatrasEnergyLossTests EnergyLossTests.cpp)
target_include_directories (FatrasEnergyLossTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (FatrasEnergyLossTests PRIVATE Acts::Core Acts::Legacy Fatras::Core)
add_test (NAME FatrasEnergyLossUnitTest COMMAND FatrasEnergyLossTests)


