Tests/Kernel/CMakeLists.txtadd_executable (FatrasScatteringTests ScatteringTests.cpp)
target_include_directories (FatrasScatteringTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (FatrasScatteringTests PRIVATE ACTSCore ACTSLegacy FatrasCore)
add_test (NAME FatrasScatteringUnitTest COMMAND FatrasScatteringTests)

add_executable (FatrasEnergyLossTests EnergyLossTests.cpp)
target_include_directories (FatrasEnergyLossTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (FatrasEnergyLossTests PRIVATE ACTSCore ACTSLegacy FatrasCore)
add_test (NAME FatrasEnergyLossUnitTest COMMAND FatrasEnergyLossTests)


