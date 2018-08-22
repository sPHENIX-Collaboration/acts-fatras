add_executable (FatrasScatteringTests ScatteringTests.cpp)
target_include_directories (FatrasScatteringTests PRIVATE ${Boost_INCLUDE_DIRS}  ../Common)
target_link_libraries (FatrasScatteringTests PRIVATE ActsCore ActsLegacy FatrasCore)
add_test (NAME FatrasScatteringUnitTest COMMAND FatrasScatteringTests)

add_executable (FatrasEnergyLossTests EnergyLossTests.cpp)
target_include_directories (FatrasEnergyLossTests PRIVATE ${Boost_INCLUDE_DIRS} ../Common)
target_link_libraries (FatrasEnergyLossTests PRIVATE ActsCore ActsLegacy FatrasCore)
add_test (NAME FatrasEnergyLossUnitTest COMMAND FatrasEnergyLossTests)

message(STATUS ${Geant4_USE_FILE})
include(${Geant4_USE_FILE})

file(GLOB_RECURSE src_files "src/*.cc" "include/*.hh")

add_executable (FatrasParametricNuclearIntTests ParametricNuclearIntTests.cpp ${src_files})
target_include_directories(FatrasParametricNuclearIntTests PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(FatrasParametricNuclearIntTests PUBLIC ${Geant4_INCLUDE_DIRS})
target_link_libraries(FatrasParametricNuclearIntTests PRIVATE Acts::Core Acts::Legacy Fatras::Core)
target_link_libraries(FatrasParametricNuclearIntTests PUBLIC ${Geant4_LIBRARIES})
add_test (NAME FatrasParametricNuclearIntUnitTest COMMAND FatrasParametricNuclearIntTests)


