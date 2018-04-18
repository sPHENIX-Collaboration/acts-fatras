add_executable (KinematicCastsTests KinematicCastsTests.cpp)
target_include_directories (KinematicCastsTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (KinematicCastsTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME KinematicCastsUnitTest COMMAND KinematicCastsTests)

add_executable (SelectorHelpersTests SelectorHelpersTests.cpp)
target_include_directories (SelectorHelpersTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (SelectorHelpersTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME SelectorHelpersUnitTest COMMAND SelectorHelpersTests)

add_executable (PdgSelectorsTests PdgSelectorsTests.cpp)
target_include_directories (PdgSelectorsTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (PdgSelectorsTests PRIVATE ACTSCore ACTSFatras)
add_test (NAME PdgSelectorsUnitTest COMMAND PdgSelectorsTests)
