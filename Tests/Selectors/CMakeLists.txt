add_executable (FatrasKinematicCastsTests KinematicCastsTests.cpp)
target_include_directories (FatrasKinematicCastsTests PRIVATE ${Boost_INCLUDE_DIRS} ../Common)
target_link_libraries (FatrasKinematicCastsTests PRIVATE Acts::Core Fatras::Core)
add_test (NAME FatrasKinematicCastsUnitTest COMMAND FatrasKinematicCastsTests)

add_executable (FatrasSelectorHelpersTests SelectorHelpersTests.cpp)
target_include_directories (FatrasSelectorHelpersTests PRIVATE ${Boost_INCLUDE_DIRS} ../Common)
target_link_libraries (FatrasSelectorHelpersTests PRIVATE Acts::Core Fatras::Core)
add_test (NAME FatrasSelectorHelpersUnitTest COMMAND FatrasSelectorHelpersTests)

add_executable (FatrasPdgSelectorsTests PdgSelectorsTests.cpp)
target_include_directories (FatrasPdgSelectorsTests PRIVATE ${Boost_INCLUDE_DIRS} ../Common)
target_link_libraries (FatrasPdgSelectorsTests PRIVATE Acts::Core Fatras::Core)
add_test (NAME FatrasPdgSelectorsUnitTest COMMAND FatrasPdgSelectorsTests)

add_executable (FatrasLimitSelectorsTests LimitSelectorsTests.cpp)
target_include_directories (FatrasLimitSelectorsTests PRIVATE ${Boost_INCLUDE_DIRS} ../Common)
target_link_libraries (FatrasLimitSelectorsTests PRIVATE Acts::Core Fatras::Core)
add_test (NAME FatrasLimitSelectorsUnitTest COMMAND FatrasLimitSelectorsTests)

