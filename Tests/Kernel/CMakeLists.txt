add_executable (ParticleTests ParticleTests.cpp)
target_include_directories (ParticleTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (ParticleTests PRIVATE ACTSCore FatrasCore)
add_test (NAME ParticleUnitTest COMMAND ParticleTests)

add_executable (FatrasPhysicsListTests PhysicsListTests.cpp)
target_include_directories (FatrasPhysicsListTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (FatrasPhysicsListTests PRIVATE ACTSCore FatrasCore)
add_test (NAME FatrasPhysicsListUnitTest COMMAND FatrasPhysicsListTests)

add_executable (FatrasSelectorListTests SelectorListTests.cpp)
target_include_directories (FatrasSelectorListTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (FatrasSelectorListTests PRIVATE ACTSCore FatrasCore)
add_test (NAME FatrasSelectorListUnitTest COMMAND FatrasSelectorListTests)

add_executable (FatrasProcessTests ProcessTests.cpp)
target_include_directories (FatrasProcessTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (FatrasProcessTests PRIVATE ACTSCore FatrasCore)
add_test (NAME FatrasProcessUnitTest COMMAND FatrasProcessTests)



