find_package(GTest CONFIG REQUIRED)

add_executable(sknano_taupol_tests
    TauPolarimeterTest.cc
    GenTauBuilderTest.cc)

target_link_libraries(sknano_taupol_tests PRIVATE
    TauPolarization
    GTest::gtest
    GTest::gtest_main)

add_test(NAME taupol_unit COMMAND sknano_taupol_tests)
