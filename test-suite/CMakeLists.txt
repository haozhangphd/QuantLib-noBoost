cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB TEST_SUITE_FILES "*.hpp" "*.cpp")

set(BENCHMARK_FILES "quantlibbenchmark.cpp" "americanoption.cpp" "asianoptions.cpp" "barrieroption.cpp"
        "basketoption.cpp" "batesmodel.cpp" "convertiblebonds.cpp" "digitaloption.cpp" "dividendoption.cpp"
        "europeanoption.cpp" "fdheston.cpp" "hestonmodel.cpp" "interpolations.cpp" "jumpdiffusion.cpp"
        "marketmodel_smm.cpp" "marketmodel_cms.cpp" "lowdiscrepancysequences.cpp" "quantooption.cpp" "riskstats.cpp"
        "shortratemodels.cpp" "utilities.cpp" "utilities.hpp" "catch.hpp" "swaptionvolstructuresutilities.hpp")

list(REMOVE_ITEM TEST_SUITE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/quantlibbenchmark.cpp)

add_executable(quantlib-test-suite ${TEST_SUITE_FILES})
target_link_libraries(quantlib-test-suite PRIVATE QuantLib)
add_executable(quantlib-benchmark ${BENCHMARK_FILES})
target_link_libraries(quantlib-benchmark PRIVATE QuantLib)

enable_testing(true)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")
include(ParseAndAddCatchTests)
ParseAndAddCatchTests(quantlib-test-suite)

install(TARGETS quantlib-test-suite quantlib-benchmark RUNTIME DESTINATION bin)
