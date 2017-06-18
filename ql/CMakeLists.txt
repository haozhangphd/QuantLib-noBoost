cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE QUANTLIB_FILES "*.hpp" "*.cpp")

add_library(QuantLib SHARED ${QUANTLIB_FILES})
#target_link_libraries(QuantLib quadmath)
