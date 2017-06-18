cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# to reference headers via <ql/foo.hpp>, we need to add the root
# directory of the project to includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(Examples)
add_subdirectory(ql)

add_subdirectory(test-suite)

target_compile_options(QuantLib PUBLIC $<$<CONFIG:Release>:-Werror -O3 -s -march=native -pipe -fstack-protector-strong>
        $<$<CONFIG:Debug>:-Werror -O0 -g -march=native -pipe -fstack-protector-strong -fvar-tracking-assignments>)

