cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/config.hpp")

file(GLOB_RECURSE QUANTLIB_FILES "*.hpp" "*.cpp")

add_library(QuantLib SHARED ${QUANTLIB_FILES})

if (MULTIPRECISION_NON_CENTRAL_CHI_SQUARED_QUADRATURE)
    target_link_libraries(QuantLib PRIVATE quadmath)
endif ()

if (QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN)
    target_compile_options(QuantLib PUBLIC -pthread)
    target_link_libraries(QuantLib PUBLIC -pthread)
endif ()

if (QL_USE_MKL)
    target_link_libraries(QuantLib PRIVATE -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group gomp)
endif ()
