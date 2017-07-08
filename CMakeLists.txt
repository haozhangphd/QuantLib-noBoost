cmake_minimum_required(VERSION 3.8)
project(QuantLib)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
################################################################
# User configuration section:                                  #
# modify the following definitions to suit your preferences.   #
################################################################
option(QL_ERROR_FUNCTIONS "if error messages should include current function information" OFF)
option(QL_ERROR_LINES "if error messages should include file and line information" OFF)
option(QL_ENABLE_TRACING "if tracing messages should be allowed" OFF)
option(QL_NEGATIVE_RATES "if negative rates should be allowed" ON)
option(QL_EXTRA_SAFETY_CHECKS "if extra safety checks should be performed (degrade performance)" OFF)
option(QL_USE_INDEXED_COUPON "if indexed coupons instead of par coupons in floating legs should be used" OFF)
option(QL_ENABLE_SESSIONS "if singletons should return different instances for different sessions" OFF)
option(QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN "if thread-safe observer pattern should be enabled (for use in environment with an async garbage collector" OFF)
option(QL_HIGH_RESOLUTION_DATE "if date resolution down to nanoseconds should be enabled" OFF)
option(QL_ENABLE_SINGLETON_THREAD_SAFE_INIT "if singleton initialization shoudl be made thread-safe" OFF)
option(QL_USE_MKL "if Intel® MKL should be used for linear algebra routines" OFF)
option(MULTIPRECISION_NON_CENTRAL_CHI_SQUARED_QUADRATURE "if multiprecision library should be used to improve the precision of the nonr-central chi-squared Gaussian quadrature" OFF)

# compiler flags for GCC
if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
    set(release_flags -O3)
    set(debug_flags -g -O0 -fvar-tracking-assignments)
    add_compile_options(-Werror -march=native -pipe -fstack-protector-strong)

# compiler flags for Clang
elseif (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    set(release_flags -O3)
    set(debug_flags -g -O0)
    add_compile_options(-Werror -march=native -pipe -fstack-protector-strong)

endif()

################################################################
# End of user configuration section                            #
################################################################

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    add_compile_options(-stdlib=libc++)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
endif()

add_compile_options("$<$<CONFIG:Release>:${release_flags}>")
add_compile_options("$<$<CONFIG:Debug>:${debug_flags}>")
if (CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -s")
endif()

if (QL_USE_MKL)
    file(GLOB_RECURSE MKLPATH /opt/intel/*/mkl.h)
    if (MKLPATH)
        string(REPLACE "/include/mkl.h" "" MKLROOT ${MKLPATH})
        include_directories(${MKLROOT}/include)
    else ()
        message(FATAL_ERROR "MKL not found, build aborted!")
    endif ()
endif ()

add_subdirectory(ql)
add_subdirectory(test-suite)
add_subdirectory(Examples)
