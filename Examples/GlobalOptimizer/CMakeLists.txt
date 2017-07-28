if (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    add_compile_options(-Wno-missing-prototypes)
endif()

add_executable(GlobalOptimizer GlobalOptimizer.cpp)
target_link_libraries(GlobalOptimizer QuantLib)
