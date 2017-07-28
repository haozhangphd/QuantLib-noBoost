if (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    add_compile_options(-Wno-missing-prototypes)
endif()

add_executable(Gaussian1dModels Gaussian1dModels.cpp)
target_link_libraries(Gaussian1dModels QuantLib)
