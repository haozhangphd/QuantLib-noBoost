if (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    add_compile_options(-Wno-missing-prototypes)
endif()

add_executable(BasketLosses BasketLosses.cpp)
target_link_libraries(BasketLosses QuantLib)
