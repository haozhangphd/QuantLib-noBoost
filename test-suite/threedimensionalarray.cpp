/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include "utilities.hpp"
#include <ql/math/threedimensionalarray.hpp>

using namespace QuantLib;


TEST_CASE("ThreeDimensionalArray_Construction", "[.]") {

    INFO("Testing 3d array constructors...");

    threeDimensionalArray arr1(3, 3, 3, 3);
    std::vector<Real> vec{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                          27};
    threeDimensionalArray arr2(3, 3, 3, vec.begin(), vec.end());

    CHECK(arr1[1][2][2] == 3);
    CHECK(arr2(0, 1, 2) == 6);
    CHECK(arr2[2][2][2] == 27);

}

TEST_CASE("ThreeDimensionalArray_Assignment", "[.]") {

    INFO("Testing assigning values to 3d array...");

    threeDimensionalArray arr1(3, 3, 3, 3);
    std::vector<Real> vec{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                          27};
    threeDimensionalArray arr2(3, 3, 3, vec.begin(), vec.end());

    arr1[2][2][2] = arr2[2][1][1];
    arr1[0][1][2] = 3.14;
    arr2(0, 0, 1) = 2.17;

    CHECK(arr1[2][2][2] == 23);
    CHECK(arr1(0, 1, 2) == 3.14);
    CHECK(arr2[0][0][1] == 2.17);

}

TEST_CASE("ThreeDimensionalArray_Iteration", "[.]") {

    INFO("Testing iterating over 3d array...");

    std::vector<Real> vec{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                          27};
    threeDimensionalArray arr(3, 3, 3, vec.begin(), vec.end());

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            for (int k = 0; k < 3; ++ k){
                CHECK(arr[i][j][k] == 1 + k + j * 3 + i * 9);
            }
        }
    }

}
