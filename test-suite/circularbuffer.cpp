/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include "utilities.hpp"
#include "circular_buffer.hpp"
#include <vector>
#include <iostream>

using QuantLib::Real;

TEST_CASE("CircularBuffer_Construction", "[.]") {


    INFO("Testing circular buffer constructors...");

    std::vector<Real> vec{1.2, 2.3, 3.4, 4.5, 5.6};

    circular_buffer test(vec.begin(), vec.end());

    CHECK(test[0] == 1.2);
    CHECK(test[2] == 3.4);


}

TEST_CASE("CircularBuffer_Assignment", "[.]") {


    INFO("Testing circular buffer assignment...");

    std::vector<Real> vec{1.2, 2.3, 3.4, 4.5, 5.6};

    circular_buffer test(vec.begin(), vec.end());

    test.push_back(3.14);
    test.push_back(2.17);
    test.push_front(3);
    test[1] = 2;

    circular_buffer expected{3, 2, 4.5, 5.6, 3.14 };

    CHECK(test == expected);

}


TEST_CASE("CircularBuffer_Iteration", "[.]") {


    INFO("Testing circular buffer assignment...");

    std::vector<Real> expected{1.2, 2.3, 3.4, 4.5, 5.6};

    circular_buffer cb(expected.begin(), expected.end());

    std::vector<Real> test;


    for (auto i: cb)
        test.push_back(i);


    CHECK(test == expected);

}
