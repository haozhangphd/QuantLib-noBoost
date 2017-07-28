/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include "utilities.hpp"
#include <ql/utilities/transformiterator.hpp>

#include <vector>
#include <numeric>

using namespace QuantLib;

TEST_CASE("PermutationIterator_Construction", "[.]") {

    INFO("Testing permutation iterator constructors...");

    std::vector<int> v(10);
    std::iota(v.begin(), v.end(), 1);

    std::vector<int> index{3, 1, 4, 1, 5, 9, 2, 6, 5, 4};

    using permIter = permutationIterator<std::vector<int>::iterator, std::vector<int>::iterator>;
    permIter begin(v.begin(), index.begin());
    permIter it = begin;
    permIter end(v.begin(), index.end());

    std::vector<int> result(begin, end);
    std::vector<int> expected{4, 2, 5, 2, 6, 10, 3, 7, 6, 5};

    CHECK(result == expected);
}

TEST_CASE("PermutationIterator_Iteration", "[.]") {
    INFO("Testing permutation iterator iteration...");

    std::vector<int> v(10);
    std::iota(v.begin(), v.end(), 1);

    std::vector<int> index{3, 1, 4, 1, 5, 9, 2, 6, 5, 4};

    using permIter = permutationIterator<std::vector<int>::iterator, std::vector<int>::iterator>;
    permIter begin(v.begin(), index.begin());
    permIter it = begin;
    permIter end(v.begin(), index.end());

    std::vector<int> ret(begin, end);

    SECTION("Elements at even i in the permutation : ") {
        it = begin;
        std::vector<int> result;
        std::vector<int> expected{4, 5, 6, 3, 6};
        for (int i = 0; i < 5; ++i, it += 2)
            result.push_back(*it);
        CHECK(result == expected);
    }

    SECTION("Permutation backwards : ") {
        it = begin + 9;
        std::vector<int> result;
        std::vector<int> expected{5, 6, 7, 3, 10, 6, 2, 5, 2, 4};
        for (int i = 0; i < 10; ++i, --it)
            result.push_back(*it);
        CHECK(result == expected);
    }

    SECTION("Iterate backward with stride 2 : ") {
        it = begin + 9;
        std::vector<int> result;
        std::vector<int> expected{5, 7, 10, 2, 2};
        for (int i = 0; i < 5; ++i, it -= 2)
            result.push_back(*it);
        CHECK(result == expected);
    }
}

TEST_CASE("PermutationIterator_Assignment", "[.]") {

    INFO("Testing assigning values to permutation iterators...");

    std::vector<int> v(10);
    std::iota(v.begin(), v.end(), 1);

    std::vector<int> index{3, 1, 4, 1, 5, 9, 2, 6, 5, 4};

    using permIter = permutationIterator<std::vector<int>::iterator, std::vector<int>::iterator>;
    permIter begin(v.begin(), index.begin());
    permIter it = begin;
    permIter end(v.begin(), index.end());

    *begin = -1;
    begin[1] = -2;
    begin[2] = -3;

    std::vector<int> result(begin, end);
    std::vector<int> expected{-1, -2, -3, -2, 6, 10, 3, 7, 6, -3};

    CHECK(result == expected);
}
