#include <ql/experimental/credit/defaultprobabilitylatentmodel.hpp>

namespace QuantLib {
    namespace detail {

        std::vector<std::vector<int>> bn{{1},
                                         {1, 1},
                                         {1, 2,  1},
                                         {1, 3,  3,   1},
                                         {1, 4,  6,   4,   1},
                                         {1, 5,  10,  10,  5,    1},
                                         {1, 6,  15,  20,  15,   6,    1},
                                         {1, 7,  21,  35,  35,   21,   7,    1},
                                         {1, 8,  28,  56,  70,   56,   28,   8,    1},
                                         {1, 9,  36,  84,  126,  126,  84,   36,   9,    1},
                                         {1, 10, 45,  120, 210,  252,  210,  120,  45,   10,   1},
                                         {1, 11, 55,  165, 330,  462,  462,  330,  165,  55,   11,   1},
                                         {1, 12, 66,  220, 495,  792,  924,  792,  495,  220,  66,   12,   1},
                                         {1, 13, 78,  286, 715,  1287, 1716, 1716, 1287, 715,  286,  78,   13,  1},
                                         {1, 14, 91,  364, 1001, 2002, 3003, 3432, 3003, 2002, 1001, 364,  91,  14,  1},
                                         {1, 15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365, 455, 105, 15, 1}};


        /* from "Algorithm 515: Generation of a Vector from the Lexicographical Index"; Buckles, B. P., and Lybanon,
        M. ACM Transactions on Mathematical Software, Vol. 3, No. 2, June 1977.
        https://stackoverflow.com/questions/561/how-to-use-combinations-of-sets-as-test-data#794  */

	/* TODO: return vector<vector<Size>> for all possible j in one function */

        std::vector<Size> combination(Size n, Size i, Size j) {

            if (i == 1) {
                return std::vector<Size>{j};
            }

            Size l, r, k = 0;
            std::vector<Size> c(i);

            for (l = 0; l < i - 1; ++l) {
                c[l] = (l != 0) ? c[l - 1] : 0;
                do {
                    c[l]++;
                    if (n - c[l] < 15) r = bn[n - c[l]][i - (l + 1)];
                    else r = QuantLib::binomialCoefficient(n - c[l], i - (l + 1));
                    k = k + r;
                } while (k < j);
                k = k - r;
            }
            c[i - 1] = c[i - 2] + j - k;

            return c;
        }
    }
}

