In this project, I port [QuantLib 1.10](https://github.com/lballabio/QuantLib) to C++17, and with **all** the Boost dependency removed. To achieve this goal, I replace Boost features with their STL counterparts as much as possible. If there is no STL replacements readily available, I either try to use alternative algorithms that gets around the particular Boost features, or re-implement the Boost feature in terms of the standard library. 

Besides removing Boost dependencies, I also try to modernize the codebase as much as I can. Deprecated language features such as `auto_ptr`, `bind1st`, `mem_fun` are removed. C arrays are replaced by std::vectors or other appropriate STL containers. I also use modern features such as lambdas, move semantics, and initializer lists as much as I can.

Both performance and correctness are extremely important. With every Boost feature removed, I make sure there is no single test failure in the test suite, except test cases that also fail in unmodified QuantLib. I also make sure the total running time of the test suite is not longer than unmodified QuantLib. In fact, the performance of QuantLib-noBoost is better than unmodified QuantLib as of now. Using the QuantLib Benchmark Suite, Quantlib-noBoost gives 3425 mflops, while unmodified QuantLib gives 1440 mflops on my Dell Inspiron 13 i7348 laptop, both running under identical conditions with identical compiler flags.

There are features in the Boost library that's impractical for me to implement on my own. Fortunately these features are either optional, or can be easily got around.

* Boost.Test in the test suite is replaced with [Catch](https://github.com/philsquared/Catch). Catch is a light-weight yet powerful unit test framework that's entirely contained inside a single header. Catch is also tightly integrated with the CMake build system. With the CTest test driver, which is a part of CMake, parallel testing essentially comes free.

* I have implemented my own linear algebra routines to replace Boost uBLAS. However there is no way my naive implementations can compete with any BLAS implementations in term of performance. For increased performance, QuantLib-noBoost can **optionally** link against [Intel® MKL](https://software.intel.com/en-us/mkl), which is freely available and widely regarded as the fastest BLAS implementation available.

* For the part of the code that **optionally** uses Boost.Multiprecision, GCC libquadmath is used as a drop-in replacement. GCC libquadmath is a part of GCC and does not require any additional installation.

## Project status:
* Porting is complete with GCC 7 on GNU/Linux. All Boost dependencies are removed and no QuantLib features are missing.

* Currently CMake is the only supported build tool.
 
* Latest commits in the original QuantLib project are regularly backported here.

* Porting to GNU Autotools and Visual Studio 2017 is planned.
