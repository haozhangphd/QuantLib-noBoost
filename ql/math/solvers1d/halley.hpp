/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*! \file halley.hpp
    \brief Halley 1-D solver
*/

#ifndef quantlib_solver1d_halley_h
#define quantlib_solver1d_halley_h

#include <ql/math/solvers1d/halleysafe.hpp>

namespace QuantLib {

    //! %Halley 1-D solver
    /*! \note This solver requires that the passed function object
              implement a method <tt>Real derivative(Real)</tt> and
	      a method <tt>Real derivative2(Real)</tt>. derivative2
	      should give the ratio between the second derivative 
	      and the first derivative divided by 2.
	      derivateive2 = f''(x)/f'(x)/2.

        \test the correctness of the returned values is tested by
              checking them against known good results.

        \ingroup solvers
    */
    class Halley : public Solver1D<Halley> {
    public:
        template<class F>
        Real solveImpl(const F &f,
                       Real xAccuracy) const {

            Real froot, dfroot, d2froot, dx0, dx;

            froot = f(root_);
            dfroot = f.derivative(root_);
            d2froot = f.derivative2(root_);
            QL_REQUIRE(dfroot != Null<Real>(),
                       "Halley requires function's derivative");
            QL_REQUIRE(d2froot != Null<Real>(),
                       "Halley requires f''(x)/f'(x) / 2");
            ++evaluationNumber_;

            while (evaluationNumber_ <= maxEvaluations_) {
                dx0 = froot / dfroot;
                dx = dx0 / std::max(0.8, std::min(1.2, (1 - dx0 * d2froot)));
                Real root_old = root_;
                root_ -= dx;
                if (lowerBoundEnforced() && root_ < xMin_)  // make sure root_ is always larger than xMin_
                    root_ = xMin_ + 0.5 * (root_old - xMin_);
                if (upperBoundEnforced() && root_ > xMax_)  // make sure root_ is always smaller than xMax_
                    root_ = xMax_ - 0.5 * (xMax_ - root_old);
                // jumped out of brackets, switch to HalleySafe
                if ((xMin_ - root_) * (root_ - xMax_) < 0.0) {
                    HalleySafe s;
                    s.setMaxEvaluations(maxEvaluations_ - evaluationNumber_);
                    return s.solve(f, xAccuracy, root_ + dx, xMin_, xMax_);
                }
                if (std::fabs(dx) < xAccuracy) {
                    f(root_);
                    ++evaluationNumber_;
                    return root_;
                }
                froot = f(root_);
                dfroot = f.derivative(root_);
                d2froot = f.derivative2(root_);
                ++evaluationNumber_;
            }

            QL_FAIL("maximum number of function evaluations ("
                            << maxEvaluations_ << ") exceeded");
        }
    };

}

#endif
