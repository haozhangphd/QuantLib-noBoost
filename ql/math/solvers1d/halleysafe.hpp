/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*! \file halleysafe.hpp
    \brief Safe (bracketed) Halley 1-D solver
*/

#ifndef quantlib_solver1d_halleysafe_h
#define quantlib_solver1d_halleysafe_h

#include <ql/math/solver1d.hpp>

namespace QuantLib {

    //! safe %Halley 1-D solver
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
    class HalleySafe : public Solver1D<HalleySafe> {
    public:
        template<class F>
        Real solveImpl(const F &f,
                       Real xAccuracy) const {


            Real froot, dfroot, d2froot, dx0, dx, dxold;
            Real xh, xl;

            // Orient the search so that f(xl) < 0
            if (fxMin_ < 0.0) {
                xl = xMin_;
                xh = xMax_;
            } else {
                xh = xMin_;
                xl = xMax_;
            }

            // the "stepsize before last"
            dxold = xMax_ - xMin_;
            // it was dxold=std::fabs(xMax_-xMin_); in Numerical Recipes
            // here (xMax_-xMin_ > 0) is verified in the constructor

            // and the last step
            dx = dxold;

            froot = f(root_);
            dfroot = f.derivative(root_);
            d2froot = f.derivative2(root_);
            QL_REQUIRE(dfroot != Null<Real>(),
                       "HalleySafe requires function's derivative");
            QL_REQUIRE(d2froot != Null<Real>(),
                       "HalleySafe requires f''(x)/f'(x) / 2");
            ++evaluationNumber_;

            while (evaluationNumber_ <= maxEvaluations_) {
                // Bisect if (out of range || not decreasing fast enough)
                if ((((root_ - xh) * dfroot - froot) *
                     ((root_ - xl) * dfroot - froot) > 0.0)
                    || (std::fabs(2.0 * froot) > std::fabs(dxold * dfroot))) {

                    dxold = dx;
                    dx = (xh - xl) / 2.0;
                    root_ = xl + dx;
                } else {
                    dxold = dx;
                    dx0 = froot / dfroot;
                    dx = dx0 / std::max(0.8, std::min(1.2, (1 - dx0 * d2froot)));
                    Real root_old = root_;
                    root_ -= dx;
                    if (lowerBoundEnforced() && root_ < xMin_)  // make sure root_ is always larger than xMin_
                        root_ = xMin_ + 0.5 * (root_old - xMin_);
                    if (upperBoundEnforced() && root_ > xMax_)  // make sure root_ is always smaller than xMax_
                        root_ = xMax_ - 0.5 * (xMax_ - root_old);
                }
                // Convergence criterion
                if (std::fabs(dx) < xAccuracy) {
                    f(root_);
                    ++evaluationNumber_;
                    return root_;
                }
                froot = f(root_);
                dfroot = f.derivative(root_);
                d2froot = f.derivative2(root_);
                ++evaluationNumber_;
                if (froot < 0.0)
                    xl = root_;
                else
                    xh = root_;
            }

            QL_FAIL("maximum number of function evaluations ("
                            << maxEvaluations_ << ") exceeded");
        }
    };

}
#endif

