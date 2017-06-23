/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*
    The implementation of the algorithm was inspired by
    "Numerical Recipes in C", 2nd edition,
    Press, Teukolsky, Vetterling, Flannery, chapter 6
*/

#include <ql/math/incompletegamma.hpp>
#include <ql/math/solvers1d/halley.hpp>

namespace QuantLib {


    Real incompleteGammaFunction(Real a, Real x, Real accuracy,
                                 Integer maxIteration) {

        QL_REQUIRE(a > 0.0, "non-positive a is not allowed");

        QL_REQUIRE(x >= 0.0, "negative x is not allowed");


        if (a > 100) {
            // Use the Gauss-Legendre quadrature representation
            return incompleteGammaFunctionQuadratureRepr(a, x);
        } else if (x < (a + 1.0)) {
            // Use the series representation
            return incompleteGammaFunctionSeriesRepr(a, x,
                                                     accuracy, maxIteration);
        } else {
            // Use the continued fraction representation
            return 1.0 - incompleteGammaFunctionContinuedFractionRepr(a, x,
                                                                      accuracy, maxIteration);
        }

    }


    Real incompleteGammaFunctionSeriesRepr(Real a, Real x, Real accuracy,
                                           Integer maxIteration) {

        if (x == 0.0) return 0.0;

        Real gln = std::lgamma(a);
        Real ap = a;
        Real del = 1.0 / a;
        Real sum = del;
        for (Integer n = 1; n <= maxIteration; n++) {
            ++ap;
            del *= x / ap;
            sum += del;
            if (std::fabs(del) < std::fabs(sum) * accuracy) {
                return sum * std::exp(-x + a * std::log(x) - gln);
            }
        }
        QL_FAIL("Incomplete gamma function cannot be calculated, accuracy not reached");
    }

    Real incompleteGammaFunctionContinuedFractionRepr(Real a, Real x,
                                                      Real accuracy,
                                                      Integer maxIteration) {

        Integer i;
        Real an, b, c, d, del, h;
        Real gln = std::lgamma(a);
        b = x + 1.0 - a;
        c = 1.0 / QL_EPSILON;
        d = 1.0 / b;
        h = d;
        for (i = 1; i <= maxIteration; i++) {
            an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (std::fabs(d) < QL_EPSILON) d = QL_EPSILON;
            c = b + an / c;
            if (std::fabs(c) < QL_EPSILON) c = QL_EPSILON;
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if (std::fabs(del - 1.0) < accuracy) {
                return std::exp(-x + a * std::log(x) - gln) * h;
            }
        }

        QL_FAIL("Incomplete gamma function cannot be calculated, accuracy not reached");
    }

    Real incompleteGammaFunctionQuadratureRepr(Real a, Real x) {

        if (x == 0.0) return 0.0;

        // the abscissas and weights for Gauss-Legendre quadradrure.
        const double y[18] = {0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
                              0.082502225484340941, 0.12007019910960293, 0.16415283300752470, 0.21442376986779355,
                              0.27051082840644336, 0.33199876341447887, 0.39843234186401943, 0.46931971407375483,
                              0.54413605556657973, 0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
                              0.87126389619061517, 0.95698180152629142};

        const double w[18] = {0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
                              0.034213810770299537, 0.040875750923643261, 0.047235083490265582, 0.053244713977759692,
                              0.058860144245324798, 0.064039797355015485, 0.068745323835736408, 0.072941885005653087,
                              0.076598410645870640, 0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
                              0.085346685739338721, 0.085983275670394821};

        Real a1 = a - 1.0;
        Real sqrta1 = std::sqrt(a1);
        Real xu;
        Real gln = std::lgamma(a);

        if (x > a1)
            xu = std::max(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
        else
            xu = std::max(0., std::min(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));

        Real sum = 0;
        for (Integer j = 0; j < 18; ++j) {
            Real t = x + (xu - x) * y[j];
            sum += w[j] * std::exp(-t + a1 * std::log(t) - gln);
        }

        sum *= (xu - x);
        return sum > 0? 1 - sum : -sum;
    }

    Real inverseIncompleteGammaFunction(Real a, Real p, Real accuracy, Integer maxIteration) {

        QL_REQUIRE(a > 0.0, "non-positive a is not allowed");

        QL_REQUIRE(p >= 0.0, "negative p non allowed");

        QL_REQUIRE(p <= 1.0, "p value greater than 1 non allowed");

        Real gln = std::lgamma(a);

        struct {

            const Real &a, p, accuracy, gln;
            const Integer &maxIteration;

            Real operator()(Real x) const { return incompleteGammaFunction(a, x, accuracy, maxIteration) - p; }

            Real derivative(Real x) const { return exp(-x - gln + log(x) * (a - 1.0)); }

            Real derivative2(Real x) const { return 0.5 * (a - 1.0) / x - 0.5; }
        } f = {a, p, accuracy, gln, maxIteration};

        // initial guess based on `Numerical Recipes, 3rd ed, Press et al` section 6.2.1
        Real guess, t, pp;

        if (a > 1.0) {
            pp = (p < 0.5) ? p : 1.0 - p;
            t = std::sqrt(-2. * std::log(pp));
            guess = (2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t;
            if (p < 0.5)
                guess = -guess;
            guess = std::max(1.0e-3, a * std::pow(1.0 - 1.0 / (9.0 * a) - guess / (3.0 * std::sqrt(a)), 3));
        } else {
            t = 1.0 - a * (0.253 + a * 0.12);
            if (p < t)
                guess = std::pow(p / t, 1.0 / a);
            else
                guess = 1.0 - std::log(1.0 - (p - t) / (1.0 - t));
        }

        if (guess < std::numeric_limits<Real>::epsilon())
            return 0;

        HalleySafe s;
        s.setLowerBound(0);
        return s.solve(f, 1e-12, guess, guess / 10, guess + 10);

    }


}
