/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2007 Klaus Spanderen

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

#include <ql/math/solvers1d/brent.hpp>
#include <ql/math/functional.hpp>
#include <ql/math/incompletegamma.hpp>
#include <ql/math/distributions/chisquaredistribution.hpp>
#include <ql/math/distributions/gammadistribution.hpp>
#include <ql/math/distributions/normaldistribution.hpp>

#if defined _MSC_VER || defined __clang__
#include <ql/math/modifiedbessel.hpp>
#endif

const long double PIl = 3.141592653589793238462643383279502884L;

namespace QuantLib {

    long double CumulativeChiSquareDistribution::operator()(long double x) const {
        return CumulativeGammaDistribution(0.5 * df_)(0.5 * x);
    }

    long double NonCentralChiSquareDistribution::operator()(long double x) const {
#if defined _MSC_VER || defined __clang__
        if (x == 0)
            return 0;
        if (df_ < 2) {
            long double v = 1 - df_ / 2;
            return 0.5 * std::exp(-0.5 * (x + ncp_)) *
                   std::pow((x / ncp_), df_ / 4 - 0.5) * (
                           modifiedBesselFunction_i_Press(v, std::sqrt(ncp_ * x)) + 2 * std::sin(PIl * v) *
                           modifiedBesselFunction_k_Press(v, std::sqrt(ncp_ * x)) / PIl);
        } 
        return 0.5 * std::exp(-0.5 * (x + ncp_)) *
               std::pow((x / ncp_), df_ / 4 - 0.5) *
               modifiedBesselFunction_i_Press(df_ / 2 - 1, std::sqrt(ncp_ * x));
#else
        if (x == 0)
            return 0;
        if (df_ < 2) {
            long double v = 1 - df_ / 2;
            return 0.5 * std::exp(-0.5 * (x + ncp_)) *
                   std::pow((x / ncp_), df_ / 4 - 0.5) * (
                   std::cyl_bessel_il(v, std::sqrt(ncp_ * x)) + 2 * std::sin(PIl * v) *
                   std::cyl_bessel_kl(v, std::sqrt(ncp_ * x)) / PIl);
        }
        return 0.5 * std::exp(-0.5 * (x + ncp_)) *
               std::pow((x / ncp_), df_ / 4 - 0.5) *
               std::cyl_bessel_il(df_ / 2 - 1, std::sqrt(ncp_ * x));
#endif
    }

    long double CumulativeNonCentralChiSquareDistribution::operator()(long double x) const {
        if (x <= 0.0)
            return 0.0;

        const long double errmax = 1e-14;
        const Size itrmax = 10000;
        long double lam = 0.5 * ncp_;
        long double f2 = 0.5 * df_;
        long double x2 = 0.5 * x;

        // for large x, algorithm described in `Computing discrete mixtures of continuous distributions:
        // noncentral chisquare, noncentral t and the distribution of the square of the sample multiple correlation coeficient.`
        // D. Benton, K. Krishnamoorthy. Computational Statistics & Data Analysis 43 (2003) 249 - 267
        // code adapted from Boost
        if (x > ncp_ + df_) {
            int k = static_cast<int>(std::round(lam));
            long double poisf = std::exp(std::log(lam) * k - lam - std::lgamma(static_cast<long double>(1 + k)));
            long double poisb = poisf * k / lam;
            long double gamf = 1 - incompleteGammaFunction(f2 + k, x2);
            long double xtermf = std::exp(std::log(x2) * (f2 + k) - x2 - std::lgamma(f2 + 1 + k));
            long double xtermb = xtermf * (f2 + k) / x2;
            long double gamb = gamf - xtermb;

            Size i;
            long double sum = -1;
            for (i = k; i - k < itrmax; ++i) {
                long double term = poisf * gamf;
                sum += term;
                poisf *= lam / (i + 1);
                gamf += xtermf;
                xtermf *= x2 / (f2 + i + 1);
                if (((sum == 0) || (fabs(term / sum) < errmax)) && (term >= poisf * gamf))
                    break;
            }
            QL_REQUIRE(i - k < itrmax, "CumulativeNonCentralChiSquareDistribution does not converge for degree of freedom = "
            << df_ <<", non-centrality = " << ncp_ <<", x = " << x );

            for (i = k; i > 0; --i) {
                long double term = poisb * gamb;
                sum += term;
                poisb *= (i - 1) / lam;
                xtermb *= (f2 + i - 1) / x2;
                gamb -= xtermb;
                if ((sum == 0) || (fabs(term / sum) < errmax))
                    break;
            }
            return -sum;
        }

        long double u = std::exp(-lam);
        long double v = u;
        long double f_x_2n = df_ - x + 2.0;

        long double t = 0.0;
        if (f2 * QL_EPSILON > 0.125 &&
            std::fabs(x2 - f2) < std::sqrt(QL_EPSILON) * f2) {
            t = std::exp((1 - t) * (2 - t / (f2 + 1))) / std::sqrt(2.0 * PIl * (f2 + 1.0));
        } else {
            t = std::exp(f2 * std::log(x2) - x2 - std::lgamma(f2 + 1));
        }

        long double ans = v * t;

        long double f_2n = df_ + 2.0;

        for (Size n = 1; n <= itrmax; ++n) {
            u *= lam / n;
            v += u;
            t *= x / f_2n;
            ans += v * t;
            f_2n += 2.0;
            f_x_2n += 2.0;
            if (f_x_2n > 0) {
                if (t * x / f_x_2n <= errmax)
                    return ans;
            }
        }
        QL_FAIL("CumulativeNonCentralChiSquareDistribution does not converge for degree of freedom = "
                << df_ <<", non-centrality = " << ncp_ <<", x = " << x );
    }

    InverseCumulativeNonCentralChiSquare::
    InverseCumulativeNonCentralChiSquare(long double df, long double ncp,
                                         Size maxEvaluations,
                                         long double accuracy)
            : nonCentralDist_(df, ncp),
              guess_(df + ncp),
              maxEvaluations_(maxEvaluations),
              accuracy_(accuracy) {
    }

    long double InverseCumulativeNonCentralChiSquare::operator()(long double x) const {

        // first find the right side of the interval
        Size evaluations = maxEvaluations_;
        long double lower, upper;
        lower = upper = guess_;
        long double step = 2;
        if (nonCentralDist_(upper) < x) {
            do {
                upper *= step;
                --evaluations;
            } while (nonCentralDist_(upper) < x && evaluations > 0);
        } else {
            do {
                lower /= step;
                --evaluations;
                step += 1.0;
            } while (nonCentralDist_(lower) > x && evaluations > 0);
        }

        // use a Brent solver for the rest
        Brent solver;
        solver.setMaxEvaluations(evaluations);
        return solver.solve(
                [this, x](long double y) { return this->nonCentralDist_(y) - x; },
                accuracy_, guess_,
                lower,
                upper);
    }
}

