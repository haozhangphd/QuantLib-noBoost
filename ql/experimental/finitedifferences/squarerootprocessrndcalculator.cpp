/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Johannes Göttker-Schnetmann
 Copyright (C) 2015 Klaus Spanderen

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


#include <ql/experimental/finitedifferences/squarerootprocessrndcalculator.hpp>
#include <ql/math/distributions/chisquaredistribution.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace QuantLib {

    SquareRootProcessRNDCalculator::SquareRootProcessRNDCalculator(
            Real v0, Real kappa, Real theta, Real sigma)
            : v0_(v0), kappa_(kappa), theta_(theta),
              d_(4 * kappa / (sigma * sigma)), df_(d_ * theta) {}


    Real SquareRootProcessRNDCalculator::pdf(Real v, Time t) const {
        const Real e = std::exp(-kappa_ * t);
        const Real k = d_ / (1 - e);
        const Real ncp = k * v0_ * e;

        return NonCentralChiSquareDistribution(df_, ncp)(v * k) * k;
    }

    Real SquareRootProcessRNDCalculator::cdf(Real v, Time t) const {
        const long double e = std::exp(-kappa_ * t);
        const long double k = d_ / (1 - e);
        const long double ncp = k * v0_ * e;

        return CumulativeNonCentralChiSquareDistribution(df_, ncp)(v * k);
    }

    Real SquareRootProcessRNDCalculator::invcdf(Real q, Time t) const {
        const Real e = std::exp(-kappa_ * t);
        const Real k = d_ / (1 - e);
        const Real ncp = k * v0_ * e;

        return InverseCumulativeNonCentralChiSquare(df_, ncp)(q) / k;
    }

    Real SquareRootProcessRNDCalculator::stationary_pdf(Real v) const {
        const Real alpha = 0.5 * df_;
        const Real beta = alpha / theta_;

        return std::pow(beta, alpha) * std::pow(v, alpha - 1)
               * std::exp(-beta * v - std::lgamma(alpha));
    }

    Real SquareRootProcessRNDCalculator::stationary_cdf(Real v) const {
        const Real alpha = 0.5 * df_;
        const Real beta = alpha / theta_;

        return boost::math::gamma_p(alpha, beta * v);
    }

    Real SquareRootProcessRNDCalculator::stationary_invcdf(Real q) const {
        const Real alpha = 0.5 * df_;
        const Real beta = alpha / theta_;

        return boost::math::gamma_p_inv(alpha, q) / beta;
    }
}
