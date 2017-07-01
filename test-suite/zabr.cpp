/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Peter Caspers

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

#include "utilities.hpp"

#include <ql/termstructures/volatility/sabrsmilesection.hpp>
#include <ql/experimental/volatility/zabrsmilesection.hpp>

using namespace QuantLib;


TEST_CASE("Zabr_Consistency", "[Zabr]") {

    Real tol = 1E-4;

    Real alpha = 0.08;
    Real beta = 0.70;
    Real nu = 0.20;
    Real rho = -0.30;
    Real tau = 5.0;
    Real forward = 0.03;

    SabrSmileSection sabr(tau, forward, std::vector<Real>{alpha, beta, nu, rho});

    ZabrSmileSection<ZabrShortMaturityLognormal> zabr0(
        tau, forward, std::vector<Real>{alpha, beta, nu, rho, 1.0});

    ZabrSmileSection<ZabrShortMaturityNormal> zabr1(
        tau, forward, std::vector<Real>{alpha, beta, nu, rho, 1.0});

    ZabrSmileSection<ZabrLocalVolatility> zabr2(
        tau, forward, std::vector<Real>{alpha, beta, nu, rho, 1.0});

    // for full finite prices reduce the number of intermediate points here
    // below the recommended value to speed up the test
    ZabrSmileSection<ZabrFullFd> zabr3(
        tau, forward, std::vector<Real>{alpha, beta, nu, rho, 1.0},
        std::vector<Real>(), 2);

    Real k = 0.0001;
    while (k <= 0.70) {
        Real c0 = sabr.optionPrice(k);
        Real z0 = zabr0.optionPrice(k);
        Real z1 = zabr1.optionPrice(k);
        Real z2 = zabr2.optionPrice(k);
        Real z3 = zabr3.optionPrice(k);
        if (std::fabs(z0 - c0) > tol)
            FAIL_CHECK("Zabr short maturity lognormal expansion price "
                          "("
                          << z0 << ") deviates from Sabr Hagan 2002 price "
                                   "by " << (z0 - c0));
        if (std::fabs(z1 - c0) > tol)
            FAIL_CHECK("Zabr short maturity normal expansion price "
                          "("
                          << z1 << ") deviates from Sabr Hagan 2002 price "
                                   "by " << (z1 - c0));
        if (std::fabs(z2 - c0) > tol)
            FAIL_CHECK("Zabr local volatility price "
                          "("
                          << z2 << ") deviates from Sabr Hagan 2002 price "
                                   "by " << (z2 - c0));
        if (std::fabs(z3 - c0) > tol)
            FAIL_CHECK("Zabr full finite difference price "
                          "("
                          << z3 << ") deviates from Sabr Hagan 2002 price "
                                   "by " << (z3 - c0));
        k += 0.0001;
    }
}
