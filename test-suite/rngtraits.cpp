/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2004 StatPro Italia srl
 Copyright (C) 2004 Walter Penschke

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
#include <ql/math/randomnumbers/rngtraits.hpp>
#include <ql/math/comparison.hpp>

using namespace QuantLib;


TEST_CASE("RngTraits_Gaussian", "[RngTraits]") {

    INFO("Testing Gaussian pseudo-random number generation...");

    PseudoRandom::rsg_type rsg =
        PseudoRandom::make_sequence_generator(100, 1234);

    const std::vector<Real>& values = rsg.nextSequence().value;
    Real sum = 0.0;
    for (Size i=0; i<values.size(); i++)
        sum += values[i];

    Real stored = 4.09916;
    Real tolerance = 1.0e-5;
    if (std::fabs(sum - stored) > tolerance)
        FAIL("the sum of the samples does not match the stored value\n"
                   << "    calculated: " << sum << "\n"
                   << "    expected:   " << stored);
}


TEST_CASE("RngTraits_DefaultPoisson", "[RngTraits]") {

    INFO("Testing Poisson pseudo-random number generation...");

    PoissonPseudoRandom::icInstance =
        std::shared_ptr<InverseCumulativePoisson>();
    PoissonPseudoRandom::rsg_type rsg =
        PoissonPseudoRandom::make_sequence_generator(100, 1234);

    const std::vector<Real>& values = rsg.nextSequence().value;
    Real sum = 0.0;
    for (Size i=0; i<values.size(); i++)
        sum += values[i];

    Real stored = 108.0;
    if (!close(sum, stored))
        FAIL("the sum of the samples does not match the stored value\n"
                   << "    calculated: " << sum << "\n"
                   << "    expected:   " << stored);
}


TEST_CASE("RngTraits_CustomPoisson", "[RngTraits]") {

    INFO("Testing custom Poisson pseudo-random number generation...");

    PoissonPseudoRandom::icInstance =
        std::make_shared<InverseCumulativePoisson>(4.0);
    PoissonPseudoRandom::rsg_type rsg =
        PoissonPseudoRandom::make_sequence_generator(100, 1234);

    const std::vector<Real>& values = rsg.nextSequence().value;
    Real sum = 0.0;
    for (Size i=0; i<values.size(); i++)
        sum += values[i];

    Real stored = 409.0;
    if (!close(sum, stored))
        FAIL("the sum of the samples does not match the stored value\n"
                   << "    calculated: " << sum << "\n"
                   << "    expected:   " << stored);
}

