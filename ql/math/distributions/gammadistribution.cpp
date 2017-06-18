/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 Sadruddin Rejeb

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

#include <ql/math/distributions/gammadistribution.hpp>

namespace QuantLib {

    long double CumulativeGammaDistribution::operator()(long double x) const {
        if (x <= 0.0) return 0.0;

        long double gln = std::lgamma(a_);

        if (x<(a_+1.0)) {
            long double ap = a_;
            long double del = 1.0/a_;
            long double sum = del;
            for (Size n=1; n<=100; n++) {
                ap += 1.0;
                del *= x/ap;
                sum += del;
                if (std::fabs(del) < std::fabs(sum)*1.0e-10)
                    return sum*std::exp(-x + a_*std::log(x) - gln);
            }
        } else {
            long double b = x + 1.0 - a_;
            long double c = QL_MAX_REAL;
            long double d = 1.0/b;
            long double h = d;
            for (Size n=1; n<=100; n++) {
                long double an = -1.0*n*(n-a_);
                b += 2.0;
                d = an*d + b;
                if (std::fabs(d) < QL_EPSILON) d = QL_EPSILON;
                c = b + an/c;
                if (std::fabs(c) < QL_EPSILON) c = QL_EPSILON;
                d = 1.0/d;
                long double del = d*c;
                h *= del;
                if (std::fabs(del - 1.0)<QL_EPSILON)
                    return 1.0-h*std::exp(-x + a_*std::log(x) - gln);
            }
        }
        QL_FAIL("too few iterations");
    }

}
