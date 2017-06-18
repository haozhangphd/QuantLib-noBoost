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

/*! \file chisquaredistribution.hpp
    \brief Chi-square (central and non-central) distributions
*/

#ifndef quantlib_chi_square_distribution_hpp
#define quantlib_chi_square_distribution_hpp

#include <ql/types.hpp>
#include <functional>

namespace QuantLib {

    class CumulativeChiSquareDistribution  {
      public:
        CumulativeChiSquareDistribution(long double df) : df_(df) {}
        long double operator()(long double x) const;
      private:
        long double df_;
    };

    class NonCentralChiSquareDistribution {
    public:
        NonCentralChiSquareDistribution(long double df, long double ncp)
                : df_(df), ncp_(ncp) {}
        long double operator()(long double x) const;
    private:
        long double df_, ncp_;
    };

    class CumulativeNonCentralChiSquareDistribution {
      public:
        CumulativeNonCentralChiSquareDistribution(long double df, long double ncp)
        : df_(df), ncp_(ncp) {}
        long double operator()(long double x) const;
      private:
        long double df_, ncp_;
    };

    class InverseCumulativeNonCentralChiSquare {
      public:
        InverseCumulativeNonCentralChiSquare(long double df, long double ncp,
                                               Size maxEvaluations=100,
                                               long double accuracy = 1e-10);
        long double operator()(long double x) const;

    private:
        CumulativeNonCentralChiSquareDistribution nonCentralDist_;
        const long double guess_;
        const Size maxEvaluations_;
        const long double accuracy_;
    };
}


#endif
