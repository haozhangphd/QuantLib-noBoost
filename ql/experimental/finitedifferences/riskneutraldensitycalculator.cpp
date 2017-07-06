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

#include <ql/math/functional.hpp>
#include <ql/math/solvers1d/brent.hpp>
#include <ql/experimental/finitedifferences/riskneutraldensitycalculator.hpp>
#include <functional>

namespace QuantLib {
    RiskNeutralDensityCalculator::InvCDFHelper::InvCDFHelper(
        const RiskNeutralDensityCalculator* calculator,
        Real guess, Real accuracy, Size maxEvaluations)
    : calculator_(calculator),
      guess_(guess),
      accuracy_(accuracy),
      maxEvaluations_(maxEvaluations) { }

    Real RiskNeutralDensityCalculator::InvCDFHelper::inverseCDF(Real p, Time t)
    const {
        const Real guessCDF = calculator_->cdf(guess_, t);

        Size evaluations = maxEvaluations_;
        Real upper = guess_, lower = guess_;

        if (guessCDF < p)
            while (calculator_->cdf(upper*=1.5, t) < p && evaluations > 0) {
                --evaluations;
            }
        else
            while (calculator_->cdf(lower*=0.75, t) > p && evaluations > 0) {
                --evaluations;
            }

        QL_REQUIRE(evaluations, "could not calculate interval");

        auto cdf = [this,t](Real x){return this->calculator_->cdf(x,t);};

        Brent solver;
        solver.setMaxEvaluations(evaluations);
        return solver.solve([cdf,p](Real x){return cdf(x) - p;},
                            accuracy_, 0.5*(lower + upper), lower, upper);
    }
}
