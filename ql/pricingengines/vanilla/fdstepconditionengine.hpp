/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2005 Joseph Wang
 Copyright (C) 2007, 2009 StatPro Italia srl

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

/*! \file fdstepconditionengine.hpp
    \brief Finite-differences step-condition engine
*/

#ifndef quantlib_fd_step_condition_engine_hpp
#define quantlib_fd_step_condition_engine_hpp

#include <ql/pricingengines/vanilla/fdvanillaengine.hpp>
#include <ql/instruments/oneassetoption.hpp>
#include <ql/methods/finitedifferences/fdtypedefs.hpp>
#include <ql/methods/finitedifferences/boundarycondition.hpp>
#include <ql/pricingengines/blackcalculator.hpp>

namespace QuantLib {

    //! Finite-differences pricing engine for American-style vanilla options
    /*! \ingroup vanillaengines */
    template <template <class> class Scheme = CrankNicolson>
    class FDStepConditionEngine :  public FDVanillaEngine {
      public:
        FDStepConditionEngine(
             const std::shared_ptr<GeneralizedBlackScholesProcess>& process,
             Size timeSteps, Size gridPoints,
             bool timeDependent = false)
        : FDVanillaEngine(process, timeSteps, gridPoints, timeDependent),
          controlBCs_(2), controlPrices_(gridPoints) {}
      protected:
        mutable std::shared_ptr<StandardStepCondition> stepCondition_;
        mutable SampledCurve prices_;
        mutable TridiagonalOperator controlOperator_;
        mutable std::vector<std::shared_ptr<bc_type> > controlBCs_;
        mutable SampledCurve controlPrices_;
        virtual void initializeStepCondition() const = 0;
        virtual void calculate(PricingEngine::results*) const;
    };


    // template definitions

    template <template <class> class Scheme>
    void FDStepConditionEngine<Scheme>::calculate(
                                            PricingEngine::results* r) const {
        OneAssetOption::results * results =
            dynamic_cast<OneAssetOption::results *>(r);
        setGridLimits();
        initializeInitialCondition();
        initializeOperator();
        initializeBoundaryConditions();
        initializeStepCondition();

        typedef FiniteDifferenceModel<ParallelEvolver<
                    Scheme<TridiagonalOperator> > > model_type;

        typename model_type::operator_type operatorSet;
        typename model_type::array_type arraySet;
        typename model_type::bc_set bcSet;
        typename model_type::condition_type conditionSet;

        prices_ = intrinsicValues_;

        controlPrices_ = intrinsicValues_;
        controlOperator_ = finiteDifferenceOperator_;
        controlBCs_[0] = BCs_[0];
        controlBCs_[1] = BCs_[1];

        operatorSet.emplace_back(finiteDifferenceOperator_);
        operatorSet.emplace_back(controlOperator_);

        arraySet.emplace_back(prices_.values());
        arraySet.emplace_back(controlPrices_.values());

        bcSet.emplace_back(BCs_);
        bcSet.emplace_back(controlBCs_);

        conditionSet.emplace_back(stepCondition_);
        conditionSet.emplace_back(std::shared_ptr<StandardStepCondition>(
                                                   new NullCondition<Array>));

        model_type model(operatorSet, bcSet);

        model.rollback(arraySet, getResidualTime(),
                       0.0, timeSteps_, conditionSet);

        prices_.values() = arraySet[0];
        controlPrices_.values() = arraySet[1];

        std::shared_ptr<StrikedTypePayoff> striked_payoff =
            std::dynamic_pointer_cast<StrikedTypePayoff>(payoff_);
        QL_REQUIRE(striked_payoff, "non-striked payoff given");

        Real variance =
            process_->blackVolatility()->blackVariance(
                                     exerciseDate_, striked_payoff->strike());
        DiscountFactor dividendDiscount =
            process_->dividendYield()->discount(exerciseDate_);
        DiscountFactor riskFreeDiscount =
            process_->riskFreeRate()->discount(exerciseDate_);
        Real spot = process_->stateVariable()->value();
        Real forwardPrice = spot * dividendDiscount / riskFreeDiscount;

        BlackCalculator black(striked_payoff, forwardPrice,
                              std::sqrt(variance), riskFreeDiscount);

        results->value = prices_.valueAtCenter()
            - controlPrices_.valueAtCenter()
            + black.value();
        results->delta = prices_.firstDerivativeAtCenter()
            - controlPrices_.firstDerivativeAtCenter()
            + black.delta(spot);
        results->gamma = prices_.secondDerivativeAtCenter()
            - controlPrices_.secondDerivativeAtCenter()
            + black.gamma(spot);
        results->additionalResults["priceCurve"] = prices_;
    }

}


#endif
