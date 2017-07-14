/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Andreas Gaida
 Copyright (C) 2008 Ralph Schreyer
 Copyright (C) 2008 Klaus Spanderen

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


#include <ql/exercise.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/utilities/fdmdividendhandler.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmsnapshotcondition.hpp>
#include <ql/methods/finitedifferences/utilities/fdminnervaluecalculator.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmamericanstepcondition.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmbermudanstepcondition.hpp>

#include <set>

namespace QuantLib {

    FdmStepConditionComposite::FdmStepConditionComposite(
        const std::list<std::vector<Time> > & stoppingTimes,
        const Conditions & conditions)
    : conditions_(conditions) {

        std::set<Real> allStoppingTimes;
        for (std::list<std::vector<Time> >::const_iterator
             iter = stoppingTimes.begin(); iter != stoppingTimes.end();
             ++iter) {
            allStoppingTimes.insert(iter->begin(), iter->end());
        }
        stoppingTimes_ = std::vector<Time>(allStoppingTimes.begin(),
                                           allStoppingTimes.end());
    }

    const FdmStepConditionComposite::Conditions&
    FdmStepConditionComposite::conditions() const {
        return conditions_;
    }

    const std::vector<Time>& FdmStepConditionComposite::stoppingTimes() const {
        return stoppingTimes_;
    }

    void FdmStepConditionComposite::applyTo(Array& a, Time t) const {
        for (Conditions::const_iterator iter = conditions_.begin();
             iter != conditions_.end(); ++iter) {
            (*iter)->applyTo(a, t);
        }
    }
    
    std::shared_ptr<FdmStepConditionComposite> 
    FdmStepConditionComposite::joinConditions(
                const std::shared_ptr<FdmSnapshotCondition>& c1,
                const std::shared_ptr<FdmStepConditionComposite>& c2) {

        std::list<std::vector<Time> > stoppingTimes;
        stoppingTimes.emplace_back(c2->stoppingTimes());
        stoppingTimes.emplace_back(std::vector<Time>(1, c1->getTime()));

        FdmStepConditionComposite::Conditions conditions;
        conditions.emplace_back(c2);
        conditions.emplace_back(c1);

        return std::make_shared<FdmStepConditionComposite>(stoppingTimes, conditions);
    }

    std::shared_ptr<FdmStepConditionComposite> 
    FdmStepConditionComposite::vanillaComposite(
                 const DividendSchedule& cashFlow,
                 const std::shared_ptr<Exercise>& exercise,
                 const std::shared_ptr<FdmMesher>& mesher,
                 const std::shared_ptr<FdmInnerValueCalculator>& calculator,
                 const Date& refDate,
                 const DayCounter& dayCounter) {
        
        std::list<std::vector<Time> > stoppingTimes;
        std::list<std::shared_ptr<StepCondition<Array> > > stepConditions;

        if(!cashFlow.empty()) {
            std::shared_ptr<FdmDividendHandler> dividendCondition = 
                std::make_shared<FdmDividendHandler>(cashFlow, mesher,
                                                     refDate, dayCounter, 0);
            stepConditions.emplace_back(dividendCondition);
            stoppingTimes.emplace_back(dividendCondition->dividendTimes());
        }

        QL_REQUIRE(   exercise->type() == Exercise::American
                   || exercise->type() == Exercise::European
                   || exercise->type() == Exercise::Bermudan,
                   "exercise type is not supported");
        if (exercise->type() == Exercise::American) {
            stepConditions.emplace_back(std::make_shared<FdmAmericanStepCondition>(mesher,calculator));
        }
        else if (exercise->type() == Exercise::Bermudan) {
            std::shared_ptr<FdmBermudanStepCondition> bermudanCondition =
                std::make_shared<FdmBermudanStepCondition>(exercise->dates(),
                                                           refDate, dayCounter,
                                                           mesher, calculator);
            stepConditions.emplace_back(bermudanCondition);
            stoppingTimes.emplace_back(bermudanCondition->exerciseTimes());
        }
        
        return std::make_shared<FdmStepConditionComposite>(stoppingTimes, stepConditions);

    }

}
