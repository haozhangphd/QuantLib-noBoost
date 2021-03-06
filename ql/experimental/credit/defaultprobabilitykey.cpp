/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2009 StatPro Italia srl
 Copyright (C) 2009 Jose Aparicio

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

#include <ql/experimental/credit/defaultprobabilitykey.hpp>
#include <algorithm>
#include <set>

namespace QuantLib {

    bool operator==(const DefaultProbKey& lhs, const DefaultProbKey& rhs) {
        if(lhs.seniority() != rhs.seniority()) return false;
        if(lhs.currency() != rhs.currency()) return false;

        Size mySize = rhs.eventTypes().size();
        if(mySize != lhs.eventTypes().size()) return false;
        // the all types must be equal in the weak sense.
        for(Size i=0; i<mySize; i++) {
            if(std::find_if(lhs.eventTypes().begin(), lhs.eventTypes().end(),
                [&rhs, i](auto j){return *j == *rhs.eventTypes()[i];}) == lhs.eventTypes().end())
                return false;
        }// naah, I bet this can be done with a double lambda
        return true;
    }

    DefaultProbKey::DefaultProbKey()
        : eventTypes_(),
          obligationCurrency_(Currency()),
          seniority_(NoSeniority) {}

    DefaultProbKey::DefaultProbKey(
        const std::vector<std::shared_ptr<DefaultType> >& eventTypes,
                   const Currency cur,
                   Seniority sen)
        : eventTypes_(eventTypes),
          obligationCurrency_(cur),
          seniority_(sen) {
        std::set<AtomicDefault::Type> buffer;
        Size numEvents = eventTypes_.size();
        for(Size i=0; i< numEvents; i++)
            buffer.insert(eventTypes_[i]->defaultType());
        QL_REQUIRE(buffer.size() == numEvents,
            "Duplicated event type in contract definition");
    }

    NorthAmericaCorpDefaultKey::NorthAmericaCorpDefaultKey(
        const Currency& currency,
        Seniority sen,
        Period graceFailureToPay,
        Real amountFailure,
        Restructuring::Type resType)
    : DefaultProbKey(std::vector<std::shared_ptr<DefaultType> >(),
                     currency, sen) {
        eventTypes_.emplace_back(std::make_shared<FailureToPay>(graceFailureToPay,
            amountFailure));
        // no specifics for Bankruptcy
        eventTypes_.emplace_back(std::make_shared<DefaultType>(AtomicDefault::Bankruptcy,
                            Restructuring::XR));
        if(resType != Restructuring::NoRestructuring)
            eventTypes_.emplace_back(std::make_shared<DefaultType>(AtomicDefault::Restructuring, resType));
    }

}
