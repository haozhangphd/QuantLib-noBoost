/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Roland Lichters

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

#include <ql/experimental/credit/riskybond.hpp>
#include <ql/experimental/credit/loss.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <ql/cashflows/cashflowvectors.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/cashflows/couponpricer.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/time/daycounters/actualactual.hpp>

using namespace std;

namespace QuantLib {

    bool RiskyBond::isExpired() const {
        return detail::simple_event(maturityDate()).hasOccurred();
    }

    void RiskyBond::setupExpired() const {
        Instrument::setupExpired();
    }

    void RiskyBond::performCalculations() const {
        NPV_ = 0;
        Date today = Settings::instance().evaluationDate();
        std::vector<std::shared_ptr<CashFlow> > cf = cashflows();
        Date d1 = effectiveDate();
        for (Size i = 0; i < cf.size(); i++) {
            Date d2 = cf[i]->date();
            if (d2 > today) {
                d1 = max(today , d1);
                Date defaultDate = d1 + (d2-d1)/2;

                Real coupon = cf[i]->amount()
                    * defaultTS_->survivalProbability(d2);
                Real recovery = notional(defaultDate) * recoveryRate_
                    * (defaultTS_->survivalProbability(d1)
                       -defaultTS_->survivalProbability(d2));
                NPV_ += coupon * yieldTS()->discount(d2);
                NPV_ += recovery * yieldTS()->discount(defaultDate);
            }
            d1 = d2;
        }
    }

    Real RiskyBond::riskfreeNPV() const {
        Date today = Settings::instance().evaluationDate();
        Real npv = 0;
        std::vector<std::shared_ptr<CashFlow> > cf = cashflows();
        for (Size i = 0; i < cf.size(); i++) {
            Date d2 = cf[i]->date();
            if (d2 > today)
                npv += cf[i]->amount() * yieldTS()->discount(d2);
        }
        return npv;
    }

    Real RiskyBond::totalFutureFlows() const {
        Date today = Settings::instance().evaluationDate();
        Real flow = 0;
        std::vector<std::shared_ptr<CashFlow> > cf = cashflows();
        for (Size i = 0; i < cf.size(); i++) {
            if (cf[i]->date() > today)
                flow += cf[i]->amount();
        }
        return flow;
    }

    std::vector<std::shared_ptr<CashFlow> > RiskyBond::expectedCashflows() {
        std::vector<std::shared_ptr<CashFlow> > expected;
        std::vector<std::shared_ptr<CashFlow> > cf = cashflows();
        Date today = Settings::instance().evaluationDate();
        Date d1 = effectiveDate();
        for (Size i = 0; i < cf.size(); i++) {
            Date d2 = cf[i]->date();
            if (d2 > today) {
                d1 = max(today , d1);
                Date defaultDate = d1 + (d2-d1)/2;

                Real coupon = cf[i]->amount()
                    * defaultTS_->survivalProbability(d2);
                Real recovery = notional(defaultDate) * recoveryRate_
                    * (defaultTS_->survivalProbability(d1)
                       -defaultTS_->survivalProbability(d2));
                std::shared_ptr<CashFlow> flow1 =
                    std::make_shared<SimpleCashFlow>(coupon, d2);
                expected.emplace_back(flow1);

                std::shared_ptr<CashFlow> flow2 =
                    std::make_shared<SimpleCashFlow>(recovery, defaultDate);
                expected.emplace_back(flow2);
            }
            d1 = d2;
        }
        return expected;
    }

    //------------------------------------------------------------------------
    RiskyFixedBond::RiskyFixedBond(
                            const std::string& name,
                            const Currency& ccy,
                            Real recoveryRate,
                            const Handle<DefaultProbabilityTermStructure>& defaultTS,
                            const Schedule& schedule,
                            Real rate,
                            const DayCounter& dayCounter,
                            BusinessDayConvention paymentConvention,
                            const std::vector<Real>& notionals,
                            const Handle<YieldTermStructure>& yieldTS)
    : RiskyBond(name, ccy, recoveryRate, defaultTS, yieldTS),
          schedule_(schedule),
          rate_(rate),
          dayCounter_(dayCounter),
          // paymentConvention_(paymentConvention),
          notionals_(notionals) {
        // FIXME: Take paymentConvention into account
        std::vector<Date> dates = schedule_.dates();
        Real previousNotional = notionals_.front();
        for (Size i = 1; i < dates.size(); i++) {
            Real currentNotional = (i < notionals_.size() ?
                             notionals_[i] :
                             notionals_.back());
            std::shared_ptr<CashFlow> interest  =
                   std::make_shared<FixedRateCoupon>(dates[i], previousNotional,
                                   rate_, dayCounter_, dates[i-1], dates[i]);
            std::shared_ptr<CashFlow> amortization =
                   std::make_shared<AmortizingPayment>(previousNotional - currentNotional, dates[i]);
            previousNotional = currentNotional;

            leg_.emplace_back(interest);
            interestLeg_.emplace_back(interest);
            if (amortization->amount() != 0){
                leg_.emplace_back(amortization);
                redemptionLeg_.emplace_back(amortization);
            }
        }

        std::shared_ptr<CashFlow> redemption =
                 std::make_shared<Redemption>(previousNotional, schedule_.dates().back());
        leg_.emplace_back(redemption);
        redemptionLeg_.emplace_back(redemption);
    }

    std::vector<std::shared_ptr<CashFlow> > RiskyFixedBond::cashflows() const{
        return leg_;
    }
    std::vector<std::shared_ptr<CashFlow> > RiskyFixedBond::interestFlows() const{
        return interestLeg_;
    }
    std::vector<std::shared_ptr<CashFlow> > RiskyFixedBond::notionalFlows() const{
        return redemptionLeg_;
    }

    Real RiskyFixedBond::notional(Date date) const {
        if (date > maturityDate())
            return 0.0;
        Real ntl = notionals_.front();
        for (Size i = 0; i < schedule_.size(); i++) {
            if (i < notionals_.size() && schedule_[i] <= date)
                ntl = notionals_[i];
            else
                break;
        }
        return ntl;
    }

    Date RiskyFixedBond::effectiveDate() const {
        return schedule_.dates().front();
    }

    Date RiskyFixedBond::maturityDate() const {
        return schedule_.dates().back();
    }

    //------------------------------------------------------------------------
    RiskyFloatingBond::RiskyFloatingBond(
                            std::string name,
                            Currency ccy,
                            Real recoveryRate,
                            Handle<DefaultProbabilityTermStructure> defaultTS,
                            Schedule schedule,
                            std::shared_ptr<IborIndex> index,
                            Integer fixingDays,
                            Real spread,
                            std::vector<Real> notionals,
                            Handle<YieldTermStructure> yieldTS)
    : RiskyBond(name, ccy, recoveryRate, defaultTS, yieldTS),
          schedule_(schedule),
          index_(index),
          fixingDays_(fixingDays),
          spread_(spread),
          notionals_(notionals) {

        // FIXME: Take paymentConvention into account
        std::vector<Date> dates = schedule_.dates();
        Real previousNotional = notionals_.front();
        for (Size i = 1; i < dates.size(); i++) {
            Real currentNotional = (i < notionals_.size() ?
                             notionals_[i] :
                             notionals_.back());
            std::shared_ptr<CashFlow> interest =
                   std::make_shared<IborCoupon>(dates[i], previousNotional, dates[i-1], dates[i],
                              fixingDays_, index_, 1.0, spread_);
            std::shared_ptr<CashFlow> amortization =
                 std::make_shared<AmortizingPayment>(previousNotional - currentNotional, dates[i]);
            previousNotional = currentNotional;

            leg_.emplace_back(interest);
            interestLeg_.emplace_back(interest);
            if (amortization->amount() != 0){
                leg_.emplace_back(amortization);
                redemptionLeg_.emplace_back(amortization);
            }
        }

        std::shared_ptr<CashFlow> redemption =
                 std::make_shared<Redemption>(previousNotional, schedule_.dates().back());
        leg_.emplace_back(redemption);
        redemptionLeg_.emplace_back(redemption);

        std::shared_ptr<IborCouponPricer> fictitiousPricer =
                std::make_shared<BlackIborCouponPricer>(Handle<OptionletVolatilityStructure>());
        setCouponPricer(leg_,fictitiousPricer);
    }

    std::vector<std::shared_ptr<CashFlow> > RiskyFloatingBond::cashflows()
        const {
        return leg_;
    }

    std::vector<std::shared_ptr<CashFlow> > RiskyFloatingBond::interestFlows()
    const {
        return interestLeg_;
    }

    std::vector<std::shared_ptr<CashFlow> > RiskyFloatingBond::notionalFlows()
    const {
        return redemptionLeg_;
    }

    Real RiskyFloatingBond::notional(Date date) const {
        if (date > maturityDate())
            return 0.0;
        Real ntl = notionals_.front();
        for (Size i = 0; i < schedule_.size(); i++) {
            if (i < notionals_.size() && schedule_[i] <= date)
                ntl = notionals_[i];
            else
                break;
        }
        return ntl;
    }

    Date RiskyFloatingBond::effectiveDate() const {
        return schedule_.dates().front();
    }

    Date RiskyFloatingBond::maturityDate() const {
        return schedule_.dates().back();
    }

}

