/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2012, 2013 Grzegorz Andruszkiewicz

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
#include <ql/types.hpp>
#include <ql/experimental/catbonds/all.hpp>
#include <ql/instruments/bonds/floatingratebond.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/calendars/unitedstates.hpp>
#include <ql/time/calendars/brazil.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <ql/time/daycounters/business252.hpp>
#include <ql/indexes/ibor/usdlibor.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/time/schedule.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/cashflows/couponpricer.hpp>
#include <ql/cashflows/cashflows.hpp>
#include <ql/pricingengines/bond/discountingbondengine.hpp>
#include <ql/pricingengines/bond/bondfunctions.hpp>
#include <ql/termstructures/yield/flatforward.hpp>

using namespace QuantLib;
using std::shared_ptr;

namespace {
    std::pair<Date, Real> data[] = {std::pair<Date, Real>(Date(1, February, 2012), 100), std::pair<Date, Real>(Date(1, July, 2013), 150), std::pair<Date, Real>(Date(5, January, 2014), 50)};
    std::shared_ptr<std::vector<std::pair<Date, Real> > > sampleEvents(new std::vector<std::pair<Date, Real> >(data, data+3));

    Date eventsStart(1, January, 2011);
    Date eventsEnd(31, December, 2014);
    
    template <typename lhs_t, typename rhs_t, typename tolerance_t>
    bool close(lhs_t lhs, rhs_t rhs, tolerance_t tolerance)
    {
	    return std::abs(rhs - lhs) <= tolerance;

    }
}

TEST_CASE( "CatBond_EventSetForWholeYears", "[CatBond]" ) {
    INFO("Testing that catastrophe events are split correctly for periods of whole years...");

	EventSet catRisk(sampleEvents, eventsStart, eventsEnd);
	std::shared_ptr<CatSimulation> simulation = catRisk.newSimulation(Date(1, January, 2015), Date(31, December, 2015));

	REQUIRE(simulation);

	std::vector<std::pair<Date, Real> > path;

	REQUIRE(simulation->nextPath(path));
	CHECK(Size(0) == path.size());

	REQUIRE(simulation->nextPath(path));
	CHECK(Size(1) == path.size());
	CHECK(Date(1, February, 2015) == path.at(0).first);
	CHECK(100 == path.at(0).second);

	REQUIRE(simulation->nextPath(path));
	CHECK(Size(1) == path.size());
	CHECK(Date(1, July, 2015) == path.at(0).first);
	CHECK(150 == path.at(0).second);

	REQUIRE(simulation->nextPath(path));
	CHECK(Size(1) == path.size());
	CHECK(Date(5, January, 2015) == path.at(0).first);
	CHECK(50 == path.at(0).second);

	REQUIRE(!simulation->nextPath(path));
}


TEST_CASE( "CatBond_EventSetForIrregularPeriods", "[CatBond]" ) {
    INFO("Testing that catastrophe events are split correctly for irregular periods...");
	
	EventSet catRisk(sampleEvents, eventsStart, eventsEnd);
	std::shared_ptr<CatSimulation> simulation = catRisk.newSimulation(Date(2, January, 2015), Date(5, January, 2016));

	REQUIRE(simulation);

	std::vector<std::pair<Date, Real> > path;

	REQUIRE(simulation->nextPath(path));
	CHECK(Size(0) == path.size());

	REQUIRE(simulation->nextPath(path));
	CHECK(Size(2) == path.size());
	CHECK(Date(1, July, 2015) == path.at(0).first);
	CHECK(150 == path.at(0).second);
	CHECK(Date(5, January, 2016) == path.at(1).first);
	CHECK(50 == path.at(1).second);

	REQUIRE(!simulation->nextPath(path));
}


TEST_CASE( "CatBond_EventSetForNoEvents", "[CatBond]" ) {
    INFO("Testing that catastrophe events are split correctly when there are no simulated events...");

    std::shared_ptr<std::vector<std::pair<Date, Real> > > emptyEvents(new std::vector<std::pair<Date, Real> >());
    EventSet catRisk(emptyEvents, eventsStart, eventsEnd);
	std::shared_ptr<CatSimulation> simulation = catRisk.newSimulation(Date(2, January, 2015), Date(5, January, 2016));

	REQUIRE(simulation);

	std::vector<std::pair<Date, Real> > path;

	REQUIRE(simulation->nextPath(path));
	CHECK(Size(0) == path.size());

    REQUIRE(simulation->nextPath(path));
	CHECK(Size(0) == path.size());

    REQUIRE(!simulation->nextPath(path));
}

TEST_CASE( "CatBond_BetaRisk", "[CatBond]" ) {
    INFO("Testing that beta risk gives correct terminal distribution...");

    const size_t PATHS = 1000000;
    BetaRisk catRisk(100.0, 100.0, 10.0, 15.0);
    std::shared_ptr<CatSimulation> simulation = catRisk.newSimulation(Date(2, January, 2015), Date(2, January, 2018));
    REQUIRE(simulation);

    std::vector<std::pair<Date, Real> > path;
    Real sum = 0.0;
    Real sumSquares = 0.0;
    Real poissonSum = 0.0;
    Real poissonSumSquares = 0.0;
    
    
    for(size_t i=0; i<PATHS; ++i)
    {
        REQUIRE(simulation->nextPath(path));
        Real processValue = 0.0;
        for(size_t j=0; j<path.size(); ++j) processValue+=path[j].second;
        sum+=processValue;
        sumSquares+=processValue*processValue;
        poissonSum+=path.size();
        poissonSumSquares+=path.size()*path.size();
    }
    Real poissonMean = poissonSum/PATHS;
    CHECK(close(Real(3.0/100.0), poissonMean, 2));
    Real poissonVar = poissonSumSquares/PATHS - poissonMean*poissonMean;
    CHECK(close(Real(3.0/100.0), poissonVar, 5));
    
    Real expectedMean = 3.0*10.0/100.0;
    Real actualMean = sum/PATHS;
    CHECK(close(expectedMean, actualMean, 1));
    
    Real expectedVar = 3.0*(15.0*15.0+10*10)/100.0;
    Real actualVar = sumSquares/PATHS - actualMean*actualMean;
    CHECK(close(expectedVar, actualVar, 1));
}

namespace {

    struct CommonVars {
        // common data
        Calendar calendar;
        Date today;
        Real faceAmount;

        // cleanup
        SavedSettings backup;

        // setup
        CommonVars() {
            calendar = TARGET();
            today = calendar.adjust(Date::todaysDate());
            Settings::instance().evaluationDate() = today;
            faceAmount = 1000000.0;
        }
    };
}

TEST_CASE( "CatBond_RiskFreeAgainstFloatingRateBond", "[CatBond]" ) {
    INFO("Testing floating-rate cat bond against risk-free floating-rate bond...");

    CommonVars vars;

    Date today(22,November,2004);
    Settings::instance().evaluationDate() = today;

    Natural settlementDays = 1;

    Handle<YieldTermStructure> riskFreeRate(flatRate(today,0.025,Actual360()));
    Handle<YieldTermStructure> discountCurve(flatRate(today,0.03,Actual360()));

    shared_ptr<IborIndex> index(new USDLibor(6*Months, riskFreeRate));
    Natural fixingDays = 1;

    Real tolerance = 1.0e-6;

    shared_ptr<IborCouponPricer> pricer(new
        BlackIborCouponPricer(Handle<OptionletVolatilityStructure>()));

    // plain

    Schedule sch(Date(30,November,2004),
                 Date(30,November,2008),
                 Period(Semiannual),
                 UnitedStates(UnitedStates::GovernmentBond),
                 ModifiedFollowing, ModifiedFollowing,
                 DateGeneration::Backward, false);

    std::shared_ptr<CatRisk> noCatRisk(new EventSet(
        std::shared_ptr<std::vector<std::pair<Date, Real> > >(new std::vector<std::pair<Date, Real> >()),
        Date(1, Jan, 2000), Date(31, Dec, 2010)));

    std::shared_ptr<EventPaymentOffset> paymentOffset(new NoOffset());
    std::shared_ptr<NotionalRisk> notionalRisk(new DigitalNotionalRisk(paymentOffset, 100));

    FloatingRateBond bond1(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           ModifiedFollowing, fixingDays,
                           std::vector<Real>(), std::vector<Spread>(),
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

	FloatingCatBond catBond1(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           notionalRisk, 
                           ModifiedFollowing, fixingDays,
                           std::vector<Real>(), std::vector<Spread>(),
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

    shared_ptr<PricingEngine> bondEngine(
                                     new DiscountingBondEngine(riskFreeRate));
    bond1.setPricingEngine(bondEngine);
    setCouponPricer(bond1.cashflows(),pricer);

    shared_ptr<PricingEngine> catBondEngine(new MonteCarloCatBondEngine(noCatRisk, riskFreeRate));
    catBond1.setPricingEngine(catBondEngine);
    setCouponPricer(catBond1.cashflows(),pricer);

    #if defined(QL_USE_INDEXED_COUPON)
    Real cachedPrice1 = 99.874645;
    #else
    Real cachedPrice1 = 99.874646;
    #endif


    Real price = bond1.cleanPrice();
    Real catPrice = catBond1.cleanPrice();
    if (std::fabs(price-cachedPrice1) > tolerance || std::fabs(catPrice-price) > tolerance) {
        FAIL("failed to reproduce floating rate bond price:\n"
                   << std::fixed
                   << "    floating bond: " << price << "\n"
                   << "    catBond bond: " << catPrice << "\n"
                   << "    expected:   " << cachedPrice1 << "\n"
                   << "    error:      " << catPrice-price);
    }

    

    // different risk-free and discount curve

    FloatingRateBond bond2(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           ModifiedFollowing, fixingDays,
                           std::vector<Rate>(), std::vector<Spread>(),
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

    FloatingCatBond catBond2(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           notionalRisk,
                           ModifiedFollowing, fixingDays,
                           std::vector<Rate>(), std::vector<Spread>(),
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

    shared_ptr<PricingEngine> bondEngine2(
                                    new DiscountingBondEngine(discountCurve));
    bond2.setPricingEngine(bondEngine2);
    setCouponPricer(bond2.cashflows(),pricer);

    shared_ptr<PricingEngine> catBondEngine2(new MonteCarloCatBondEngine(noCatRisk, discountCurve));
    catBond2.setPricingEngine(catBondEngine2);
    setCouponPricer(catBond2.cashflows(),pricer);

    #if defined(QL_USE_INDEXED_COUPON)
    Real cachedPrice2 = 97.955904;
    #else
    Real cachedPrice2 = 97.955904;
    #endif

    price = bond2.cleanPrice();
    catPrice = catBond2.cleanPrice();
    if (std::fabs(price-cachedPrice2) > tolerance || std::fabs(catPrice-price) > tolerance) {
        FAIL("failed to reproduce floating rate bond price:\n"
                   << std::fixed
                   << "    floating bond: " << price << "\n"
                   << "    catBond bond: " << catPrice << "\n"
                   << "    expected:   " << cachedPrice2 << "\n"
                   << "    error:      " << catPrice-price);
    }

    // varying spread

    std::vector<Rate> spreads(4);
    spreads[0] = 0.001;
    spreads[1] = 0.0012;
    spreads[2] = 0.0014;
    spreads[3] = 0.0016;

    FloatingRateBond bond3(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           ModifiedFollowing, fixingDays,
                           std::vector<Real>(), spreads,
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

    FloatingCatBond catBond3(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           notionalRisk,
                           ModifiedFollowing, fixingDays,
                           std::vector<Real>(), spreads,
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

    bond3.setPricingEngine(bondEngine2);
    setCouponPricer(bond3.cashflows(),pricer);

    catBond3.setPricingEngine(catBondEngine2);
    setCouponPricer(catBond3.cashflows(),pricer);

    #if defined(QL_USE_INDEXED_COUPON)
    Real cachedPrice3 = 98.495458;
    #else
    Real cachedPrice3 = 98.495459;
    #endif

    price = bond3.cleanPrice();
    catPrice = catBond3.cleanPrice();
    if (std::fabs(price-cachedPrice3) > tolerance || std::fabs(catPrice-price) > tolerance) {
        FAIL("failed to reproduce floating rate bond price:\n"
                   << std::fixed
                   << "    floating bond: " << price << "\n"
                   << "    catBond bond: " << catPrice << "\n"
                   << "    expected:   " << cachedPrice2 << "\n"
                   << "    error:      " << catPrice-price);
    }
}



TEST_CASE( "CatBond_CatBondInDoomScenario", "[CatBond]" ) {
    INFO("Testing floating-rate cat bond in a doom scenario (certain default)...");

    CommonVars vars;

    Date today(22,November,2004);
    Settings::instance().evaluationDate() = today;

    Natural settlementDays = 1;

    Handle<YieldTermStructure> riskFreeRate(flatRate(today,0.025,Actual360()));
    Handle<YieldTermStructure> discountCurve(flatRate(today,0.03,Actual360()));

    shared_ptr<IborIndex> index(new USDLibor(6*Months, riskFreeRate));
    Natural fixingDays = 1;

    Real tolerance = 1.0e-6;

    shared_ptr<IborCouponPricer> pricer(new
        BlackIborCouponPricer(Handle<OptionletVolatilityStructure>()));

    Schedule sch(Date(30,November,2004),
                 Date(30,November,2008),
                 Period(Semiannual),
                 UnitedStates(UnitedStates::GovernmentBond),
                 ModifiedFollowing, ModifiedFollowing,
                 DateGeneration::Backward, false);

    std::shared_ptr<std::vector<std::pair<Date, Real> > > events(new std::vector<std::pair<Date, Real> >());
    events->emplace_back(std::pair<Date, Real>(Date(30,November,2004), 1000));
    std::shared_ptr<CatRisk> doomCatRisk(new EventSet(
        events, 
        Date(30,November,2004), Date(30,November,2008)));

    std::shared_ptr<EventPaymentOffset> paymentOffset(new NoOffset());
    std::shared_ptr<NotionalRisk> notionalRisk(new DigitalNotionalRisk(paymentOffset, 100));

    FloatingCatBond catBond(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           notionalRisk,
                           ModifiedFollowing, fixingDays,
                           std::vector<Rate>(), std::vector<Spread>(),
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

    shared_ptr<PricingEngine> catBondEngine(new MonteCarloCatBondEngine(doomCatRisk, discountCurve));
    catBond.setPricingEngine(catBondEngine);
    setCouponPricer(catBond.cashflows(),pricer);

    Real price = catBond.cleanPrice();
    CHECK(0 == price);

    Real lossProbability = catBond.lossProbability();
    Real exhaustionProbability = catBond.exhaustionProbability();
    Real expectedLoss = catBond.expectedLoss();

    CHECK(close(Real(1.0), lossProbability, tolerance));
    CHECK(close(Real(1.0), exhaustionProbability, tolerance));
    CHECK(close(Real(1.0), expectedLoss, tolerance));
}


TEST_CASE( "CatBond_CatBondWithDoomOnceInTenYears", "[CatBond]" ) {
    INFO("Testing floating-rate cat bond in a doom once in 10 years scenario...");

    CommonVars vars;

    Date today(22,November,2004);
    Settings::instance().evaluationDate() = today;

    Natural settlementDays = 1;

    Handle<YieldTermStructure> riskFreeRate(flatRate(today,0.025,Actual360()));
    Handle<YieldTermStructure> discountCurve(flatRate(today,0.03,Actual360()));

    shared_ptr<IborIndex> index(new USDLibor(6*Months, riskFreeRate));
    Natural fixingDays = 1;

    Real tolerance = 1.0e-6;

    shared_ptr<IborCouponPricer> pricer(new
        BlackIborCouponPricer(Handle<OptionletVolatilityStructure>()));

    Schedule sch(Date(30,November,2004),
                 Date(30,November,2008),
                 Period(Semiannual),
                 UnitedStates(UnitedStates::GovernmentBond),
                 ModifiedFollowing, ModifiedFollowing,
                 DateGeneration::Backward, false);

    std::shared_ptr<std::vector<std::pair<Date, Real> > > events(new std::vector<std::pair<Date, Real> >());
    events->emplace_back(std::pair<Date, Real>(Date(30,November,2008), 1000));
    std::shared_ptr<CatRisk> doomCatRisk(new EventSet(
        events, 
        Date(30,November,2004), Date(30,November,2044)));

    std::shared_ptr<CatRisk> noCatRisk(new EventSet(
        std::shared_ptr<std::vector<std::pair<Date, Real> > >(new std::vector<std::pair<Date, Real> >()),
        Date(1, Jan, 2000), Date(31, Dec, 2010)));

    std::shared_ptr<EventPaymentOffset> paymentOffset(new NoOffset());
    std::shared_ptr<NotionalRisk> notionalRisk(new DigitalNotionalRisk(paymentOffset, 100));

    FloatingCatBond catBond(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           notionalRisk,
                           ModifiedFollowing, fixingDays,
                           std::vector<Rate>(), std::vector<Spread>(),
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

    shared_ptr<PricingEngine> catBondEngine(new MonteCarloCatBondEngine(doomCatRisk, discountCurve));
    catBond.setPricingEngine(catBondEngine);
    setCouponPricer(catBond.cashflows(),pricer);

    Real price = catBond.cleanPrice();
    Real yield = catBond.yield(ActualActual(ActualActual::ISMA), Simple, Annual);
    Real lossProbability = catBond.lossProbability();
    Real exhaustionProbability = catBond.exhaustionProbability();
    Real expectedLoss = catBond.expectedLoss();

    CHECK(close(Real(0.1), lossProbability, tolerance));
    CHECK(close(Real(0.1), exhaustionProbability, tolerance));
    CHECK(close(Real(0.1), expectedLoss, tolerance));

    shared_ptr<PricingEngine> catBondEngineRF(new MonteCarloCatBondEngine(noCatRisk, discountCurve));
    catBond.setPricingEngine(catBondEngineRF);

    Real riskFreePrice = catBond.cleanPrice();
    Real riskFreeYield = catBond.yield(ActualActual(ActualActual::ISMA), Simple, Annual);
    Real riskFreeLossProbability = catBond.lossProbability();
    Real riskFreeExhaustionProbability = catBond.exhaustionProbability();
    Real riskFreeExpectedLoss = catBond.expectedLoss();
    
    CHECK(close(Real(0.0), riskFreeLossProbability, tolerance));
    CHECK(close(Real(0.0), riskFreeExhaustionProbability, tolerance));
    CHECK(std::abs(riskFreeExpectedLoss) < tolerance);
    
    CHECK(close(riskFreePrice*0.9, price, tolerance));
    CHECK(riskFreeYield < yield);
}

TEST_CASE( "CatBond_CatBondWithDoomOnceInTenYearsProportional", "[CatBond]" ) {
    INFO("Testing floating-rate cat bond in a doom once in 10 years scenario with proportional notional reduction...");

    CommonVars vars;

    Date today(22,November,2004);
    Settings::instance().evaluationDate() = today;

    Natural settlementDays = 1;

    Handle<YieldTermStructure> riskFreeRate(flatRate(today,0.025,Actual360()));
    Handle<YieldTermStructure> discountCurve(flatRate(today,0.03,Actual360()));

    shared_ptr<IborIndex> index(new USDLibor(6*Months, riskFreeRate));
    Natural fixingDays = 1;

    Real tolerance = 1.0e-6;

    shared_ptr<IborCouponPricer> pricer(new
        BlackIborCouponPricer(Handle<OptionletVolatilityStructure>()));

    Schedule sch(Date(30,November,2004),
                 Date(30,November,2008),
                 Period(Semiannual),
                 UnitedStates(UnitedStates::GovernmentBond),
                 ModifiedFollowing, ModifiedFollowing,
                 DateGeneration::Backward, false);

    std::shared_ptr<std::vector<std::pair<Date, Real> > > events(new std::vector<std::pair<Date, Real> >());
    events->emplace_back(std::pair<Date, Real>(Date(30,November,2008), 1000));
    std::shared_ptr<CatRisk> doomCatRisk(new EventSet(
        events, 
        Date(30,November,2004), Date(30,November,2044)));

    std::shared_ptr<CatRisk> noCatRisk(new EventSet(
        std::shared_ptr<std::vector<std::pair<Date, Real> > >(new std::vector<std::pair<Date, Real> >()),
        Date(1, Jan, 2000), Date(31, Dec, 2010)));

    std::shared_ptr<EventPaymentOffset> paymentOffset(new NoOffset());
    std::shared_ptr<NotionalRisk> notionalRisk(new ProportionalNotionalRisk(paymentOffset, 500, 1500));

    FloatingCatBond catBond(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           notionalRisk,
                           ModifiedFollowing, fixingDays,
                           std::vector<Rate>(), std::vector<Spread>(),
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

    shared_ptr<PricingEngine> catBondEngine(new MonteCarloCatBondEngine(doomCatRisk, discountCurve));
    catBond.setPricingEngine(catBondEngine);
    setCouponPricer(catBond.cashflows(),pricer);

    Real price = catBond.cleanPrice();
    Real yield = catBond.yield(ActualActual(ActualActual::ISMA), Simple, Annual);
    Real lossProbability = catBond.lossProbability();
    Real exhaustionProbability = catBond.exhaustionProbability();
    Real expectedLoss = catBond.expectedLoss();

    CHECK(close(Real(0.1), lossProbability, tolerance));
    CHECK(close(Real(0.0), exhaustionProbability, tolerance));
    CHECK(close(Real(0.05), expectedLoss, tolerance));

    shared_ptr<PricingEngine> catBondEngineRF(new MonteCarloCatBondEngine(noCatRisk, discountCurve));
    catBond.setPricingEngine(catBondEngineRF);

    Real riskFreePrice = catBond.cleanPrice();
    Real riskFreeYield = catBond.yield(ActualActual(ActualActual::ISMA), Simple, Annual);
    Real riskFreeLossProbability = catBond.lossProbability();
    Real riskFreeExpectedLoss = catBond.expectedLoss();
    
    CHECK(close(Real(0.0), riskFreeLossProbability, tolerance));
    CHECK(std::abs(riskFreeExpectedLoss) < tolerance);
    
    CHECK(close(riskFreePrice*0.95, price, tolerance));
    CHECK(riskFreeYield < yield);
}


TEST_CASE( "CatBond_CatBondWithGeneratedEventsProportional", "[CatBond]" ) {
    INFO("Testing floating-rate cat bond in a generated scenario with proportional notional reduction...");

    CommonVars vars;

    Date today(22,November,2004);
    Settings::instance().evaluationDate() = today;

    Natural settlementDays = 1;

    Handle<YieldTermStructure> riskFreeRate(flatRate(today,0.025,Actual360()));
    Handle<YieldTermStructure> discountCurve(flatRate(today,0.03,Actual360()));

    shared_ptr<IborIndex> index(new USDLibor(6*Months, riskFreeRate));
    Natural fixingDays = 1;

    Real tolerance = 1.0e-6;

    shared_ptr<IborCouponPricer> pricer(new
        BlackIborCouponPricer(Handle<OptionletVolatilityStructure>()));

    Schedule sch(Date(30,November,2004),
                 Date(30,November,2008),
                 Period(Semiannual),
                 UnitedStates(UnitedStates::GovernmentBond),
                 ModifiedFollowing, ModifiedFollowing,
                 DateGeneration::Backward, false);

	std::shared_ptr<CatRisk> betaCatRisk(new BetaRisk(5000, 50, 500, 500));

    std::shared_ptr<CatRisk> noCatRisk(new EventSet(
        std::shared_ptr<std::vector<std::pair<Date, Real> > >(new std::vector<std::pair<Date, Real> >()),
        Date(1, Jan, 2000), Date(31, Dec, 2010)));

    std::shared_ptr<EventPaymentOffset> paymentOffset(new NoOffset());
    std::shared_ptr<NotionalRisk> notionalRisk(new ProportionalNotionalRisk(paymentOffset, 500, 1500));

    FloatingCatBond catBond(settlementDays, vars.faceAmount, sch,
                           index, ActualActual(ActualActual::ISMA),
                           notionalRisk,
                           ModifiedFollowing, fixingDays,
                           std::vector<Rate>(), std::vector<Spread>(),
                           std::vector<Rate>(), std::vector<Rate>(),
                           false,
                           100.0, Date(30,November,2004));

    shared_ptr<PricingEngine> catBondEngine(new MonteCarloCatBondEngine(betaCatRisk, discountCurve));
    catBond.setPricingEngine(catBondEngine);
    setCouponPricer(catBond.cashflows(),pricer);

    Real price = catBond.cleanPrice();
    Real yield = catBond.yield(ActualActual(ActualActual::ISMA), Simple, Annual);
    Real lossProbability = catBond.lossProbability();
    Real exhaustionProbability = catBond.exhaustionProbability();
    Real expectedLoss = catBond.expectedLoss();

    CHECK(lossProbability<1.0);
    CHECK(lossProbability>0.0);
    CHECK(exhaustionProbability<1.0);
    CHECK(exhaustionProbability>0.0);
    CHECK(expectedLoss>0.0);

    shared_ptr<PricingEngine> catBondEngineRF(new MonteCarloCatBondEngine(noCatRisk, discountCurve));
    catBond.setPricingEngine(catBondEngineRF);

    Real riskFreePrice = catBond.cleanPrice();
    Real riskFreeYield = catBond.yield(ActualActual(ActualActual::ISMA), Simple, Annual);
    Real riskFreeLossProbability = catBond.lossProbability();
    Real riskFreeExpectedLoss = catBond.expectedLoss();
    
    CHECK(close(Real(0.0), riskFreeLossProbability, tolerance));
    CHECK(std::abs(riskFreeExpectedLoss) < tolerance);
    
    CHECK(riskFreePrice > price);
    CHECK(riskFreeYield < yield);
}

