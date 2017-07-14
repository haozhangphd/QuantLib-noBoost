/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2005, 2008 Klaus Spanderen
 Copyright (C) 2007 StatPro Italia srl

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
#include <ql/time/calendars/target.hpp>
#include <ql/processes/batesprocess.hpp>
#include <ql/processes/merton76process.hpp>
#include <ql/instruments/europeanoption.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/pricingengines/vanilla/batesengine.hpp>
#include <ql/pricingengines/vanilla/jumpdiffusionengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/mceuropeanhestonengine.hpp>
#include <ql/pricingengines/vanilla/fdbatesvanillaengine.hpp>
#include <ql/models/equity/batesmodel.hpp>
#include <ql/models/equity/hestonmodelhelper.hpp>
#include <ql/time/period.hpp>
#include <ql/quotes/simplequote.hpp>

using namespace QuantLib;


namespace {

    Real getCalibrationError(
               std::vector<std::shared_ptr<CalibrationHelper> > & options) {
        Real sse = 0;
        for (Size i = 0; i < options.size(); ++i) {
            const Real diff = options[i]->calibrationError()*100.0;
            sse += diff*diff;
        }
        return sse;
    }

}


TEST_CASE("BatesModel_AnalyticVsBlack", "[BatesModel]") {

    INFO("Testing analytic Bates engine against Black formula...");

    SavedSettings backup;

    Date settlementDate = Date::todaysDate();
    Settings::instance().evaluationDate() = settlementDate;

    DayCounter dayCounter = ActualActual();
    Date exerciseDate = settlementDate + 6*Months;

    std::shared_ptr<StrikedTypePayoff> payoff =
                                     std::make_shared<PlainVanillaPayoff>(Option::Put, 30);
    std::shared_ptr<Exercise> exercise = std::make_shared<EuropeanExercise>(exerciseDate);

    Handle<YieldTermStructure> riskFreeTS(flatRate(0.1, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0.04, dayCounter));
    Handle<Quote> s0(std::make_shared<SimpleQuote>(32.0));

    Real yearFraction = dayCounter.yearFraction(settlementDate, exerciseDate);
    Real forwardPrice = s0->value()*std::exp((0.1-0.04)*yearFraction);
    Real expected = blackFormula(payoff->optionType(), payoff->strike(),
        forwardPrice, std::sqrt(0.05*yearFraction)) *
                                            std::exp(-0.1*yearFraction);
    const Real v0 = 0.05;
    const Real kappa = 5.0;
    const Real theta = 0.05;
    const Real sigma = 1.0e-4;
    const Real rho = 0.0;
    const Real lambda = 0.0001;
    const Real nu = 0.0; 
    const Real delta = 0.0001;

    VanillaOption option(payoff, exercise);

    std::shared_ptr<BatesProcess> process =
        std::make_shared<BatesProcess>(riskFreeTS, dividendTS, s0, v0, 
                         kappa, theta, sigma, rho, lambda, nu, delta);

    std::shared_ptr<PricingEngine> engine = std::make_shared<BatesEngine>(
        std::make_shared<BatesModel>(process), 64);

    option.setPricingEngine(engine);
    Real calculated = option.NPV();

    Real tolerance = 2.0e-7;
    Real error = std::fabs(calculated - expected);
    if (error > tolerance) {
        FAIL_CHECK("failed to reproduce Black price with BatesEngine"
                    << QL_FIXED
                    << "\n    calculated: " << calculated
                    << "\n    expected:   " << expected
                    << QL_SCIENTIFIC
                    << "\n    error:      " << error);
    }

    engine = std::make_shared<BatesDetJumpEngine>(
        std::make_shared<BatesDetJumpModel>(process, 1.0, 0.0001), 64);

    option.setPricingEngine(engine);
    calculated = option.NPV();

    error = std::fabs(calculated - expected);
    if (error > tolerance) {
        FAIL_CHECK("failed to reproduce Black price with " \
                    "BatesDetJumpEngine"
                    << QL_FIXED
                    << "\n    calculated: " << calculated
                    << "\n    expected:   " << expected
                    << QL_SCIENTIFIC
                    << "\n    error:      " << error);
    }

    engine = std::make_shared<BatesDoubleExpEngine>(
        std::make_shared<BatesDoubleExpModel>(process, 0.0001, 0.0001, 0.0001), 64);

    option.setPricingEngine(engine);
    calculated = option.NPV();

    error = std::fabs(calculated - expected);
    if (error > tolerance) {
        FAIL_CHECK("failed to reproduce Black price with BatesDoubleExpEngine"
                    << QL_FIXED
                    << "\n    calculated: " << calculated
                    << "\n    expected:   " << expected
                    << QL_SCIENTIFIC
                    << "\n    error:      " << error);
    }

    engine = std::make_shared<BatesDoubleExpDetJumpEngine>(
        std::make_shared<BatesDoubleExpDetJumpModel>(
                process, 0.0001, 0.0001, 0.0001, 0.5, 1.0, 0.0001), 64);

    option.setPricingEngine(engine);
    calculated = option.NPV();

    error = std::fabs(calculated - expected);
    if (error > tolerance) {
        FAIL_CHECK("failed to reproduce Black price with " \
                    "BatesDoubleExpDetJumpEngine"
                    << QL_FIXED
                    << "\n    calculated: " << calculated
                    << "\n    expected:   " << expected
                    << QL_SCIENTIFIC
                    << "\n    error:      " << error);
    }
}


TEST_CASE("BatesModel_AnalyticAndMcVsJumpDiffusion", "[BatesModel]") {

    INFO("Testing analytic Bates engine against Merton-76 engine...");

    SavedSettings backup;

    Date settlementDate = Date::todaysDate();
    Settings::instance().evaluationDate() = settlementDate;

    DayCounter dayCounter = ActualActual();

    std::shared_ptr<StrikedTypePayoff> payoff =
                                     std::make_shared<PlainVanillaPayoff>(Option::Put, 95);

    Handle<YieldTermStructure> riskFreeTS(flatRate(0.1, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0.04, dayCounter));
    Handle<Quote> s0(std::make_shared<SimpleQuote>(100));

    Real v0 = 0.0433;
    // FLOATING_POINT_EXCEPTION
    std::shared_ptr<SimpleQuote> vol = std::make_shared<SimpleQuote>(std::sqrt(v0));
    std::shared_ptr<BlackVolTermStructure> volTS =
        flatVol(settlementDate, vol, dayCounter);

    const Real kappa = 0.5;
    const Real theta = v0;
    const Real sigma = 1.0e-4;
    const Real rho = 0.0;

    std::shared_ptr<SimpleQuote> jumpIntensity = std::make_shared<SimpleQuote>(2);
    std::shared_ptr<SimpleQuote> meanLogJump = std::make_shared<SimpleQuote>(-0.2);
    std::shared_ptr<SimpleQuote> jumpVol = std::make_shared<SimpleQuote>(0.2);

    std::shared_ptr<BatesProcess> batesProcess =
        std::make_shared<BatesProcess>(riskFreeTS, dividendTS, s0, v0, kappa, theta, sigma, rho,
        jumpIntensity->value(), meanLogJump->value(), jumpVol->value());

    std::shared_ptr<Merton76Process> mertonProcess =
        std::make_shared<Merton76Process>(s0, dividendTS, riskFreeTS,
                            Handle<BlackVolTermStructure>(volTS),
                            Handle<Quote>(jumpIntensity),
                            Handle<Quote>(meanLogJump),
                            Handle<Quote>(jumpVol));

    std::shared_ptr<PricingEngine> batesEngine = std::make_shared<BatesEngine>(
        std::make_shared<BatesModel>(batesProcess), 160);

    const Real mcTol = 0.1;
    std::shared_ptr<PricingEngine> mcBatesEngine =
        MakeMCEuropeanHestonEngine<PseudoRandom>(batesProcess)
            .withStepsPerYear(2)
            .withAntitheticVariate()
            .withAbsoluteTolerance(mcTol)
            .withSeed(1234);

    std::shared_ptr<PricingEngine> mertonEngine =
        std::make_shared<JumpDiffusionEngine>(mertonProcess, 1e-10, 1000);

    for (Integer i=1; i<=5; i+=2) {
        Date exerciseDate = settlementDate + i*Years;
        std::shared_ptr<Exercise> exercise =
            std::make_shared<EuropeanExercise>(exerciseDate);

        VanillaOption batesOption(payoff, exercise);

        batesOption.setPricingEngine(batesEngine);
        Real calculated = batesOption.NPV();

        batesOption.setPricingEngine(mcBatesEngine);
        Real mcCalculated = batesOption.NPV();

        EuropeanOption mertonOption(payoff, exercise);
        mertonOption.setPricingEngine(mertonEngine);
        Real expected = mertonOption.NPV();

        Real tolerance = 2e-8;
        Real relError = std::fabs(calculated - expected)/expected;
        if (relError > tolerance) {
            FAIL("failed to reproduce Merton76 price with semi "
                       "analytic BatesEngine"
                       << QL_FIXED << std::setprecision(8)
                       << "\n    calculated: " << calculated
                       << "\n    expected:   " << expected
                       << "\n    rel. error: " << relError
                       << "\n    tolerance:  " << tolerance);
        }

        Real mcError = std::fabs(expected - mcCalculated);
        if (mcError > 3*mcTol) {
            FAIL("failed to reproduce Merton76 price with Monte-Carlo "
                       "BatesEngine"
                       << QL_FIXED << std::setprecision(8)
                       << "\n    calculated: " << mcCalculated
                       << "\n    expected:   " << expected
                       << "\n    error: "      << mcError
                       << "\n    tolerance:  " << mcTol);
        }
    }
}

namespace {
    struct HestonModelData {
        const char* const name;
        Real v0;
        Real kappa;
        Real theta;
        Real sigma;
        Real rho;
        Real r;
        Real q;
    };
    
    HestonModelData hestonModels[] = {
        // ADI finite difference schemes for option pricing in the 
        // Heston model with correlation, K.J. in t'Hout and S. Foulon,
        {"'t Hout case 1", 0.04, 1.5, 0.04, 0.3, -0.9, 0.025, 0.0},
        // Efficient numerical methods for pricing American options under 
        // stochastic volatility, Samuli Ikonen and Jari Toivanen,
        {"Ikonen-Toivanen", 0.0625, 5, 0.16, 0.9, 0.1, 0.1, 0.0},
        // Not-so-complex logarithms in the Heston model, 
        // Christian Kahl and Peter JГ¤ckel
        {"Kahl-Jaeckel", 0.16, 1.0, 0.16, 2.0, -0.8, 0.0, 0.0},
        // self defined test cases
        {"Equity case", 0.07, 2.0, 0.04, 0.55, -0.8, 0.03, 0.035 },
    };
}

TEST_CASE("BatesModel_AnalyticVsMCPricing", "[BatesModel]") {
    INFO("Testing analytic Bates engine against Monte-Carlo "
                       "engine...");

    SavedSettings backup;

    Date settlementDate(30, March, 2007);
    Settings::instance().evaluationDate() = settlementDate;

    DayCounter dayCounter = ActualActual();
    Date exerciseDate(30, March, 2012);

    std::shared_ptr<StrikedTypePayoff> payoff =
                                   std::make_shared<PlainVanillaPayoff>(Option::Put, 100);
    std::shared_ptr<Exercise> exercise = std::make_shared<EuropeanExercise>(exerciseDate);

    
    for (Size i=0; i < LENGTH(hestonModels); ++i) { 
        Handle<YieldTermStructure> riskFreeTS(flatRate(hestonModels[i].r, 
                                                       dayCounter));
        Handle<YieldTermStructure> dividendTS(flatRate(hestonModels[i].q, 
                                                       dayCounter));
        Handle<Quote> s0(std::make_shared<SimpleQuote>(100));

        std::shared_ptr<BatesProcess> batesProcess =
                       std::make_shared<BatesProcess>(riskFreeTS, dividendTS, s0,
                       hestonModels[i].v0, 
                       hestonModels[i].kappa, 
                       hestonModels[i].theta, 
                       hestonModels[i].sigma, 
                       hestonModels[i].rho, 2.0, -0.2, 0.1);
    
        const Real mcTolerance = 0.5;
        std::shared_ptr<PricingEngine> mcEngine =
                MakeMCEuropeanHestonEngine<PseudoRandom>(batesProcess)
                .withStepsPerYear(20)
                .withAntitheticVariate()
                .withAbsoluteTolerance(mcTolerance)
                .withSeed(1234);
    
        std::shared_ptr<BatesModel> batesModel = std::make_shared<BatesModel>(batesProcess);
        
        std::shared_ptr<PricingEngine> fdEngine =
                            std::make_shared<FdBatesVanillaEngine>(batesModel, 50, 100, 30);
    
        std::shared_ptr<PricingEngine> analyticEngine =
                                             std::make_shared<BatesEngine>(batesModel, 160);
    
        VanillaOption option(payoff, exercise);
    
        option.setPricingEngine(mcEngine);
        const Real calculated = option.NPV();
    
        option.setPricingEngine(analyticEngine);
        const Real expected = option.NPV();
    
        option.setPricingEngine(fdEngine);
        const Real fdCalculated = option.NPV();
        
        const Real mcError = std::fabs(calculated - expected);
        if (mcError > 3*mcTolerance) {
            FAIL("failed to reproduce Monte-Carlo price for BatesEngine"
                       << "\n    parameter:  " << hestonModels[i].name
                       << QL_FIXED << std::setprecision(8)
                       << "\n    calculated: " << calculated
                       << "\n    expected:   " << expected
                       << "\n    error: "      << mcError
                       << "\n    tolerance:  " << mcTolerance);
        }
        const Real fdTolerance = 0.2;
        const Real fdError = std::fabs(fdCalculated - expected);
        if (fdError > fdTolerance) {
            FAIL("failed to reproduce PIDE price for BatesEngine"
                       << "\n    parameter:  " << hestonModels[i].name
                       << QL_FIXED << std::setprecision(8)
                       << "\n    calculated: " << fdCalculated
                       << "\n    expected:   " << expected
                       << "\n    error: "      << fdError
                       << "\n    tolerance:  " << fdTolerance);
        }
    }
}

TEST_CASE("BatesModel_DAXCalibration", "[BatesModel]") {
    /* this example is taken from A. Sepp
       Pricing European-Style Options under Jump Diffusion Processes
       with Stochstic Volatility: Applications of Fourier Transform
       http://math.ut.ee/~spartak/papers/stochjumpvols.pdf
    */

    INFO(
             "Testing Bates model calibration using DAX volatility data...");

    SavedSettings backup;

    Date settlementDate(5, July, 2002);
    Settings::instance().evaluationDate() = settlementDate;

    DayCounter dayCounter = Actual365Fixed();
    Calendar calendar = TARGET();

    Integer t[] = { 13, 41, 75, 165, 256, 345, 524, 703 };
    Rate r[] = { 0.0357,0.0349,0.0341,0.0355,0.0359,0.0368,0.0386,0.0401 };

    std::vector<Date> dates;
    std::vector<Rate> rates;
    dates.emplace_back(settlementDate);
    rates.emplace_back(0.0357);
    for (Size i = 0; i < 8; ++i) {
        dates.emplace_back(settlementDate + t[i]);
        rates.emplace_back(r[i]);
    }
     // FLOATING_POINT_EXCEPTION
    Handle<YieldTermStructure> riskFreeTS(
                       std::make_shared<ZeroCurve>(dates, rates, dayCounter));

    Handle<YieldTermStructure> dividendTS(
                                   flatRate(settlementDate, 0.0, dayCounter));

    Volatility v[] =
      { 0.6625,0.4875,0.4204,0.3667,0.3431,0.3267,0.3121,0.3121,
        0.6007,0.4543,0.3967,0.3511,0.3279,0.3154,0.2984,0.2921,
        0.5084,0.4221,0.3718,0.3327,0.3155,0.3027,0.2919,0.2889,
        0.4541,0.3869,0.3492,0.3149,0.2963,0.2926,0.2819,0.2800,
        0.4060,0.3607,0.3330,0.2999,0.2887,0.2811,0.2751,0.2775,
        0.3726,0.3396,0.3108,0.2781,0.2788,0.2722,0.2661,0.2686,
        0.3550,0.3277,0.3012,0.2781,0.2781,0.2661,0.2661,0.2681,
        0.3428,0.3209,0.2958,0.2740,0.2688,0.2627,0.2580,0.2620,
        0.3302,0.3062,0.2799,0.2631,0.2573,0.2533,0.2504,0.2544,
        0.3343,0.2959,0.2705,0.2540,0.2504,0.2464,0.2448,0.2462,
        0.3460,0.2845,0.2624,0.2463,0.2425,0.2385,0.2373,0.2422,
        0.3857,0.2860,0.2578,0.2399,0.2357,0.2327,0.2312,0.2351,
        0.3976,0.2860,0.2607,0.2356,0.2297,0.2268,0.2241,0.2320 };

    Handle<Quote> s0(std::make_shared<SimpleQuote>(4468.17));
    Real strike[] = { 3400,3600,3800,4000,4200,4400,
                      4500,4600,4800,5000,5200,5400,5600 };


    Real v0 = 0.0433;
    std::shared_ptr<SimpleQuote> vol = std::make_shared<SimpleQuote>(std::sqrt(v0));
    std::shared_ptr<BlackVolTermStructure> volTS =
        flatVol(settlementDate, vol, dayCounter);

    const Real kappa = 1.0;
    const Real theta = v0;
    const Real sigma = 1.0;
    const Real rho = 0.0;
    const Real lambda = 1.1098;
    const Real nu = -0.1285;
    const Real delta = 0.1702;

    std::shared_ptr<BatesProcess> process =
        std::make_shared<BatesProcess>(riskFreeTS, dividendTS, s0, v0, 
                         kappa, theta, sigma, rho, lambda, nu, delta);

    std::shared_ptr<BatesModel> batesModel = std::make_shared<BatesModel>(process);

    std::shared_ptr<PricingEngine> batesEngine =
                                            std::make_shared<BatesEngine>(batesModel, 64);

    std::vector<std::shared_ptr<CalibrationHelper> > options;

    for (Size s = 0; s < 13; ++s) {
        for (Size m = 0; m < 8; ++m) {
            Handle<Quote> vol(std::make_shared<SimpleQuote>(v[s*8+m]));

            Period maturity((int)((t[m]+3)/7.), Weeks); // round to weeks

            // this is the calibration helper for the bates models
            // FLOATING_POINT_EXCEPTION
            options.emplace_back(std::make_shared<HestonModelHelper>(maturity, calendar,
                                          s0->value(), strike[s], vol,
                                          riskFreeTS, dividendTS, 
                                          CalibrationHelper::ImpliedVolError));
            options.back()->setPricingEngine(batesEngine);
        }
    }

    // check calibration engine
    LevenbergMarquardt om;
    batesModel->calibrate(options, om,
                          EndCriteria(400, 40, 1.0e-8, 1.0e-8, 1.0e-8));

    Real expected = 36.6;
    Real calculated = getCalibrationError(options);

    if (std::fabs(calculated - expected) > 2.5)
        FAIL_CHECK("failed to calibrate the bates model"
                    << "\n    calculated: " << calculated
                    << "\n    expected:   " << expected);

    //check pricing of derived engines
    std::vector<std::shared_ptr<PricingEngine> > pricingEngines;
    
    process = std::make_shared<BatesProcess>(riskFreeTS, dividendTS, s0, v0, 
                         kappa, theta, sigma, rho, 1.0, -0.1, 0.1);

    pricingEngines.emplace_back(std::make_shared<BatesDetJumpEngine>(
            std::make_shared<BatesDetJumpModel>(process), 64)) ;

    std::shared_ptr<HestonProcess> hestonProcess =
        std::make_shared<HestonProcess>(riskFreeTS, dividendTS, s0, v0, kappa, theta, sigma, rho);

    pricingEngines.emplace_back(std::make_shared<BatesDoubleExpEngine>(
            std::make_shared<BatesDoubleExpModel>(
                         hestonProcess, 1.0), 64));

    pricingEngines.emplace_back(std::make_shared<BatesDoubleExpDetJumpEngine>(
            std::make_shared<BatesDoubleExpDetJumpModel>(hestonProcess, 1.0), 64));

    Real expectedValues[] = { 5896.37,
                              5499.29,
                              6497.89};

    Real tolerance=0.1;
    for (Size i = 0; i < pricingEngines.size(); ++i) {
        for (Size j = 0; j < options.size(); ++j) {
            options[j]->setPricingEngine(pricingEngines[i]);
        }

        Real calculated = std::fabs(getCalibrationError(options));
        if (std::fabs(calculated - expectedValues[i]) > tolerance)
            FAIL_CHECK("failed to calculated prices for derived Bates models"
                        << "\n    calculated: " << calculated
                        << "\n    expected:   " << expectedValues[i]);
    }
}
