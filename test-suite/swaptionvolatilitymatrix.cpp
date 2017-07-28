/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2008 Ferdinando Ametrano
 Copyright (C) 2006 François du Vignaud
 Copyright (C) 2007 Cristina Duminuco

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

#include "swaptionvolstructuresutilities.hpp"
#include "utilities.hpp"
#include <ql/utilities/dataformatters.hpp>
#include <ql/indexes/swap/euriborswap.hpp>
#include <ql/instruments/makeswaption.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/math/comparison.hpp>
#include <string>

using namespace QuantLib;


namespace {

    struct CommonVars {
        // global data
        Date referenceDate_;
        SwaptionMarketConventions conventions_;
        AtmVolatility atm_;
        RelinkableHandle<YieldTermStructure> termStructure_;
        RelinkableHandle<SwaptionVolatilityStructure> atmVolMatrix_;
        Real tolerance_;

        // cleanup
        SavedSettings backup_;

        // setup
        CommonVars() {
            conventions_.setConventions();
            atm_.setMarketData();
            Settings::instance().evaluationDate() =
                conventions_.calendar.adjust(Date::todaysDate());
            atmVolMatrix_ = RelinkableHandle<SwaptionVolatilityStructure>(
                std::make_shared<SwaptionVolatilityMatrix>(conventions_.calendar,
                                             conventions_.optionBdc,
                                             atm_.tenors.options,
                                             atm_.tenors.swaps,
                                             atm_.volsHandle,
                                             conventions_.dayCounter));
            termStructure_.linkTo(
                std::make_shared<FlatForward>(0, conventions_.calendar,
                                0.05, Actual365Fixed()));
        }

        // utilities
        void makeObservabilityTest(
                const std::string& description,
                const std::shared_ptr<SwaptionVolatilityStructure>& vol,
                bool mktDataFloating,
                bool referenceDateFloating) {
            Rate dummyStrike = .02;
            Date referenceDate = Settings::instance().evaluationDate();
            Volatility initialVol = vol->volatility(
                        referenceDate + atm_.tenors.options[0],
                        atm_.tenors.swaps[0], dummyStrike, false);
            // testing evaluation date change ...
            Settings::instance().evaluationDate() =
                referenceDate - Period(1, Years);
            Volatility newVol =  vol->volatility(
                        referenceDate + atm_.tenors.options[0],
                        atm_.tenors.swaps[0], dummyStrike, false);
            Settings::instance().evaluationDate() = referenceDate;
            if (referenceDateFloating && (initialVol == newVol))
                FAIL_CHECK(description <<
                            " the volatility should change when the reference date is changed !");
            if (!referenceDateFloating && (initialVol != newVol))
                FAIL_CHECK(description <<
                            " the volatility should not change when the reference date is changed !");

            // test market data change...
            if (mktDataFloating){
                Volatility initialVolatility = atm_.volsHandle[0][0]->value();
                std::dynamic_pointer_cast<SimpleQuote>(
                              atm_.volsHandle[0][0].currentLink())->setValue(10);
                newVol = vol->volatility(
                    referenceDate + atm_.tenors.options[0],
                    atm_.tenors.swaps[0], dummyStrike, false);
                std::dynamic_pointer_cast<SimpleQuote>(
                    atm_.volsHandle[0][0].currentLink())
                    ->setValue(initialVolatility);
                if (initialVol == newVol)
                    FAIL_CHECK(description << " the volatility should change when"
                                " the market data is changed !");
            }
        }

        void makeCoherenceTest(
                const std::string& description,
                const std::shared_ptr<SwaptionVolatilityDiscrete>& vol) {

            for (Size i=0; i<atm_.tenors.options.size(); ++i) {
                Date optionDate =
                    vol->optionDateFromTenor(atm_.tenors.options[i]);
                if (optionDate!=vol->optionDates()[i])
                    FAIL(
                         "optionDateFromTenor failure for " <<
                         description << ":"
                         "\n       option tenor: " << atm_.tenors.options[i] <<
                         "\nactual option date : " << optionDate <<
                         "\n  exp. option date : " << vol->optionDates()[i]);
                Time optionTime = vol->timeFromReference(optionDate);
                if (!close(optionTime,vol->optionTimes()[i]))
                    FAIL(
                         "timeFromReference failure for " <<
                         description << ":"
                         "\n       option tenor: " << atm_.tenors.options[i] <<
                         "\n       option date : " << optionDate <<
                         "\nactual option time : " << optionTime <<
                         "\n  exp. option time : " << vol->optionTimes()[i]);
            }

            std::shared_ptr<BlackSwaptionEngine> engine =
                std::make_shared<BlackSwaptionEngine>(termStructure_,
                                    Handle<SwaptionVolatilityStructure>(vol));

            for (Size j=0; j<atm_.tenors.swaps.size(); j++) {
                Time swapLength = vol->swapLength(atm_.tenors.swaps[j]);
                if (!close(swapLength,years(atm_.tenors.swaps[j])))
                    FAIL("convertSwapTenor failure for " <<
                               description << ":"
                               "\n        swap tenor : " << atm_.tenors.swaps[j] <<
                               "\n actual swap length: " << swapLength <<
                               "\n   exp. swap length: " << years(atm_.tenors.swaps[j]));

                std::shared_ptr<SwapIndex> swapIndex =
                    std::make_shared<EuriborSwapIsdaFixA>(atm_.tenors.swaps[j], termStructure_);

                for (Size i=0; i<atm_.tenors.options.size(); ++i) {
                    Real error, tolerance = 1.0e-16;
                    Volatility actVol, expVol = atm_.vols[i][j];

                    actVol = vol->volatility(atm_.tenors.options[i],
                                             atm_.tenors.swaps[j], 0.05, true);
                    error = std::abs(expVol-actVol);
                    if (error>tolerance)
                        FAIL(
                              "recovery of atm vols failed for " <<
                              description << ":"
                              "\noption tenor = " << atm_.tenors.options[i] <<
                              "\n swap length = " << atm_.tenors.swaps[j] <<
                              "\nexpected vol = " << io::volatility(expVol) <<
                              "\n  actual vol = " << io::volatility(actVol) <<
                              "\n       error = " << io::volatility(error) <<
                              "\n   tolerance = " << tolerance);

                    Date optionDate =
                        vol->optionDateFromTenor(atm_.tenors.options[i]);
                    actVol = vol->volatility(optionDate,
                                             atm_.tenors.swaps[j], 0.05, true);
                    error = std::abs(expVol-actVol);
                    if (error>tolerance)
                        FAIL(
                             "recovery of atm vols failed for " <<
                             description << ":"
                             "\noption tenor: " << atm_.tenors.options[i] <<
                             "\noption date : " << optionDate <<
                             "\n  swap tenor: " << atm_.tenors.swaps[j] <<
                             "\n   exp. vol: " << io::volatility(expVol) <<
                             "\n actual vol: " << io::volatility(actVol) <<
                             "\n      error: " << io::volatility(error) <<
                             "\n  tolerance: " << tolerance);

                    Time optionTime = vol->timeFromReference(optionDate);
                    actVol = vol->volatility(optionTime, swapLength,
                                             0.05, true);
                    error = std::abs(expVol-actVol);
                    if (error>tolerance)
                        FAIL(
                             "recovery of atm vols failed for " <<
                             description << ":"
                             "\noption tenor: " << atm_.tenors.options[i] <<
                             "\noption time : " << optionTime <<
                             "\n  swap tenor: " << atm_.tenors.swaps[j] <<
                             "\n swap length: " << swapLength <<
                             "\n    exp. vol: " << io::volatility(expVol) <<
                             "\n  actual vol: " << io::volatility(actVol) <<
                             "\n       error: " << io::volatility(error) <<
                             "\n   tolerance: " << tolerance);

                    // ATM swaption
                    Swaption swaption =
                        MakeSwaption(swapIndex, atm_.tenors.options[i])
                        .withPricingEngine(engine);

                    Date exerciseDate = swaption.exercise()->dates().front();
                    if (exerciseDate!=vol->optionDates()[i])
                        FAIL(
                             "\noptionDateFromTenor mismatch for " <<
                             description << ":"
                             "\n      option tenor: " << atm_.tenors.options[i] <<
                             "\nactual option date: " << exerciseDate <<
                             "\n  exp. option date: " << vol->optionDates()[i]);

                    Date start = swaption.underlyingSwap()->startDate();
                    Date end = swaption.underlyingSwap()->maturityDate();
                    Time swapLength2 = vol->swapLength(start, end);
                    if (!close(swapLength2,swapLength))
                        FAIL("\nswapLength failure for " <<
                                   description << ":"
                                   "\n   exp. swap length: " << swapLength <<
                                   "\n actual swap length: " << swapLength2 <<
                                   "\n        swap tenor : " << atm_.tenors.swaps[j] <<
                                   "\n  swap index tenor : " << swapIndex->tenor() <<
                                   "\n        option date: " << exerciseDate <<
                                   "\n         start date: " << start <<
                                   "\n      maturity date: " << end
                                   );

                    Real npv = swaption.NPV();
                    actVol = swaption.impliedVolatility(npv, termStructure_,
                                                        expVol*0.98, 1e-6,
                                                        100, 10.0e-7, 4.0,
                                                        ShiftedLognormal, 0.0);
                    error = std::abs(expVol-actVol);
                    Real tolerance2 = 0.000001;
                    if (error>tolerance2)
                        FAIL(
                             "recovery of atm vols through BlackSwaptionEngine failed for " <<
                             description << ":"
                             "\noption tenor: " << atm_.tenors.options[i] <<
                             "\noption time : " << optionTime <<
                             "\n  swap tenor: " << atm_.tenors.swaps[j] <<
                             "\n swap length: " << swapLength <<
                             "\n   exp. vol: " << io::volatility(expVol) <<
                             "\n actual vol: " << io::volatility(actVol) <<
                             "\n      error: " << io::volatility(error) <<
                             "\n  tolerance: " << tolerance2);
                }
            }
        }
    };

}


TEST_CASE("SwaptionVolatilityMatrix_SwaptionVolMatrixObservability", "[SwaptionVolatilityMatrix]") {

    INFO("Testing swaption volatility matrix observability...");

    CommonVars vars;

    std::shared_ptr<SwaptionVolatilityMatrix> vol;
    std::string description;

    //floating reference date, floating market data
    description = "floating reference date, floating market data";
    vol = std::make_shared<SwaptionVolatilityMatrix>(vars.conventions_.calendar,
                                 vars.conventions_.optionBdc,
                                 vars.atm_.tenors.options,
                                 vars.atm_.tenors.swaps,
                                 vars.atm_.volsHandle,
                                 vars.conventions_.dayCounter);
    vars.makeObservabilityTest(description, vol, true, true);

    //fixed reference date, floating market data
    description = "fixed reference date, floating market data";
    vol = std::make_shared<SwaptionVolatilityMatrix>(Settings::instance().evaluationDate(),
                                 vars.conventions_.calendar,
                                 vars.conventions_.optionBdc,
                                 vars.atm_.tenors.options,
                                 vars.atm_.tenors.swaps,
                                 vars.atm_.volsHandle,
                                 vars.conventions_.dayCounter);
    vars.makeObservabilityTest(description, vol, true, false);

    // floating reference date, fixed market data
    description = "floating reference date, fixed market data";
    vol = std::make_shared<SwaptionVolatilityMatrix>(vars.conventions_.calendar,
                                 vars.conventions_.optionBdc,
                                 vars.atm_.tenors.options,
                                 vars.atm_.tenors.swaps,
                                 vars.atm_.volsHandle,
                                 vars.conventions_.dayCounter);
    vars.makeObservabilityTest(description, vol, false, true);

    // fixed reference date, fixed market data
    description = "fixed reference date, fixed market data";
    vol = std::make_shared<SwaptionVolatilityMatrix>(Settings::instance().evaluationDate(),
                                 vars.conventions_.calendar,
                                 vars.conventions_.optionBdc,
                                 vars.atm_.tenors.options,
                                 vars.atm_.tenors.swaps,
                                 vars.atm_.volsHandle,
                                 vars.conventions_.dayCounter);
    vars.makeObservabilityTest(description, vol, false, false);

   // fixed reference date and fixed market data, option dates
        //SwaptionVolatilityMatrix(const Date& referenceDate,
        //                         const std::vector<Date>& exerciseDates,
        //                         const std::vector<Period>& swapTenors,
        //                         const Matrix& volatilities,
        //                         const DayCounter& dayCounter);
}


TEST_CASE("SwaptionVolatilityMatrix_SwaptionVolMatrixCoherence", "[SwaptionVolatilityMatrix]") {

    INFO("Testing swaption volatility matrix...");

    CommonVars vars;

    std::shared_ptr<SwaptionVolatilityMatrix> vol;
    std::string description;

    //floating reference date, floating market data
    description = "floating reference date, floating market data";
    vol = std::make_shared<SwaptionVolatilityMatrix>(vars.conventions_.calendar,
                                 vars.conventions_.optionBdc,
                                 vars.atm_.tenors.options,
                                 vars.atm_.tenors.swaps,
                                 vars.atm_.volsHandle,
                                 vars.conventions_.dayCounter);
    vars.makeCoherenceTest(description, vol);

    //fixed reference date, floating market data
    description = "fixed reference date, floating market data";
    vol = std::make_shared<SwaptionVolatilityMatrix>(Settings::instance().evaluationDate(),
                                 vars.conventions_.calendar,
                                 vars.conventions_.optionBdc,
                                 vars.atm_.tenors.options,
                                 vars.atm_.tenors.swaps,
                                 vars.atm_.volsHandle,
                                 vars.conventions_.dayCounter);
    vars.makeCoherenceTest(description, vol);

    // floating reference date, fixed market data
    description = "floating reference date, fixed market data";
    vol = std::make_shared<SwaptionVolatilityMatrix>(vars.conventions_.calendar,
                                 vars.conventions_.optionBdc,
                                 vars.atm_.tenors.options,
                                 vars.atm_.tenors.swaps,
                                 vars.atm_.volsHandle,
                                 vars.conventions_.dayCounter);
    vars.makeCoherenceTest(description, vol);

    // fixed reference date, fixed market data
    description = "fixed reference date, fixed market data";
    vol = std::make_shared<SwaptionVolatilityMatrix>(Settings::instance().evaluationDate(),
                                 vars.conventions_.calendar,
                                 vars.conventions_.optionBdc,
                                 vars.atm_.tenors.options,
                                 vars.atm_.tenors.swaps,
                                 vars.atm_.volsHandle,
                                 vars.conventions_.dayCounter);
    vars.makeCoherenceTest(description, vol);
}
