/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Jose Aparicio

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

#include <ql/qldefines.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/instruments/makevanillaswap.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/pricingengines/swap/cvaswapengine.hpp>
#include <ql/termstructures/yield/piecewiseyieldcurve.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/termstructures/credit/interpolatedhazardratecurve.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <ql/time/daycounters/actual360.hpp>

#include <chrono>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace QuantLib;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

    Integer sessionId() { return 0; }

}
#endif

/*
  This example reproduces Table 2 on page 11 of 
  A Formula for Interest Rate Swaps Valuation under
  Counterparty Risk in presence of Netting Agreements
  
  Damiano Brigo and Massimo Masetti; May 4, 2005
 */

int main(int, char *[]) {

    try {

        std::chrono::time_point<std::chrono::steady_clock> startT = std::chrono::steady_clock::now();
        std::cout << std::endl;

        Calendar calendar = TARGET();
        Date todaysDate(10, March, 2004);
        // must be a business day
        todaysDate = calendar.adjust(todaysDate);

        Settings::instance().evaluationDate() = todaysDate;

        std::shared_ptr < IborIndex > yieldIndx = std::make_shared<Euribor3M>();
        Size tenorsSwapMkt[] = {5, 10, 15, 20, 25, 30};

        // rates ignoring counterparty risk:
        Rate ratesSwapmkt[] = {.03249, .04074, .04463, .04675, .04775, .04811};

        vector<std::shared_ptr<RateHelper> > swapHelpers;
        for (Size i = 0; i < sizeof(tenorsSwapMkt) / sizeof(Size); i++)
            swapHelpers.emplace_back(std::make_shared<SwapRateHelper>(
                    Handle<Quote>(std::make_shared<SimpleQuote>(ratesSwapmkt[i])),
                    tenorsSwapMkt[i] * Years,
                    TARGET(),
                    Quarterly,
                    ModifiedFollowing,
                    ActualActual(ActualActual::ISDA),
                    yieldIndx));

        std::shared_ptr < YieldTermStructure > swapTS =
                std::make_shared<PiecewiseYieldCurve<Discount, LogLinear>>(
                        2, TARGET(), swapHelpers, ActualActual(ActualActual::ISDA),
                        1.0e-12);
        swapTS->enableExtrapolation();

        std::shared_ptr < PricingEngine > riskFreeEngine(
                std::make_shared<DiscountingSwapEngine>(
                        Handle<YieldTermStructure>(swapTS)));

        std::vector<Handle<DefaultProbabilityTermStructure> >
                defaultIntensityTS;

        Size defaultTenors[] = {0, 12, 36, 60, 84, 120, 180, 240, 300,
                                360};// months
        // Three risk levels:
        Real intensitiesLow[] = {0.0036, 0.0036, 0.0065, 0.0099, 0.0111,
                                 0.0177, 0.0177, 0.0177, 0.0177, 0.0177,
                                 0.0177};
        Real intensitiesMedium[] = {0.0202, 0.0202, 0.0231, 0.0266, 0.0278,
                                    0.0349, 0.0349, 0.0349, 0.0349, 0.0349,
                                    0.0349};
        Real intensitiesHigh[] = {0.0534, 0.0534, 0.0564, 0.06, 0.0614, 0.0696,
                                  0.0696, 0.0696, 0.0696, 0.0696, 0.0696};
        // Recovery rates:
        Real ctptyRRLow = 0.4, ctptyRRMedium = 0.35, ctptyRRHigh = 0.3;

        std::vector<Date> defaultTSDates;
        std::vector<Real> intesitiesVLow, intesitiesVMedium, intesitiesVHigh;

        for (Size i = 0; i < sizeof(defaultTenors) / sizeof(Size); i++) {
            defaultTSDates.emplace_back(TARGET().advance(todaysDate,
                                                         Period(defaultTenors[i], Months)));
            intesitiesVLow.emplace_back(intensitiesLow[i]);
            intesitiesVMedium.emplace_back(intensitiesMedium[i]);
            intesitiesVHigh.emplace_back(intensitiesHigh[i]);
        }

        defaultIntensityTS.emplace_back(Handle<DefaultProbabilityTermStructure>(
                std::make_shared<InterpolatedHazardRateCurve<BackwardFlat>>(
                        defaultTSDates,
                        intesitiesVLow,
                        Actual360(),
                        TARGET())));
        defaultIntensityTS.emplace_back(Handle<DefaultProbabilityTermStructure>(
                std::make_shared<InterpolatedHazardRateCurve<BackwardFlat>>(
                        defaultTSDates,
                        intesitiesVMedium,
                        Actual360(),
                        TARGET())));
        defaultIntensityTS.emplace_back(Handle<DefaultProbabilityTermStructure>(
                std::make_shared<InterpolatedHazardRateCurve<BackwardFlat>>(
                        defaultTSDates,
                        intesitiesVHigh,
                        Actual360(),
                        TARGET())));

        Volatility blackVol = 0.15;
        std::shared_ptr < PricingEngine > ctptySwapCvaLow =
                std::make_shared<CounterpartyAdjSwapEngine>(
                        Handle<YieldTermStructure>(swapTS),
                        blackVol,
                        defaultIntensityTS[0],
                        ctptyRRLow
                );

        std::shared_ptr < PricingEngine > ctptySwapCvaMedium =
                std::make_shared<CounterpartyAdjSwapEngine>(
                        Handle<YieldTermStructure>(swapTS),
                        blackVol,
                        defaultIntensityTS[1],
                        ctptyRRMedium);
        std::shared_ptr < PricingEngine > ctptySwapCvaHigh =
                std::make_shared<CounterpartyAdjSwapEngine>(
                        Handle<YieldTermStructure>(swapTS),
                        blackVol,
                        defaultIntensityTS[2],
                        ctptyRRHigh);

        defaultIntensityTS[0]->enableExtrapolation();
        defaultIntensityTS[1]->enableExtrapolation();
        defaultIntensityTS[2]->enableExtrapolation();


        /// SWAP RISKY REPRICE----------------------------------------------

        // fixed leg
        Frequency fixedLegFrequency = Quarterly;
        BusinessDayConvention fixedLegConvention = ModifiedFollowing;
        DayCounter fixedLegDayCounter = ActualActual(ActualActual::ISDA);
        DayCounter floatingLegDayCounter = ActualActual(ActualActual::ISDA);

        VanillaSwap::Type swapType =
                //VanillaSwap::Receiver ;
                VanillaSwap::Payer;
        std::shared_ptr < IborIndex > yieldIndxS =
                std::make_shared<Euribor3M>(Handle<YieldTermStructure>(swapTS));
        std::vector<VanillaSwap> riskySwaps;
        for (Size i = 0; i < sizeof(tenorsSwapMkt) / sizeof(Size); i++)
            riskySwaps.emplace_back(MakeVanillaSwap(tenorsSwapMkt[i] * Years,
                                                    yieldIndxS,
                                                    ratesSwapmkt[i],
                                                    0 * Days)
                                            .withSettlementDays(2)
                                            .withFixedLegDayCount(fixedLegDayCounter)
                                            .withFixedLegTenor(Period(fixedLegFrequency))
                                            .withFixedLegConvention(fixedLegConvention)
                                            .withFixedLegTerminationDateConvention(fixedLegConvention)
                                            .withFixedLegCalendar(calendar)
                                            .withFloatingLegCalendar(calendar)
                                            .withNominal(100.)
                                            .withType(swapType));

        cout << "-- Correction in the contract fix rate in bp --" << endl;
        /* The paper plots correction to be substracted, here is printed
           with its sign 
        */
        for (Size i = 0; i < riskySwaps.size(); i++) {
            cout << fixed << setprecision(3);
            cout << setw(4);
            riskySwaps[i].setPricingEngine(riskFreeEngine);
            // should recover the input here:
            Real nonRiskyFair = riskySwaps[i].fairRate();
            cout << tenorsSwapMkt[i];
            cout << setw(5);

            cout << " | " << io::rate(nonRiskyFair);
            cout << fixed << setprecision(2);
            cout << setw(5);
            // Low Risk:
            riskySwaps[i].setPricingEngine(ctptySwapCvaLow);
            cout << " | " << setw(6)
                 << 10000. * (riskySwaps[i].fairRate() - nonRiskyFair);
            //cout << " | " << setw(6) << riskySwaps[i].NPV() ;

            // Medium Risk:
            riskySwaps[i].setPricingEngine(ctptySwapCvaMedium);
            cout << " | " << setw(6)
                 << 10000. * (riskySwaps[i].fairRate() - nonRiskyFair);
            //cout << " | " << setw(6) << riskySwaps[i].NPV() ;

            riskySwaps[i].setPricingEngine(ctptySwapCvaHigh);
            cout << " | " << setw(6)
                 << 10000. * (riskySwaps[i].fairRate() - nonRiskyFair);
            //cout << " | " << setw(6) << riskySwaps[i].NPV() ;

            cout << endl;
        }

        cout << endl;

        std::chrono::time_point<std::chrono::steady_clock> endT = std::chrono::steady_clock::now();
        Real seconds = static_cast<double>((endT - startT).count()) / 1.0e9;
        Integer hours = Integer(seconds / 3600);
        seconds -= hours * 3600;
        Integer minutes = Integer(seconds / 60);
        seconds -= minutes * 60;
        cout << "Run completed in ";
        if (hours > 0)
            cout << hours << " h ";
        if (hours > 0 || minutes > 0)
            cout << minutes << " m ";
        cout << fixed << setprecision(0)
             << seconds << " s" << endl;

        return 0;
    } catch (exception &e) {
        cerr << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "unknown error" << endl;
        return 1;
    }
}

