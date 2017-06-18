/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007 Giorgio Facchinetti

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
#include <ql/indexes/swap/euriborswap.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/cashflows/rangeaccrual.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/termstructures/volatility/interpolatedsmilesection.hpp>
#include <ql/termstructures/volatility/flatsmilesection.hpp>

#include <ql/cashflows/conundrumpricer.hpp>
#include <ql/cashflows/cashflowvectors.hpp>
#include <ql/cashflows/cashflows.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolmatrix.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolcube2.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolcube1.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/time/schedule.hpp>

using namespace QuantLib;


namespace {

    struct CommonVars {
        // General settings
        Date referenceDate, today, settlement;
        Calendar calendar;

        // Volatility Stuctures
        std::vector<Handle<SwaptionVolatilityStructure> > swaptionVolatilityStructures;
        Handle<SwaptionVolatilityStructure> atmVol;
        Handle<SwaptionVolatilityStructure> flatSwaptionVolatilityCube1;
        Handle<SwaptionVolatilityStructure> flatSwaptionVolatilityCube2;
        Handle<SwaptionVolatilityStructure> swaptionVolatilityCubeBySabr;

        std::vector<Period> atmOptionTenors, optionTenors;
        std::vector<Period> atmSwapTenors, swapTenors;
        std::vector<Spread> strikeSpreads;

        Matrix atmVolMatrix, volSpreadsMatrix;
        std::vector<std::vector<Handle<Quote> > > volSpreads;

        DayCounter dayCounter;
        BusinessDayConvention optionBDC;
        Natural swapSettlementDays;
        bool vegaWeightedSmileFit;

        // Range Accrual valuation
        Rate infiniteLowerStrike, infiniteUpperStrike;
        Real gearing, correlation;
        Spread spread;
        Date startDate;
        Date endDate;
        Date paymentDate;
        Natural fixingDays;
        DayCounter rangeCouponDayCount;
        std::shared_ptr<Schedule> observationSchedule;
        // Observation Schedule conventions
        Frequency observationsFrequency;
        BusinessDayConvention observationsConvention;

        // Term Structure
        RelinkableHandle<YieldTermStructure> termStructure;

        // indices and index conventions
        Frequency fixedLegFrequency;
        BusinessDayConvention fixedLegConvention;
        DayCounter fixedLegDayCounter;
        std::shared_ptr<IborIndex> iborIndex;
        std::shared_ptr<SwapIndex> swapIndexBase;
        std::shared_ptr<SwapIndex> shortSwapIndexBase;

        // Range accrual pricers properties
        std::vector<bool> byCallSpread;
        Real flatVol;
        std::vector<std::shared_ptr<SmileSection> > smilesOnExpiry;
        std::vector<std::shared_ptr<SmileSection> > smilesOnPayment;

        //test parameters
        Real rateTolerance;
        Real priceTolerance;


        // cleanup
        SavedSettings backup;

        void createYieldCurve() {

            // Yield Curve
            std::vector<Date> dates;
            dates.emplace_back(Date(39147));    dates.emplace_back(Date(39148));    dates.emplace_back(Date(39151));
            dates.emplace_back(Date(39153));    dates.emplace_back(Date(39159));    dates.emplace_back(Date(39166));
            dates.emplace_back(Date(39183));    dates.emplace_back(Date(39294));    dates.emplace_back(Date(39384));
            dates.emplace_back(Date(39474));    dates.emplace_back(Date(39567));    dates.emplace_back(Date(39658));
            dates.emplace_back(Date(39748));    dates.emplace_back(Date(39839));    dates.emplace_back(Date(39931));
            dates.emplace_back(Date(40250));    dates.emplace_back(Date(40614));    dates.emplace_back(Date(40978));
            dates.emplace_back(Date(41344));    dates.emplace_back(Date(41709));    dates.emplace_back(Date(42074));
            dates.emplace_back(Date(42441));    dates.emplace_back(Date(42805));    dates.emplace_back(Date(43170));
            dates.emplace_back(Date(43535));    dates.emplace_back(Date(43900));    dates.emplace_back(Date(44268));
            dates.emplace_back(Date(44632));    dates.emplace_back(Date(44996));    dates.emplace_back(Date(45361));
            dates.emplace_back(Date(45727));    dates.emplace_back(Date(46092));    dates.emplace_back(Date(46459));
            dates.emplace_back(Date(46823));    dates.emplace_back(Date(47188));    dates.emplace_back(Date(47553));
            dates.emplace_back(Date(47918));    dates.emplace_back(Date(48283));    dates.emplace_back(Date(48650));
            dates.emplace_back(Date(49014));    dates.emplace_back(Date(49379));    dates.emplace_back(Date(49744));
            dates.emplace_back(Date(50110));    dates.emplace_back(Date(53762));    dates.emplace_back(Date(57415));
            dates.emplace_back(Date(61068));

            std::vector<Rate> zeroRates;
            zeroRates.emplace_back(0.02676568527);    zeroRates.emplace_back(0.02676568527);
            zeroRates.emplace_back(0.02676333038);    zeroRates.emplace_back(0.02682286201);
            zeroRates.emplace_back(0.02682038347);    zeroRates.emplace_back(0.02683030208);
            zeroRates.emplace_back(0.02700136766);    zeroRates.emplace_back(0.02932526033);
            zeroRates.emplace_back(0.03085568949);    zeroRates.emplace_back(0.03216370631);
            zeroRates.emplace_back(0.03321234116);    zeroRates.emplace_back(0.03404978072);
            zeroRates.emplace_back(0.03471117149);    zeroRates.emplace_back(0.03527141916);
            zeroRates.emplace_back(0.03574660393);    zeroRates.emplace_back(0.03691715582);
            zeroRates.emplace_back(0.03796468718);    zeroRates.emplace_back(0.03876457629);
            zeroRates.emplace_back(0.03942029708);    zeroRates.emplace_back(0.03999925325);
            zeroRates.emplace_back(0.04056663618);    zeroRates.emplace_back(0.04108743922);
            zeroRates.emplace_back(0.04156156761);    zeroRates.emplace_back(0.0419979179);
            zeroRates.emplace_back(0.04239486483);    zeroRates.emplace_back(0.04273799032);
            zeroRates.emplace_back(0.04305531203);    zeroRates.emplace_back(0.04336417578);
            zeroRates.emplace_back(0.04364017665);    zeroRates.emplace_back(0.04388153459);
            zeroRates.emplace_back(0.04408005012);    zeroRates.emplace_back(0.04424764425);
            zeroRates.emplace_back(0.04437504759);    zeroRates.emplace_back(0.04447696334);
            zeroRates.emplace_back(0.04456212318);    zeroRates.emplace_back(0.04464090072);
            zeroRates.emplace_back(0.0447068707);     zeroRates.emplace_back(0.04475921774);
            zeroRates.emplace_back(0.04477418345);    zeroRates.emplace_back(0.04477880755);
            zeroRates.emplace_back(0.04476692489);    zeroRates.emplace_back(0.04473779454);
            zeroRates.emplace_back(0.04468646066);    zeroRates.emplace_back(0.04430951558);
            zeroRates.emplace_back(0.04363922313);    zeroRates.emplace_back(0.04363601992);

            termStructure.linkTo( std::shared_ptr<YieldTermStructure>(
                new ZeroCurve(dates, zeroRates, Actual365Fixed())));
        }

        void createVolatilityStructures() {

            // ATM swaptionvol matrix
            optionBDC = Following;

            atmOptionTenors = std::vector<Period>();
            atmOptionTenors.emplace_back(Period(1, Months));
            atmOptionTenors.emplace_back(Period(6, Months));
            atmOptionTenors.emplace_back(Period(1, Years));
            atmOptionTenors.emplace_back(Period(5, Years));
            atmOptionTenors.emplace_back(Period(10, Years));
            atmOptionTenors.emplace_back(Period(30, Years));

            atmSwapTenors = std::vector<Period>();
            atmSwapTenors.emplace_back(Period(1, Years));
            atmSwapTenors.emplace_back(Period(5, Years));
            atmSwapTenors.emplace_back(Period(10, Years));
            atmSwapTenors.emplace_back(Period(30, Years));

            atmVolMatrix = Matrix(atmOptionTenors.size(), atmSwapTenors.size());
            //atmVolMatrix[0][0]=0.1300; atmVolMatrix[0][1]=0.1560; atmVolMatrix[0][2]=0.1390; atmVolMatrix[0][3]=0.1220;
            //atmVolMatrix[1][0]=0.1440; atmVolMatrix[1][1]=0.1580; atmVolMatrix[1][2]=0.1460; atmVolMatrix[1][3]=0.1260;
            //atmVolMatrix[2][0]=0.1600; atmVolMatrix[2][1]=0.1590; atmVolMatrix[2][2]=0.1470; atmVolMatrix[2][3]=0.1290;
            //atmVolMatrix[3][0]=0.1640; atmVolMatrix[3][1]=0.1470; atmVolMatrix[3][2]=0.1370; atmVolMatrix[3][3]=0.1220;
            //atmVolMatrix[4][0]=0.1400; atmVolMatrix[4][1]=0.1300; atmVolMatrix[4][2]=0.1250; atmVolMatrix[4][3]=0.1100;
            //atmVolMatrix[5][0]=0.1130; atmVolMatrix[5][1]=0.1090; atmVolMatrix[5][2]=0.1070; atmVolMatrix[5][3]=0.0930;

            atmVolMatrix[0][0]=flatVol; atmVolMatrix[0][1]=flatVol; atmVolMatrix[0][2]=flatVol; atmVolMatrix[0][3]=flatVol;
            atmVolMatrix[1][0]=flatVol; atmVolMatrix[1][1]=flatVol; atmVolMatrix[1][2]=flatVol; atmVolMatrix[1][3]=flatVol;
            atmVolMatrix[2][0]=flatVol; atmVolMatrix[2][1]=flatVol; atmVolMatrix[2][2]=flatVol; atmVolMatrix[2][3]=flatVol;
            atmVolMatrix[3][0]=flatVol; atmVolMatrix[3][1]=flatVol; atmVolMatrix[3][2]=flatVol; atmVolMatrix[3][3]=flatVol;
            atmVolMatrix[4][0]=flatVol; atmVolMatrix[4][1]=flatVol; atmVolMatrix[4][2]=flatVol; atmVolMatrix[4][3]=flatVol;
            atmVolMatrix[5][0]=flatVol; atmVolMatrix[5][1]=flatVol; atmVolMatrix[5][2]=flatVol; atmVolMatrix[5][3]=flatVol;

            Size nRowsAtmVols = atmVolMatrix.rows();
            Size nColsAtmVols = atmVolMatrix.columns();


            //swaptionvolcube
            optionTenors = std::vector<Period>();
            optionTenors.emplace_back(Period(1, Years));
            optionTenors.emplace_back(Period(10, Years));
            optionTenors.emplace_back(Period(30, Years));

            swapTenors = std::vector<Period>();
            swapTenors.emplace_back(Period(2, Years));
            swapTenors.emplace_back(Period(10, Years));
            swapTenors.emplace_back(Period(30, Years));

            strikeSpreads = std::vector<Rate>();
            strikeSpreads.emplace_back(-0.020);
            strikeSpreads.emplace_back(-0.005);
            strikeSpreads.emplace_back(+0.000);
            strikeSpreads.emplace_back(+0.005);
            strikeSpreads.emplace_back(+0.020);

            Size nRows = optionTenors.size()*swapTenors.size();
            Size nCols = strikeSpreads.size();
            volSpreadsMatrix = Matrix(nRows, nCols);
            volSpreadsMatrix[0][0]=0.0599; volSpreadsMatrix[0][1]=0.0049;
            volSpreadsMatrix[0][2]=0.0000;
            volSpreadsMatrix[0][3]=-0.0001; volSpreadsMatrix[0][4]=0.0127;

            volSpreadsMatrix[1][0]=0.0729; volSpreadsMatrix[1][1]=0.0086;
            volSpreadsMatrix[1][2]=0.0000;
            volSpreadsMatrix[1][3]=-0.0024; volSpreadsMatrix[1][4]=0.0098;

            volSpreadsMatrix[2][0]=0.0738; volSpreadsMatrix[2][1]=0.0102;
            volSpreadsMatrix[2][2]=0.0000;
            volSpreadsMatrix[2][3]=-0.0039; volSpreadsMatrix[2][4]=0.0065;

            volSpreadsMatrix[3][0]=0.0465; volSpreadsMatrix[3][1]=0.0063;
            volSpreadsMatrix[3][2]=0.0000;
            volSpreadsMatrix[3][3]=-0.0032; volSpreadsMatrix[3][4]=-0.0010;

            volSpreadsMatrix[4][0]=0.0558; volSpreadsMatrix[4][1]=0.0084;
            volSpreadsMatrix[4][2]=0.0000;
            volSpreadsMatrix[4][3]=-0.0050; volSpreadsMatrix[4][4]=-0.0057;

            volSpreadsMatrix[5][0]=0.0576; volSpreadsMatrix[5][1]=0.0083;
            volSpreadsMatrix[5][2]=0.0000;
            volSpreadsMatrix[5][3]=-0.0043; volSpreadsMatrix[5][4]=-0.0014;

            volSpreadsMatrix[6][0]=0.0437; volSpreadsMatrix[6][1]=0.0059;
            volSpreadsMatrix[6][2]=0.0000;
            volSpreadsMatrix[6][3]=-0.0030; volSpreadsMatrix[6][4]=-0.0006;

            volSpreadsMatrix[7][0]=0.0533; volSpreadsMatrix[7][1]=0.0078;
            volSpreadsMatrix[7][2]=0.0000;
            volSpreadsMatrix[7][3]=-0.0045; volSpreadsMatrix[7][4]=-0.0046;

            volSpreadsMatrix[8][0]=0.0545; volSpreadsMatrix[8][1]=0.0079;
            volSpreadsMatrix[8][2]=0.0000;
            volSpreadsMatrix[8][3]=-0.0042; volSpreadsMatrix[8][4]=-0.0020;


            swapSettlementDays = 2;
            fixedLegFrequency = Annual;
            fixedLegConvention = Unadjusted;
            fixedLegDayCounter = Thirty360();
            std::shared_ptr<SwapIndex> swapIndexBase(new
                EuriborSwapIsdaFixA(2*Years, termStructure));

            std::shared_ptr<SwapIndex> shortSwapIndexBase(new
                EuriborSwapIsdaFixA(1*Years, termStructure));

            vegaWeightedSmileFit = false;

            // ATM Volatility structure
            std::vector<std::vector<Handle<Quote> > > atmVolsHandle;
            atmVolsHandle =
                std::vector<std::vector<Handle<Quote> > >(nRowsAtmVols);
            Size i;
            for (i=0; i<nRowsAtmVols; i++){
                atmVolsHandle[i] = std::vector<Handle<Quote> >(nColsAtmVols);
                for (Size j=0; j<nColsAtmVols; j++) {
                    atmVolsHandle[i][j] =
                        Handle<Quote>(std::shared_ptr<Quote>(new
                            SimpleQuote(atmVolMatrix[i][j])));
                }
            }

            dayCounter = Actual365Fixed();

            atmVol = Handle<SwaptionVolatilityStructure>(
                std::shared_ptr<SwaptionVolatilityStructure>(new
                    SwaptionVolatilityMatrix(calendar,
                                             optionBDC,
                                             atmOptionTenors,
                                             atmSwapTenors,
                                             atmVolsHandle,
                                             dayCounter)));

            // Volatility Cube without smile
            std::vector<std::vector<Handle<Quote> > > parametersGuess(
                                      optionTenors.size()*swapTenors.size());
            for (i=0; i<optionTenors.size()*swapTenors.size(); i++) {
                parametersGuess[i] = std::vector<Handle<Quote> >(4);
                parametersGuess[i][0] =
                    Handle<Quote>(std::shared_ptr<Quote>(
                                                       new SimpleQuote(0.2)));
                parametersGuess[i][1] =
                    Handle<Quote>(std::shared_ptr<Quote>(
                                                       new SimpleQuote(0.5)));
                parametersGuess[i][2] =
                    Handle<Quote>(std::shared_ptr<Quote>(
                                                       new SimpleQuote(0.4)));
                parametersGuess[i][3] =
                    Handle<Quote>(std::shared_ptr<Quote>(
                                                       new SimpleQuote(0.0)));
            }
            std::vector<bool> isParameterFixed(4, false);
            isParameterFixed[1]=true;

            std::vector<std::vector<Handle<Quote> > > nullVolSpreads(nRows);
            for (i=0; i<optionTenors.size()*swapTenors.size(); i++){
                nullVolSpreads[i] = std::vector<Handle<Quote> >(nCols);
                for (Size j=0; j<strikeSpreads.size(); j++) {
                    nullVolSpreads[i][j] =
                        Handle<Quote>(std::shared_ptr<Quote>(
                                                        new SimpleQuote(0.)));
                }
            }

            std::shared_ptr<SwaptionVolCube1>
            flatSwaptionVolatilityCube1ptr(new SwaptionVolCube1(
                atmVol,
                optionTenors,
                swapTenors,
                strikeSpreads,
                nullVolSpreads,
                swapIndexBase,
                shortSwapIndexBase,
                vegaWeightedSmileFit,
                parametersGuess,
                isParameterFixed,
                false));
            flatSwaptionVolatilityCube1 = Handle<SwaptionVolatilityStructure>(
                std::shared_ptr<SwaptionVolatilityStructure>(
                                             flatSwaptionVolatilityCube1ptr));
            flatSwaptionVolatilityCube1->enableExtrapolation();

            std::shared_ptr<SwaptionVolCube2>
            flatSwaptionVolatilityCube2ptr(new SwaptionVolCube2(atmVol,
                                                             optionTenors,
                                                             swapTenors,
                                                             strikeSpreads,
                                                             nullVolSpreads,
                                                             swapIndexBase,
                                                             shortSwapIndexBase,
                                                             vegaWeightedSmileFit));
            flatSwaptionVolatilityCube2 = Handle<SwaptionVolatilityStructure>(
                std::shared_ptr<SwaptionVolatilityStructure>(flatSwaptionVolatilityCube2ptr));
            flatSwaptionVolatilityCube2->enableExtrapolation();


            // Volatility Cube with smile
            volSpreads = std::vector<std::vector<Handle<Quote> > >(nRows);
            for (i=0; i<optionTenors.size()*swapTenors.size(); i++){
                volSpreads[i] = std::vector<Handle<Quote> >(nCols);
                for (Size j=0; j<strikeSpreads.size(); j++) {
                    volSpreads[i][j] =
                        Handle<Quote>(std::shared_ptr<Quote>(new
                            SimpleQuote(volSpreadsMatrix[i][j])));
                }
            }

            std::shared_ptr<SwaptionVolCube1>
            swaptionVolatilityCubeBySabrPtr(new SwaptionVolCube1(
                atmVol,
                optionTenors,
                swapTenors,
                strikeSpreads,
                volSpreads,
                swapIndexBase,
                shortSwapIndexBase,
                vegaWeightedSmileFit,
                parametersGuess,
                isParameterFixed,
                false));
            swaptionVolatilityCubeBySabr = Handle<SwaptionVolatilityStructure>(
            std::shared_ptr<SwaptionVolatilityStructure>(swaptionVolatilityCubeBySabrPtr));
            swaptionVolatilityCubeBySabr->enableExtrapolation();

            swaptionVolatilityStructures = std::vector<Handle<SwaptionVolatilityStructure> >();
            //swaptionVolatilityStructures.emplace_back(atmVol);
            //swaptionVolatilityStructures.emplace_back(flatSwaptionVolatilityCube1);
            swaptionVolatilityStructures.emplace_back(flatSwaptionVolatilityCube2);
            swaptionVolatilityStructures.emplace_back(swaptionVolatilityCubeBySabr);
        }


        void createSmileSections() {
            std::vector<Rate> strikes, stdDevsOnExpiry, stdDevsOnPayment;
            strikes.emplace_back(0.003);   stdDevsOnExpiry.emplace_back(2.45489828353233);    stdDevsOnPayment.emplace_back(1.66175264544155);
            strikes.emplace_back(0.004);   stdDevsOnExpiry.emplace_back(2.10748097295326);    stdDevsOnPayment.emplace_back(1.46691241671427);
            strikes.emplace_back(0.005);   stdDevsOnExpiry.emplace_back(1.87317517200074);    stdDevsOnPayment.emplace_back(1.32415790098009);
            strikes.emplace_back(0.006);   stdDevsOnExpiry.emplace_back(1.69808302023488);    stdDevsOnPayment.emplace_back(1.21209617319357);
            strikes.emplace_back(0.007);   stdDevsOnExpiry.emplace_back(1.55911989073644);    stdDevsOnPayment.emplace_back(1.12016686638666);
            strikes.emplace_back(0.008);   stdDevsOnExpiry.emplace_back(1.44436083444893);    stdDevsOnPayment.emplace_back(1.04242066059821);
            strikes.emplace_back(0.009);   stdDevsOnExpiry.emplace_back(1.34687413874126);    stdDevsOnPayment.emplace_back(0.975173254741177);
            strikes.emplace_back(0.01);    stdDevsOnExpiry.emplace_back(1.26228953588707);    stdDevsOnPayment.emplace_back(0.916013813275761);
            strikes.emplace_back(0.011);   stdDevsOnExpiry.emplace_back(1.18769456816136);    stdDevsOnPayment.emplace_back(0.863267064731419);
            strikes.emplace_back(0.012);   stdDevsOnExpiry.emplace_back(1.12104324191799);    stdDevsOnPayment.emplace_back(0.815743793189994);
            strikes.emplace_back(0.013);   stdDevsOnExpiry.emplace_back(1.06085561121201);    stdDevsOnPayment.emplace_back(0.772552896805455);
            strikes.emplace_back(0.014);   stdDevsOnExpiry.emplace_back(1.00603120341767);    stdDevsOnPayment.emplace_back(0.733033340026564);
            strikes.emplace_back(0.015);   stdDevsOnExpiry.emplace_back(0.955725690399709);   stdDevsOnPayment.emplace_back(0.696673144338147);
            strikes.emplace_back(0.016);   stdDevsOnExpiry.emplace_back(0.909281318404816);   stdDevsOnPayment.emplace_back(0.663070503816902);
            strikes.emplace_back(0.017);   stdDevsOnExpiry.emplace_back(0.866185798452041);   stdDevsOnPayment.emplace_back(0.631911102538957);
            strikes.emplace_back(0.018);   stdDevsOnExpiry.emplace_back(0.826018547612582);   stdDevsOnPayment.emplace_back(0.602948672357772);
            strikes.emplace_back(0.019);   stdDevsOnExpiry.emplace_back(0.788447526732122);   stdDevsOnPayment.emplace_back(0.575982310311697);
            strikes.emplace_back(0.02);    stdDevsOnExpiry.emplace_back(0.753200779931885);   stdDevsOnPayment.emplace_back(0.550849997883271);
            strikes.emplace_back(0.021);   stdDevsOnExpiry.emplace_back(0.720053785498);      stdDevsOnPayment.emplace_back(0.527428600999225);
            strikes.emplace_back(0.022);   stdDevsOnExpiry.emplace_back(0.688823131326177);   stdDevsOnPayment.emplace_back(0.505604706697337);
            strikes.emplace_back(0.023);   stdDevsOnExpiry.emplace_back(0.659357028088728);   stdDevsOnPayment.emplace_back(0.485294065348527);
            strikes.emplace_back(0.024);   stdDevsOnExpiry.emplace_back(0.631532146956907);   stdDevsOnPayment.emplace_back(0.466418908064414);
            strikes.emplace_back(0.025);   stdDevsOnExpiry.emplace_back(0.605247295045587);   stdDevsOnPayment.emplace_back(0.448904706326966);
            strikes.emplace_back(0.026);   stdDevsOnExpiry.emplace_back(0.580413928580285);   stdDevsOnPayment.emplace_back(0.432686652729201);
            strikes.emplace_back(0.027);   stdDevsOnExpiry.emplace_back(0.556962477452476);   stdDevsOnPayment.emplace_back(0.417699939864133);
            strikes.emplace_back(0.028);   stdDevsOnExpiry.emplace_back(0.534829696108958);   stdDevsOnPayment.emplace_back(0.403876519954429);
            strikes.emplace_back(0.029);   stdDevsOnExpiry.emplace_back(0.513968150384827);   stdDevsOnPayment.emplace_back(0.391145104852406);
            strikes.emplace_back(0.03);    stdDevsOnExpiry.emplace_back(0.494330406115181);   stdDevsOnPayment.emplace_back(0.379434406410383);
            strikes.emplace_back(0.031);   stdDevsOnExpiry.emplace_back(0.475869029135118);   stdDevsOnPayment.emplace_back(0.368669896110328);
            strikes.emplace_back(0.032);   stdDevsOnExpiry.emplace_back(0.458549234390376);   stdDevsOnPayment.emplace_back(0.358777045434208);
            strikes.emplace_back(0.033);   stdDevsOnExpiry.emplace_back(0.442329912271372);   stdDevsOnPayment.emplace_back(0.349678085493644);
            strikes.emplace_back(0.034);   stdDevsOnExpiry.emplace_back(0.427163628613205);   stdDevsOnPayment.emplace_back(0.341304968511301);
            strikes.emplace_back(0.035);   stdDevsOnExpiry.emplace_back(0.413009273806291);   stdDevsOnPayment.emplace_back(0.333586406339497);
            strikes.emplace_back(0.036);   stdDevsOnExpiry.emplace_back(0.399819413685729);   stdDevsOnPayment.emplace_back(0.326457591571248);
            strikes.emplace_back(0.037);   stdDevsOnExpiry.emplace_back(0.387546614086615);   stdDevsOnPayment.emplace_back(0.31985630909585);
            strikes.emplace_back(0.038);   stdDevsOnExpiry.emplace_back(0.376137116288728);   stdDevsOnPayment.emplace_back(0.313728768765505);
            strikes.emplace_back(0.039);   stdDevsOnExpiry.emplace_back(0.365540323849504);   stdDevsOnPayment.emplace_back(0.308024420802767);
            strikes.emplace_back(0.04);    stdDevsOnExpiry.emplace_back(0.35570564032638);    stdDevsOnPayment.emplace_back(0.30269822405978);
            strikes.emplace_back(0.041);   stdDevsOnExpiry.emplace_back(0.346572982443814);   stdDevsOnPayment.emplace_back(0.297710321981251);
            strikes.emplace_back(0.042);   stdDevsOnExpiry.emplace_back(0.338091753759242);   stdDevsOnPayment.emplace_back(0.293025394530372);
            strikes.emplace_back(0.043);   stdDevsOnExpiry.emplace_back(0.330211357830103);   stdDevsOnPayment.emplace_back(0.288612334151791);
            strikes.emplace_back(0.044);   stdDevsOnExpiry.emplace_back(0.322881198213832);   stdDevsOnPayment.emplace_back(0.284443273660505);
            strikes.emplace_back(0.045);   stdDevsOnExpiry.emplace_back(0.316056686795423);   stdDevsOnPayment.emplace_back(0.280494558352965);
            strikes.emplace_back(0.046);   stdDevsOnExpiry.emplace_back(0.309691654321036);   stdDevsOnPayment.emplace_back(0.276744153710797);
            strikes.emplace_back(0.047);   stdDevsOnExpiry.emplace_back(0.303745307408855);   stdDevsOnPayment.emplace_back(0.273174237697079);
            strikes.emplace_back(0.048);   stdDevsOnExpiry.emplace_back(0.298180014954725);   stdDevsOnPayment.emplace_back(0.269767960385995);
            strikes.emplace_back(0.049);   stdDevsOnExpiry.emplace_back(0.292961308132149);   stdDevsOnPayment.emplace_back(0.266511064148011);
            strikes.emplace_back(0.05);    stdDevsOnExpiry.emplace_back(0.288057880392292);   stdDevsOnPayment.emplace_back(0.263391235575797);
            strikes.emplace_back(0.051);   stdDevsOnExpiry.emplace_back(0.283441587463978);   stdDevsOnPayment.emplace_back(0.260399077595342);
            strikes.emplace_back(0.052);   stdDevsOnExpiry.emplace_back(0.279088079809224);   stdDevsOnPayment.emplace_back(0.257518712391935);
            strikes.emplace_back(0.053);   stdDevsOnExpiry.emplace_back(0.274968896929089);   stdDevsOnPayment.emplace_back(0.254747223632261);
            strikes.emplace_back(0.054);   stdDevsOnExpiry.emplace_back(0.271067594979739);   stdDevsOnPayment.emplace_back(0.252074566168237);
            strikes.emplace_back(0.055);   stdDevsOnExpiry.emplace_back(0.267364567839682);   stdDevsOnPayment.emplace_back(0.249494259259166);
            strikes.emplace_back(0.056);   stdDevsOnExpiry.emplace_back(0.263842422981787);   stdDevsOnPayment.emplace_back(0.246999498127314);
            strikes.emplace_back(0.057);   stdDevsOnExpiry.emplace_back(0.26048629770105);    stdDevsOnPayment.emplace_back(0.244584774143087);
            strikes.emplace_back(0.058);   stdDevsOnExpiry.emplace_back(0.257282594203533);   stdDevsOnPayment.emplace_back(0.242244902713927);
            strikes.emplace_back(0.059);   stdDevsOnExpiry.emplace_back(0.254218979606362);   stdDevsOnPayment.emplace_back(0.23997567135838);
            strikes.emplace_back(0.06);    stdDevsOnExpiry.emplace_back(0.251284385937726);   stdDevsOnPayment.emplace_back(0.237772543557956);
            strikes.emplace_back(0.061);   stdDevsOnExpiry.emplace_back(0.248469326364644);   stdDevsOnPayment.emplace_back(0.235632278942307);
            strikes.emplace_back(0.062);   stdDevsOnExpiry.emplace_back(0.245764630281902);   stdDevsOnPayment.emplace_back(0.233550665029978);
            strikes.emplace_back(0.063);   stdDevsOnExpiry.emplace_back(0.243162391995349);   stdDevsOnPayment.emplace_back(0.231525109524691);
            strikes.emplace_back(0.064);   stdDevsOnExpiry.emplace_back(0.240655338266368);   stdDevsOnPayment.emplace_back(0.22955269609313);
            strikes.emplace_back(0.065);   stdDevsOnExpiry.emplace_back(0.238237144539637);   stdDevsOnPayment.emplace_back(0.227630508401982);
            strikes.emplace_back(0.066);   stdDevsOnExpiry.emplace_back(0.235901802487603);   stdDevsOnPayment.emplace_back(0.225756278192003);
            strikes.emplace_back(0.067);   stdDevsOnExpiry.emplace_back(0.233643936238243);   stdDevsOnPayment.emplace_back(0.223927413166912);
            strikes.emplace_back(0.068);   stdDevsOnExpiry.emplace_back(0.2314584861473);     stdDevsOnPayment.emplace_back(0.222142617178571);
            strikes.emplace_back(0.069);   stdDevsOnExpiry.emplace_back(0.229341341253818);   stdDevsOnPayment.emplace_back(0.220398973893664);
            strikes.emplace_back(0.07);    stdDevsOnExpiry.emplace_back(0.22728807436907);    stdDevsOnPayment.emplace_back(0.218695187164053);
            strikes.emplace_back(0.071);   stdDevsOnExpiry.emplace_back(0.225295206987632);   stdDevsOnPayment.emplace_back(0.217029636804562);
            strikes.emplace_back(0.072);   stdDevsOnExpiry.emplace_back(0.223359576831843);   stdDevsOnPayment.emplace_back(0.215400702630017);
            strikes.emplace_back(0.073);   stdDevsOnExpiry.emplace_back(0.221477389168511);   stdDevsOnPayment.emplace_back(0.213806764455244);
            strikes.emplace_back(0.074);   stdDevsOnExpiry.emplace_back(0.219646430403273);   stdDevsOnPayment.emplace_back(0.212246202095067);
            strikes.emplace_back(0.075);   stdDevsOnExpiry.emplace_back(0.21786353825847);    stdDevsOnPayment.emplace_back(0.210718367475417);
            strikes.emplace_back(0.076);   stdDevsOnExpiry.emplace_back(0.21612649913974);    stdDevsOnPayment.emplace_back(0.20922164041112);
            strikes.emplace_back(0.077);   stdDevsOnExpiry.emplace_back(0.214433415680486);   stdDevsOnPayment.emplace_back(0.20775504879107);
            strikes.emplace_back(0.078);   stdDevsOnExpiry.emplace_back(0.212781441830814);   stdDevsOnPayment.emplace_back(0.206317296467129);
            strikes.emplace_back(0.079);   stdDevsOnExpiry.emplace_back(0.21116931267966);    stdDevsOnPayment.emplace_back(0.20490741132819);
            strikes.emplace_back(0.08);    stdDevsOnExpiry.emplace_back(0.209594814632662);   stdDevsOnPayment.emplace_back(0.203524745300185);
            strikes.emplace_back(0.081);   stdDevsOnExpiry.emplace_back(0.20805636655099);    stdDevsOnPayment.emplace_back(0.202168002234973);
            strikes.emplace_back(0.082);   stdDevsOnExpiry.emplace_back(0.20655270352358);    stdDevsOnPayment.emplace_back(0.20083621002145);
            strikes.emplace_back(0.083);   stdDevsOnExpiry.emplace_back(0.20508161195607);    stdDevsOnPayment.emplace_back(0.199529044622581);
            strikes.emplace_back(0.084);   stdDevsOnExpiry.emplace_back(0.203642775620693);   stdDevsOnPayment.emplace_back(0.198245209890227);
            strikes.emplace_back(0.085);   stdDevsOnExpiry.emplace_back(0.202233980923088);   stdDevsOnPayment.emplace_back(0.196984381787351);
            strikes.emplace_back(0.086);   stdDevsOnExpiry.emplace_back(0.200854279179957);   stdDevsOnPayment.emplace_back(0.195745912239886);
            strikes.emplace_back(0.087);   stdDevsOnExpiry.emplace_back(0.199503037935767);   stdDevsOnPayment.emplace_back(0.19452850509969);
            strikes.emplace_back(0.088);   stdDevsOnExpiry.emplace_back(0.198178676051688);   stdDevsOnPayment.emplace_back(0.193332160366764);
            strikes.emplace_back(0.089);   stdDevsOnExpiry.emplace_back(0.196880244844423);   stdDevsOnPayment.emplace_back(0.192155905930003);
            strikes.emplace_back(0.09);    stdDevsOnExpiry.emplace_back(0.195606795630673);   stdDevsOnPayment.emplace_back(0.190999417752372);
            strikes.emplace_back(0.091);   stdDevsOnExpiry.emplace_back(0.194357695954907);   stdDevsOnPayment.emplace_back(0.189861723722766);
            strikes.emplace_back(0.092);   stdDevsOnExpiry.emplace_back(0.19313168090606);    stdDevsOnPayment.emplace_back(0.188742823841186);
            strikes.emplace_back(0.093);   stdDevsOnExpiry.emplace_back(0.191928434256365);   stdDevsOnPayment.emplace_back(0.187641745996527);
            strikes.emplace_back(0.094);   stdDevsOnExpiry.emplace_back(0.190746691094761);   stdDevsOnPayment.emplace_back(0.186558166151753);
            strikes.emplace_back(0.095);   stdDevsOnExpiry.emplace_back(0.189586451421245);   stdDevsOnPayment.emplace_back(0.185491436232795);
            strikes.emplace_back(0.096);   stdDevsOnExpiry.emplace_back(0.188446134096988);   stdDevsOnPayment.emplace_back(0.184441556239653);
            strikes.emplace_back(0.097);   stdDevsOnExpiry.emplace_back(0.18732573912199);    stdDevsOnPayment.emplace_back(0.183407878098257);
            strikes.emplace_back(0.098);   stdDevsOnExpiry.emplace_back(0.186224317812954);   stdDevsOnPayment.emplace_back(0.182390725845642);
            strikes.emplace_back(0.099);   stdDevsOnExpiry.emplace_back(0.185141553942112);   stdDevsOnPayment.emplace_back(0.181386859111458);
            strikes.emplace_back(0.1);     stdDevsOnExpiry.emplace_back(0.184076498826167);   stdDevsOnPayment.emplace_back(0.180399194229021);
            strikes.emplace_back(0.101);   stdDevsOnExpiry.emplace_back(0.18302915246512);    stdDevsOnPayment.emplace_back(0.17942643505019);
            strikes.emplace_back(0.102);   stdDevsOnExpiry.emplace_back(0.181999514858969);   stdDevsOnPayment.emplace_back(0.178466637352756);
            strikes.emplace_back(0.103);   stdDevsOnExpiry.emplace_back(0.180984739957821);   stdDevsOnPayment.emplace_back(0.177521421321893);
            strikes.emplace_back(0.104);   stdDevsOnExpiry.emplace_back(0.179986725128272);   stdDevsOnPayment.emplace_back(0.176590462920567);
            strikes.emplace_back(0.105);   stdDevsOnExpiry.emplace_back(0.179004521687023);   stdDevsOnPayment.emplace_back(0.175677650593196);
            strikes.emplace_back(0.106);   stdDevsOnExpiry.emplace_back(0.178041924367268);   stdDevsOnPayment.emplace_back(0.17476516230286);
            strikes.emplace_back(0.107);   stdDevsOnExpiry.emplace_back(0.177083754236237);   stdDevsOnPayment.emplace_back(0.173873088345724);
            strikes.emplace_back(0.108);   stdDevsOnExpiry.emplace_back(0.176145822682231);   stdDevsOnPayment.emplace_back(0.173000456610684);
            strikes.emplace_back(0.109);   stdDevsOnExpiry.emplace_back(0.175227181021952);   stdDevsOnPayment.emplace_back(0.172122316246049);
            strikes.emplace_back(0.11);    stdDevsOnExpiry.emplace_back(0.174309488044971);   stdDevsOnPayment.emplace_back(0.171266858473859);
            strikes.emplace_back(0.111);   stdDevsOnExpiry.emplace_back(0.173412982328314);   stdDevsOnPayment.emplace_back(0.170434407331149);
            strikes.emplace_back(0.112);   stdDevsOnExpiry.emplace_back(0.172536715188681);   stdDevsOnPayment.emplace_back(0.169585106262623);
            strikes.emplace_back(0.113);   stdDevsOnExpiry.emplace_back(0.171706301075121);   stdDevsOnPayment.emplace_back(0.168765292564274);
            strikes.emplace_back(0.114);   stdDevsOnExpiry.emplace_back(0.17079651379229);    stdDevsOnPayment.emplace_back(0.167976586421278);
            strikes.emplace_back(0.115);   stdDevsOnExpiry.emplace_back(0.169963569856602);   stdDevsOnPayment.emplace_back(0.167267917425907);
            strikes.emplace_back(0.116);   stdDevsOnExpiry.emplace_back(0.169192922790819);   stdDevsOnPayment.emplace_back(0.166364178135514);
            strikes.emplace_back(0.117);   stdDevsOnExpiry.emplace_back(0.168289776291075);   stdDevsOnPayment.emplace_back(0.165629586177349);
            strikes.emplace_back(0.118);   stdDevsOnExpiry.emplace_back(0.167505847659119);   stdDevsOnPayment.emplace_back(0.165014239848036);
            strikes.emplace_back(0.119);   stdDevsOnExpiry.emplace_back(0.166813308851542);   stdDevsOnPayment.emplace_back(0.164618590628398);
            strikes.emplace_back(0.12);    stdDevsOnExpiry.emplace_back(0.166305130831553);   stdDevsOnPayment.emplace_back(0.164530452554899);
            strikes.emplace_back(0.121);   stdDevsOnExpiry.emplace_back(0.166077130612255);   stdDevsOnPayment.emplace_back(0.162925173083904);
            strikes.emplace_back(0.122);   stdDevsOnExpiry.emplace_back(0.164586116695486);   stdDevsOnPayment.emplace_back(0.162717141307485);
            strikes.emplace_back(0.123);   stdDevsOnExpiry.emplace_back(0.164242693341591);   stdDevsOnPayment.emplace_back(0.162840275380755);
            strikes.emplace_back(0.124);   stdDevsOnExpiry.emplace_back(0.164213284159352);   stdDevsOnPayment.emplace_back(0.163289714748189);
            strikes.emplace_back(0.125);   stdDevsOnExpiry.emplace_back(0.164516546586962);   stdDevsOnPayment.emplace_back(0.16401944615083);
            strikes.emplace_back(0.126);   stdDevsOnExpiry.emplace_back(0.165118644253458);   stdDevsOnPayment.emplace_back(0.164961421811344);
            strikes.emplace_back(0.127);   stdDevsOnExpiry.emplace_back(0.165959810111063);   stdDevsOnPayment.emplace_back(0.166058935248619);
            strikes.emplace_back(0.128);   stdDevsOnExpiry.emplace_back(0.166976798606573);   stdDevsOnPayment.emplace_back(0.16725625209265);
            strikes.emplace_back(0.129);   stdDevsOnExpiry.emplace_back(0.168115851019766);   stdDevsOnPayment.emplace_back(0.16851675615849);
            strikes.emplace_back(0.13);    stdDevsOnExpiry.emplace_back(0.169332063007866);   stdDevsOnPayment.emplace_back(0.16981808889073);
            strikes.emplace_back(0.131);   stdDevsOnExpiry.emplace_back(0.170600136349594);   stdDevsOnPayment.emplace_back(0.171139511919136);
            strikes.emplace_back(0.132);   stdDevsOnExpiry.emplace_back(0.171891926773773);   stdDevsOnPayment.emplace_back(0.172468711836379);
            strikes.emplace_back(0.133);   stdDevsOnExpiry.emplace_back(0.173201742180614);   stdDevsOnPayment.emplace_back(0.173801476161007);
            strikes.emplace_back(0.134);   stdDevsOnExpiry.emplace_back(0.17451282249852);    stdDevsOnPayment.emplace_back(0.175129703967145);
            strikes.emplace_back(0.135);   stdDevsOnExpiry.emplace_back(0.175823902816426);   stdDevsOnPayment.emplace_back(0.17645371929183);
            strikes.emplace_back(0.136);   stdDevsOnExpiry.emplace_back(0.177132453312204);   stdDevsOnPayment.emplace_back(0.177767365431397);
            strikes.emplace_back(0.137);   stdDevsOnExpiry.emplace_back(0.178433098113831);   stdDevsOnPayment.emplace_back(0.179076475052476);
            strikes.emplace_back(0.138);   stdDevsOnExpiry.emplace_back(0.17972646967684);    stdDevsOnPayment.emplace_back(0.180372947229192);
            strikes.emplace_back(0.139);   stdDevsOnExpiry.emplace_back(0.181011935545698);   stdDevsOnPayment.emplace_back(0.181660994443001);
            strikes.emplace_back(0.14);    stdDevsOnExpiry.emplace_back(0.182286965898278);   stdDevsOnPayment.emplace_back(0.182938996508727);
            strikes.emplace_back(0.141);   stdDevsOnExpiry.emplace_back(0.18355314187341);    stdDevsOnPayment.emplace_back(0.18420889764858);
            strikes.emplace_back(0.142);   stdDevsOnExpiry.emplace_back(0.184810147243326);   stdDevsOnPayment.emplace_back(0.185468105566281);
            strikes.emplace_back(0.143);   stdDevsOnExpiry.emplace_back(0.186056717096965);   stdDevsOnPayment.emplace_back(0.186718888521073);
            strikes.emplace_back(0.144);   stdDevsOnExpiry.emplace_back(0.187295381256453);   stdDevsOnPayment.emplace_back(0.187958006142609);
            strikes.emplace_back(0.145);   stdDevsOnExpiry.emplace_back(0.188523609899662);   stdDevsOnPayment.emplace_back(0.189190318986411);
            strikes.emplace_back(0.146);   stdDevsOnExpiry.emplace_back(0.189745197759785);   stdDevsOnPayment.emplace_back(0.190412586682131);
            strikes.emplace_back(0.147);   stdDevsOnExpiry.emplace_back(0.190955085192566);   stdDevsOnPayment.emplace_back(0.191624809229768);
            strikes.emplace_back(0.148);   stdDevsOnExpiry.emplace_back(0.186502914474815);   stdDevsOnPayment.emplace_back(0.192830226999672);
            strikes.emplace_back(0.149);   stdDevsOnExpiry.emplace_back(0.187658094504074);   stdDevsOnPayment.emplace_back(0.194024951547423);
            strikes.emplace_back(0.15);    stdDevsOnExpiry.emplace_back(0.188817069266526);   stdDevsOnPayment.emplace_back(0.195212547280407);
            strikes.emplace_back(0.151);   stdDevsOnExpiry.emplace_back(0.189958019046315);   stdDevsOnPayment.emplace_back(0.196391394013447);
            strikes.emplace_back(0.152);   stdDevsOnExpiry.emplace_back(0.191090746904187);   stdDevsOnPayment.emplace_back(0.197560195598405);
            strikes.emplace_back(0.153);   stdDevsOnExpiry.emplace_back(0.192215885295675);   stdDevsOnPayment.emplace_back(0.19871895203528);
            strikes.emplace_back(0.154);   stdDevsOnExpiry.emplace_back(0.193335331587374);   stdDevsOnPayment.emplace_back(0.199872523879597);
            strikes.emplace_back(0.155);   stdDevsOnExpiry.emplace_back(0.194446555957158);   stdDevsOnPayment.emplace_back(0.195112095799581);
            strikes.emplace_back(0.156);   stdDevsOnExpiry.emplace_back(0.195547028582896);   stdDevsOnPayment.emplace_back(0.196220302459009);
            strikes.emplace_back(0.157);   stdDevsOnExpiry.emplace_back(0.196646236297571);   stdDevsOnPayment.emplace_back(0.197317167822215);
            strikes.emplace_back(0.158);   stdDevsOnExpiry.emplace_back(0.197736589634797);   stdDevsOnPayment.emplace_back(0.198405608222512);
            strikes.emplace_back(0.159);   stdDevsOnExpiry.emplace_back(0.198811131583722);   stdDevsOnPayment.emplace_back(0.19949340054874);
            strikes.emplace_back(0.16);    stdDevsOnExpiry.emplace_back(0.199887570899243);   stdDevsOnPayment.emplace_back(0.200565963134326);
            strikes.emplace_back(0.161);   stdDevsOnExpiry.emplace_back(0.20095167733189);    stdDevsOnPayment.emplace_back(0.201636905534738);
            strikes.emplace_back(0.162);   stdDevsOnExpiry.emplace_back(0.20200756184262);    stdDevsOnPayment.emplace_back(0.202695534527823);
            strikes.emplace_back(0.163);   stdDevsOnExpiry.emplace_back(0.203061232758988);   stdDevsOnPayment.emplace_back(0.203753839483873);
            strikes.emplace_back(0.164);   stdDevsOnExpiry.emplace_back(0.204112690080994);   stdDevsOnPayment.emplace_back(0.204791730106723);
            strikes.emplace_back(0.165);   stdDevsOnExpiry.emplace_back(0.205146754875869);   stdDevsOnPayment.emplace_back(0.205839341840621);
            strikes.emplace_back(0.166);   stdDevsOnExpiry.emplace_back(0.206178289848616);   stdDevsOnPayment.emplace_back(0.206869779611668);
            strikes.emplace_back(0.167);   stdDevsOnExpiry.emplace_back(0.207207294999235);   stdDevsOnPayment.emplace_back(0.207893412604981);
            strikes.emplace_back(0.168);   stdDevsOnExpiry.emplace_back(0.208224599722511);   stdDevsOnPayment.emplace_back(0.208916397524225);
            strikes.emplace_back(0.169);   stdDevsOnExpiry.emplace_back(0.209234947434935);   stdDevsOnPayment.emplace_back(0.209924476739862);
            strikes.emplace_back(0.17);    stdDevsOnExpiry.emplace_back(0.210235175858846);   stdDevsOnPayment.emplace_back(0.210934824214744);
            strikes.emplace_back(0.171);   stdDevsOnExpiry.emplace_back(0.211231609549565);   stdDevsOnPayment.emplace_back(0.211933506356369);
            strikes.emplace_back(0.172);   stdDevsOnExpiry.emplace_back(0.212231205517945);   stdDevsOnPayment.emplace_back(0.212931216386889);
            strikes.emplace_back(0.173);   stdDevsOnExpiry.emplace_back(0.213219101058981);   stdDevsOnPayment.emplace_back(0.213916613010082);
            strikes.emplace_back(0.174);   stdDevsOnExpiry.emplace_back(0.214192133895015);   stdDevsOnPayment.emplace_back(0.2149003894481);
            strikes.emplace_back(0.175);   stdDevsOnExpiry.emplace_back(0.215167064097645);   stdDevsOnPayment.emplace_back(0.215876064960245);
            strikes.emplace_back(0.176);   stdDevsOnExpiry.emplace_back(0.216146105261233);   stdDevsOnPayment.emplace_back(0.216845259731692);
            strikes.emplace_back(0.177);   stdDevsOnExpiry.emplace_back(0.217099215748008);   stdDevsOnPayment.emplace_back(0.217806029540231);
            strikes.emplace_back(0.178);   stdDevsOnExpiry.emplace_back(0.218056437195741);   stdDevsOnPayment.emplace_back(0.218763558978421);
            strikes.emplace_back(0.179);   stdDevsOnExpiry.emplace_back(0.219005120493791);   stdDevsOnPayment.emplace_back(0.219722060527715);
            strikes.emplace_back(0.18);    stdDevsOnExpiry.emplace_back(0.219951273969714);   stdDevsOnPayment.emplace_back(0.220671489040032);
            strikes.emplace_back(0.181);   stdDevsOnExpiry.emplace_back(0.220885410790527);   stdDevsOnPayment.emplace_back(0.221608280107987);
            strikes.emplace_back(0.182);   stdDevsOnExpiry.emplace_back(0.221831248038684);   stdDevsOnPayment.emplace_back(0.222542154842628);
            strikes.emplace_back(0.183);   stdDevsOnExpiry.emplace_back(0.222757162937581);   stdDevsOnPayment.emplace_back(0.223469224799535);
            strikes.emplace_back(0.184);   stdDevsOnExpiry.emplace_back(0.223673907231264);   stdDevsOnPayment.emplace_back(0.224396942830512);
            strikes.emplace_back(0.185);   stdDevsOnExpiry.emplace_back(0.224599189674629);   stdDevsOnPayment.emplace_back(0.225315911861546);
            strikes.emplace_back(0.186);   stdDevsOnExpiry.emplace_back(0.225503601085437);   stdDevsOnPayment.emplace_back(0.226230992448161);
            strikes.emplace_back(0.187);   stdDevsOnExpiry.emplace_back(0.226412755912736);   stdDevsOnPayment.emplace_back(0.227133759627449);
            strikes.emplace_back(0.188);   stdDevsOnExpiry.emplace_back(0.227313372590352);   stdDevsOnPayment.emplace_back(0.228049488288135);
            strikes.emplace_back(0.189);   stdDevsOnExpiry.emplace_back(0.228216519090096);   stdDevsOnPayment.emplace_back(0.228941886282305);
            strikes.emplace_back(0.19);    stdDevsOnExpiry.emplace_back(0.229108597618029);   stdDevsOnPayment.emplace_back(0.229833312165371);
            strikes.emplace_back(0.191);   stdDevsOnExpiry.emplace_back(0.229988343263088);   stdDevsOnPayment.emplace_back(0.230724738048437);
            strikes.emplace_back(0.192);   stdDevsOnExpiry.emplace_back(0.230883267840916);   stdDevsOnPayment.emplace_back(0.231599638042722);
            strikes.emplace_back(0.193);   stdDevsOnExpiry.emplace_back(0.231748467008738);   stdDevsOnPayment.emplace_back(0.232481018777706);
            strikes.emplace_back(0.194);   stdDevsOnExpiry.emplace_back(0.232617460909752);   stdDevsOnPayment.emplace_back(0.233360131253445);
            strikes.emplace_back(0.195);   stdDevsOnExpiry.emplace_back(0.233487087266298);   stdDevsOnPayment.emplace_back(0.234224662062612);
            strikes.emplace_back(0.196);   stdDevsOnExpiry.emplace_back(0.234348491700928);   stdDevsOnPayment.emplace_back(0.235093405353234);
            strikes.emplace_back(0.197);   stdDevsOnExpiry.emplace_back(0.2352057851746);     stdDevsOnPayment.emplace_back(0.235946918903214);
            strikes.emplace_back(0.198);   stdDevsOnExpiry.emplace_back(0.236064976014868);   stdDevsOnPayment.emplace_back(0.236808209342033);
            strikes.emplace_back(0.199);   stdDevsOnExpiry.emplace_back(0.236907723011302);   stdDevsOnPayment.emplace_back(0.237655242151315);
            strikes.emplace_back(0.2);     stdDevsOnExpiry.emplace_back(0.237747940185609);   stdDevsOnPayment.emplace_back(0.238496766331003);
            strikes.emplace_back(0.201);   stdDevsOnExpiry.emplace_back(0.238590687182044);   stdDevsOnPayment.emplace_back(0.239337318399586);
            strikes.emplace_back(0.202);   stdDevsOnExpiry.emplace_back(0.239419203929008);   stdDevsOnPayment.emplace_back(0.24017916661631);
            strikes.emplace_back(0.203);   stdDevsOnExpiry.emplace_back(0.240250882953632);   stdDevsOnPayment.emplace_back(0.241003840870182);
            strikes.emplace_back(0.204);   stdDevsOnExpiry.emplace_back(0.241075604967404);   stdDevsOnPayment.emplace_back(0.241832727605508);
            strikes.emplace_back(0.205);   stdDevsOnExpiry.emplace_back(0.24190317303107);    stdDevsOnPayment.emplace_back(0.242654809563101);
            strikes.emplace_back(0.206);   stdDevsOnExpiry.emplace_back(0.24272251917282);    stdDevsOnPayment.emplace_back(0.243478511705869);
            strikes.emplace_back(0.207);   stdDevsOnExpiry.emplace_back(0.243530164887227);   stdDevsOnPayment.emplace_back(0.244289576404275);
            strikes.emplace_back(0.208);   stdDevsOnExpiry.emplace_back(0.24434287024589);    stdDevsOnPayment.emplace_back(0.245096752658261);
            strikes.emplace_back(0.209);   stdDevsOnExpiry.emplace_back(0.245145772543807);   stdDevsOnPayment.emplace_back(0.245903280838178);
            strikes.emplace_back(0.21);    stdDevsOnExpiry.emplace_back(0.245951204663852);   stdDevsOnPayment.emplace_back(0.246714021499549);
            strikes.emplace_back(0.211);   stdDevsOnExpiry.emplace_back(0.246737030662404);   stdDevsOnPayment.emplace_back(0.247507912235104);
            strikes.emplace_back(0.212);   stdDevsOnExpiry.emplace_back(0.247526967621914);   stdDevsOnPayment.emplace_back(0.248306663526183);
            strikes.emplace_back(0.213);   stdDevsOnExpiry.emplace_back(0.24831563967036);    stdDevsOnPayment.emplace_back(0.2490927773729);
            strikes.emplace_back(0.214);   stdDevsOnExpiry.emplace_back(0.249111584957424);   stdDevsOnPayment.emplace_back(0.249878243145547);
            strikes.emplace_back(0.215);   stdDevsOnExpiry.emplace_back(0.249880650884377);   stdDevsOnPayment.emplace_back(0.250649127251622);
            strikes.emplace_back(0.216);   stdDevsOnExpiry.emplace_back(0.250656357594417);   stdDevsOnPayment.emplace_back(0.251433620913165);
            strikes.emplace_back(0.217);   stdDevsOnExpiry.emplace_back(0.251434594126584);   stdDevsOnPayment.emplace_back(0.25220580116738);
            strikes.emplace_back(0.218);   stdDevsOnExpiry.emplace_back(0.252199549092579);   stdDevsOnPayment.emplace_back(0.252970852606827);
            strikes.emplace_back(0.219);   stdDevsOnExpiry.emplace_back(0.252961025553147);   stdDevsOnPayment.emplace_back(0.253737200194414);
            strikes.emplace_back(0.22);    stdDevsOnExpiry.emplace_back(0.253727877885738);   stdDevsOnPayment.emplace_back(0.254501279522756);
            strikes.emplace_back(0.221);   stdDevsOnExpiry.emplace_back(0.254480499968858);   stdDevsOnPayment.emplace_back(0.255269571332552);
            strikes.emplace_back(0.222);   stdDevsOnExpiry.emplace_back(0.25523533564634);    stdDevsOnPayment.emplace_back(0.256016476698044);
            strikes.emplace_back(0.223);   stdDevsOnExpiry.emplace_back(0.255984162996268);   stdDevsOnPayment.emplace_back(0.256770834915338);
            strikes.emplace_back(0.224);   stdDevsOnExpiry.emplace_back(0.25673583639609);    stdDevsOnPayment.emplace_back(0.257510611466062);
            strikes.emplace_back(0.225);   stdDevsOnExpiry.emplace_back(0.257477706735166);   stdDevsOnPayment.emplace_back(0.258255572609344);
            strikes.emplace_back(0.226);   stdDevsOnExpiry.emplace_back(0.258220525757539);   stdDevsOnPayment.emplace_back(0.25900280201187);
            strikes.emplace_back(0.227);   stdDevsOnExpiry.emplace_back(0.258953541719166);   stdDevsOnPayment.emplace_back(0.259739986266314);
            strikes.emplace_back(0.228);   stdDevsOnExpiry.emplace_back(0.259691301097284);   stdDevsOnPayment.emplace_back(0.260471337854129);
            strikes.emplace_back(0.229);   stdDevsOnExpiry.emplace_back(0.260414197770398);   stdDevsOnPayment.emplace_back(0.26120171733084);
            strikes.emplace_back(0.23);    stdDevsOnExpiry.emplace_back(0.261138359354577);   stdDevsOnPayment.emplace_back(0.26193598525197);
            strikes.emplace_back(0.231);   stdDevsOnExpiry.emplace_back(0.261857461294499);   stdDevsOnPayment.emplace_back(0.262664096469436);
            strikes.emplace_back(0.232);   stdDevsOnExpiry.emplace_back(0.262581939106444);   stdDevsOnPayment.emplace_back(0.263373413538876);
            strikes.emplace_back(0.233);   stdDevsOnExpiry.emplace_back(0.263292502896683);   stdDevsOnPayment.emplace_back(0.264094071904539);
            strikes.emplace_back(0.234);   stdDevsOnExpiry.emplace_back(0.264006545192349);   stdDevsOnPayment.emplace_back(0.26480533319619);
            strikes.emplace_back(0.235);   stdDevsOnExpiry.emplace_back(0.264707305921843);   stdDevsOnPayment.emplace_back(0.265520482932259);
            strikes.emplace_back(0.236);   stdDevsOnExpiry.emplace_back(0.265417869712082);   stdDevsOnPayment.emplace_back(0.266215866409198);
            strikes.emplace_back(0.237);   stdDevsOnExpiry.emplace_back(0.266128433502322);   stdDevsOnPayment.emplace_back(0.266921619071255);
            strikes.emplace_back(0.238);   stdDevsOnExpiry.emplace_back(0.266818442487771);   stdDevsOnPayment.emplace_back(0.267621215029648);
            strikes.emplace_back(0.239);   stdDevsOnExpiry.emplace_back(0.267506237878858);   stdDevsOnPayment.emplace_back(0.268319190802866);
            strikes.emplace_back(0.24);    stdDevsOnExpiry.emplace_back(0.268213955619203);   stdDevsOnPayment.emplace_back(0.269024295390853);
            strikes.emplace_back(0.241);   stdDevsOnExpiry.emplace_back(0.268901118554758);   stdDevsOnPayment.emplace_back(0.269714494275234);
            strikes.emplace_back(0.242);   stdDevsOnExpiry.emplace_back(0.269581956934992);   stdDevsOnPayment.emplace_back(0.270383630752344);
            strikes.emplace_back(0.243);   stdDevsOnExpiry.emplace_back(0.270257103215438);   stdDevsOnPayment.emplace_back(0.271079662303353);
            strikes.emplace_back(0.244);   stdDevsOnExpiry.emplace_back(0.270943317467695);   stdDevsOnPayment.emplace_back(0.271764028521105);
            strikes.emplace_back(0.245);   stdDevsOnExpiry.emplace_back(0.271623207164631);   stdDevsOnPayment.emplace_back(0.272445154368508);
            strikes.emplace_back(0.246);   stdDevsOnExpiry.emplace_back(0.272295191167417);   stdDevsOnPayment.emplace_back(0.273122067734456);
            strikes.emplace_back(0.247);   stdDevsOnExpiry.emplace_back(0.272981721647439);   stdDevsOnPayment.emplace_back(0.273801573396684);
            strikes.emplace_back(0.248);   stdDevsOnExpiry.emplace_back(0.2736334670732);     stdDevsOnPayment.emplace_back(0.274467145466411);
            strikes.emplace_back(0.249);   stdDevsOnExpiry.emplace_back(0.274298494065133);   stdDevsOnPayment.emplace_back(0.275112951277007);
            strikes.emplace_back(0.25);    stdDevsOnExpiry.emplace_back(0.274975221484409);   stdDevsOnPayment.emplace_back(0.275792456939235);
            strikes.emplace_back(0.251);   stdDevsOnExpiry.emplace_back(0.275627599365702);   stdDevsOnPayment.emplace_back(0.276450252120124);
            strikes.emplace_back(0.252);   stdDevsOnExpiry.emplace_back(0.276287250485613);   stdDevsOnPayment.emplace_back(0.277106427115838);
            strikes.emplace_back(0.253);   stdDevsOnExpiry.emplace_back(0.27693614986148);    stdDevsOnPayment.emplace_back(0.277760009815272);
            strikes.emplace_back(0.254);   stdDevsOnExpiry.emplace_back(0.277595168525859);   stdDevsOnPayment.emplace_back(0.278409055996218);
            strikes.emplace_back(0.255);   stdDevsOnExpiry.emplace_back(0.278247230179386);   stdDevsOnPayment.emplace_back(0.279076248251119);
            strikes.emplace_back(0.256);   stdDevsOnExpiry.emplace_back(0.27889233482206);    stdDevsOnPayment.emplace_back(0.27972043387654);
            strikes.emplace_back(0.257);   stdDevsOnExpiry.emplace_back(0.279533328503776);   stdDevsOnPayment.emplace_back(0.280366887761207);
            strikes.emplace_back(0.258);   stdDevsOnExpiry.emplace_back(0.280163886669214);   stdDevsOnPayment.emplace_back(0.281011073386628);
            strikes.emplace_back(0.259);   stdDevsOnExpiry.emplace_back(0.280801401845504);   stdDevsOnPayment.emplace_back(0.281636464864025);
            strikes.emplace_back(0.26);    stdDevsOnExpiry.emplace_back(0.281444609121582);   stdDevsOnPayment.emplace_back(0.282278058193167);
            strikes.emplace_back(0.261);   stdDevsOnExpiry.emplace_back(0.282062518176379);   stdDevsOnPayment.emplace_back(0.2829177073001);
            strikes.emplace_back(0.262);   stdDevsOnExpiry.emplace_back(0.282703195630329);   stdDevsOnPayment.emplace_back(0.28354990355523);
            strikes.emplace_back(0.263);   stdDevsOnExpiry.emplace_back(0.283321420912892);   stdDevsOnPayment.emplace_back(0.284173026773382);
            strikes.emplace_back(0.264);   stdDevsOnExpiry.emplace_back(0.283962098366842);   stdDevsOnPayment.emplace_back(0.284814296065489);
            strikes.emplace_back(0.265);   stdDevsOnExpiry.emplace_back(0.284589810482385);   stdDevsOnPayment.emplace_back(0.285413116506022);
            strikes.emplace_back(0.266);   stdDevsOnExpiry.emplace_back(0.285200130070798);   stdDevsOnPayment.emplace_back(0.286036563761209);
            strikes.emplace_back(0.267);   stdDevsOnExpiry.emplace_back(0.285821833858787);   stdDevsOnPayment.emplace_back(0.286665843683024);
            strikes.emplace_back(0.268);   stdDevsOnExpiry.emplace_back(0.286418239425495);   stdDevsOnPayment.emplace_back(0.287277625604954);
            strikes.emplace_back(0.269);   stdDevsOnExpiry.emplace_back(0.287055438374019);   stdDevsOnPayment.emplace_back(0.28788454697136);
            strikes.emplace_back(0.27);    stdDevsOnExpiry.emplace_back(0.287650895257428);   stdDevsOnPayment.emplace_back(0.288504105782127);
            strikes.emplace_back(0.271);   stdDevsOnExpiry.emplace_back(0.288259001251479);   stdDevsOnPayment.emplace_back(0.289124960741035);
            strikes.emplace_back(0.272);   stdDevsOnExpiry.emplace_back(0.288866474789997);   stdDevsOnPayment.emplace_back(0.289720216774184);
            strikes.emplace_back(0.273);   stdDevsOnExpiry.emplace_back(0.289457504484683);   stdDevsOnPayment.emplace_back(0.29032389777024);
            strikes.emplace_back(0.274);   stdDevsOnExpiry.emplace_back(0.290065294250967);   stdDevsOnPayment.emplace_back(0.290925958581123);
            strikes.emplace_back(0.275);   stdDevsOnExpiry.emplace_back(0.290661383589909);   stdDevsOnPayment.emplace_back(0.291521538651306);
            strikes.emplace_back(0.276);   stdDevsOnExpiry.emplace_back(0.291270122039491);   stdDevsOnPayment.emplace_back(0.292117442758525);
            strikes.emplace_back(0.277);   stdDevsOnExpiry.emplace_back(0.291862732873007);   stdDevsOnPayment.emplace_back(0.292714643013883);
            strikes.emplace_back(0.278);   stdDevsOnExpiry.emplace_back(0.292448070467904);   stdDevsOnPayment.emplace_back(0.293307306750752);
            strikes.emplace_back(0.279);   stdDevsOnExpiry.emplace_back(0.293031826923971);   stdDevsOnPayment.emplace_back(0.293893489746923);
            strikes.emplace_back(0.28);    stdDevsOnExpiry.emplace_back(0.293630129857275);   stdDevsOnPayment.emplace_back(0.294490365965247);
            strikes.emplace_back(0.281);   stdDevsOnExpiry.emplace_back(0.294205348163659);   stdDevsOnPayment.emplace_back(0.295054514443043);
            strikes.emplace_back(0.282);   stdDevsOnExpiry.emplace_back(0.294780250242278);   stdDevsOnPayment.emplace_back(0.295640697439214);
            strikes.emplace_back(0.283);   stdDevsOnExpiry.emplace_back(0.295369066342601);   stdDevsOnPayment.emplace_back(0.296229796768699);
            strikes.emplace_back(0.284);   stdDevsOnExpiry.emplace_back(0.295941122371326);   stdDevsOnPayment.emplace_back(0.296821812431499);
            strikes.emplace_back(0.285);   stdDevsOnExpiry.emplace_back(0.296506537616964);   stdDevsOnPayment.emplace_back(0.297392117612959);
            strikes.emplace_back(0.286);   stdDevsOnExpiry.emplace_back(0.297095986172819);   stdDevsOnPayment.emplace_back(0.297958858387035);
            strikes.emplace_back(0.287);   stdDevsOnExpiry.emplace_back(0.297663931240585);   stdDevsOnPayment.emplace_back(0.298539532753612);
            strikes.emplace_back(0.288);   stdDevsOnExpiry.emplace_back(0.298226500436329);   stdDevsOnPayment.emplace_back(0.299106597564723);
            strikes.emplace_back(0.289);   stdDevsOnExpiry.emplace_back(0.298813735397823);   stdDevsOnPayment.emplace_back(0.299686299820195);
            strikes.emplace_back(0.29);    stdDevsOnExpiry.emplace_back(0.299351322600051);   stdDevsOnPayment.emplace_back(0.3002384589277);
            strikes.emplace_back(0.291);   stdDevsOnExpiry.emplace_back(0.299931284322926);   stdDevsOnPayment.emplace_back(0.300800663183287);
            strikes.emplace_back(0.292);   stdDevsOnExpiry.emplace_back(0.300491639924308);   stdDevsOnPayment.emplace_back(0.301358978994454);
            strikes.emplace_back(0.293);   stdDevsOnExpiry.emplace_back(0.30104408983154);    stdDevsOnPayment.emplace_back(0.30192507169446);
            strikes.emplace_back(0.294);   stdDevsOnExpiry.emplace_back(0.30161646208803);    stdDevsOnPayment.emplace_back(0.30247496254272);
            strikes.emplace_back(0.295);   stdDevsOnExpiry.emplace_back(0.302157527795685);   stdDevsOnPayment.emplace_back(0.303038462946447);
            strikes.emplace_back(0.296);   stdDevsOnExpiry.emplace_back(0.302717883397067);   stdDevsOnPayment.emplace_back(0.303587057646567);
            strikes.emplace_back(0.297);   stdDevsOnExpiry.emplace_back(0.303255786827061);   stdDevsOnPayment.emplace_back(0.304149261902154);
            strikes.emplace_back(0.298);   stdDevsOnExpiry.emplace_back(0.303781989829713);   stdDevsOnPayment.emplace_back(0.304700124861519);
            strikes.emplace_back(0.299);   stdDevsOnExpiry.emplace_back(0.304330645003752);   stdDevsOnPayment.emplace_back(0.305239646524661);
            strikes.emplace_back(0.3);     stdDevsOnExpiry.emplace_back(0.304881197544388);   stdDevsOnPayment.emplace_back(0.305780140298908);
            strikes.emplace_back(0.301);   stdDevsOnExpiry.emplace_back(0.305453569800878);   stdDevsOnPayment.emplace_back(0.306327114813854);
            strikes.emplace_back(0.302);   stdDevsOnExpiry.emplace_back(0.305970602198316);   stdDevsOnPayment.emplace_back(0.306869876847346);
            strikes.emplace_back(0.303);   stdDevsOnExpiry.emplace_back(0.306495540289904);   stdDevsOnPayment.emplace_back(0.307408102362348);
            strikes.emplace_back(0.304);   stdDevsOnExpiry.emplace_back(0.307054947207988);   stdDevsOnPayment.emplace_back(0.307955400914329);
            strikes.emplace_back(0.305);   stdDevsOnExpiry.emplace_back(0.307578936616277);   stdDevsOnPayment.emplace_back(0.308478072651655);
            strikes.emplace_back(0.306);   stdDevsOnExpiry.emplace_back(0.308105139618929);   stdDevsOnPayment.emplace_back(0.30902828753695);
            strikes.emplace_back(0.307);   stdDevsOnExpiry.emplace_back(0.308639564543498);   stdDevsOnPayment.emplace_back(0.309550959274277);
            strikes.emplace_back(0.308);   stdDevsOnExpiry.emplace_back(0.30916892982381);    stdDevsOnPayment.emplace_back(0.310072982937534);
            strikes.emplace_back(0.309);   stdDevsOnExpiry.emplace_back(0.309708414392635);   stdDevsOnPayment.emplace_back(0.310612504600676);
            strikes.emplace_back(0.31);    stdDevsOnExpiry.emplace_back(0.310226711701136);   stdDevsOnPayment.emplace_back(0.311149433967539);
            strikes.emplace_back(0.311);   stdDevsOnExpiry.emplace_back(0.310757025664747);   stdDevsOnPayment.emplace_back(0.311663356704923);
            strikes.emplace_back(0.312);   stdDevsOnExpiry.emplace_back(0.311267733506864);   stdDevsOnPayment.emplace_back(0.312172742923818);
            strikes.emplace_back(0.313);   stdDevsOnExpiry.emplace_back(0.31179520142058);    stdDevsOnPayment.emplace_back(0.312693470438935);
            strikes.emplace_back(0.314);   stdDevsOnExpiry.emplace_back(0.312293892607588);   stdDevsOnPayment.emplace_back(0.313225215213239);
            strikes.emplace_back(0.315);   stdDevsOnExpiry.emplace_back(0.312807130271834);   stdDevsOnPayment.emplace_back(0.313724880321086);
            strikes.emplace_back(0.316);   stdDevsOnExpiry.emplace_back(0.313339025374274);   stdDevsOnPayment.emplace_back(0.314244311688064);
            strikes.emplace_back(0.317);   stdDevsOnExpiry.emplace_back(0.313838981472347);   stdDevsOnPayment.emplace_back(0.31478739775859);
            strikes.emplace_back(0.318);   stdDevsOnExpiry.emplace_back(0.314357595008614);   stdDevsOnPayment.emplace_back(0.315282202310914);
            strikes.emplace_back(0.319);   stdDevsOnExpiry.emplace_back(0.314857234878921);   stdDevsOnPayment.emplace_back(0.315793208714983);
            strikes.emplace_back(0.32);    stdDevsOnExpiry.emplace_back(0.315366361582208);   stdDevsOnPayment.emplace_back(0.316289633452481);
            strikes.emplace_back(0.321);   stdDevsOnExpiry.emplace_back(0.315885291346242);   stdDevsOnPayment.emplace_back(0.31680258407876);
            strikes.emplace_back(0.322);   stdDevsOnExpiry.emplace_back(0.316385879899846);   stdDevsOnPayment.emplace_back(0.317304193408817);
            strikes.emplace_back(0.323);   stdDevsOnExpiry.emplace_back(0.316888682047813);   stdDevsOnPayment.emplace_back(0.317796729701896);
            strikes.emplace_back(0.324);   stdDevsOnExpiry.emplace_back(0.317366185974499);   stdDevsOnPayment.emplace_back(0.318296394809743);
            strikes.emplace_back(0.325);   stdDevsOnExpiry.emplace_back(0.317897448621407);   stdDevsOnPayment.emplace_back(0.318796059917591);
            strikes.emplace_back(0.326);   stdDevsOnExpiry.emplace_back(0.318374952548092);   stdDevsOnPayment.emplace_back(0.319321972025266);
            strikes.emplace_back(0.327);   stdDevsOnExpiry.emplace_back(0.318880916973719);   stdDevsOnPayment.emplace_back(0.319809971799856);
            strikes.emplace_back(0.328);   stdDevsOnExpiry.emplace_back(0.319367907733385);   stdDevsOnPayment.emplace_back(0.320308988833634);
            strikes.emplace_back(0.329);   stdDevsOnExpiry.emplace_back(0.319854898493051);   stdDevsOnPayment.emplace_back(0.320818051015494);
            strikes.emplace_back(0.33);    stdDevsOnExpiry.emplace_back(0.320354538363358);   stdDevsOnPayment.emplace_back(0.321299570049386);
            strikes.emplace_back(0.331);   stdDevsOnExpiry.emplace_back(0.320847853678344);   stdDevsOnPayment.emplace_back(0.321782061194382);
            strikes.emplace_back(0.332);   stdDevsOnExpiry.emplace_back(0.321319033049709);   stdDevsOnPayment.emplace_back(0.322270709043042);
            strikes.emplace_back(0.333);   stdDevsOnExpiry.emplace_back(0.321799699254055);   stdDevsOnPayment.emplace_back(0.322769726076819);
            strikes.emplace_back(0.334);   stdDevsOnExpiry.emplace_back(0.322302501402021);   stdDevsOnPayment.emplace_back(0.323246384555187);
            strikes.emplace_back(0.335);   stdDevsOnExpiry.emplace_back(0.322783167606367);   stdDevsOnPayment.emplace_back(0.323742809292685);
            strikes.emplace_back(0.336);   stdDevsOnExpiry.emplace_back(0.323279645199013);   stdDevsOnPayment.emplace_back(0.324182851586107);
            strikes.emplace_back(0.337);   stdDevsOnExpiry.emplace_back(0.323735013182078);   stdDevsOnPayment.emplace_back(0.324681868619885);
            strikes.emplace_back(0.338);   stdDevsOnExpiry.emplace_back(0.324212517108763);   stdDevsOnPayment.emplace_back(0.325174404912964);
            strikes.emplace_back(0.339);   stdDevsOnExpiry.emplace_back(0.324693183313109);   stdDevsOnPayment.emplace_back(0.325647498983947);
            strikes.emplace_back(0.34);    stdDevsOnExpiry.emplace_back(0.325170687239794);   stdDevsOnPayment.emplace_back(0.326114112314233);
            strikes.emplace_back(0.341);   stdDevsOnExpiry.emplace_back(0.325638704333499);   stdDevsOnPayment.emplace_back(0.326590446755566);
            strikes.emplace_back(0.342);   stdDevsOnExpiry.emplace_back(0.326132019648485);   stdDevsOnPayment.emplace_back(0.327066781196899);
            strikes.emplace_back(0.343);   stdDevsOnExpiry.emplace_back(0.326571576243249);   stdDevsOnPayment.emplace_back(0.327536634897533);
            strikes.emplace_back(0.344);   stdDevsOnExpiry.emplace_back(0.327049080169934);   stdDevsOnPayment.emplace_back(0.328019450079565);
            strikes.emplace_back(0.345);   stdDevsOnExpiry.emplace_back(0.327507610430659);   stdDevsOnPayment.emplace_back(0.32848606340985);
            strikes.emplace_back(0.346);   stdDevsOnExpiry.emplace_back(0.327997763467985);   stdDevsOnPayment.emplace_back(0.328949436369786);
            strikes.emplace_back(0.347);   stdDevsOnExpiry.emplace_back(0.328443644618068);   stdDevsOnPayment.emplace_back(0.329409568959373);
            strikes.emplace_back(0.348);   stdDevsOnExpiry.emplace_back(0.328933797655394);   stdDevsOnPayment.emplace_back(0.329876182289659);
            strikes.emplace_back(0.349);   stdDevsOnExpiry.emplace_back(0.329363867417177);   stdDevsOnPayment.emplace_back(0.330346035990293);
            strikes.emplace_back(0.35);    stdDevsOnExpiry.emplace_back(0.329841371343863);   stdDevsOnPayment.emplace_back(0.330799687839182);
            strikes.emplace_back(0.351);   stdDevsOnExpiry.emplace_back(0.330284090216286);   stdDevsOnPayment.emplace_back(0.331279262650864);
            strikes.emplace_back(0.352);   stdDevsOnExpiry.emplace_back(0.330755269587652);   stdDevsOnPayment.emplace_back(0.331723193388705);
            strikes.emplace_back(0.353);   stdDevsOnExpiry.emplace_back(0.331204313015395);   stdDevsOnPayment.emplace_back(0.332163883756196);
            strikes.emplace_back(0.354);   stdDevsOnExpiry.emplace_back(0.331647031887819);   stdDevsOnPayment.emplace_back(0.332624016345783);
            strikes.emplace_back(0.355);   stdDevsOnExpiry.emplace_back(0.332140347202805);   stdDevsOnPayment.emplace_back(0.333093870046418);
            strikes.emplace_back(0.356);   stdDevsOnExpiry.emplace_back(0.332579903797569);   stdDevsOnPayment.emplace_back(0.333560483376703);
            strikes.emplace_back(0.357);   stdDevsOnExpiry.emplace_back(0.333022622669992);   stdDevsOnPayment.emplace_back(0.334001173744195);
            strikes.emplace_back(0.358);   stdDevsOnExpiry.emplace_back(0.333474828375396);   stdDevsOnPayment.emplace_back(0.334464546704131);
            strikes.emplace_back(0.359);   stdDevsOnExpiry.emplace_back(0.33392387180314);    stdDevsOnPayment.emplace_back(0.334889035219877);
            strikes.emplace_back(0.36);    stdDevsOnExpiry.emplace_back(0.334341292454282);   stdDevsOnPayment.emplace_back(0.335342687068766);
            strikes.emplace_back(0.361);   stdDevsOnExpiry.emplace_back(0.334799822715007);   stdDevsOnPayment.emplace_back(0.335793098547305);
            strikes.emplace_back(0.362);   stdDevsOnExpiry.emplace_back(0.335271002086372);   stdDevsOnPayment.emplace_back(0.336246750396193);
            strikes.emplace_back(0.363);   stdDevsOnExpiry.emplace_back(0.335701071848155);   stdDevsOnPayment.emplace_back(0.33666799854159);
            strikes.emplace_back(0.364);   stdDevsOnExpiry.emplace_back(0.336159602108879);   stdDevsOnPayment.emplace_back(0.337144332982923);
            strikes.emplace_back(0.365);   stdDevsOnExpiry.emplace_back(0.336583347315342);   stdDevsOnPayment.emplace_back(0.337555860017272);
            strikes.emplace_back(0.366);   stdDevsOnExpiry.emplace_back(0.336994443411164);   stdDevsOnPayment.emplace_back(0.337993310014414);
            strikes.emplace_back(0.367);   stdDevsOnExpiry.emplace_back(0.337459298227208);   stdDevsOnPayment.emplace_back(0.338440481122604);
            strikes.emplace_back(0.368);   stdDevsOnExpiry.emplace_back(0.337892530266652);   stdDevsOnPayment.emplace_back(0.33886496963835);
            strikes.emplace_back(0.369);   stdDevsOnExpiry.emplace_back(0.338335249139075);   stdDevsOnPayment.emplace_back(0.339331582968636);
            strikes.emplace_back(0.37);    stdDevsOnExpiry.emplace_back(0.338768481178518);   stdDevsOnPayment.emplace_back(0.339749590743683);
            strikes.emplace_back(0.371);   stdDevsOnExpiry.emplace_back(0.339163765886039);   stdDevsOnPayment.emplace_back(0.340180560000127);
            strikes.emplace_back(0.372);   stdDevsOnExpiry.emplace_back(0.339634945257404);   stdDevsOnPayment.emplace_back(0.340601808145523);
            strikes.emplace_back(0.373);   stdDevsOnExpiry.emplace_back(0.340039716797906);   stdDevsOnPayment.emplace_back(0.341058700364761);
            strikes.emplace_back(0.374);   stdDevsOnExpiry.emplace_back(0.34048559794799);    stdDevsOnPayment.emplace_back(0.341473467769459);
            strikes.emplace_back(0.375);   stdDevsOnExpiry.emplace_back(0.340896694043811);   stdDevsOnPayment.emplace_back(0.341897956285205);
            strikes.emplace_back(0.376);   stdDevsOnExpiry.emplace_back(0.341336250638575);   stdDevsOnPayment.emplace_back(0.342348367763744);
            strikes.emplace_back(0.377);   stdDevsOnExpiry.emplace_back(0.341753671289717);   stdDevsOnPayment.emplace_back(0.34277285627949);
            strikes.emplace_back(0.378);   stdDevsOnExpiry.emplace_back(0.34217741649618);    stdDevsOnPayment.emplace_back(0.343174662202791);
            strikes.emplace_back(0.379);   stdDevsOnExpiry.emplace_back(0.342604323980302);   stdDevsOnPayment.emplace_back(0.343608871829585);
            strikes.emplace_back(0.38);    stdDevsOnExpiry.emplace_back(0.343021744631445);   stdDevsOnPayment.emplace_back(0.344000956641838);
            strikes.emplace_back(0.381);   stdDevsOnExpiry.emplace_back(0.343404380228325);   stdDevsOnPayment.emplace_back(0.344428685527933);
            strikes.emplace_back(0.382);   stdDevsOnExpiry.emplace_back(0.343843936823088);   stdDevsOnPayment.emplace_back(0.344875856636123);
            strikes.emplace_back(0.383);   stdDevsOnExpiry.emplace_back(0.344277168862531);   stdDevsOnPayment.emplace_back(0.34525822033733);
            strikes.emplace_back(0.384);   stdDevsOnExpiry.emplace_back(0.344700914068994);   stdDevsOnPayment.emplace_back(0.34570215107517);
            strikes.emplace_back(0.385);   stdDevsOnExpiry.emplace_back(0.345121496997796);   stdDevsOnPayment.emplace_back(0.346139601072313);
            strikes.emplace_back(0.386);   stdDevsOnExpiry.emplace_back(0.345538917648939);   stdDevsOnPayment.emplace_back(0.346502522551424);
            strikes.emplace_back(0.387);   stdDevsOnExpiry.emplace_back(0.345912066412839);   stdDevsOnPayment.emplace_back(0.346956174400312);
            strikes.emplace_back(0.388);   stdDevsOnExpiry.emplace_back(0.346354785285262);   stdDevsOnPayment.emplace_back(0.34737418217536);
            strikes.emplace_back(0.389);   stdDevsOnExpiry.emplace_back(0.346746907715123);   stdDevsOnPayment.emplace_back(0.347756545876566);
            strikes.emplace_back(0.39);    stdDevsOnExpiry.emplace_back(0.347151679255625);   stdDevsOnPayment.emplace_back(0.348171313281264);
            strikes.emplace_back(0.391);   stdDevsOnExpiry.emplace_back(0.347562775351446);   stdDevsOnPayment.emplace_back(0.348576359574914);
            strikes.emplace_back(0.392);   stdDevsOnExpiry.emplace_back(0.347986520557909);   stdDevsOnPayment.emplace_back(0.348968444387168);
            strikes.emplace_back(0.393);   stdDevsOnExpiry.emplace_back(0.348407103486711);   stdDevsOnPayment.emplace_back(0.349415615495358);
            strikes.emplace_back(0.394);   stdDevsOnExpiry.emplace_back(0.348773927695291);   stdDevsOnPayment.emplace_back(0.349820661789009);
            strikes.emplace_back(0.395);   stdDevsOnExpiry.emplace_back(0.349203997457074);   stdDevsOnPayment.emplace_back(0.350232188823358);
            strikes.emplace_back(0.396);   stdDevsOnExpiry.emplace_back(0.349624580385876);   stdDevsOnPayment.emplace_back(0.350637235117008);
            strikes.emplace_back(0.397);   stdDevsOnExpiry.emplace_back(0.349981917761475);   stdDevsOnPayment.emplace_back(0.351029319929262);
            strikes.emplace_back(0.398);   stdDevsOnExpiry.emplace_back(0.350396176134957);   stdDevsOnPayment.emplace_back(0.351434366222912);
            strikes.emplace_back(0.399);   stdDevsOnExpiry.emplace_back(0.350797785397799);   stdDevsOnPayment.emplace_back(0.351803768442721);
            strikes.emplace_back(0.4);     stdDevsOnExpiry.emplace_back(0.351174096439359);   stdDevsOnPayment.emplace_back(0.352250939550912);
            strikes.emplace_back(0.401);   stdDevsOnExpiry.emplace_back(0.351597841645821);   stdDevsOnPayment.emplace_back(0.352620341770721);
            strikes.emplace_back(0.402);   stdDevsOnExpiry.emplace_back(0.351986801798022);   stdDevsOnPayment.emplace_back(0.35303186880507);
            strikes.emplace_back(0.403);   stdDevsOnExpiry.emplace_back(0.352356788284262);   stdDevsOnPayment.emplace_back(0.353407751765577);
            strikes.emplace_back(0.404);   stdDevsOnExpiry.emplace_back(0.352793182601365);   stdDevsOnPayment.emplace_back(0.353799836577831);
            strikes.emplace_back(0.405);   stdDevsOnExpiry.emplace_back(0.353194791864206);   stdDevsOnPayment.emplace_back(0.354208123241831);
            strikes.emplace_back(0.406);   stdDevsOnExpiry.emplace_back(0.353548966962145);   stdDevsOnPayment.emplace_back(0.354587246572688);
            strikes.emplace_back(0.407);   stdDevsOnExpiry.emplace_back(0.353953738502647);   stdDevsOnPayment.emplace_back(0.354995533236687);
            strikes.emplace_back(0.408);   stdDevsOnExpiry.emplace_back(0.354333211821867);   stdDevsOnPayment.emplace_back(0.355374656567544);
            strikes.emplace_back(0.409);   stdDevsOnExpiry.emplace_back(0.354725334251728);   stdDevsOnPayment.emplace_back(0.355776462490846);
            strikes.emplace_back(0.41);    stdDevsOnExpiry.emplace_back(0.355104807570948);   stdDevsOnPayment.emplace_back(0.356129662858909);
            strikes.emplace_back(0.411);   stdDevsOnExpiry.emplace_back(0.355503254556129);   stdDevsOnPayment.emplace_back(0.356554151374654);
            strikes.emplace_back(0.412);   stdDevsOnExpiry.emplace_back(0.355863754209388);   stdDevsOnPayment.emplace_back(0.356910592113067);
            strikes.emplace_back(0.413);   stdDevsOnExpiry.emplace_back(0.35627168802755);    stdDevsOnPayment.emplace_back(0.357279994332876);
            strikes.emplace_back(0.414);   stdDevsOnExpiry.emplace_back(0.356660648179751);   stdDevsOnPayment.emplace_back(0.357691521367225);
            strikes.emplace_back(0.415);   stdDevsOnExpiry.emplace_back(0.357036959221311);   stdDevsOnPayment.emplace_back(0.358080365809129);
            strikes.emplace_back(0.416);   stdDevsOnExpiry.emplace_back(0.35740062115223);    stdDevsOnPayment.emplace_back(0.358472450621383);
            strikes.emplace_back(0.417);   stdDevsOnExpiry.emplace_back(0.357795905859751);   stdDevsOnPayment.emplace_back(0.358838612470843);
            strikes.emplace_back(0.418);   stdDevsOnExpiry.emplace_back(0.35815008095769);    stdDevsOnPayment.emplace_back(0.359211255061002);
            strikes.emplace_back(0.419);   stdDevsOnExpiry.emplace_back(0.358561177053512);   stdDevsOnPayment.emplace_back(0.359635743576747);
            strikes.emplace_back(0.42);    stdDevsOnExpiry.emplace_back(0.358915352151451);   stdDevsOnPayment.emplace_back(0.35999218431516);
            strikes.emplace_back(0.421);   stdDevsOnExpiry.emplace_back(0.359313799136632);   stdDevsOnPayment.emplace_back(0.360361586534969);
            strikes.emplace_back(0.422);   stdDevsOnExpiry.emplace_back(0.359674298789891);   stdDevsOnPayment.emplace_back(0.360756911717572);
            strikes.emplace_back(0.423);   stdDevsOnExpiry.emplace_back(0.36002847388783);    stdDevsOnPayment.emplace_back(0.361116592826334);
            strikes.emplace_back(0.424);   stdDevsOnExpiry.emplace_back(0.360455381371953);   stdDevsOnPayment.emplace_back(0.361502196897889);
            strikes.emplace_back(0.425);   stdDevsOnExpiry.emplace_back(0.360809556469892);   stdDevsOnPayment.emplace_back(0.361891041339793);
            strikes.emplace_back(0.426);   stdDevsOnExpiry.emplace_back(0.361141595624209);   stdDevsOnPayment.emplace_back(0.362260443559603);
            strikes.emplace_back(0.427);   stdDevsOnExpiry.emplace_back(0.361562178553012);   stdDevsOnPayment.emplace_back(0.362594201705571);
            strikes.emplace_back(0.428);   stdDevsOnExpiry.emplace_back(0.361922678206271);   stdDevsOnPayment.emplace_back(0.363021930591666);
            strikes.emplace_back(0.429);   stdDevsOnExpiry.emplace_back(0.36227369102655);    stdDevsOnPayment.emplace_back(0.363358929107983);
            strikes.emplace_back(0.43);    stdDevsOnExpiry.emplace_back(0.362637352957469);   stdDevsOnPayment.emplace_back(0.363741292809189);
            strikes.emplace_back(0.431);   stdDevsOnExpiry.emplace_back(0.362994690333068);   stdDevsOnPayment.emplace_back(0.364075050955157);
            strikes.emplace_back(0.432);   stdDevsOnExpiry.emplace_back(0.363342540875687);   stdDevsOnPayment.emplace_back(0.364473616508109);
            strikes.emplace_back(0.433);   stdDevsOnExpiry.emplace_back(0.363737825583208);   stdDevsOnPayment.emplace_back(0.364820336135474);
            strikes.emplace_back(0.434);   stdDevsOnExpiry.emplace_back(0.364095162958807);   stdDevsOnPayment.emplace_back(0.365186497984934);
            strikes.emplace_back(0.435);   stdDevsOnExpiry.emplace_back(0.364455662612066);   stdDevsOnPayment.emplace_back(0.365562380945441);
            strikes.emplace_back(0.436);   stdDevsOnExpiry.emplace_back(0.364866758707888);   stdDevsOnPayment.emplace_back(0.365889658350711);
            strikes.emplace_back(0.437);   stdDevsOnExpiry.emplace_back(0.365208284695186);   stdDevsOnPayment.emplace_back(0.366291464274012);
            strikes.emplace_back(0.438);   stdDevsOnExpiry.emplace_back(0.365584595736746);   stdDevsOnPayment.emplace_back(0.366657626123472);
            strikes.emplace_back(0.439);   stdDevsOnExpiry.emplace_back(0.365907148058083);   stdDevsOnPayment.emplace_back(0.367001105380488);
            strikes.emplace_back(0.44);    stdDevsOnExpiry.emplace_back(0.366270809989003);   stdDevsOnPayment.emplace_back(0.367380228711345);
            strikes.emplace_back(0.441);   stdDevsOnExpiry.emplace_back(0.366612335976301);   stdDevsOnPayment.emplace_back(0.367710746486964);
            strikes.emplace_back(0.442);   stdDevsOnExpiry.emplace_back(0.366998133850841);   stdDevsOnPayment.emplace_back(0.368083389077122);
            strikes.emplace_back(0.443);   stdDevsOnExpiry.emplace_back(0.367333335282819);   stdDevsOnPayment.emplace_back(0.368452791296931);
            strikes.emplace_back(0.444);   stdDevsOnExpiry.emplace_back(0.367693834936078);   stdDevsOnPayment.emplace_back(0.368793030183598);
            strikes.emplace_back(0.445);   stdDevsOnExpiry.emplace_back(0.368073308255299);   stdDevsOnPayment.emplace_back(0.369162432403407);
            strikes.emplace_back(0.446);   stdDevsOnExpiry.emplace_back(0.368440132463878);   stdDevsOnPayment.emplace_back(0.369492950179026);
            strikes.emplace_back(0.447);   stdDevsOnExpiry.emplace_back(0.368765847062875);   stdDevsOnPayment.emplace_back(0.369878554250581);
            strikes.emplace_back(0.448);   stdDevsOnExpiry.emplace_back(0.369104210772513);   stdDevsOnPayment.emplace_back(0.370196110544803);
            strikes.emplace_back(0.449);   stdDevsOnExpiry.emplace_back(0.369490008647054);   stdDevsOnPayment.emplace_back(0.370559032023914);
            strikes.emplace_back(0.45);    stdDevsOnExpiry.emplace_back(0.369825210079032);   stdDevsOnPayment.emplace_back(0.370896030540231);
            strikes.emplace_back(0.451);   stdDevsOnExpiry.emplace_back(0.370204683398252);   stdDevsOnPayment.emplace_back(0.371278394241437);
            strikes.emplace_back(0.452);   stdDevsOnExpiry.emplace_back(0.370514586608948);   stdDevsOnPayment.emplace_back(0.37159271016531);
            strikes.emplace_back(0.453);   stdDevsOnExpiry.emplace_back(0.370856112596247);   stdDevsOnPayment.emplace_back(0.371968593125818);
            strikes.emplace_back(0.454);   stdDevsOnExpiry.emplace_back(0.371197638583545);   stdDevsOnPayment.emplace_back(0.37232503386423);
            strikes.emplace_back(0.455);   stdDevsOnExpiry.emplace_back(0.371586598735746);   stdDevsOnPayment.emplace_back(0.372665272750896);
            strikes.emplace_back(0.456);   stdDevsOnExpiry.emplace_back(0.371874366002821);   stdDevsOnPayment.emplace_back(0.373002271267214);
            strikes.emplace_back(0.457);   stdDevsOnExpiry.emplace_back(0.37222854110076);    stdDevsOnPayment.emplace_back(0.373358712005626);
            strikes.emplace_back(0.458);   stdDevsOnExpiry.emplace_back(0.372579553921038);   stdDevsOnPayment.emplace_back(0.373685989410896);
            strikes.emplace_back(0.459);   stdDevsOnExpiry.emplace_back(0.372917917630676);   stdDevsOnPayment.emplace_back(0.374045670519657);
            strikes.emplace_back(0.46);    stdDevsOnExpiry.emplace_back(0.373332176004159);   stdDevsOnPayment.emplace_back(0.374382669035975);
            strikes.emplace_back(0.461);   stdDevsOnExpiry.emplace_back(0.373597807327613);   stdDevsOnPayment.emplace_back(0.374748830885435);
            strikes.emplace_back(0.462);   stdDevsOnExpiry.emplace_back(0.373977280646833);   stdDevsOnPayment.emplace_back(0.375095550512799);
            strikes.emplace_back(0.463);   stdDevsOnExpiry.emplace_back(0.374287183857529);   stdDevsOnPayment.emplace_back(0.375406626066323);
            strikes.emplace_back(0.464);   stdDevsOnExpiry.emplace_back(0.374616060734187);   stdDevsOnPayment.emplace_back(0.375727422730894);
            strikes.emplace_back(0.465);   stdDevsOnExpiry.emplace_back(0.375001858608727);   stdDevsOnPayment.emplace_back(0.376096824950703);
            strikes.emplace_back(0.466);   stdDevsOnExpiry.emplace_back(0.375311761819424);   stdDevsOnPayment.emplace_back(0.376420861985624);
            strikes.emplace_back(0.467);   stdDevsOnExpiry.emplace_back(0.375599529086499);   stdDevsOnPayment.emplace_back(0.376754620131592);
            strikes.emplace_back(0.468);   stdDevsOnExpiry.emplace_back(0.376004300627001);   stdDevsOnPayment.emplace_back(0.377098099388607);
            strikes.emplace_back(0.469);   stdDevsOnExpiry.emplace_back(0.376292067894076);   stdDevsOnPayment.emplace_back(0.377425376793877);
            strikes.emplace_back(0.47);    stdDevsOnExpiry.emplace_back(0.376646242992015);   stdDevsOnPayment.emplace_back(0.377755894569496);
            strikes.emplace_back(0.471);   stdDevsOnExpiry.emplace_back(0.376952983925051);   stdDevsOnPayment.emplace_back(0.378079931604416);
            strikes.emplace_back(0.472);   stdDevsOnExpiry.emplace_back(0.377332457244272);   stdDevsOnPayment.emplace_back(0.378449333824225);
            strikes.emplace_back(0.473);   stdDevsOnExpiry.emplace_back(0.377661334120929);   stdDevsOnPayment.emplace_back(0.378776611229495);
            strikes.emplace_back(0.474);   stdDevsOnExpiry.emplace_back(0.377999697830567);   stdDevsOnPayment.emplace_back(0.379113609745812);
            strikes.emplace_back(0.475);   stdDevsOnExpiry.emplace_back(0.378277978264662);   stdDevsOnPayment.emplace_back(0.379463569743526);
            strikes.emplace_back(0.476);   stdDevsOnExpiry.emplace_back(0.378635315640261);   stdDevsOnPayment.emplace_back(0.379761683815653);
            strikes.emplace_back(0.477);   stdDevsOnExpiry.emplace_back(0.378948381128618);   stdDevsOnPayment.emplace_back(0.380095441961621);
            strikes.emplace_back(0.478);   stdDevsOnExpiry.emplace_back(0.379264608894634);   stdDevsOnPayment.emplace_back(0.380471324922129);
            strikes.emplace_back(0.479);   stdDevsOnExpiry.emplace_back(0.379574512105331);   stdDevsOnPayment.emplace_back(0.380779160105303);
            strikes.emplace_back(0.48);    stdDevsOnExpiry.emplace_back(0.379900226704328);   stdDevsOnPayment.emplace_back(0.381122639362319);
            strikes.emplace_back(0.481);   stdDevsOnExpiry.emplace_back(0.380276537745888);   stdDevsOnPayment.emplace_back(0.381414272693747);
            strikes.emplace_back(0.482);   stdDevsOnExpiry.emplace_back(0.380605414622546);   stdDevsOnPayment.emplace_back(0.381773953802509);
            strikes.emplace_back(0.483);   stdDevsOnExpiry.emplace_back(0.380896344167281);   stdDevsOnPayment.emplace_back(0.382068827504286);
            strikes.emplace_back(0.484);   stdDevsOnExpiry.emplace_back(0.381234707876919);   stdDevsOnPayment.emplace_back(0.382405826020603);
            strikes.emplace_back(0.485);   stdDevsOnExpiry.emplace_back(0.381595207530179);   stdDevsOnPayment.emplace_back(0.382716901574127);
            strikes.emplace_back(0.486);   stdDevsOnExpiry.emplace_back(0.381835540632351);   stdDevsOnPayment.emplace_back(0.383031217498);
            strikes.emplace_back(0.487);   stdDevsOnExpiry.emplace_back(0.382139119287727);   stdDevsOnPayment.emplace_back(0.38331637008873);
            strikes.emplace_back(0.488);   stdDevsOnExpiry.emplace_back(0.382486969830346);   stdDevsOnPayment.emplace_back(0.383685772308539);
            strikes.emplace_back(0.489);   stdDevsOnExpiry.emplace_back(0.382847469483605);   stdDevsOnPayment.emplace_back(0.383964444158571);
            strikes.emplace_back(0.49);    stdDevsOnExpiry.emplace_back(0.38313207447302);    stdDevsOnPayment.emplace_back(0.384301442674888);
            strikes.emplace_back(0.491);   stdDevsOnExpiry.emplace_back(0.383479925015639);   stdDevsOnPayment.emplace_back(0.384661123783649);
            strikes.emplace_back(0.492);   stdDevsOnExpiry.emplace_back(0.383796152781656);   stdDevsOnPayment.emplace_back(0.384917113041236);
            strikes.emplace_back(0.493);   stdDevsOnExpiry.emplace_back(0.384106055992352);   stdDevsOnPayment.emplace_back(0.385280034520347);
            strikes.emplace_back(0.494);   stdDevsOnExpiry.emplace_back(0.384457068812631);   stdDevsOnPayment.emplace_back(0.38559435044422);
            strikes.emplace_back(0.495);   stdDevsOnExpiry.emplace_back(0.384744836079706);   stdDevsOnPayment.emplace_back(0.385931348960537);
            strikes.emplace_back(0.496);   stdDevsOnExpiry.emplace_back(0.385048414735082);   stdDevsOnPayment.emplace_back(0.386222982291966);
            strikes.emplace_back(0.497);   stdDevsOnExpiry.emplace_back(0.38537412933408);    stdDevsOnPayment.emplace_back(0.386559980808283);
            strikes.emplace_back(0.498);   stdDevsOnExpiry.emplace_back(0.385709330766058);   stdDevsOnPayment.emplace_back(0.386848373769362);
            strikes.emplace_back(0.499);   stdDevsOnExpiry.emplace_back(0.386038207642715);   stdDevsOnPayment.emplace_back(0.387146487841489);
            strikes.emplace_back(0.5);     stdDevsOnExpiry.emplace_back(0.386338624020431);   stdDevsOnPayment.emplace_back(0.387460803765362);
            strikes.emplace_back(0.501);   stdDevsOnExpiry.emplace_back(0.386623229009846);   stdDevsOnPayment.emplace_back(0.387765398578187);
            strikes.emplace_back(0.502);   stdDevsOnExpiry.emplace_back(0.386952105886504);   stdDevsOnPayment.emplace_back(0.388108877835202);
            strikes.emplace_back(0.503);   stdDevsOnExpiry.emplace_back(0.38724935998656);    stdDevsOnPayment.emplace_back(0.388374588203837);
            strikes.emplace_back(0.504);   stdDevsOnExpiry.emplace_back(0.387584561418537);   stdDevsOnPayment.emplace_back(0.388718067460853);
            strikes.emplace_back(0.505);   stdDevsOnExpiry.emplace_back(0.387856517297312);   stdDevsOnPayment.emplace_back(0.388980537459138);
            strikes.emplace_back(0.506);   stdDevsOnExpiry.emplace_back(0.388213854672911);   stdDevsOnPayment.emplace_back(0.3893402185679);
            strikes.emplace_back(0.507);   stdDevsOnExpiry.emplace_back(0.388514271050627);   stdDevsOnPayment.emplace_back(0.389605928936535);
            strikes.emplace_back(0.508);   stdDevsOnExpiry.emplace_back(0.388748279597479);   stdDevsOnPayment.emplace_back(0.389994773378439);
            strikes.emplace_back(0.509);   stdDevsOnExpiry.emplace_back(0.389073994196477);   stdDevsOnPayment.emplace_back(0.390247522265677);
            strikes.emplace_back(0.51);    stdDevsOnExpiry.emplace_back(0.389431331572076);   stdDevsOnPayment.emplace_back(0.39056183818955);
            strikes.emplace_back(0.511);   stdDevsOnExpiry.emplace_back(0.389643204175307);   stdDevsOnPayment.emplace_back(0.390882634854121);
            strikes.emplace_back(0.512);   stdDevsOnExpiry.emplace_back(0.390025839772187);   stdDevsOnPayment.emplace_back(0.391193710407644);
            strikes.emplace_back(0.513);   stdDevsOnExpiry.emplace_back(0.390351554371185);   stdDevsOnPayment.emplace_back(0.391491824479771);
            strikes.emplace_back(0.514);   stdDevsOnExpiry.emplace_back(0.390607698861658);   stdDevsOnPayment.emplace_back(0.391828822996088);
            strikes.emplace_back(0.515);   stdDevsOnExpiry.emplace_back(0.390901790684054);   stdDevsOnPayment.emplace_back(0.392081571883326);
            strikes.emplace_back(0.516);   stdDevsOnExpiry.emplace_back(0.391218018450071);   stdDevsOnPayment.emplace_back(0.392399128177548);
            strikes.emplace_back(0.517);   stdDevsOnExpiry.emplace_back(0.391524759383107);   stdDevsOnPayment.emplace_back(0.392726405582818);
            strikes.emplace_back(0.518);   stdDevsOnExpiry.emplace_back(0.391809364372522);   stdDevsOnPayment.emplace_back(0.393034240765992);
            strikes.emplace_back(0.519);   stdDevsOnExpiry.emplace_back(0.392093969361938);   stdDevsOnPayment.emplace_back(0.393261066690437);
            strikes.emplace_back(0.52);    stdDevsOnExpiry.emplace_back(0.392419683960935);   stdDevsOnPayment.emplace_back(0.393617507428849);
            strikes.emplace_back(0.521);   stdDevsOnExpiry.emplace_back(0.392723262616311);   stdDevsOnPayment.emplace_back(0.393948025204468);
            strikes.emplace_back(0.522);   stdDevsOnExpiry.emplace_back(0.393001543050406);   stdDevsOnPayment.emplace_back(0.39422345668415);
            strikes.emplace_back(0.523);   stdDevsOnExpiry.emplace_back(0.393336744482384);   stdDevsOnPayment.emplace_back(0.394476205571388);
            strikes.emplace_back(0.524);   stdDevsOnExpiry.emplace_back(0.393611862638818);   stdDevsOnPayment.emplace_back(0.394813204087705);
            strikes.emplace_back(0.525);   stdDevsOnExpiry.emplace_back(0.393940739515476);   stdDevsOnPayment.emplace_back(0.395075674085991);
            strikes.emplace_back(0.526);   stdDevsOnExpiry.emplace_back(0.394168423507008);   stdDevsOnPayment.emplace_back(0.395435355194752);
            strikes.emplace_back(0.527);   stdDevsOnExpiry.emplace_back(0.394459353051744);   stdDevsOnPayment.emplace_back(0.39568810408199);
            strikes.emplace_back(0.528);   stdDevsOnExpiry.emplace_back(0.394800879039042);   stdDevsOnPayment.emplace_back(0.396012141116911);
            strikes.emplace_back(0.529);   stdDevsOnExpiry.emplace_back(0.395025400752914);   stdDevsOnPayment.emplace_back(0.396274611115196);
            strikes.emplace_back(0.53);    stdDevsOnExpiry.emplace_back(0.395347953074251);   stdDevsOnPayment.emplace_back(0.396572725187323);
            strikes.emplace_back(0.531);   stdDevsOnExpiry.emplace_back(0.395613584397705);   stdDevsOnPayment.emplace_back(0.396893521851894);
            strikes.emplace_back(0.532);   stdDevsOnExpiry.emplace_back(0.395951948107343);   stdDevsOnPayment.emplace_back(0.397162472590878);
            strikes.emplace_back(0.533);   stdDevsOnExpiry.emplace_back(0.39626817587336);    stdDevsOnPayment.emplace_back(0.397483269255449);
            strikes.emplace_back(0.534);   stdDevsOnExpiry.emplace_back(0.396505346697872);   stdDevsOnPayment.emplace_back(0.397774902586878);
            strikes.emplace_back(0.535);   stdDevsOnExpiry.emplace_back(0.39683106129687);    stdDevsOnPayment.emplace_back(0.398034132214814);
            strikes.emplace_back(0.536);   stdDevsOnExpiry.emplace_back(0.397137802229906);   stdDevsOnPayment.emplace_back(0.398329005916591);
            strikes.emplace_back(0.537);   stdDevsOnExpiry.emplace_back(0.39739078444272);    stdDevsOnPayment.emplace_back(0.398656283321861);
            strikes.emplace_back(0.538);   stdDevsOnExpiry.emplace_back(0.397672227154475);   stdDevsOnPayment.emplace_back(0.39890579183875);
            strikes.emplace_back(0.539);   stdDevsOnExpiry.emplace_back(0.397931533922608);   stdDevsOnPayment.emplace_back(0.399184463688781);
            strikes.emplace_back(0.54);    stdDevsOnExpiry.emplace_back(0.398269897632246);   stdDevsOnPayment.emplace_back(0.399463135538813);
            strikes.emplace_back(0.541);   stdDevsOnExpiry.emplace_back(0.398484932513138);   stdDevsOnPayment.emplace_back(0.399770970721987);
            strikes.emplace_back(0.542);   stdDevsOnExpiry.emplace_back(0.398826458500436);   stdDevsOnPayment.emplace_back(0.400023719609225);
            strikes.emplace_back(0.543);   stdDevsOnExpiry.emplace_back(0.399126874878152);   stdDevsOnPayment.emplace_back(0.400334795162749);
            strikes.emplace_back(0.544);   stdDevsOnExpiry.emplace_back(0.399398830756927);   stdDevsOnPayment.emplace_back(0.400603745901733);
            strikes.emplace_back(0.545);   stdDevsOnExpiry.emplace_back(0.399680273468682);   stdDevsOnPayment.emplace_back(0.400914821455256);
            strikes.emplace_back(0.546);   stdDevsOnExpiry.emplace_back(0.399917444293194);   stdDevsOnPayment.emplace_back(0.40122589700878);
            strikes.emplace_back(0.547);   stdDevsOnExpiry.emplace_back(0.400274781668793);   stdDevsOnPayment.emplace_back(0.401465684414621);
            strikes.emplace_back(0.548);   stdDevsOnExpiry.emplace_back(0.400543575269907);   stdDevsOnPayment.emplace_back(0.401812404041986);
            strikes.emplace_back(0.549);   stdDevsOnExpiry.emplace_back(0.400828180259323);   stdDevsOnPayment.emplace_back(0.402016547373986);
            strikes.emplace_back(0.55);    stdDevsOnExpiry.emplace_back(0.401103298415757);   stdDevsOnPayment.emplace_back(0.402347065149604);
            strikes.emplace_back(0.551);   stdDevsOnExpiry.emplace_back(0.401356280628571);   stdDevsOnPayment.emplace_back(0.402612775518239);
            strikes.emplace_back(0.552);   stdDevsOnExpiry.emplace_back(0.401666183839267);   stdDevsOnPayment.emplace_back(0.402875245516525);
            strikes.emplace_back(0.553);   stdDevsOnExpiry.emplace_back(0.401874894164838);   stdDevsOnPayment.emplace_back(0.403186321070048);
            strikes.emplace_back(0.554);   stdDevsOnExpiry.emplace_back(0.402238556095758);   stdDevsOnPayment.emplace_back(0.403468233290429);
            strikes.emplace_back(0.555);   stdDevsOnExpiry.emplace_back(0.40246624008729);    stdDevsOnPayment.emplace_back(0.403714501436968);
            strikes.emplace_back(0.556);   stdDevsOnExpiry.emplace_back(0.402747682799045);   stdDevsOnPayment.emplace_back(0.403989932916651);
            strikes.emplace_back(0.557);   stdDevsOnExpiry.emplace_back(0.403070235120382);   stdDevsOnPayment.emplace_back(0.404297768099825);
            strikes.emplace_back(0.558);   stdDevsOnExpiry.emplace_back(0.403326379610856);   stdDevsOnPayment.emplace_back(0.404589401431253);
            strikes.emplace_back(0.559);   stdDevsOnExpiry.emplace_back(0.403535089936427);   stdDevsOnPayment.emplace_back(0.404884275133031);
            strikes.emplace_back(0.56);    stdDevsOnExpiry.emplace_back(0.403797558982221);   stdDevsOnPayment.emplace_back(0.405104620316777);
            strikes.emplace_back(0.561);   stdDevsOnExpiry.emplace_back(0.404126435858878);   stdDevsOnPayment.emplace_back(0.405393013277856);
            strikes.emplace_back(0.562);   stdDevsOnExpiry.emplace_back(0.404423689958934);   stdDevsOnPayment.emplace_back(0.405668444757538);
            strikes.emplace_back(0.563);   stdDevsOnExpiry.emplace_back(0.404733593169631);   stdDevsOnPayment.emplace_back(0.405901751422681);
            strikes.emplace_back(0.564);   stdDevsOnExpiry.emplace_back(0.404977088549464);   stdDevsOnPayment.emplace_back(0.406203105865157);
            strikes.emplace_back(0.565);   stdDevsOnExpiry.emplace_back(0.405220583929296);   stdDevsOnPayment.emplace_back(0.406510941048331);
            strikes.emplace_back(0.566);   stdDevsOnExpiry.emplace_back(0.405533649417653);   stdDevsOnPayment.emplace_back(0.406741007343125);
            strikes.emplace_back(0.567);   stdDevsOnExpiry.emplace_back(0.405742359743224);   stdDevsOnPayment.emplace_back(0.407042361785601);
            strikes.emplace_back(0.568);   stdDevsOnExpiry.emplace_back(0.406083885730522);   stdDevsOnPayment.emplace_back(0.407295110672839);
            strikes.emplace_back(0.569);   stdDevsOnExpiry.emplace_back(0.406295758333754);   stdDevsOnPayment.emplace_back(0.407593224744966);
            strikes.emplace_back(0.57);    stdDevsOnExpiry.emplace_back(0.406643608876372);   stdDevsOnPayment.emplace_back(0.40786217548395);
            strikes.emplace_back(0.571);   stdDevsOnExpiry.emplace_back(0.406814371870021);   stdDevsOnPayment.emplace_back(0.408082520667695);
            strikes.emplace_back(0.572);   stdDevsOnExpiry.emplace_back(0.407114788247737);   stdDevsOnPayment.emplace_back(0.408416278813663);
            strikes.emplace_back(0.573);   stdDevsOnExpiry.emplace_back(0.407434178291414);   stdDevsOnPayment.emplace_back(0.408636623997409);
            strikes.emplace_back(0.574);   stdDevsOnExpiry.emplace_back(0.407623914951024);   stdDevsOnPayment.emplace_back(0.40891853621779);
            strikes.emplace_back(0.575);   stdDevsOnExpiry.emplace_back(0.407936980439381);   stdDevsOnPayment.emplace_back(0.409252294363758);
            strikes.emplace_back(0.576);   stdDevsOnExpiry.emplace_back(0.408265857316039);   stdDevsOnPayment.emplace_back(0.409449956955059);
            strikes.emplace_back(0.577);   stdDevsOnExpiry.emplace_back(0.408496703585231);   stdDevsOnPayment.emplace_back(0.409748071027186);
            strikes.emplace_back(0.578);   stdDevsOnExpiry.emplace_back(0.408774984019326);   stdDevsOnPayment.emplace_back(0.410000819914424);
            strikes.emplace_back(0.579);   stdDevsOnExpiry.emplace_back(0.409024803954479);   stdDevsOnPayment.emplace_back(0.410295693616202);
            strikes.emplace_back(0.58);    stdDevsOnExpiry.emplace_back(0.409261974778992);   stdDevsOnPayment.emplace_back(0.410567884725535);
            strikes.emplace_back(0.581);   stdDevsOnExpiry.emplace_back(0.409518119269465);   stdDevsOnPayment.emplace_back(0.410755826205789);
            strikes.emplace_back(0.582);   stdDevsOnExpiry.emplace_back(0.409831184757822);   stdDevsOnPayment.emplace_back(0.411125228425598);
            strikes.emplace_back(0.583);   stdDevsOnExpiry.emplace_back(0.410109465191917);   stdDevsOnPayment.emplace_back(0.411394179164582);
            strikes.emplace_back(0.584);   stdDevsOnExpiry.emplace_back(0.41034979829409);    stdDevsOnPayment.emplace_back(0.411614524348328);
            strikes.emplace_back(0.585);   stdDevsOnExpiry.emplace_back(0.410583806840942);   stdDevsOnPayment.emplace_back(0.411867273235566);
            strikes.emplace_back(0.586);   stdDevsOnExpiry.emplace_back(0.410931657383561);   stdDevsOnPayment.emplace_back(0.412171868048391);
            strikes.emplace_back(0.587);   stdDevsOnExpiry.emplace_back(0.41109925809955);    stdDevsOnPayment.emplace_back(0.41237925175074);
            strikes.emplace_back(0.588);   stdDevsOnExpiry.emplace_back(0.411374376255984);   stdDevsOnPayment.emplace_back(0.412670885082168);
            strikes.emplace_back(0.589);   stdDevsOnExpiry.emplace_back(0.411690604022001);   stdDevsOnPayment.emplace_back(0.41294631656185);
            strikes.emplace_back(0.59);    stdDevsOnExpiry.emplace_back(0.411943586234814);   stdDevsOnPayment.emplace_back(0.413215267300834);
            strikes.emplace_back(0.591);   stdDevsOnExpiry.emplace_back(0.412215542113589);   stdDevsOnPayment.emplace_back(0.413458295077025);
            strikes.emplace_back(0.592);   stdDevsOnExpiry.emplace_back(0.412348357775316);   stdDevsOnPayment.emplace_back(0.413707803593913);
            strikes.emplace_back(0.593);   stdDevsOnExpiry.emplace_back(0.412677234651974);   stdDevsOnPayment.emplace_back(0.413931389148008);
            strikes.emplace_back(0.594);   stdDevsOnExpiry.emplace_back(0.412908080921166);   stdDevsOnPayment.emplace_back(0.4141679361835);
            strikes.emplace_back(0.595);   stdDevsOnExpiry.emplace_back(0.41318003679994);    stdDevsOnPayment.emplace_back(0.414479011737024);
            strikes.emplace_back(0.596);   stdDevsOnExpiry.emplace_back(0.413382422570191);   stdDevsOnPayment.emplace_back(0.414747962476008);
            strikes.emplace_back(0.597);   stdDevsOnExpiry.emplace_back(0.413701812613868);   stdDevsOnPayment.emplace_back(0.415000711363246);
            strikes.emplace_back(0.598);   stdDevsOnExpiry.emplace_back(0.413954794826682);   stdDevsOnPayment.emplace_back(0.415282623583626);
            strikes.emplace_back(0.599);   stdDevsOnExpiry.emplace_back(0.414229912983116);   stdDevsOnPayment.emplace_back(0.415538612841214);
            strikes.emplace_back(0.6);     stdDevsOnExpiry.emplace_back(0.414508193417211);   stdDevsOnPayment.emplace_back(0.415778400247055);
            strikes.emplace_back(0.601);   stdDevsOnExpiry.emplace_back(0.414748526519384);   stdDevsOnPayment.emplace_back(0.416005226171499);
            strikes.emplace_back(0.602);   stdDevsOnExpiry.emplace_back(0.415020482398158);   stdDevsOnPayment.emplace_back(0.41632602283607);
            strikes.emplace_back(0.603);   stdDevsOnExpiry.emplace_back(0.415276626888632);   stdDevsOnPayment.emplace_back(0.416513964316324);
            strikes.emplace_back(0.604);   stdDevsOnExpiry.emplace_back(0.415583367821668);   stdDevsOnPayment.emplace_back(0.416795876536705);
            strikes.emplace_back(0.605);   stdDevsOnExpiry.emplace_back(0.415731994871696);   stdDevsOnPayment.emplace_back(0.417119913571625);
            strikes.emplace_back(0.606);   stdDevsOnExpiry.emplace_back(0.416060871748354);   stdDevsOnPayment.emplace_back(0.417217124682101);
            strikes.emplace_back(0.607);   stdDevsOnExpiry.emplace_back(0.416335989904788);   stdDevsOnPayment.emplace_back(0.417550882828069);
            strikes.emplace_back(0.608);   stdDevsOnExpiry.emplace_back(0.416532051119719);   stdDevsOnPayment.emplace_back(0.417861958381593);
            strikes.emplace_back(0.609);   stdDevsOnExpiry.emplace_back(0.416772384221892);   stdDevsOnPayment.emplace_back(0.41811794763918);
            strikes.emplace_back(0.61);    stdDevsOnExpiry.emplace_back(0.417028528712365);   stdDevsOnPayment.emplace_back(0.418286446897339);
            strikes.emplace_back(0.611);   stdDevsOnExpiry.emplace_back(0.417281510925179);   stdDevsOnPayment.emplace_back(0.418587801339814);
            strikes.emplace_back(0.612);   stdDevsOnExpiry.emplace_back(0.417455436196488);   stdDevsOnPayment.emplace_back(0.418817867634608);
            strikes.emplace_back(0.613);   stdDevsOnExpiry.emplace_back(0.417721067519942);   stdDevsOnPayment.emplace_back(0.419047933929401);
            strikes.emplace_back(0.614);   stdDevsOnExpiry.emplace_back(0.418068918062561);   stdDevsOnPayment.emplace_back(0.419313644298036);
            strikes.emplace_back(0.615);   stdDevsOnExpiry.emplace_back(0.418204896001948);   stdDevsOnPayment.emplace_back(0.41962471985156);
            strikes.emplace_back(0.616);   stdDevsOnExpiry.emplace_back(0.418533772878605);   stdDevsOnPayment.emplace_back(0.419812661331814);
            strikes.emplace_back(0.617);   stdDevsOnExpiry.emplace_back(0.418770943703118);   stdDevsOnPayment.emplace_back(0.420036246885909);
            strikes.emplace_back(0.618);   stdDevsOnExpiry.emplace_back(0.419042899581892);   stdDevsOnPayment.emplace_back(0.420282515032448);
            strikes.emplace_back(0.619);   stdDevsOnExpiry.emplace_back(0.419229473963842);   stdDevsOnPayment.emplace_back(0.420596830956321);
            strikes.emplace_back(0.62);    stdDevsOnExpiry.emplace_back(0.419529890341558);   stdDevsOnPayment.emplace_back(0.420807455029019);
            strikes.emplace_back(0.621);   stdDevsOnExpiry.emplace_back(0.41974492522245);    stdDevsOnPayment.emplace_back(0.421037521323813);
            strikes.emplace_back(0.622);   stdDevsOnExpiry.emplace_back(0.420026367934205);   stdDevsOnPayment.emplace_back(0.421335635395939);
            strikes.emplace_back(0.623);   stdDevsOnExpiry.emplace_back(0.420282512424678);   stdDevsOnPayment.emplace_back(0.421523576876193);
            strikes.emplace_back(0.624);   stdDevsOnExpiry.emplace_back(0.420427977197046);   stdDevsOnPayment.emplace_back(0.421857335022161);
            strikes.emplace_back(0.625);   stdDevsOnExpiry.emplace_back(0.420658823466239);   stdDevsOnPayment.emplace_back(0.422061478354161);
            strikes.emplace_back(0.626);   stdDevsOnExpiry.emplace_back(0.421019323119498);   stdDevsOnPayment.emplace_back(0.422259140945463);
            strikes.emplace_back(0.627);   stdDevsOnExpiry.emplace_back(0.421193248390807);   stdDevsOnPayment.emplace_back(0.422589658721081);
            strikes.emplace_back(0.628);   stdDevsOnExpiry.emplace_back(0.421395634161058);   stdDevsOnPayment.emplace_back(0.422774359830986);
            strikes.emplace_back(0.629);   stdDevsOnExpiry.emplace_back(0.421702375094094);   stdDevsOnPayment.emplace_back(0.423036829829271);
            strikes.emplace_back(0.63);    stdDevsOnExpiry.emplace_back(0.421882624920724);   stdDevsOnPayment.emplace_back(0.423273376864763);
            strikes.emplace_back(0.631);   stdDevsOnExpiry.emplace_back(0.422176716743119);   stdDevsOnPayment.emplace_back(0.423493722048509);
            strikes.emplace_back(0.632);   stdDevsOnExpiry.emplace_back(0.422401238456991);   stdDevsOnPayment.emplace_back(0.423720547972953);
            strikes.emplace_back(0.633);   stdDevsOnExpiry.emplace_back(0.422651058392145);   stdDevsOnPayment.emplace_back(0.423931172045652);
            strikes.emplace_back(0.634);   stdDevsOnExpiry.emplace_back(0.434664551223124);   stdDevsOnPayment.emplace_back(0.424248728339874);
            strikes.emplace_back(0.635);   stdDevsOnExpiry.emplace_back(0.434920695713598);   stdDevsOnPayment.emplace_back(0.424478794634667);
            strikes.emplace_back(0.636);   stdDevsOnExpiry.emplace_back(0.423403680475265);   stdDevsOnPayment.emplace_back(0.424715341670159);
            strikes.emplace_back(0.637);   stdDevsOnExpiry.emplace_back(0.435432984694545);   stdDevsOnPayment.emplace_back(0.424990773149841);
            strikes.emplace_back(0.638);   stdDevsOnExpiry.emplace_back(0.423846399347688);   stdDevsOnPayment.emplace_back(0.425172233889397);
            strikes.emplace_back(0.639);   stdDevsOnExpiry.emplace_back(0.435910488621231);   stdDevsOnPayment.emplace_back(0.425428223146984);
            strikes.emplace_back(0.64);    stdDevsOnExpiry.emplace_back(0.436163470834044);   stdDevsOnPayment.emplace_back(0.425645327960381);
            strikes.emplace_back(0.641);   stdDevsOnExpiry.emplace_back(0.436413290769197);   stdDevsOnPayment.emplace_back(0.425894836477269);
            strikes.emplace_back(0.642);   stdDevsOnExpiry.emplace_back(0.43664413703839);    stdDevsOnPayment.emplace_back(0.42617674869765);
            strikes.emplace_back(0.643);   stdDevsOnExpiry.emplace_back(0.436897119251203);   stdDevsOnPayment.emplace_back(0.426361449807555);
            strikes.emplace_back(0.644);   stdDevsOnExpiry.emplace_back(0.437159588296997);   stdDevsOnPayment.emplace_back(0.426607717954094);
            strikes.emplace_back(0.645);   stdDevsOnExpiry.emplace_back(0.437358811789588);   stdDevsOnPayment.emplace_back(0.42678917869365);
            strikes.emplace_back(0.646);   stdDevsOnExpiry.emplace_back(0.437602307169421);   stdDevsOnPayment.emplace_back(0.427165061654157);
            strikes.emplace_back(0.647);   stdDevsOnExpiry.emplace_back(0.437852127104574);   stdDevsOnPayment.emplace_back(0.427340041653014);
            strikes.emplace_back(0.648);   stdDevsOnExpiry.emplace_back(0.438076648818446);   stdDevsOnPayment.emplace_back(0.427608992391998);
            strikes.emplace_back(0.649);   stdDevsOnExpiry.emplace_back(0.43833279330892);    stdDevsOnPayment.emplace_back(0.427780732020506);
            strikes.emplace_back(0.65);    stdDevsOnExpiry.emplace_back(0.438569964133432);   stdDevsOnPayment.emplace_back(0.428033480907744);
            strikes.emplace_back(0.651);   stdDevsOnExpiry.emplace_back(0.438810297235605);   stdDevsOnPayment.emplace_back(0.428308912387426);
            strikes.emplace_back(0.652);   stdDevsOnExpiry.emplace_back(0.439056954893098);   stdDevsOnPayment.emplace_back(0.428548699793267);
            strikes.emplace_back(0.653);   stdDevsOnExpiry.emplace_back(0.43927831432931);    stdDevsOnPayment.emplace_back(0.428713958681077);
            strikes.emplace_back(0.654);   stdDevsOnExpiry.emplace_back(0.439534458819784);   stdDevsOnPayment.emplace_back(0.429012072753204);
            strikes.emplace_back(0.655);   stdDevsOnExpiry.emplace_back(0.439765305088976);   stdDevsOnPayment.emplace_back(0.429193533492759);
            strikes.emplace_back(0.656);   stdDevsOnExpiry.emplace_back(0.440015125024129);   stdDevsOnPayment.emplace_back(0.429443042009648);
            strikes.emplace_back(0.657);   stdDevsOnExpiry.emplace_back(0.440264944959282);   stdDevsOnPayment.emplace_back(0.429686069785838);
            strikes.emplace_back(0.658);   stdDevsOnExpiry.emplace_back(0.440454681618893);   stdDevsOnPayment.emplace_back(0.429857809414346);
            strikes.emplace_back(0.659);   stdDevsOnExpiry.emplace_back(0.440723475220007);   stdDevsOnPayment.emplace_back(0.430185086819615);
            strikes.emplace_back(0.66);    stdDevsOnExpiry.emplace_back(0.440957483766859);   stdDevsOnPayment.emplace_back(0.430402191633012);
            strikes.emplace_back(0.661);   stdDevsOnExpiry.emplace_back(0.441210465979673);   stdDevsOnPayment.emplace_back(0.43057393126152);
            strikes.emplace_back(0.662);   stdDevsOnExpiry.emplace_back(0.441403364916943);   stdDevsOnPayment.emplace_back(0.430758632371424);
            strikes.emplace_back(0.663);   stdDevsOnExpiry.emplace_back(0.441681645351038);   stdDevsOnPayment.emplace_back(0.431076188665647);
            strikes.emplace_back(0.664);   stdDevsOnExpiry.emplace_back(0.441858732900007);   stdDevsOnPayment.emplace_back(0.431290053108694);
            strikes.emplace_back(0.665);   stdDevsOnExpiry.emplace_back(0.44209590372452);    stdDevsOnPayment.emplace_back(0.431484475329646);
            strikes.emplace_back(0.666);   stdDevsOnExpiry.emplace_back(0.442374184158615);   stdDevsOnPayment.emplace_back(0.431721022365138);
            strikes.emplace_back(0.667);   stdDevsOnExpiry.emplace_back(0.442592381317166);   stdDevsOnPayment.emplace_back(0.432022376807614);
            strikes.emplace_back(0.668);   stdDevsOnExpiry.emplace_back(0.442813740753378);   stdDevsOnPayment.emplace_back(0.432194116436122);
            strikes.emplace_back(0.669);   stdDevsOnExpiry.emplace_back(0.443076209799172);   stdDevsOnPayment.emplace_back(0.432417701990217);
            strikes.emplace_back(0.67);    stdDevsOnExpiry.emplace_back(0.443275433291763);   stdDevsOnPayment.emplace_back(0.432644527914661);
            strikes.emplace_back(0.671);   stdDevsOnExpiry.emplace_back(0.443531577782236);   stdDevsOnPayment.emplace_back(0.43293616124609);
            strikes.emplace_back(0.672);   stdDevsOnExpiry.emplace_back(0.443705503053546);   stdDevsOnPayment.emplace_back(0.433030131986216);
            strikes.emplace_back(0.673);   stdDevsOnExpiry.emplace_back(0.443964809821679);   stdDevsOnPayment.emplace_back(0.433376851613581);
            strikes.emplace_back(0.674);   stdDevsOnExpiry.emplace_back(0.444183006980231);   stdDevsOnPayment.emplace_back(0.433603677538025);
            strikes.emplace_back(0.675);   stdDevsOnExpiry.emplace_back(0.444439151470705);   stdDevsOnPayment.emplace_back(0.433772176796184);
            strikes.emplace_back(0.676);   stdDevsOnExpiry.emplace_back(0.444676322295217);   stdDevsOnPayment.emplace_back(0.434083252349708);
            strikes.emplace_back(0.677);   stdDevsOnExpiry.emplace_back(0.44491665539739);    stdDevsOnPayment.emplace_back(0.434290636052057);
            strikes.emplace_back(0.678);   stdDevsOnExpiry.emplace_back(0.44510322977934);    stdDevsOnPayment.emplace_back(0.43440080864393);
            strikes.emplace_back(0.679);   stdDevsOnExpiry.emplace_back(0.445324589215552);   stdDevsOnPayment.emplace_back(0.43464383642012);
            strikes.emplace_back(0.68);    stdDevsOnExpiry.emplace_back(0.445602869649647);   stdDevsOnPayment.emplace_back(0.434909546788755);
            strikes.emplace_back(0.681);   stdDevsOnExpiry.emplace_back(0.445757821254995);   stdDevsOnPayment.emplace_back(0.435084526787612);
            strikes.emplace_back(0.682);   stdDevsOnExpiry.emplace_back(0.446029777133769);   stdDevsOnPayment.emplace_back(0.43537616011904);
            strikes.emplace_back(0.683);   stdDevsOnExpiry.emplace_back(0.446282759346583);   stdDevsOnPayment.emplace_back(0.435567341969643);
            strikes.emplace_back(0.684);   stdDevsOnExpiry.emplace_back(0.446478820561513);   stdDevsOnPayment.emplace_back(0.435719639376056);
            strikes.emplace_back(0.685);   stdDevsOnExpiry.emplace_back(0.446655908110483);   stdDevsOnPayment.emplace_back(0.435969147892944);
            strikes.emplace_back(0.686);   stdDevsOnExpiry.emplace_back(0.446940513099898);   stdDevsOnPayment.emplace_back(0.43627050233542);
            strikes.emplace_back(0.687);   stdDevsOnExpiry.emplace_back(0.44716187253611);    stdDevsOnPayment.emplace_back(0.448587150032745);
            strikes.emplace_back(0.688);   stdDevsOnExpiry.emplace_back(0.447424341581904);   stdDevsOnPayment.emplace_back(0.436633423814531);
            strikes.emplace_back(0.689);   stdDevsOnExpiry.emplace_back(0.447604591408533);   stdDevsOnPayment.emplace_back(0.436827846035483);
            strikes.emplace_back(0.69);    stdDevsOnExpiry.emplace_back(0.447800652623464);   stdDevsOnPayment.emplace_back(0.44929679113922);
            strikes.emplace_back(0.691);   stdDevsOnExpiry.emplace_back(0.448085257612879);   stdDevsOnPayment.emplace_back(0.449510655582268);
            strikes.emplace_back(0.692);   stdDevsOnExpiry.emplace_back(0.44830029249377);    stdDevsOnPayment.emplace_back(0.449750442988109);
            strikes.emplace_back(0.693);   stdDevsOnExpiry.emplace_back(0.448505840541681);   stdDevsOnPayment.emplace_back(0.449986990023601);
            strikes.emplace_back(0.694);   stdDevsOnExpiry.emplace_back(0.448708226311932);   stdDevsOnPayment.emplace_back(0.450226777429442);
            strikes.emplace_back(0.695);   stdDevsOnExpiry.emplace_back(0.448935910303464);   stdDevsOnPayment.emplace_back(0.450447122613188);
            strikes.emplace_back(0.696);   stdDevsOnExpiry.emplace_back(0.449141458351375);   stdDevsOnPayment.emplace_back(0.450638304463791);
            strikes.emplace_back(0.697);   stdDevsOnExpiry.emplace_back(0.449365980065247);   stdDevsOnPayment.emplace_back(0.450878091869632);
            strikes.emplace_back(0.698);   stdDevsOnExpiry.emplace_back(0.449565203557838);   stdDevsOnPayment.emplace_back(0.451091956312679);
            strikes.emplace_back(0.699);   stdDevsOnExpiry.emplace_back(0.449871944490874);   stdDevsOnPayment.emplace_back(0.451296099644679);
            strikes.emplace_back(0.7);     stdDevsOnExpiry.emplace_back(0.450061681150484);   stdDevsOnPayment.emplace_back(0.451519685198774);
            strikes.emplace_back(0.701);   stdDevsOnExpiry.emplace_back(0.450229281866473);   stdDevsOnPayment.emplace_back(0.45174003038252);
            strikes.emplace_back(0.702);   stdDevsOnExpiry.emplace_back(0.450472777246306);   stdDevsOnPayment.emplace_back(0.451979817788361);
            strikes.emplace_back(0.703);   stdDevsOnExpiry.emplace_back(0.45073208401444);    stdDevsOnPayment.emplace_back(0.452222845564551);
            strikes.emplace_back(0.704);   stdDevsOnExpiry.emplace_back(0.450940794340011);   stdDevsOnPayment.emplace_back(0.452436710007599);
            strikes.emplace_back(0.705);   stdDevsOnExpiry.emplace_back(0.451162153776223);   stdDevsOnPayment.emplace_back(0.452653814820995);
            strikes.emplace_back(0.706);   stdDevsOnExpiry.emplace_back(0.451374026379454);   stdDevsOnPayment.emplace_back(0.452848237041948);
            strikes.emplace_back(0.707);   stdDevsOnExpiry.emplace_back(0.451570087594384);   stdDevsOnPayment.emplace_back(0.453094505188487);
            strikes.emplace_back(0.708);   stdDevsOnExpiry.emplace_back(0.451813582974217);   stdDevsOnPayment.emplace_back(0.453292167779789);
            strikes.emplace_back(0.709);   stdDevsOnExpiry.emplace_back(0.452028617855109);   stdDevsOnPayment.emplace_back(0.453525474444931);
            strikes.emplace_back(0.71);    stdDevsOnExpiry.emplace_back(0.45223416590302);    stdDevsOnPayment.emplace_back(0.453729617776931);
            strikes.emplace_back(0.711);   stdDevsOnExpiry.emplace_back(0.45243338939561);    stdDevsOnPayment.emplace_back(0.453959684071725);
            strikes.emplace_back(0.712);   stdDevsOnExpiry.emplace_back(0.452635775165861);   stdDevsOnPayment.emplace_back(0.454212432958962);
            strikes.emplace_back(0.713);   stdDevsOnExpiry.emplace_back(0.452869783712714);   stdDevsOnPayment.emplace_back(0.454390653328169);
            strikes.emplace_back(0.714);   stdDevsOnExpiry.emplace_back(0.453056358094664);   stdDevsOnPayment.emplace_back(0.454617479252613);
            strikes.emplace_back(0.715);   stdDevsOnExpiry.emplace_back(0.453274555253215);   stdDevsOnPayment.emplace_back(0.454837824436359);
            strikes.emplace_back(0.716);   stdDevsOnExpiry.emplace_back(0.453578133908591);   stdDevsOnPayment.emplace_back(0.455019285175914);
            strikes.emplace_back(0.717);   stdDevsOnExpiry.emplace_back(0.453761546012881);   stdDevsOnPayment.emplace_back(0.455252591841057);
            strikes.emplace_back(0.718);   stdDevsOnExpiry.emplace_back(0.453970256338452);   stdDevsOnPayment.emplace_back(0.455479417765501);
            strikes.emplace_back(0.719);   stdDevsOnExpiry.emplace_back(0.454156830720402);   stdDevsOnPayment.emplace_back(0.455657638134707);
            strikes.emplace_back(0.72);    stdDevsOnExpiry.emplace_back(0.454346567380012);   stdDevsOnPayment.emplace_back(0.45589094479985);
            strikes.emplace_back(0.721);   stdDevsOnExpiry.emplace_back(0.454590062759845);   stdDevsOnPayment.emplace_back(0.456098328502199);
            strikes.emplace_back(0.722);   stdDevsOnExpiry.emplace_back(0.454795610807756);   stdDevsOnPayment.emplace_back(0.45634135627839);
            strikes.emplace_back(0.723);   stdDevsOnExpiry.emplace_back(0.455032781632269);   stdDevsOnPayment.emplace_back(0.45654225924004);
            strikes.emplace_back(0.724);   stdDevsOnExpiry.emplace_back(0.455257303346141);   stdDevsOnPayment.emplace_back(0.45674640257204);
            strikes.emplace_back(0.725);   stdDevsOnExpiry.emplace_back(0.455469175949372);   stdDevsOnPayment.emplace_back(0.456960267015087);
            strikes.emplace_back(0.726);   stdDevsOnExpiry.emplace_back(0.455684210830264);   stdDevsOnPayment.emplace_back(0.457200054420929);
            strikes.emplace_back(0.727);   stdDevsOnExpiry.emplace_back(0.455908732544135);   stdDevsOnPayment.emplace_back(0.457378274790135);
            strikes.emplace_back(0.728);   stdDevsOnExpiry.emplace_back(0.456054197316503);   stdDevsOnPayment.emplace_back(0.457621302566325);
            strikes.emplace_back(0.729);   stdDevsOnExpiry.emplace_back(0.456281881308035);   stdDevsOnPayment.emplace_back(0.457835167009372);
            strikes.emplace_back(0.73);    stdDevsOnExpiry.emplace_back(0.456525376687868);   stdDevsOnPayment.emplace_back(0.458039310341372);
            strikes.emplace_back(0.731);   stdDevsOnExpiry.emplace_back(0.456683490570877);   stdDevsOnPayment.emplace_back(0.458233732562325);
            strikes.emplace_back(0.732);   stdDevsOnExpiry.emplace_back(0.45693647278369);    stdDevsOnPayment.emplace_back(0.458467039227467);
            strikes.emplace_back(0.733);   stdDevsOnExpiry.emplace_back(0.45712937172096);    stdDevsOnPayment.emplace_back(0.458671182559467);
            strikes.emplace_back(0.734);   stdDevsOnExpiry.emplace_back(0.457322270658231);   stdDevsOnPayment.emplace_back(0.458907729594959);
            strikes.emplace_back(0.735);   stdDevsOnExpiry.emplace_back(0.457521494150821);   stdDevsOnPayment.emplace_back(0.459079469223467);
            strikes.emplace_back(0.736);   stdDevsOnExpiry.emplace_back(0.457736529031713);   stdDevsOnPayment.emplace_back(0.459299814407213);
            strikes.emplace_back(0.737);   stdDevsOnExpiry.emplace_back(0.457970537578565);   stdDevsOnPayment.emplace_back(0.459481275146768);
            strikes.emplace_back(0.738);   stdDevsOnExpiry.emplace_back(0.458141300572214);   stdDevsOnPayment.emplace_back(0.459701620330514);
            strikes.emplace_back(0.739);   stdDevsOnExpiry.emplace_back(0.458406931895669);   stdDevsOnPayment.emplace_back(0.459931686625307);
            strikes.emplace_back(0.74);    stdDevsOnExpiry.emplace_back(0.458587181722298);   stdDevsOnPayment.emplace_back(0.460116387735212);
            strikes.emplace_back(0.741);   stdDevsOnExpiry.emplace_back(0.458846488490432);   stdDevsOnPayment.emplace_back(0.460317290696863);
            strikes.emplace_back(0.742);   stdDevsOnExpiry.emplace_back(0.459045711983022);   stdDevsOnPayment.emplace_back(0.460534395510259);
            strikes.emplace_back(0.743);   stdDevsOnExpiry.emplace_back(0.459216474976672);   stdDevsOnPayment.emplace_back(0.460757981064355);
            strikes.emplace_back(0.744);   stdDevsOnExpiry.emplace_back(0.459412536191602);   stdDevsOnPayment.emplace_back(0.46098156661845);
            strikes.emplace_back(0.745);   stdDevsOnExpiry.emplace_back(0.459656031571435);   stdDevsOnPayment.emplace_back(0.461156546617307);
            strikes.emplace_back(0.746);   stdDevsOnExpiry.emplace_back(0.459829956842744);   stdDevsOnPayment.emplace_back(0.461370411060354);
            strikes.emplace_back(0.747);   stdDevsOnExpiry.emplace_back(0.460054478556616);   stdDevsOnPayment.emplace_back(0.461545391059211);
            strikes.emplace_back(0.748);   stdDevsOnExpiry.emplace_back(0.460244215216226);   stdDevsOnPayment.emplace_back(0.46179165920575);
            strikes.emplace_back(0.749);   stdDevsOnExpiry.emplace_back(0.460396004543914);   stdDevsOnPayment.emplace_back(0.461966639204608);
            strikes.emplace_back(0.75);    stdDevsOnExpiry.emplace_back(0.460671122700349);   stdDevsOnPayment.emplace_back(0.462212907351147);
            strikes.emplace_back(0.751);   stdDevsOnExpiry.emplace_back(0.460857697082299);   stdDevsOnPayment.emplace_back(0.462378166238956);
            strikes.emplace_back(0.752);   stdDevsOnExpiry.emplace_back(0.461028460075948);   stdDevsOnPayment.emplace_back(0.462598511422702);
            strikes.emplace_back(0.753);   stdDevsOnExpiry.emplace_back(0.461215034457898);   stdDevsOnPayment.emplace_back(0.462792933643655);
            strikes.emplace_back(0.754);   stdDevsOnExpiry.emplace_back(0.461458529837731);   stdDevsOnPayment.emplace_back(0.4630132788274);
            strikes.emplace_back(0.755);   stdDevsOnExpiry.emplace_back(0.461657753330322);   stdDevsOnPayment.emplace_back(0.4632174221594);
            strikes.emplace_back(0.756);   stdDevsOnExpiry.emplace_back(0.461869625933553);   stdDevsOnPayment.emplace_back(0.463411844380352);
            strikes.emplace_back(0.757);   stdDevsOnExpiry.emplace_back(0.462122608146366);   stdDevsOnPayment.emplace_back(0.463632189564098);
            strikes.emplace_back(0.758);   stdDevsOnExpiry.emplace_back(0.462264910641074);   stdDevsOnPayment.emplace_back(0.463820131044352);
            strikes.emplace_back(0.759);   stdDevsOnExpiry.emplace_back(0.462470458688985);   stdDevsOnPayment.emplace_back(0.464014553265304);
            strikes.emplace_back(0.76);    stdDevsOnExpiry.emplace_back(0.462676006736896);   stdDevsOnPayment.emplace_back(0.464254340671146);
            strikes.emplace_back(0.761);   stdDevsOnExpiry.emplace_back(0.462910015283748);   stdDevsOnPayment.emplace_back(0.464432561040352);
            strikes.emplace_back(0.762);   stdDevsOnExpiry.emplace_back(0.463087102832718);   stdDevsOnPayment.emplace_back(0.464646425483399);
            strikes.emplace_back(0.763);   stdDevsOnExpiry.emplace_back(0.463245216715726);   stdDevsOnPayment.emplace_back(0.464850568815399);
            strikes.emplace_back(0.764);   stdDevsOnExpiry.emplace_back(0.463482387540239);   stdDevsOnPayment.emplace_back(0.464999625851462);
            strikes.emplace_back(0.765);   stdDevsOnExpiry.emplace_back(0.46369109786581);    stdDevsOnPayment.emplace_back(0.465281538071843);
            strikes.emplace_back(0.766);   stdDevsOnExpiry.emplace_back(0.463909295024361);   stdDevsOnPayment.emplace_back(0.465443556589303);
            strikes.emplace_back(0.767);   stdDevsOnExpiry.emplace_back(0.46407373346269);    stdDevsOnPayment.emplace_back(0.465650940291652);
            strikes.emplace_back(0.768);   stdDevsOnExpiry.emplace_back(0.464285606065922);   stdDevsOnPayment.emplace_back(0.465861564364351);
            strikes.emplace_back(0.769);   stdDevsOnExpiry.emplace_back(0.464475342725532);   stdDevsOnPayment.emplace_back(0.466036544363208);
            strikes.emplace_back(0.77);    stdDevsOnExpiry.emplace_back(0.464668241662802);   stdDevsOnPayment.emplace_back(0.466218005102763);
            strikes.emplace_back(0.771);   stdDevsOnExpiry.emplace_back(0.464845329211771);   stdDevsOnPayment.emplace_back(0.466483715471398);
            strikes.emplace_back(0.772);   stdDevsOnExpiry.emplace_back(0.465031903593721);   stdDevsOnPayment.emplace_back(0.466632772507461);
            strikes.emplace_back(0.773);   stdDevsOnExpiry.emplace_back(0.465256425307593);   stdDevsOnPayment.emplace_back(0.466807752506318);
            strikes.emplace_back(0.774);   stdDevsOnExpiry.emplace_back(0.465427188301242);   stdDevsOnPayment.emplace_back(0.467060501393556);
            strikes.emplace_back(0.775);   stdDevsOnExpiry.emplace_back(0.465645385459794);   stdDevsOnPayment.emplace_back(0.467229000651715);
            strikes.emplace_back(0.776);   stdDevsOnExpiry.emplace_back(0.465806661620462);   stdDevsOnPayment.emplace_back(0.467416942131968);
            strikes.emplace_back(0.777);   stdDevsOnExpiry.emplace_back(0.466021696501354);   stdDevsOnPayment.emplace_back(0.467591922130826);
            strikes.emplace_back(0.778);   stdDevsOnExpiry.emplace_back(0.466255705048206);   stdDevsOnPayment.emplace_back(0.467815507684921);
            strikes.emplace_back(0.779);   stdDevsOnExpiry.emplace_back(0.466489713595059);   stdDevsOnPayment.emplace_back(0.467996968424476);
            strikes.emplace_back(0.78);    stdDevsOnExpiry.emplace_back(0.466676287977009);   stdDevsOnPayment.emplace_back(0.468188150275079);
            strikes.emplace_back(0.781);   stdDevsOnExpiry.emplace_back(0.466843888692998);   stdDevsOnPayment.emplace_back(0.468395533977428);
            strikes.emplace_back(0.782);   stdDevsOnExpiry.emplace_back(0.466992515743026);   stdDevsOnPayment.emplace_back(0.468576994716984);
            strikes.emplace_back(0.783);   stdDevsOnExpiry.emplace_back(0.467213875179237);   stdDevsOnPayment.emplace_back(0.468823262863523);
            strikes.emplace_back(0.784);   stdDevsOnExpiry.emplace_back(0.467470019669711);   stdDevsOnPayment.emplace_back(0.469017685084475);
            strikes.emplace_back(0.785);   stdDevsOnExpiry.emplace_back(0.467656594051661);   stdDevsOnPayment.emplace_back(0.469228309157174);
            strikes.emplace_back(0.786);   stdDevsOnExpiry.emplace_back(0.467811545657009);   stdDevsOnPayment.emplace_back(0.46944541397057);
            strikes.emplace_back(0.787);   stdDevsOnExpiry.emplace_back(0.467938036763416);   stdDevsOnPayment.emplace_back(0.46956854804384);
            strikes.emplace_back(0.788);   stdDevsOnExpiry.emplace_back(0.468172045310268);   stdDevsOnPayment.emplace_back(0.469788893227586);
            strikes.emplace_back(0.789);   stdDevsOnExpiry.emplace_back(0.468374431080519);   stdDevsOnPayment.emplace_back(0.470035161374125);
            strikes.emplace_back(0.79);    stdDevsOnExpiry.emplace_back(0.468595790516731);   stdDevsOnPayment.emplace_back(0.470200420261935);
            strikes.emplace_back(0.791);   stdDevsOnExpiry.emplace_back(0.468801338564642);   stdDevsOnPayment.emplace_back(0.470407803964284);
            strikes.emplace_back(0.792);   stdDevsOnExpiry.emplace_back(0.468972101558291);   stdDevsOnPayment.emplace_back(0.470514736185807);
            strikes.emplace_back(0.793);   stdDevsOnExpiry.emplace_back(0.46913970227428);    stdDevsOnPayment.emplace_back(0.470770725443395);
            strikes.emplace_back(0.794);   stdDevsOnExpiry.emplace_back(0.469307302990269);   stdDevsOnPayment.emplace_back(0.470974868775394);
            strikes.emplace_back(0.795);   stdDevsOnExpiry.emplace_back(0.469550798370102);   stdDevsOnPayment.emplace_back(0.471127166181807);
            strikes.emplace_back(0.796);   stdDevsOnExpiry.emplace_back(0.46970891225311);    stdDevsOnPayment.emplace_back(0.471337790254505);
            strikes.emplace_back(0.797);   stdDevsOnExpiry.emplace_back(0.469952407632943);   stdDevsOnPayment.emplace_back(0.47156461617895);
            strikes.emplace_back(0.798);   stdDevsOnExpiry.emplace_back(0.470151631125534);   stdDevsOnPayment.emplace_back(0.471733115437108);
            strikes.emplace_back(0.799);   stdDevsOnExpiry.emplace_back(0.47026231084364);    stdDevsOnPayment.emplace_back(0.471937258769108);
            strikes.emplace_back(0.8);     stdDevsOnExpiry.emplace_back(0.470499481668152);   stdDevsOnPayment.emplace_back(0.47205067173133);
            strikes.emplace_back(0.801);   stdDevsOnExpiry.emplace_back(0.470670244661801);   stdDevsOnPayment.emplace_back(0.472283978396473);
            strikes.emplace_back(0.802);   stdDevsOnExpiry.emplace_back(0.470878954987373);   stdDevsOnPayment.emplace_back(0.472442756543584);
            strikes.emplace_back(0.803);   stdDevsOnExpiry.emplace_back(0.4710275820374);     stdDevsOnPayment.emplace_back(0.472737630245361);
            strikes.emplace_back(0.804);   stdDevsOnExpiry.emplace_back(0.471220480974671);   stdDevsOnPayment.emplace_back(0.472880206540726);
            strikes.emplace_back(0.805);   stdDevsOnExpiry.emplace_back(0.471432353577902);   stdDevsOnPayment.emplace_back(0.473081109502377);
            strikes.emplace_back(0.806);   stdDevsOnExpiry.emplace_back(0.471615765682192);   stdDevsOnPayment.emplace_back(0.473246368390186);
            strikes.emplace_back(0.807);   stdDevsOnExpiry.emplace_back(0.471837125118404);   stdDevsOnPayment.emplace_back(0.473450511722186);
            strikes.emplace_back(0.808);   stdDevsOnExpiry.emplace_back(0.472026861778014);   stdDevsOnPayment.emplace_back(0.473625491721043);
            strikes.emplace_back(0.809);   stdDevsOnExpiry.emplace_back(0.472184975661022);   stdDevsOnPayment.emplace_back(0.4738004717199);
            strikes.emplace_back(0.81);    stdDevsOnExpiry.emplace_back(0.472346251821691);   stdDevsOnPayment.emplace_back(0.473991653570503);
            strikes.emplace_back(0.811);   stdDevsOnExpiry.emplace_back(0.472573935813223);   stdDevsOnPayment.emplace_back(0.474179595050757);
            strikes.emplace_back(0.812);   stdDevsOnExpiry.emplace_back(0.47271307603027);    stdDevsOnPayment.emplace_back(0.474377257642059);
            strikes.emplace_back(0.813);   stdDevsOnExpiry.emplace_back(0.472940760021802);   stdDevsOnPayment.emplace_back(0.474513353196725);
            strikes.emplace_back(0.814);   stdDevsOnExpiry.emplace_back(0.473168444013335);   stdDevsOnPayment.emplace_back(0.474730458010122);
            strikes.emplace_back(0.815);   stdDevsOnExpiry.emplace_back(0.473310746508042);   stdDevsOnPayment.emplace_back(0.474931360971772);
            strikes.emplace_back(0.816);   stdDevsOnExpiry.emplace_back(0.473528943666594);   stdDevsOnPayment.emplace_back(0.475148465785169);
            strikes.emplace_back(0.817);   stdDevsOnExpiry.emplace_back(0.473658597050661);   stdDevsOnPayment.emplace_back(0.475359089857867);
            strikes.emplace_back(0.818);   stdDevsOnExpiry.emplace_back(0.473857820543251);   stdDevsOnPayment.emplace_back(0.475485464301486);
            strikes.emplace_back(0.819);   stdDevsOnExpiry.emplace_back(0.474012772148599);   stdDevsOnPayment.emplace_back(0.475702569114883);
            strikes.emplace_back(0.82);    stdDevsOnExpiry.emplace_back(0.474291052582694);   stdDevsOnPayment.emplace_back(0.475906712446883);
            strikes.emplace_back(0.821);   stdDevsOnExpiry.emplace_back(0.474436517355062);   stdDevsOnPayment.emplace_back(0.476046048371899);
            strikes.emplace_back(0.822);   stdDevsOnExpiry.emplace_back(0.474661039068934);   stdDevsOnPayment.emplace_back(0.476269633925994);
            strikes.emplace_back(0.823);   stdDevsOnExpiry.emplace_back(0.474765394231719);   stdDevsOnPayment.emplace_back(0.476399248739962);
            strikes.emplace_back(0.824);   stdDevsOnExpiry.emplace_back(0.475027863277513);   stdDevsOnPayment.emplace_back(0.476626074664406);
            strikes.emplace_back(0.825);   stdDevsOnExpiry.emplace_back(0.475198626271162);   stdDevsOnPayment.emplace_back(0.476807535403962);
            strikes.emplace_back(0.826);   stdDevsOnExpiry.emplace_back(0.475321955099909);   stdDevsOnPayment.emplace_back(0.476953352069676);
            strikes.emplace_back(0.827);   stdDevsOnExpiry.emplace_back(0.475552801369101);   stdDevsOnPayment.emplace_back(0.477144533920279);
            strikes.emplace_back(0.828);   stdDevsOnExpiry.emplace_back(0.475685617030829);   stdDevsOnPayment.emplace_back(0.477355157992977);
            strikes.emplace_back(0.829);   stdDevsOnExpiry.emplace_back(0.47590065191172);    stdDevsOnPayment.emplace_back(0.477539859102882);
            strikes.emplace_back(0.83);    stdDevsOnExpiry.emplace_back(0.476074577183029);   stdDevsOnPayment.emplace_back(0.477718079472088);
            strikes.emplace_back(0.831);   stdDevsOnExpiry.emplace_back(0.476201068289436);   stdDevsOnPayment.emplace_back(0.477876857619199);
            strikes.emplace_back(0.832);   stdDevsOnExpiry.emplace_back(0.47645721277991);    stdDevsOnPayment.emplace_back(0.478097202802945);
            strikes.emplace_back(0.833);   stdDevsOnExpiry.emplace_back(0.476593190719297);   stdDevsOnPayment.emplace_back(0.478281903912849);
            strikes.emplace_back(0.834);   stdDevsOnExpiry.emplace_back(0.476861984320411);   stdDevsOnPayment.emplace_back(0.478489287615198);
            strikes.emplace_back(0.835);   stdDevsOnExpiry.emplace_back(0.476975826316177);   stdDevsOnPayment.emplace_back(0.478589739096024);
            strikes.emplace_back(0.836);   stdDevsOnExpiry.emplace_back(0.477203510307709);   stdDevsOnPayment.emplace_back(0.478803603539071);
            strikes.emplace_back(0.837);   stdDevsOnExpiry.emplace_back(0.477323676858796);   stdDevsOnPayment.emplace_back(0.478991545019325);
            strikes.emplace_back(0.838);   stdDevsOnExpiry.emplace_back(0.477538711739687);   stdDevsOnPayment.emplace_back(0.479166525018182);
            strikes.emplace_back(0.839);   stdDevsOnExpiry.emplace_back(0.477722123843977);   stdDevsOnPayment.emplace_back(0.479399831683325);
            strikes.emplace_back(0.84);    stdDevsOnExpiry.emplace_back(0.477946645557849);   stdDevsOnPayment.emplace_back(0.479503523534499);
            strikes.emplace_back(0.841);   stdDevsOnExpiry.emplace_back(0.478098434885537);   stdDevsOnPayment.emplace_back(0.479688224644404);
            strikes.emplace_back(0.842);   stdDevsOnExpiry.emplace_back(0.478291333822807);   stdDevsOnPayment.emplace_back(0.479979857975832);
            strikes.emplace_back(0.843);   stdDevsOnExpiry.emplace_back(0.478370390764311);   stdDevsOnPayment.emplace_back(0.480073828715959);
            strikes.emplace_back(0.844);   stdDevsOnExpiry.emplace_back(0.478550640590941);   stdDevsOnPayment.emplace_back(0.48027473167761);
            strikes.emplace_back(0.845);   stdDevsOnExpiry.emplace_back(0.478797298248434);   stdDevsOnPayment.emplace_back(0.480427029084022);
            strikes.emplace_back(0.846);   stdDevsOnExpiry.emplace_back(0.478987034908044);   stdDevsOnPayment.emplace_back(0.480624691675324);
            strikes.emplace_back(0.847);   stdDevsOnExpiry.emplace_back(0.479119850569771);   stdDevsOnPayment.emplace_back(0.48080291204453);
            strikes.emplace_back(0.848);   stdDevsOnExpiry.emplace_back(0.479338047728323);   stdDevsOnPayment.emplace_back(0.48096493056199);
            strikes.emplace_back(0.849);   stdDevsOnExpiry.emplace_back(0.479568893997515);   stdDevsOnPayment.emplace_back(0.481162593153292);
            strikes.emplace_back(0.85);    stdDevsOnExpiry.emplace_back(0.479749143824145);   stdDevsOnPayment.emplace_back(0.48137321722599);
            strikes.emplace_back(0.851);   stdDevsOnExpiry.emplace_back(0.47985033670927);    stdDevsOnPayment.emplace_back(0.481509312780656);
            strikes.emplace_back(0.852);   stdDevsOnExpiry.emplace_back(0.480043235646541);   stdDevsOnPayment.emplace_back(0.481719936853355);
            strikes.emplace_back(0.853);   stdDevsOnExpiry.emplace_back(0.48022348547317);    stdDevsOnPayment.emplace_back(0.481875474630116);
            strikes.emplace_back(0.854);   stdDevsOnExpiry.emplace_back(0.480324678358296);   stdDevsOnPayment.emplace_back(0.48210554092491);
            strikes.emplace_back(0.855);   stdDevsOnExpiry.emplace_back(0.48060295879239);    stdDevsOnPayment.emplace_back(0.48226755944237);
            strikes.emplace_back(0.856);   stdDevsOnExpiry.emplace_back(0.480716800788156);   stdDevsOnPayment.emplace_back(0.482432818330179);
            strikes.emplace_back(0.857);   stdDevsOnExpiry.emplace_back(0.480988756666931);   stdDevsOnPayment.emplace_back(0.482607798329037);
            strikes.emplace_back(0.858);   stdDevsOnExpiry.emplace_back(0.481096274107377);   stdDevsOnPayment.emplace_back(0.482773057216846);
            strikes.emplace_back(0.859);   stdDevsOnExpiry.emplace_back(0.481238576602084);   stdDevsOnPayment.emplace_back(0.482977200548846);
            strikes.emplace_back(0.86);    stdDevsOnExpiry.emplace_back(0.481456773760636);   stdDevsOnPayment.emplace_back(0.48316190165875);
            strikes.emplace_back(0.861);   stdDevsOnExpiry.emplace_back(0.481573778034062);   stdDevsOnPayment.emplace_back(0.483278554991322);
            strikes.emplace_back(0.862);   stdDevsOnExpiry.emplace_back(0.481833084802196);   stdDevsOnPayment.emplace_back(0.483495659804718);
            strikes.emplace_back(0.863);   stdDevsOnExpiry.emplace_back(0.481978549574564);   stdDevsOnPayment.emplace_back(0.483647957211131);
            strikes.emplace_back(0.864);   stdDevsOnExpiry.emplace_back(0.482136663457572);   stdDevsOnPayment.emplace_back(0.483819696839639);
            strikes.emplace_back(0.865);   stdDevsOnExpiry.emplace_back(0.482345373783143);   stdDevsOnPayment.emplace_back(0.484001157579194);
            strikes.emplace_back(0.866);   stdDevsOnExpiry.emplace_back(0.482462378056569);   stdDevsOnPayment.emplace_back(0.484208541281543);
            strikes.emplace_back(0.867);   stdDevsOnExpiry.emplace_back(0.482705873436402);   stdDevsOnPayment.emplace_back(0.484289550540273);
            strikes.emplace_back(0.868);   stdDevsOnExpiry.emplace_back(0.482908259206653);   stdDevsOnPayment.emplace_back(0.48454878016821);
            strikes.emplace_back(0.869);   stdDevsOnExpiry.emplace_back(0.483072697644982);   stdDevsOnPayment.emplace_back(0.484769125351956);
            strikes.emplace_back(0.87);    stdDevsOnExpiry.emplace_back(0.483180215085428);   stdDevsOnPayment.emplace_back(0.484911701647321);
            strikes.emplace_back(0.871);   stdDevsOnExpiry.emplace_back(0.483461657797183);   stdDevsOnPayment.emplace_back(0.485086681646178);
            strikes.emplace_back(0.872);   stdDevsOnExpiry.emplace_back(0.483537552461027);   stdDevsOnPayment.emplace_back(0.485219536830495);
            strikes.emplace_back(0.873);   stdDevsOnExpiry.emplace_back(0.483692504066375);   stdDevsOnPayment.emplace_back(0.485446362754939);
            strikes.emplace_back(0.874);   stdDevsOnExpiry.emplace_back(0.483831644283422);   stdDevsOnPayment.emplace_back(0.485559775717161);
            strikes.emplace_back(0.875);   stdDevsOnExpiry.emplace_back(0.484037192331333);   stdDevsOnPayment.emplace_back(0.485760678678812);
            strikes.emplace_back(0.876);   stdDevsOnExpiry.emplace_back(0.484242740379244);   stdDevsOnPayment.emplace_back(0.485987504603256);
            strikes.emplace_back(0.877);   stdDevsOnExpiry.emplace_back(0.484366069207991);   stdDevsOnPayment.emplace_back(0.486104157935828);
            strikes.emplace_back(0.878);   stdDevsOnExpiry.emplace_back(0.484596915477183);   stdDevsOnPayment.emplace_back(0.486282378305034);
            strikes.emplace_back(0.879);   stdDevsOnExpiry.emplace_back(0.484707595195289);   stdDevsOnPayment.emplace_back(0.486489762007383);
            strikes.emplace_back(0.88);    stdDevsOnExpiry.emplace_back(0.484849897689996);   stdDevsOnPayment.emplace_back(0.486612896080653);
            strikes.emplace_back(0.881);   stdDevsOnExpiry.emplace_back(0.48509655534749);    stdDevsOnPayment.emplace_back(0.486758712746367);
            strikes.emplace_back(0.882);   stdDevsOnExpiry.emplace_back(0.485254669230498);   stdDevsOnPayment.emplace_back(0.486956375337668);
            strikes.emplace_back(0.883);   stdDevsOnExpiry.emplace_back(0.485498164610331);   stdDevsOnPayment.emplace_back(0.48715403792897);
            strikes.emplace_back(0.884);   stdDevsOnExpiry.emplace_back(0.485574059274175);   stdDevsOnPayment.emplace_back(0.487290133483636);
            strikes.emplace_back(0.885);   stdDevsOnExpiry.emplace_back(0.485691063547601);   stdDevsOnPayment.emplace_back(0.487416507927255);
            strikes.emplace_back(0.886);   stdDevsOnExpiry.emplace_back(0.485899773873172);   stdDevsOnPayment.emplace_back(0.487610930148207);
            strikes.emplace_back(0.887);   stdDevsOnExpiry.emplace_back(0.486171729751947);   stdDevsOnPayment.emplace_back(0.487782669776715);
            strikes.emplace_back(0.888);   stdDevsOnExpiry.emplace_back(0.486307707691334);   stdDevsOnPayment.emplace_back(0.488009495701159);
            strikes.emplace_back(0.889);   stdDevsOnExpiry.emplace_back(0.486519580294565);   stdDevsOnPayment.emplace_back(0.488165033477921);
            strikes.emplace_back(0.89);    stdDevsOnExpiry.emplace_back(0.486627097735011);   stdDevsOnPayment.emplace_back(0.488314090513985);
            strikes.emplace_back(0.891);   stdDevsOnExpiry.emplace_back(0.486769400229719);   stdDevsOnPayment.emplace_back(0.488476109031445);
            strikes.emplace_back(0.892);   stdDevsOnExpiry.emplace_back(0.486943325501028);   stdDevsOnPayment.emplace_back(0.488631646808207);
            strikes.emplace_back(0.893);   stdDevsOnExpiry.emplace_back(0.487107763939357);   stdDevsOnPayment.emplace_back(0.488829309399508);
            strikes.emplace_back(0.894);   stdDevsOnExpiry.emplace_back(0.487297500598967);   stdDevsOnPayment.emplace_back(0.488962164583825);
            strikes.emplace_back(0.895);   stdDevsOnExpiry.emplace_back(0.487420829427713);   stdDevsOnPayment.emplace_back(0.489172788656524);
            strikes.emplace_back(0.896);   stdDevsOnExpiry.emplace_back(0.487642188863925);   stdDevsOnPayment.emplace_back(0.489318605322238);
            strikes.emplace_back(0.897);   stdDevsOnExpiry.emplace_back(0.487809789579914);   stdDevsOnPayment.emplace_back(0.489493585321095);
            strikes.emplace_back(0.898);   stdDevsOnExpiry.emplace_back(0.487952092074622);   stdDevsOnPayment.emplace_back(0.489655603838555);
            strikes.emplace_back(0.899);   stdDevsOnExpiry.emplace_back(0.48810388140231);    stdDevsOnPayment.emplace_back(0.48987918939265);
            strikes.emplace_back(0.9);     stdDevsOnExpiry.emplace_back(0.488334727671502);   stdDevsOnPayment.emplace_back(0.49000232346592);
            strikes.emplace_back(0.901);   stdDevsOnExpiry.emplace_back(0.48847703016621);    stdDevsOnPayment.emplace_back(0.490232389760713);
            strikes.emplace_back(0.902);   stdDevsOnExpiry.emplace_back(0.488660442270499);   stdDevsOnPayment.emplace_back(0.490349043093285);
            strikes.emplace_back(0.903);   stdDevsOnExpiry.emplace_back(0.488847016652449);   stdDevsOnPayment.emplace_back(0.490475417536904);
            strikes.emplace_back(0.904);   stdDevsOnExpiry.emplace_back(0.488938722704594);   stdDevsOnPayment.emplace_back(0.490699003090999);
            strikes.emplace_back(0.905);   stdDevsOnExpiry.emplace_back(0.48904624014504);    stdDevsOnPayment.emplace_back(0.490838339016015);
            strikes.emplace_back(0.906);   stdDevsOnExpiry.emplace_back(0.489264437303591);   stdDevsOnPayment.emplace_back(0.491026280496268);
            strikes.emplace_back(0.907);   stdDevsOnExpiry.emplace_back(0.489409902075959);   stdDevsOnPayment.emplace_back(0.491175337532332);
            strikes.emplace_back(0.908);   stdDevsOnExpiry.emplace_back(0.489647072900472);   stdDevsOnPayment.emplace_back(0.491392442345728);
            strikes.emplace_back(0.909);   stdDevsOnExpiry.emplace_back(0.489760914896238);   stdDevsOnPayment.emplace_back(0.491515576418998);
            strikes.emplace_back(0.91);    stdDevsOnExpiry.emplace_back(0.48999492344309);    stdDevsOnPayment.emplace_back(0.491658152714363);
            strikes.emplace_back(0.911);   stdDevsOnExpiry.emplace_back(0.49018149782504);    stdDevsOnPayment.emplace_back(0.491872017157411);
            strikes.emplace_back(0.912);   stdDevsOnExpiry.emplace_back(0.490190984658021);   stdDevsOnPayment.emplace_back(0.491975709008585);
            strikes.emplace_back(0.913);   stdDevsOnExpiry.emplace_back(0.490466102814455);   stdDevsOnPayment.emplace_back(0.492170131229537);
            strikes.emplace_back(0.914);   stdDevsOnExpiry.emplace_back(0.490611567586823);   stdDevsOnPayment.emplace_back(0.49236131308014);
            strikes.emplace_back(0.915);   stdDevsOnExpiry.emplace_back(0.490785492858132);   stdDevsOnPayment.emplace_back(0.492529812338299);
            strikes.emplace_back(0.916);   stdDevsOnExpiry.emplace_back(0.490908821686879);   stdDevsOnPayment.emplace_back(0.492701551966807);
            strikes.emplace_back(0.917);   stdDevsOnExpiry.emplace_back(0.491032150515626);   stdDevsOnPayment.emplace_back(0.492831166780775);
            strikes.emplace_back(0.918);   stdDevsOnExpiry.emplace_back(0.491244023118857);   stdDevsOnPayment.emplace_back(0.492983464187188);
            strikes.emplace_back(0.919);   stdDevsOnExpiry.emplace_back(0.491373676502924);   stdDevsOnPayment.emplace_back(0.493194088259886);
            strikes.emplace_back(0.92);    stdDevsOnExpiry.emplace_back(0.491534952663592);   stdDevsOnPayment.emplace_back(0.493245934185473);
            strikes.emplace_back(0.921);   stdDevsOnExpiry.emplace_back(0.491737338433843);   stdDevsOnPayment.emplace_back(0.493557009738997);
            strikes.emplace_back(0.922);   stdDevsOnExpiry.emplace_back(0.491952373314735);   stdDevsOnPayment.emplace_back(0.493641259368076);
            strikes.emplace_back(0.923);   stdDevsOnExpiry.emplace_back(0.49205672847752);    stdDevsOnPayment.emplace_back(0.493845402700076);
            strikes.emplace_back(0.924);   stdDevsOnExpiry.emplace_back(0.492227491471169);   stdDevsOnPayment.emplace_back(0.493978257884393);
            strikes.emplace_back(0.925);   stdDevsOnExpiry.emplace_back(0.492376118521197);   stdDevsOnPayment.emplace_back(0.494130555290806);
            strikes.emplace_back(0.926);   stdDevsOnExpiry.emplace_back(0.492569017458467);   stdDevsOnPayment.emplace_back(0.49431849677106);
            strikes.emplace_back(0.927);   stdDevsOnExpiry.emplace_back(0.492682859454234);   stdDevsOnPayment.emplace_back(0.494457832696075);
            strikes.emplace_back(0.928);   stdDevsOnExpiry.emplace_back(0.492882082946824);   stdDevsOnPayment.emplace_back(0.494665216398424);
            strikes.emplace_back(0.929);   stdDevsOnExpiry.emplace_back(0.49298643810961);    stdDevsOnPayment.emplace_back(0.494749466027504);
            strikes.emplace_back(0.93);    stdDevsOnExpiry.emplace_back(0.493214122101142);   stdDevsOnPayment.emplace_back(0.494986013062996);
            strikes.emplace_back(0.931);   stdDevsOnExpiry.emplace_back(0.493350100040529);   stdDevsOnPayment.emplace_back(0.495105906765916);
            strikes.emplace_back(0.932);   stdDevsOnExpiry.emplace_back(0.493524025311838);   stdDevsOnPayment.emplace_back(0.495329492320011);
            strikes.emplace_back(0.933);   stdDevsOnExpiry.emplace_back(0.493713761971448);   stdDevsOnPayment.emplace_back(0.495446145652583);
            strikes.emplace_back(0.934);   stdDevsOnExpiry.emplace_back(0.493887687242758);   stdDevsOnPayment.emplace_back(0.495562798985154);
            strikes.emplace_back(0.935);   stdDevsOnExpiry.emplace_back(0.494052125681086);   stdDevsOnPayment.emplace_back(0.495763701946804);
            strikes.emplace_back(0.936);   stdDevsOnExpiry.emplace_back(0.494143831733231);   stdDevsOnPayment.emplace_back(0.495948403056709);
            strikes.emplace_back(0.937);   stdDevsOnExpiry.emplace_back(0.494355704336463);   stdDevsOnPayment.emplace_back(0.496084498611376);
            strikes.emplace_back(0.938);   stdDevsOnExpiry.emplace_back(0.494539116440752);   stdDevsOnPayment.emplace_back(0.496282161202677);
            strikes.emplace_back(0.939);   stdDevsOnExpiry.emplace_back(0.494615011104596);   stdDevsOnPayment.emplace_back(0.496437698979439);
            strikes.emplace_back(0.94);    stdDevsOnExpiry.emplace_back(0.494782611820585);   stdDevsOnPayment.emplace_back(0.496599717496899);
            strikes.emplace_back(0.941);   stdDevsOnExpiry.emplace_back(0.494886966983371);   stdDevsOnPayment.emplace_back(0.496790899347502);
            strikes.emplace_back(0.942);   stdDevsOnExpiry.emplace_back(0.495174734250446);   stdDevsOnPayment.emplace_back(0.496878389346931);
            strikes.emplace_back(0.943);   stdDevsOnExpiry.emplace_back(0.495323361300474);   stdDevsOnPayment.emplace_back(0.497037167494042);
            strikes.emplace_back(0.944);   stdDevsOnExpiry.emplace_back(0.495459339239861);   stdDevsOnPayment.emplace_back(0.497176503419057);
            strikes.emplace_back(0.945);   stdDevsOnExpiry.emplace_back(0.495519422515405);   stdDevsOnPayment.emplace_back(0.497338521936518);
            strikes.emplace_back(0.946);   stdDevsOnExpiry.emplace_back(0.495721808285655);   stdDevsOnPayment.emplace_back(0.497562107490613);
            strikes.emplace_back(0.947);   stdDevsOnExpiry.emplace_back(0.495921031778246);   stdDevsOnPayment.emplace_back(0.497701443415628);
            strikes.emplace_back(0.948);   stdDevsOnExpiry.emplace_back(0.496066496550614);   stdDevsOnPayment.emplace_back(0.497837538970295);
            strikes.emplace_back(0.949);   stdDevsOnExpiry.emplace_back(0.496218285878302);   stdDevsOnPayment.emplace_back(0.49801899970985);
            strikes.emplace_back(0.95);    stdDevsOnExpiry.emplace_back(0.496385886594291);   stdDevsOnPayment.emplace_back(0.498219902671501);
            strikes.emplace_back(0.951);   stdDevsOnExpiry.emplace_back(0.496525026811338);   stdDevsOnPayment.emplace_back(0.498268508226739);
            strikes.emplace_back(0.952);   stdDevsOnExpiry.emplace_back(0.496676816139026);   stdDevsOnPayment.emplace_back(0.498498574521533);
            strikes.emplace_back(0.953);   stdDevsOnExpiry.emplace_back(0.496793820412452);   stdDevsOnPayment.emplace_back(0.498628189335501);
            strikes.emplace_back(0.954);   stdDevsOnExpiry.emplace_back(0.497049964902926);   stdDevsOnPayment.emplace_back(0.498812890445405);
            strikes.emplace_back(0.955);   stdDevsOnExpiry.emplace_back(0.497189105119973);   stdDevsOnPayment.emplace_back(0.499013793407056);
            strikes.emplace_back(0.956);   stdDevsOnExpiry.emplace_back(0.497347219002982);   stdDevsOnPayment.emplace_back(0.499172571554167);
            strikes.emplace_back(0.957);   stdDevsOnExpiry.emplace_back(0.497476872387049);   stdDevsOnPayment.emplace_back(0.49924385970185);
            strikes.emplace_back(0.958);   stdDevsOnExpiry.emplace_back(0.4976792581573);     stdDevsOnPayment.emplace_back(0.499503089329786);
            strikes.emplace_back(0.959);   stdDevsOnExpiry.emplace_back(0.497770964209444);   stdDevsOnPayment.emplace_back(0.499593819699564);
            strikes.emplace_back(0.96);    stdDevsOnExpiry.emplace_back(0.497998648200977);   stdDevsOnPayment.emplace_back(0.499726674883881);
            strikes.emplace_back(0.961);   stdDevsOnExpiry.emplace_back(0.498042920088219);   stdDevsOnPayment.emplace_back(0.499937298956579);
            strikes.emplace_back(0.962);   stdDevsOnExpiry.emplace_back(0.498204196248888);   stdDevsOnPayment.emplace_back(0.500031269696706);
            strikes.emplace_back(0.963);   stdDevsOnExpiry.emplace_back(0.49842871796276);    stdDevsOnPayment.emplace_back(0.500203009325214);
            strikes.emplace_back(0.964);   stdDevsOnExpiry.emplace_back(0.498517261737244);   stdDevsOnPayment.emplace_back(0.500345585620579);
            strikes.emplace_back(0.965);   stdDevsOnExpiry.emplace_back(0.498672213342592);   stdDevsOnPayment.emplace_back(0.500436315990357);
            strikes.emplace_back(0.966);   stdDevsOnExpiry.emplace_back(0.49881767811496);    stdDevsOnPayment.emplace_back(0.500689064877594);
            strikes.emplace_back(0.967);   stdDevsOnExpiry.emplace_back(0.499076984883094);   stdDevsOnPayment.emplace_back(0.500792756728769);
            strikes.emplace_back(0.968);   stdDevsOnExpiry.emplace_back(0.499200313711841);   stdDevsOnPayment.emplace_back(0.501006621171816);
            strikes.emplace_back(0.969);   stdDevsOnExpiry.emplace_back(0.499415348592732);   stdDevsOnPayment.emplace_back(0.50107466894915);
            strikes.emplace_back(0.97);    stdDevsOnExpiry.emplace_back(0.499513379200197);   stdDevsOnPayment.emplace_back(0.501343619688134);
            strikes.emplace_back(0.971);   stdDevsOnExpiry.emplace_back(0.499665168527885);   stdDevsOnPayment.emplace_back(0.501431109687562);
            strikes.emplace_back(0.972);   stdDevsOnExpiry.emplace_back(0.499769523690671);   stdDevsOnPayment.emplace_back(0.50160284931607);
            strikes.emplace_back(0.973);   stdDevsOnExpiry.emplace_back(0.499911826185378);   stdDevsOnPayment.emplace_back(0.50172598338934);
            strikes.emplace_back(0.974);   stdDevsOnExpiry.emplace_back(0.500098400567328);   stdDevsOnPayment.emplace_back(0.501894482647498);
            strikes.emplace_back(0.975);   stdDevsOnExpiry.emplace_back(0.500281812671618);   stdDevsOnPayment.emplace_back(0.502037058942863);
            strikes.emplace_back(0.976);   stdDevsOnExpiry.emplace_back(0.500370356446103);   stdDevsOnPayment.emplace_back(0.502263884867308);
            strikes.emplace_back(0.977);   stdDevsOnExpiry.emplace_back(0.50050000983017);    stdDevsOnPayment.emplace_back(0.502344894126038);
            strikes.emplace_back(0.978);   stdDevsOnExpiry.emplace_back(0.500702395600421);   stdDevsOnPayment.emplace_back(0.502506912643498);
            strikes.emplace_back(0.979);   stdDevsOnExpiry.emplace_back(0.500844698095128);   stdDevsOnPayment.emplace_back(0.502568479680133);
            strikes.emplace_back(0.98);    stdDevsOnExpiry.emplace_back(0.501091355752621);   stdDevsOnPayment.emplace_back(0.502827709308069);
            strikes.emplace_back(0.981);   stdDevsOnExpiry.emplace_back(0.501129303084543);   stdDevsOnPayment.emplace_back(0.502980006714482);
            strikes.emplace_back(0.982);   stdDevsOnExpiry.emplace_back(0.501338013410114);   stdDevsOnPayment.emplace_back(0.503096660047053);
            strikes.emplace_back(0.983);   stdDevsOnExpiry.emplace_back(0.501540399180365);   stdDevsOnPayment.emplace_back(0.503375331897085);
            strikes.emplace_back(0.984);   stdDevsOnExpiry.emplace_back(0.501644754343151);   stdDevsOnPayment.emplace_back(0.503423937452323);
            strikes.emplace_back(0.985);   stdDevsOnExpiry.emplace_back(0.501720649006995);   stdDevsOnPayment.emplace_back(0.503621600043624);
            strikes.emplace_back(0.986);   stdDevsOnExpiry.emplace_back(0.501954657553847);   stdDevsOnPayment.emplace_back(0.503806301153529);
            strikes.emplace_back(0.987);   stdDevsOnExpiry.emplace_back(0.502096960048555);   stdDevsOnPayment.emplace_back(0.503929435226798);
            strikes.emplace_back(0.988);   stdDevsOnExpiry.emplace_back(0.502210802044321);   stdDevsOnPayment.emplace_back(0.504078492262862);
            strikes.emplace_back(0.989);   stdDevsOnExpiry.emplace_back(0.502406863259251);   stdDevsOnPayment.emplace_back(0.504175703373338);
            strikes.emplace_back(0.99);    stdDevsOnExpiry.emplace_back(0.502596599918861);   stdDevsOnPayment.emplace_back(0.50437012559429);
            strikes.emplace_back(0.991);   stdDevsOnExpiry.emplace_back(0.502659845472065);   stdDevsOnPayment.emplace_back(0.504447894482671);
            strikes.emplace_back(0.992);   stdDevsOnExpiry.emplace_back(0.502811634799753);   stdDevsOnPayment.emplace_back(0.504713604851306);
            strikes.emplace_back(0.993);   stdDevsOnExpiry.emplace_back(0.502903340851898);   stdDevsOnPayment.emplace_back(0.504713604851306);
            strikes.emplace_back(0.994);   stdDevsOnExpiry.emplace_back(0.503083590678527);   stdDevsOnPayment.emplace_back(0.504924228924004);
            strikes.emplace_back(0.995);   stdDevsOnExpiry.emplace_back(0.503304950114739);   stdDevsOnPayment.emplace_back(0.505095968552512);
            strikes.emplace_back(0.996);   stdDevsOnExpiry.emplace_back(0.503384007056243);   stdDevsOnPayment.emplace_back(0.505189939292639);
            strikes.emplace_back(0.997);   stdDevsOnExpiry.emplace_back(0.503529471828611);   stdDevsOnPayment.emplace_back(0.505390842254289);
            strikes.emplace_back(0.998);   stdDevsOnExpiry.emplace_back(0.503712883932901);   stdDevsOnPayment.emplace_back(0.505611187438035);
            strikes.emplace_back(0.999);   stdDevsOnExpiry.emplace_back(0.503858348705269);   stdDevsOnPayment.emplace_back(0.505695437067115);
            strikes.emplace_back(1);       stdDevsOnExpiry.emplace_back(0.504029111698918);   stdDevsOnPayment.emplace_back(0.505818571140384);
            strikes.emplace_back(1.001);   stdDevsOnExpiry.emplace_back(0.504127142306383);   stdDevsOnPayment.emplace_back(0.505964387806099);
            strikes.emplace_back(1.002);   stdDevsOnExpiry.emplace_back(0.504301067577692);   stdDevsOnPayment.emplace_back(0.506139367804955);

            //Create smiles on Expiry Date
            smilesOnExpiry = std::vector<std::shared_ptr<SmileSection> >();
            smilesOnExpiry.emplace_back(std::shared_ptr<SmileSection>(
                new FlatSmileSection(startDate, flatVol, rangeCouponDayCount)));
            Real dummyAtmLevel = 0;
            smilesOnExpiry.emplace_back(std::shared_ptr<SmileSection>(new
                InterpolatedSmileSection<Linear>(startDate,
                                                 strikes,
                                                 stdDevsOnExpiry,
                                                 dummyAtmLevel,
                                                 rangeCouponDayCount)));
            //smilesOnExpiry.emplace_back(
            //    swaptionVolatilityStructures_[0]->smileSection(startDate,
            //                                                   Period(6, Months)));
            //Create smiles on Payment Date
            smilesOnPayment = std::vector<std::shared_ptr<SmileSection> >();
            smilesOnPayment.emplace_back(std::shared_ptr<SmileSection>(
                new FlatSmileSection(endDate, flatVol, rangeCouponDayCount)));
            smilesOnPayment.emplace_back(std::shared_ptr<SmileSection>(new
                InterpolatedSmileSection<Linear>(endDate,
                                                 strikes,
                                                 stdDevsOnPayment,
                                                 dummyAtmLevel,
                                                 rangeCouponDayCount)));
            //vars.smilesOnPayment.emplace_back(
            //    swaptionVolatilityStructures_[0]->smileSection(vars.endDate,
            //                                                   Period(6, Months)));

            QL_REQUIRE(smilesOnExpiry.size()==smilesOnPayment.size(),
                       "smilesOnExpiry.size()!=smilesOnPayment.size()");
        }

        CommonVars() {

            //General Settings
            calendar = TARGET();
            today = Date(39147); // 6 Mar 2007
            Settings::instance().evaluationDate() = today;
            settlement = today;
            //create Yield Curve
            createYieldCurve();
            referenceDate = termStructure->referenceDate();
            // Ibor index
            iborIndex =
                std::shared_ptr<IborIndex>(new Euribor6M(termStructure));

            // create Volatility Structures
            flatVol = 0.1;
            createVolatilityStructures();

            // Range Accrual valuation
            gearing = 1.0;
            spread = 0.0;
            infiniteLowerStrike = 1.e-9;
            infiniteUpperStrike = 1.0;
            correlation = 1.0;

            startDate = Date(42800); //6 Mar 2017
            endDate = Date(42984);   //6 Sep 2017
            paymentDate = endDate;   //6 Sep 2017
            fixingDays = 2;
            rangeCouponDayCount = iborIndex->dayCounter();

            // observations schedule
            observationsConvention = ModifiedFollowing;
            observationsFrequency = Daily;//Monthly;
            observationSchedule = std::shared_ptr<Schedule>(new
                                Schedule(startDate,endDate,
                                         Period(observationsFrequency),
                                         calendar,observationsConvention,
                                         observationsConvention,
                                         DateGeneration::Forward, false));
            // Range accrual pricers properties
            byCallSpread = std::vector<bool>();
            byCallSpread.emplace_back(true);
            byCallSpread.emplace_back(false);

            std::vector<Rate> strikes;

            //Create smiles sections
            createSmileSections();

            //test parameters
            rateTolerance = 2.0e-8;
            priceTolerance = 2.0e-4;
        }
    };

}


//******************************************************************************************//
//******************************************************************************************//
TEST_CASE( "RangeAccrual_InfiniteRange", "[RangeAccrual]" )  {

    INFO("Testing infinite range accrual floaters...");

    CommonVars vars;

    //Coupon
    RangeAccrualFloatersCoupon coupon(vars.paymentDate,
                                      1.0,
                                      vars.iborIndex,
                                      vars.startDate,
                                      vars.endDate,
                                      vars.fixingDays,
                                      vars.rangeCouponDayCount,
                                      vars.gearing, vars.spread,
                                      vars.startDate, vars.endDate,
                                      vars.observationSchedule,
                                      vars.infiniteLowerStrike,
                                      vars.infiniteUpperStrike);

    Date fixingDate = coupon.fixingDate();

    for (Size z = 0; z < vars.smilesOnPayment.size(); z++) {
        for (Size i = 0; i < vars.byCallSpread.size(); i++){
            std::shared_ptr<RangeAccrualPricer> bgmPricer(new
                RangeAccrualPricerByBgm(vars.correlation,
                                        vars.smilesOnExpiry[z],
                                        vars.smilesOnPayment[z],
                                        true,
                                        vars.byCallSpread[i]));

                coupon.setPricer(bgmPricer);

                //Computation
                const Rate rate = coupon.rate();
                const Rate indexfixing = vars.iborIndex->fixing(fixingDate);
                const Rate difference =  rate-indexfixing;

                if (std::fabs(difference) > vars.rateTolerance) {
                    FAIL_CHECK("\n" <<
                                "i:\t" << i << "\n"
                                "fixingDate:\t" << fixingDate << "\n"
                                "startDate:\t" << vars.startDate << "\n"
                                "range accrual rate:\t" << io::rate(rate) << "\n"
                                "index fixing:\t" << io::rate(indexfixing) << "\n"
                                "difference:\t" << io::rate(difference) << "\n"
                                "tolerance: \t" << io::rate(vars.rateTolerance));
                }
        }
    }
}

TEST_CASE( "RangeAccrual_PriceMonotonicityWithRespectToLowerStrike", "[RangeAccrual]" ) {

    INFO(
            "Testing price monotonicity with respect to the lower strike...");

    CommonVars vars;

    for (Size z = 0; z < vars.smilesOnPayment.size(); z++) {
        for (Size i = 0; i < vars.byCallSpread.size(); i++){
            std::shared_ptr<RangeAccrualPricer> bgmPricer(new
                RangeAccrualPricerByBgm(vars.correlation,
                                        vars.smilesOnExpiry[z],
                                        vars.smilesOnPayment[z],
                                        true,
                                        vars.byCallSpread[i]));

            Real effectiveLowerStrike;
            Real previousPrice = 100.;

            for (Size k = 1; k < 100; k++){
                effectiveLowerStrike = 0.005 + k*0.001;
                RangeAccrualFloatersCoupon coupon(
                                            vars.paymentDate,
                                            1.,
                                            vars.iborIndex,
                                            vars.startDate,
                                            vars.endDate,
                                            vars.fixingDays,
                                            vars.rangeCouponDayCount,
                                            vars.gearing, vars.spread,
                                            vars.startDate, vars.endDate,
                                            vars.observationSchedule,
                                            effectiveLowerStrike,
                                            vars.infiniteUpperStrike);

                coupon.setPricer(bgmPricer);

                //Computation
                const Rate price = coupon.price(vars.termStructure);

                if (previousPrice <= price) {
                    FAIL_CHECK("\n" <<
                                "i:\t" << i << "\n"
                                "k:\t" << k << "\n"
                                "Price at lower strike\t" << effectiveLowerStrike-0.001 <<
                                ": \t" << previousPrice << "\n"
                                "Price at lower strike\t" << effectiveLowerStrike <<
                                ": \t" << price << "\n");
                }
                previousPrice = price;
            }
        }
    }
}


TEST_CASE( "RangeAccrual_PriceMonotonicityWithRespectToUpperStrike", "[RangeAccrual]" ) {

    INFO(
            "Testing price monotonicity with respect to the upper strike...");

    CommonVars vars;

    for (Size z = 0; z < vars.smilesOnPayment.size(); z++) {
        for (Size i = 0; i < vars.byCallSpread.size(); i++){
            std::shared_ptr<RangeAccrualPricer> bgmPricer(new
                RangeAccrualPricerByBgm(vars.correlation,
                                        vars.smilesOnExpiry[z],
                                        vars.smilesOnPayment[z],
                                        true,
                                        vars.byCallSpread[i]));

            Real effectiveUpperStrike;
            Real previousPrice = 0.;

            for (Size k = 1; k < 95; k++){
                effectiveUpperStrike = 0.006 + k*0.001;
                RangeAccrualFloatersCoupon coupon(
                                            vars.paymentDate,
                                            1.,
                                            vars.iborIndex,
                                            vars.startDate,
                                            vars.endDate,
                                            vars.fixingDays,
                                            vars.rangeCouponDayCount,
                                            vars.gearing, vars.spread,
                                            vars.startDate, vars.endDate,
                                            vars.observationSchedule,
                                            .004,
                                            effectiveUpperStrike);

                coupon.setPricer(bgmPricer);

                //Computation
                const Rate price = coupon.price(vars.termStructure);

                if (previousPrice > price) {
                    FAIL_CHECK("\n" <<
                                "i:\t" << i << "\n"
                                "k:\t" << k << "\n"
                                "Price at upper strike\t" << effectiveUpperStrike-0.001 <<
                                ": \t" << previousPrice << "\n"
                                "Price at upper strike\t" << effectiveUpperStrike <<
                                ": \t" << price << "\n");
                }
                previousPrice = price;
            }
        }
    }
}
