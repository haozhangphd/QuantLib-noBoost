/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Jose Aparicio

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
#include <ql/experimental/credit/randomdefaultlatentmodel.hpp>
#include <ql/termstructures/credit/flathazardrate.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <chrono>
#include <memory>
#include <algorithm>
#include <functional>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace QuantLib;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

    Integer sessionId() { return 0; }

}
#endif


/* This sample code shows basic usage of a Latent variable model.
   The data and correlation problem presented is the same as in:
     'Modelling Dependent Defaults: Asset Correlations Are Not Enough!'
     Frey R., A. J. McNeil and M. A. Nyfeler RiskLab publications March 2001
*/
int main(int, char* []) {

    try {

        std::chrono::time_point<std::chrono::steady_clock> startT =  std::chrono::steady_clock::now();
	std::cout << std::endl;

        Calendar calendar = TARGET();
        Date todaysDate(19, March, 2014);
        // must be a business day
        todaysDate = calendar.adjust(todaysDate);

        Settings::instance().evaluationDate() = todaysDate;


        /* --------------------------------------------------------------
                        SET UP BASKET PORTFOLIO
        -------------------------------------------------------------- */
        // build curves and issuers into a basket of three names
        std::vector<Real> hazardRates(3, -std::log(1.-0.01));
        std::vector<std::string> names;
        for(Size i=0; i<hazardRates.size(); i++)
            names.emplace_back(std::string("Acme") + 
                to_string(i));
        std::vector<Handle<DefaultProbabilityTermStructure> > defTS;
        for(Size i=0; i<hazardRates.size(); i++)
            defTS.emplace_back(Handle<DefaultProbabilityTermStructure>(
                std::make_shared<FlatHazardRate>(0, TARGET(), hazardRates[i], 
                    Actual365Fixed())));
        std::vector<Issuer> issuers;
        for(Size i=0; i<hazardRates.size(); i++) {
            std::vector<QuantLib::Issuer::key_curve_pair> curves(1,
                std::make_pair(NorthAmericaCorpDefaultKey(
                    EURCurrency(), QuantLib::SeniorSec,
                    Period(), 1. // amount threshold
                    ), defTS[i]));
            issuers.emplace_back(Issuer(curves));
        }

        std::shared_ptr<Pool> thePool = std::make_shared<Pool>();
        for(Size i=0; i<hazardRates.size(); i++)
            thePool->add(names[i], issuers[i], NorthAmericaCorpDefaultKey(
                    EURCurrency(), QuantLib::SeniorSec, Period(), 1.));

        std::vector<DefaultProbKey> defaultKeys(hazardRates.size(), 
            NorthAmericaCorpDefaultKey(EURCurrency(), SeniorSec, Period(), 1.));
        // Recoveries are irrelevant in this example but must be given as the 
        //   lib stands.
        std::vector<std::shared_ptr<RecoveryRateModel> > rrModels(
            hazardRates.size(), std::make_shared<ConstantRecoveryModel>(
            ConstantRecoveryModel(0.5, SeniorSec)));
        std::shared_ptr<Basket> theBskt = std::make_shared<Basket>(
            todaysDate, names, std::vector<Real>(hazardRates.size(), 100.), 
            thePool);
        /* --------------------------------------------------------------
                        SET UP JOINT DEFAULT EVENT LATENT MODELS
        -------------------------------------------------------------- */
        // Latent model factors, corresponds to the first entry in Table1 of the
        //   publication mentioned. It is a single factor model
        std::vector<std::vector<Real> > fctrsWeights(hazardRates.size(), 
            std::vector<Real>(1, std::sqrt(0.1)));
        // --- Default Latent models -------------------------------------
        // Gaussian integrable joint default model:
        std::shared_ptr<GaussianDefProbLM> lmG(new 
            GaussianDefProbLM(fctrsWeights, 
            LatentModelIntegrationType::GaussianQuadrature,
			GaussianCopulaPolicy::initTraits() // otherwise gcc screams
			));
        // Define StudentT copula
        // this is as far as we can be from the Gaussian, 2 T_3 factors:
        std::vector<Integer> ordersT(2, 3);
        TCopulaPolicy::initTraits iniT;
        iniT.tOrders = ordersT;
        // StudentT integrable joint default model:
        std::shared_ptr<TDefProbLM> lmT(new TDefProbLM(fctrsWeights, 
            // LatentModelIntegrationType::GaussianQuadrature,
            LatentModelIntegrationType::Trapezoid,
            iniT));

        // --- Default Loss models ----------------------------------------
        // Gaussian random joint default model:
        Size numSimulations = 100000;
        // Size numCoresUsed = 4;
        // Sobol, many cores
        std::shared_ptr<DefaultLossModel> rdlmG(
            std::make_shared<RandomDefaultLM<GaussianCopulaPolicy> >(lmG, 
                std::vector<Real>(), numSimulations, 1.e-6, 2863311530));
        // StudentT random joint default model:
        std::shared_ptr<DefaultLossModel> rdlmT(
            std::make_shared<RandomDefaultLM<TCopulaPolicy> >(lmT, 
            std::vector<Real>(), numSimulations, 1.e-6, 2863311530));

        /* --------------------------------------------------------------
                        DUMP SOME RESULTS
        -------------------------------------------------------------- */
        /* Default correlations in a T copula should be below those of the 
        gaussian for the same factors.
        The calculations on the MC show dispersion on both copulas (thats
        ok) and too large values with very large dispersions on the T case.
        Computations are ok, within the dispersion, for the gaussian; compare
        with the direct integration in both cases.
        However the T does converge to the gaussian value for large value of
        the parameters.
        */
        Date calcDate(TARGET().advance(Settings::instance().evaluationDate(), 
            Period(120, Months)));
        std::vector<Probability> probEventsTLatent, probEventsTRandLoss;
        std::vector<Probability> probEventsGLatent, probEventsGRandLoss;
        //
        lmT->resetBasket(theBskt);
        for(Size numEvts=0; numEvts <=theBskt->size(); numEvts++) {
            probEventsTLatent.emplace_back(lmT->probAtLeastNEvents(numEvts, 
                calcDate));
         }
        //
        theBskt->setLossModel(rdlmT);
        for(Size numEvts=0; numEvts <=theBskt->size(); numEvts++) {
            probEventsTRandLoss.emplace_back(theBskt->probAtLeastNEvents(numEvts, 
                calcDate));
         }
        //
        lmG->resetBasket(theBskt);
        for(Size numEvts=0; numEvts <=theBskt->size(); numEvts++) {
            probEventsGLatent.emplace_back(lmG->probAtLeastNEvents(numEvts, 
                calcDate));
         }
        //
        theBskt->setLossModel(rdlmG);
        for(Size numEvts=0; numEvts <=theBskt->size(); numEvts++) {
            probEventsGRandLoss.emplace_back(theBskt->probAtLeastNEvents(numEvts, 
                calcDate));
         }

        Date correlDate = TARGET().advance(
            Settings::instance().evaluationDate(), Period(12, Months));
        std::vector<std::vector<Real> > correlsGlm, correlsTlm, correlsGrand, 
            correlsTrand;
        //
        lmT->resetBasket(theBskt);
        for(Size iName1=0; iName1 <theBskt->size(); iName1++) {
            std::vector<Real> tmp;
            for(Size iName2=0; iName2 <theBskt->size(); iName2++)
                tmp.emplace_back(lmT->defaultCorrelation(correlDate, 
                    iName1, iName2));
            correlsTlm.emplace_back(tmp);
        }
        //
        theBskt->setLossModel(rdlmT);
        for(Size iName1=0; iName1 <theBskt->size(); iName1++) {
            std::vector<Real> tmp;
            for(Size iName2=0; iName2 <theBskt->size(); iName2++)
                tmp.emplace_back(theBskt->defaultCorrelation(correlDate, 
                    iName1, iName2));
            correlsTrand.emplace_back(tmp);
        }
        //
        lmG->resetBasket(theBskt);
        for(Size iName1=0; iName1 <theBskt->size(); iName1++) {
            std::vector<Real> tmp;
            for(Size iName2=0; iName2 <theBskt->size(); iName2++)
                tmp.emplace_back(lmG->defaultCorrelation(correlDate, 
                    iName1, iName2));
            correlsGlm.emplace_back(tmp);
        }
        //
        theBskt->setLossModel(rdlmG);
        for(Size iName1=0; iName1 <theBskt->size(); iName1++) {
            std::vector<Real> tmp;
            for(Size iName2=0; iName2 <theBskt->size(); iName2++)
                tmp.emplace_back(theBskt->defaultCorrelation(correlDate, 
                    iName1, iName2));
            correlsGrand.emplace_back(tmp);
        }

        std::cout << 
            " Gaussian versus T prob of extreme event (random and integrable)-" 
            << std::endl;
        for(Size numEvts=0; numEvts <=theBskt->size(); numEvts++) {
            std::cout << "-Prob of " << numEvts << " events... " <<
                probEventsGLatent[numEvts] << " ** " <<
                probEventsTLatent[numEvts] << " ** " << 
                probEventsGRandLoss[numEvts]<< " ** " <<
                probEventsTRandLoss[numEvts] 
            << std::endl;
        }

        cout << endl;
        cout << "-- Default correlations G,T,GRand,TRand--" << endl;
        cout << "-----------------------------------------" << endl;
        for(Size iName1=0; iName1 <theBskt->size(); iName1++) {
            for(Size iName2=0; iName2 <theBskt->size(); iName2++)
                cout << 
                    correlsGlm[iName1][iName2] << " , ";
            cout << endl;
        }
        cout << endl;
        for(Size iName1=0; iName1 <theBskt->size(); iName1++) {
            for(Size iName2=0; iName2 <theBskt->size(); iName2++)
                cout << 
                    correlsTlm[iName1][iName2] << " , ";
            ;
                cout << endl;
        }
        cout << endl;
        for(Size iName1=0; iName1 <theBskt->size(); iName1++) {
            for(Size iName2=0; iName2 <theBskt->size(); iName2++)
                cout << 
                    correlsGrand[iName1][iName2] << " , ";
            cout << endl;
        }
        cout << endl;
        for(Size iName1=0; iName1 <theBskt->size(); iName1++) {
            for(Size iName2=0; iName2 <theBskt->size(); iName2++)
                cout << 
                    correlsTrand[iName1][iName2] << " , ";
            ;
                cout << endl;
        }



        std::chrono::time_point<std::chrono::steady_clock> endT = std::chrono::steady_clock::now();
        Real seconds = static_cast<double>((endT - startT).count()) / 1.0e9;
        Integer hours = Integer(seconds/3600);
        seconds -= hours * 3600;
        Integer minutes = Integer(seconds/60);
        seconds -= minutes * 60;
        cout << "Run completed in ";
        if (hours > 0)
            cout << hours << " h ";
        if (hours > 0 || minutes > 0)
            cout << minutes << " m ";
        cout << fixed << setprecision(0)
             << seconds << " s" << endl;

        return 0;
    } catch (exception& e) {
        cerr << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "unknown error" << endl;
        return 1;
    }
}

