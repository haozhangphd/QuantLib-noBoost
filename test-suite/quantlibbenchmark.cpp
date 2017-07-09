/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2008, 2010 Klaus Spanderen

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


/*
 QuantLib Benchmark Suite

 Measures the performance of a preselected set of numerically intensive
 test cases. The overall QuantLib Benchmark Index is given by the average
 performance in mflops.

 The number of floating point operations of a given test case was measured
 using the perfex library, http://user.it.uu.se/~mikpe/linux/perfctr
 and PAPI, http://icl.cs.utk.edu/papi

 Example results: 1. i7 4702HQ@2.2GHz       :6524.9 mflops
 	 	 	 	  2. i7 870@2.93GHz         :4759.2 mflops
                  3. Core2 Q9300@2.5Ghz     :2272.6 mflops
                  4. Core2 Q6600@2.4Ghz     :1984.0 mflops
                  5. i3 540@3.1Ghz          :1755.3 mflops
                  6. Core2 Dual@2.0Ghz      : 835.9 mflops
                  7. Athlon 64 X2 4400+     : 824.2 mflops
                  8. Core2 Dual@2.0Ghz      : 754.1 mflops
                  9. Pentium4 Dual@2.8Ghz   : 423.8 mflops
                 10. Raspberry Pi3@1.2GHz   : 309.2 mflops
                 11. Pentium4@3.0Ghz        : 266.3 mflops
                 12. PentiumIII@1.1Ghz      : 146.2 mflops
                 13. Alpha 2xEV68@833Mhz    : 184.6 mflops
                 14. Raspberry Pi ARM@700Mhz:  28.3 mflops
                 15. Strong ARM@206Mhz      :   1.4 mflops

 Remarks: OS: Linux, static libs
  1. g++-4.8.1 -O3 -ffast-math -march=core-avx2
      Remark: eight processes
  2. gcc-4.6.3, -O3 -ffast-math -mfpmath=sse,387 -march=corei7
      Remark: eight processes
  3. icc-11.0,  -gcc-version=420 -fast -fp-model fast=2 -ipo-jobs2
      Remark: four processes
  4. icc-11.0,  -gcc-version=420 -fast -fp-model fast=2 -ipo-jobs2
      Remark: four processes
  5. gcc-4.4.5, -O3 -ffast-math -mfpmath=sse,387 -msse4.2 -march=core2
      Remark: four processes
  6. icc-11.0,  -gcc-version=420 -fast -fp-model fast=2 -ipo-jobs2
      Remark: two processes
  7. icc-11.0,  -gcc-version=420 -xSSSE3 -O3 -ipo -no-prec-div -static
                -fp-model fast=2 -ipo-jobs2, Remark: two processes
  8. gcc-4.2.1, -O3 -ffast-math -mfpmath=sse,387 -msse3 -funroll-all-loops
      Remark: two processes
  9. gcc-4.0.1, -O3 -march=pentium4 -ffast-math
      -mfpmath=sse,387 -msse2 -funroll-all-loops, Remark: two processes
 10. gcc-4.9.2  -O2, Remark: four processes
 11. gcc-4.0.1, -O3 -march=pentium4 -ffast-math
                -mfpmath=sse,387 -msse2 -funroll-all-loops
 12. gcc-4.1.1, -O3 -march=pentium3 -ffast-math
                -mfpmath=sse,387 -msse -funroll-all-loops
 13. gcc-3.3.5, -O3 -mcpu=e67 -funroll-all-loops, Remark: two processes
 14. gcc-4.6.3, -O3
 15. gcc-3.4.3, -O2 -g on a Zaurus PDA

  This benchmark is derived from quantlibtestsuite.cpp. Please see the
  copyrights therein.
*/
#define CATCH_CONFIG_RUNNER

#include "catch.hpp"

#include <ql/types.hpp>
#include <ql/version.hpp>
#include <iostream>
#include <iomanip>
#include <list>
#include <string>
#include <chrono>

/* PAPI code
#include <stdio.h
#include <papi.h>
*/

/* uncomment the following lines to unmask floating-point exceptions.
   See http://www.wilmott.com/messageview.cfm?catid=10&threadid=9481
*/
// #  include <float.h>
//   namespace { unsigned int u = _controlfp(_EM_INEXACT, _MCW_EM); }

#include "utilities.hpp"

namespace {

    class Benchmark {
    public:
        Benchmark(const std::string &name, double mflop)
                : name_(name), mflop_(mflop) {
        }

        double getMflop() const {
            return mflop_;
        }

        std::string getName() const {
            return name_;
        }

    private:
        const std::string name_;
        const double mflop_; // total number of mega floating
        // point operations (not per sec!)
    };

    std::list<double> runTimes;
    std::list<Benchmark> bm;
    std::chrono::time_point<std::chrono::steady_clock> startT;
    std::chrono::time_point<std::chrono::steady_clock> endT;

    /* PAPI code
    float real_time, proc_time, mflops;
    long_long lflop, flop=0;
    */

    void startTimer() {
        startT = std::chrono::steady_clock::now();

        /* PAPI code
        lflop = flop;
        PAPI_flops(&real_time, &proc_time, &flop, &mflops);
        */
    }

    void stopTimer() {
        endT = std::chrono::steady_clock::now();
        runTimes.emplace_back(static_cast<double>((endT - startT).count()) / 1.0e9);

        /* PAPI code
        PAPI_flops(&real_time, &proc_time, &flop, &mflops);
        printf("Real_time: %f Proc_time: %f Total mflop: %f\n",
               real_time, proc_time, (flop-lflop)/1e6);
        */
    }

    void printResults() {
        std::string header = "Benchmark Suite "
                "QuantLib " QL_VERSION;

        std::cout << std::endl
                  << std::string(75, '-') << std::endl;
        std::cout << header << std::endl;
        std::cout << std::string(75, '-')
                  << std::endl << std::endl;

        double sum = 0;
        std::list<double>::const_iterator iterT = runTimes.begin();
        std::list<Benchmark>::const_iterator iterBM = bm.begin();

        while (iterT != runTimes.end()) {
            const double mflopsPerSec = iterBM->getMflop() / (*iterT);
            std::cout << iterBM->getName()
                      << std::string(59 - iterBM->getName().length(), ' ') << ":"
                      << std::fixed << std::setw(8) << std::setprecision(1)
                      << mflopsPerSec
                      << " mflops" << std::endl;

            sum += mflopsPerSec;
            ++iterT;
            ++iterBM;
        }
        std::cout << std::string(75, '-') << std::endl
                  << "QuantLib Benchmark Index                                   :"
                  << std::fixed << std::setw(8) << std::setprecision(1)
                  << sum / runTimes.size()
                  << " mflops" << std::endl;
    }
}

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {
    Integer sessionId() { return 0; }
}
#endif

int main() {

    bm.emplace_back(Benchmark("AmericanOption_FdAmericanGreeks", 518.31));
    bm.emplace_back(Benchmark("AmericanOption_FdShoutGreeks", 546.58));
    bm.emplace_back(Benchmark("AsianOption_MCDiscreteArithmeticAveragePrice", 5186.13));
    bm.emplace_back(Benchmark("BarrierOption_BabsiriValues", 880.8));
    bm.emplace_back(Benchmark("BasketOption_EuroTwoValues", 340.04));
    bm.emplace_back(Benchmark("BasketOption_TavellaValues", 933.80));
    bm.emplace_back(Benchmark("BasketOption_OddSamples", 642.46));
    bm.emplace_back(Benchmark("BatesModel_DAXCalibration", 1993.35));
    bm.emplace_back(Benchmark("ConvertibleBond_Bond", 159.85));
    bm.emplace_back(Benchmark("DigitalOption_MCCashAtHit", 995.87));
    bm.emplace_back(Benchmark("DividendOption_FdEuropeanGreeks", 949.52));
    bm.emplace_back(Benchmark("DividendOption_FdAmericanGreeks", 1113.74));
    bm.emplace_back(Benchmark("EuropeanOption_McEngines", 1988.63));
    bm.emplace_back(Benchmark("EuropeanOption_ImpliedVol", 131.51));
    bm.emplace_back(Benchmark("EuropeanOption_FdEngines", 148.43));
    bm.emplace_back(Benchmark("EuropeanOption_PriceCurve", 414.76));
    bm.emplace_back(Benchmark("FdHeston_FdmHestonAmerican", 234.21));
    bm.emplace_back(Benchmark("HestonModel_DAXCalibration", 555.19));
    bm.emplace_back(Benchmark("Interpolation_SabrInterpolation", 2266.06));
    bm.emplace_back(Benchmark("JumpDiffusion_Greeks", 433.77));
    bm.emplace_back(Benchmark("MarketModelCms_MultiStepCmSwapsAndSwaptions",
                              11497.73));
    bm.emplace_back(Benchmark("MarketModelSmm_MultiStepCoterminalSwapsAndSwaptions",
                              11244.95));
    bm.emplace_back(Benchmark("QuantoOption_ForwardGreeks", 90.98));
    bm.emplace_back(Benchmark("LowDiscrepancy_MersenneTwisterDiscrepancy", 951.98));
    bm.emplace_back(Benchmark("RiskStatistics_Results", 300.28));
    bm.emplace_back(Benchmark("ShortRateModel_Swaps", 454.73));


    Catch::Session session; // There must be exactly one instance
    Catch::ConfigData configData;

    for (std::list<Benchmark>::const_iterator iter = bm.begin();
         iter != bm.end(); ++iter) {
        configData.testsOrTags = std::vector<std::string>{iter->getName()};
        session.useConfigData(configData);
#ifndef _MSC_VER
        session.configData().outputFilename= std::string{"/dev/null"};
#else
        session.configData().outputFilename = std::string{ "NUL" };
#endif
        startTimer();
        session.run();
        stopTimer();
    }

    printResults();


}
