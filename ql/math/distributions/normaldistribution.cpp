/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2002, 2003 Ferdinando Ametrano
 Copyright (C) 2008 StatPro Italia srl
 Copyright (C) 2010 Kakhkhor Abdijalilov

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

#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/comparison.hpp>
#include <cmath>

#include <random>

namespace QuantLib {

    long double CumulativeNormalDistribution::operator()(long double z) const {
        //QL_REQUIRE(!(z >= average_ && 2.0*average_-z > average_),
        //           "not a real number. ");

        z = (z - average_) / sigma_;

        long double result = 0.5 * ( 1.0 + std::erf( z * M_SQRT1_2l ));

        if (result<=1e-8) { //todo: investigate the threshold level
            // Asymptotic expansion for very negative z following (26.2.12)
            // on page 408 in M. Abramowitz and A. Stegun,
            // Pocketbook of Mathematical Functions, ISBN 3-87144818-4.
            long double sum=1.0, zsqr=z*z, i=1.0, g=1.0, x, y,
                    a=QL_MAX_REAL, lasta;
            do {
                lasta=a;
                x = (4.0*i-3.0)/zsqr;
                y = x*((4.0*i-1)/zsqr);
                a = g*(x-y);
                sum -= a;
                g *= y;
                ++i;
                a = std::fabs(a);
           } while (lasta>a && a>=std::fabs(sum*QL_EPSILON));
            result = -gaussian_(z)/z*sum;
        }
        return result;
    }

    #if !defined(QL_PATCH_SOLARIS)
    const CumulativeNormalDistribution InverseCumulativeNormal::f_;
    #endif

    // Coefficients for the rational approximation.
    const long double InverseCumulativeNormal::a1_ = -3.969683028665376e+01;
    const long double InverseCumulativeNormal::a2_ =  2.209460984245205e+02;
    const long double InverseCumulativeNormal::a3_ = -2.759285104469687e+02;
    const long double InverseCumulativeNormal::a4_ =  1.383577518672690e+02;
    const long double InverseCumulativeNormal::a5_ = -3.066479806614716e+01;
    const long double InverseCumulativeNormal::a6_ =  2.506628277459239e+00;

    const long double InverseCumulativeNormal::b1_ = -5.447609879822406e+01;
    const long double InverseCumulativeNormal::b2_ =  1.615858368580409e+02;
    const long double InverseCumulativeNormal::b3_ = -1.556989798598866e+02;
    const long double InverseCumulativeNormal::b4_ =  6.680131188771972e+01;
    const long double InverseCumulativeNormal::b5_ = -1.328068155288572e+01;

    const long double InverseCumulativeNormal::c1_ = -7.784894002430293e-03;
    const long double InverseCumulativeNormal::c2_ = -3.223964580411365e-01;
    const long double InverseCumulativeNormal::c3_ = -2.400758277161838e+00;
    const long double InverseCumulativeNormal::c4_ = -2.549732539343734e+00;
    const long double InverseCumulativeNormal::c5_ =  4.374664141464968e+00;
    const long double InverseCumulativeNormal::c6_ =  2.938163982698783e+00;

    const long double InverseCumulativeNormal::d1_ =  7.784695709041462e-03;
    const long double InverseCumulativeNormal::d2_ =  3.224671290700398e-01;
    const long double InverseCumulativeNormal::d3_ =  2.445134137142996e+00;
    const long double InverseCumulativeNormal::d4_ =  3.754408661907416e+00;

    // Limits of the approximation regions
    const long double InverseCumulativeNormal::x_low_ = 0.02425;
    const long double InverseCumulativeNormal::x_high_= 1.0 - x_low_;

    long double InverseCumulativeNormal::tail_value(long double x) {
        if (x <= 0.0 || x >= 1.0) {
            // try to recover if due to numerical error
            if (close_enough(x, 1.0)) {
                return QL_MAX_REAL; // largest value available
            } else if (std::fabs(x) < QL_EPSILON) {
                return QL_MIN_REAL; // largest negative value available
            } else {
                QL_FAIL("InverseCumulativeNormal(" << x
                        << ") undefined: must be 0 < x < 1");
            }
        }

        long double z;
        if (x < x_low_) {
            // Rational approximation for the lower region 0<x<u_low
            z = std::sqrt(-2.0*std::log(x));
            z = (((((c1_*z+c2_)*z+c3_)*z+c4_)*z+c5_)*z+c6_) /
                ((((d1_*z+d2_)*z+d3_)*z+d4_)*z+1.0);
        } else {
            // Rational approximation for the upper region u_high<x<1
            z = std::sqrt(-2.0*std::log(1.0-x));
            z = -(((((c1_*z+c2_)*z+c3_)*z+c4_)*z+c5_)*z+c6_) /
                ((((d1_*z+d2_)*z+d3_)*z+d4_)*z+1.0);
        }

        return z;
    }

    const long double MoroInverseCumulativeNormal::a0_ =  2.50662823884;
    const long double MoroInverseCumulativeNormal::a1_ =-18.61500062529;
    const long double MoroInverseCumulativeNormal::a2_ = 41.39119773534;
    const long double MoroInverseCumulativeNormal::a3_ =-25.44106049637;

    const long double MoroInverseCumulativeNormal::b0_ = -8.47351093090;
    const long double MoroInverseCumulativeNormal::b1_ = 23.08336743743;
    const long double MoroInverseCumulativeNormal::b2_ =-21.06224101826;
    const long double MoroInverseCumulativeNormal::b3_ =  3.13082909833;

    const long double MoroInverseCumulativeNormal::c0_ = 0.3374754822726147;
    const long double MoroInverseCumulativeNormal::c1_ = 0.9761690190917186;
    const long double MoroInverseCumulativeNormal::c2_ = 0.1607979714918209;
    const long double MoroInverseCumulativeNormal::c3_ = 0.0276438810333863;
    const long double MoroInverseCumulativeNormal::c4_ = 0.0038405729373609;
    const long double MoroInverseCumulativeNormal::c5_ = 0.0003951896511919;
    const long double MoroInverseCumulativeNormal::c6_ = 0.0000321767881768;
    const long double MoroInverseCumulativeNormal::c7_ = 0.0000002888167364;
    const long double MoroInverseCumulativeNormal::c8_ = 0.0000003960315187;

    long double MoroInverseCumulativeNormal::operator()(long double x) const {
        QL_REQUIRE(x > 0.0 && x < 1.0,
                   "MoroInverseCumulativeNormal(" << x
                   << ") undefined: must be 0<x<1");

        long double result;
        long double temp=x-0.5;

        if (std::fabs(temp) < 0.42) {
            // Beasley and Springer, 1977
            result=temp*temp;
            result=temp*
                (((a3_*result+a2_)*result+a1_)*result+a0_) /
                ((((b3_*result+b2_)*result+b1_)*result+b0_)*result+1.0);
        } else {
            // improved approximation for the tail (Moro 1995)
            if (x<0.5)
                result = x;
            else
                result=1.0-x;
            result = std::log(-std::log(result));
            result = c0_+result*(c1_+result*(c2_+result*(c3_+result*
                                   (c4_+result*(c5_+result*(c6_+result*
                                                       (c7_+result*c8_)))))));
            if (x<0.5)
                result=-result;
        }

        return average_ + result*sigma_;
    }

    MaddockInverseCumulativeNormal::MaddockInverseCumulativeNormal(
        long double average, long double sigma)
    : average_(average), sigma_(sigma) {}

    long double MaddockInverseCumulativeNormal::operator()(long double x) const {

        x *= 2;

        if (x <= 0.5) {
            const long double Y = 0.0891314744949340820313;
            const std::array<long double, 8> P = {
                     -0.000508781949658280665617,
                     -0.00836874819741736770379,
                     0.0334806625409744615033,
                     -0.0126926147662974029034,
                     -0.0365637971411762664006,
                     0.0219878681111168899165,
                     0.00822687874676915743155,
                     -0.00538772965071242932965 };
            const std::array<long double, 10> Q = {
                     1.0,
                     -0.970005043303290640362,
                     -1.56574558234175846809,
                     1.56221558398423026363,
                     0.662328840472002992063,
                     -0.71228902341542847553,
                     -0.0527396382340099713954,
                     0.0795283687341571680018,
                     -0.00233393759374190016776,
                     0.000886216390456424707504};
            long double g = x * (x + 10);
            long double r =
                    (((((((P[7] * x + P[6]) * x + P[5]) * x + P[4]) * x + P[3]) * x + P[2]) * x + P[1]) * x + P[0]) /
                    (((((((((Q[9] * x + Q[8]) * x + Q[7]) * x + Q[6]) * x + Q[5]) * x + Q[4]) * x + Q[3]) * x + Q[2]) * x + Q[1]) * x + Q[0]);
            return sigma_ * std::sqrt(2.0) * (g * Y + g * r) + average_;
        } else if (x <= 0.75) {
            static const float Y = 2.249481201171875f;
            const std::array<long double, 9> P = {
                     -0.202433508355938759655,
                     0.105264680699391713268,
                     8.37050328343119927838,
                     17.6447298408374015486,
                     -18.8510648058714251895,
                     -44.6382324441786960818,
                     17.445385985570866523,
                     21.1294655448340526258,
                     -3.67192254707729348546 };
            const std::array<long double, 9> Q = {
                     1.0,
                     6.24264124854247537712,
                     3.9713437953343869095,
                     -28.6608180499800029974,
                     -20.1432634680485188801,
                     48.5609213108739935468,
                     10.8268667355460159008,
                     -22.6436933413139721736,
                     1.72114765761200282724 };

            long double g = sqrt(-2 * log(1 - x));
            long double xs = 1 - x - 0.25;
            long double r =
                    ((((((((P[8] * xs + P[7]) * xs + P[6]) * xs + P[5]) * xs + P[4]) * xs + P[3]) * xs + P[2]) * xs + P[1]) * xs + P[0]) /
                    ((((((((Q[8] * xs + Q[7]) * xs + Q[6]) * xs + Q[5]) * xs + Q[4]) * xs + Q[3]) * xs + Q[2]) * xs + Q[1]) * xs + Q[0]);
            return  sigma_ * std::sqrt(2.0) * g / (Y + r) + average_;
        } else {
            long double z = sqrt(-log(1 - x));
            if (z < 3) {
                // Max error found: 1.089051e-20
                static const float Y = 0.807220458984375f;
                const std::array<long double, 11> P = {
                        -0.131102781679951906451,
                        -0.163794047193317060787,
                        0.117030156341995252019,
                        0.387079738972604337464,
                        0.337785538912035898924,
                        0.142869534408157156766,
                        0.0290157910005329060432,
                        0.00214558995388805277169,
                        -0.679465575181126350155e-6,
                        0.285225331782217055858e-7,
                        -0.681149956853776992068e-9};
                const std::array<long double, 8> Q = {
                         1.0,
                         3.46625407242567245975,
                         5.38168345707006855425,
                         4.77846592945843778382,
                         2.59301921623620271374,
                         0.848854343457902036425,
                         0.152264338295331783612,
                         0.01105924229346489121 };

                long double xs = z - 1.125;
                long double R =
                   ((((((((((P[10] * xs + P[9]) * xs + P[8]) * xs + P[7]) * xs + P[6]) * xs + P[5]) * xs + P[4]) * xs + P[3]) * xs + P[2]) * xs + P[1]) * xs + P[0]) /
                   (((((((Q[7] * xs + Q[6]) * xs + Q[5]) * xs + Q[4]) * xs + Q[3]) * xs + Q[2]) * xs + Q[1]) * xs + Q[0]);
                return sigma_ * std::sqrt(2.0) * (Y * z + R * z) + average_;
            } else if (z < 6) {
                // Max error found: 8.389174e-21
                static const float Y = 0.93995571136474609375f;
                const std::array<long double, 9> P = {
                         -0.0350353787183177984712,
                         -0.00222426529213447927281,
                         0.0185573306514231072324,
                         0.00950804701325919603619,
                         0.00187123492819559223345,
                         0.000157544617424960554631,
                         0.460469890584317994083e-5,
                         -0.230404776911882601748e-9,
                         0.266339227425782031962e-11 };
                const std::array<long double, 7> Q = {
                         1.0,
                         1.3653349817554063097,
                         0.762059164553623404043,
                         0.220091105764131249824,
                         0.0341589143670947727934,
                         0.00263861676657015992959,
                         0.764675292302794483503e-4 };

                long double xs = z - 3;
                long double R =
                   ((((((((P[8] * xs + P[7]) * xs + P[6]) * xs + P[5]) * xs + P[4]) * xs + P[3]) * xs + P[2]) * xs + P[1]) * xs + P[0]) /
                   ((((((Q[6] * xs + Q[5]) * xs + Q[4]) * xs + Q[3]) * xs + Q[2]) * xs + Q[1]) * xs + Q[0]);
                return sigma_ * std::sqrt(2.0) * (Y * z + R * z) + average_;
            } else if (z < 18) {
                // Max error found: 1.481312e-19
                static const float Y = 0.98362827301025390625f;
                const std::array<long double, 9> P = {
                         -0.0167431005076633737133,
                         -0.00112951438745580278863,
                         0.00105628862152492910091,
                         0.000209386317487588078668,
                         0.149624783758342370182e-4,
                         0.449696789927706453732e-6,
                         0.462596163522878599135e-8,
                         -0.281128735628831791805e-13,
                         0.99055709973310326855e-16 };

                const std::array<long double, 7> Q = {
                         1.0,
                         0.591429344886417493481,
                         0.138151865749083321638,
                         0.0160746087093676504695,
                         0.000964011807005165528527,
                         0.275335474764726041141e-4,
                         0.282243172016108031869e-6};

                long double xs = z - 6;
                long double R =
                   ((((((((P[8] * xs + P[7]) * xs + P[6]) * xs + P[5]) * xs + P[4]) * xs + P[3]) * xs + P[2]) * xs + P[1]) * xs + P[0]) /
                   ((((((Q[6] * xs + Q[5]) * xs + Q[4]) * xs + Q[3]) * xs + Q[2]) * xs + Q[1]) * xs + Q[0]);
                return sigma_ * std::sqrt(2.0) * (Y * z + R * z) + average_;
            } else if (z < 44) {
                // Max error found: 5.697761e-20
                static const float Y = 0.99714565277099609375f;
                const std::array<long double, 8> P = {
                         -0.0024978212791898131227,
                         -0.779190719229053954292e-5,
                         0.254723037413027451751e-4,
                         0.162397777342510920873e-5,
                         0.396341011304801168516e-7,
                         0.411632831190944208473e-9,
                         0.145596286718675035587e-11,
                         -0.116765012397184275695e-17 };
                const std::array<long double, 7> Q = {
                         1.0,
                         0.207123112214422517181,
                         0.0169410838120975906478,
                         0.000690538265622684595676,
                         0.145007359818232637924e-4,
                         0.144437756628144157666e-6,
                         0.509761276599778486139e-9};

                long double xs = z - 18;
                long double R =
                   (((((((P[7] * xs + P[6]) * xs + P[5]) * xs + P[4]) * xs + P[3]) * xs + P[2]) * xs + P[1]) * xs + P[0]) /
                   ((((((Q[6] * xs + Q[5]) * xs + Q[4]) * xs + Q[3]) * xs + Q[2]) * xs + Q[1]) * xs + Q[0]);
                return sigma_ * std::sqrt(2.0) * (Y * z + R * z) + average_;
            } else {
                // Max error found: 1.279746e-20
                static const float Y = 0.99941349029541015625f;
                const std::array<long double, 8> P = {
                         -0.000539042911019078575891,
                         -0.28398759004727721098e-6,
                         0.899465114892291446442e-6,
                         0.229345859265920864296e-7,
                         0.225561444863500149219e-9,
                         0.947846627503022684216e-12,
                         0.135880130108924861008e-14,
                         -0.348890393399948882918e-21 };
                const std::array<long double, 7> Q = {
                        1.0,
                        0.0845746234001899436914,
                        0.00282092984726264681981,
                        0.468292921940894236786e-4,
                        0.399968812193862100054e-6,
                        0.161809290887904476097e-8,
                        0.231558608310259605225e-11};
                long double xs = z - 44;
                long double R =
                   (((((((P[7] * xs + P[6]) * xs + P[5]) * xs + P[4]) * xs + P[3]) * xs + P[2]) * xs + P[1]) * xs + P[0]) /
                   ((((((Q[6] * xs + Q[5]) * xs + Q[4]) * xs + Q[3]) * xs + Q[2]) * xs + Q[1]) * xs + Q[0]);
                return sigma_ * std::sqrt(2.0) * (Y * z + R * z) + average_;
            }
        }
    }


    MaddockCumulativeNormal::MaddockCumulativeNormal(
        long double average, long double sigma)
    : average_(average), sigma_(sigma) {}

    long double MaddockCumulativeNormal::operator()(long double x) const {
        long double diff = (x - average_) / (sigma_ * std::sqrt(2.0));
        return std::erfc(diff) / 2;
    }
}
