/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 Ferdinando Ametrano
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003 StatPro Italia srl
 Copyright (C) 2005 Gary Kennedy

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

#include <ql/math/distributions/bivariatenormaldistribution.hpp>
#include <ql/math/integrals/gaussianquadratures.hpp>

namespace QuantLib {

    // Drezner 1978

    const long double BivariateCumulativeNormalDistributionDr78::x_[] = {
        0.24840615,
        0.39233107,
        0.21141819,
        0.03324666,
        0.00082485334
    };

    const long double BivariateCumulativeNormalDistributionDr78::y_[] = {
        0.10024215,
        0.48281397,
        1.06094980,
        1.77972940,
        2.66976040000
    };

    BivariateCumulativeNormalDistributionDr78::
    BivariateCumulativeNormalDistributionDr78(long double rho)
    : rho_(rho), rho2_(rho*rho) {

        QL_REQUIRE(rho>=-1.0,
                   "rho must be >= -1.0 (" << rho << " not allowed)");
        QL_REQUIRE(rho<=1.0,
                   "rho must be <= 1.0 (" << rho << " not allowed)");
    }

    long double BivariateCumulativeNormalDistributionDr78::operator()(long double a,
                                                               long double b) const {

        CumulativeNormalDistribution cumNormalDist;
        long double CumNormDistA = cumNormalDist(a);
        long double CumNormDistB = cumNormalDist(b);
        long double MaxCumNormDistAB = std::max(CumNormDistA, CumNormDistB);
        long double MinCumNormDistAB = std::min(CumNormDistA, CumNormDistB);

        if (1.0-MaxCumNormDistAB<1e-15)
            return MinCumNormDistAB;

        if (MinCumNormDistAB<1e-15)
            return MinCumNormDistAB;

        long double a1 = a / std::sqrt(2.0 * (1.0 - rho2_));
        long double b1 = b / std::sqrt(2.0 * (1.0 - rho2_));

        long double result=-1.0;

        if (a<=0.0 && b<=0 && rho_<=0) {
            long double sum=0.0;
            for (Size i=0; i<5; i++) {
                for (Size j=0;j<5; j++) {
                    sum += x_[i]*x_[j]*
                        std::exp(a1*(2.0*y_[i]-a1)+b1*(2.0*y_[j]-b1)
                                 +2.0*rho_*(y_[i]-a1)*(y_[j]-b1));
                }
            }
            result = std::sqrt(1.0 - rho2_)/M_PI*sum;
        } else if (a<=0 && b>=0 && rho_>=0) {
            BivariateCumulativeNormalDistributionDr78 bivCumNormalDist(-rho_);
            result= CumNormDistA - bivCumNormalDist(a, -b);
        } else if (a>=0.0 && b<=0.0 && rho_>=0.0) {
            BivariateCumulativeNormalDistributionDr78 bivCumNormalDist(-rho_);
            result= CumNormDistB - bivCumNormalDist(-a, b);
        } else if (a>=0.0 && b>=0.0 && rho_<=0.0) {
            result= CumNormDistA + CumNormDistB -1.0 + (*this)(-a, -b);
        } else if (a*b*rho_>0.0) {
            long double rho1 = (rho_*a-b)*(a>0.0 ? 1.0: -1.0)/
                std::sqrt(a*a-2.0*rho_*a*b+b*b);
            BivariateCumulativeNormalDistributionDr78 bivCumNormalDist(rho1);

            long double rho2 = (rho_*b-a)*(b>0.0 ? 1.0: -1.0)/
                std::sqrt(a*a-2.0*rho_*a*b+b*b);
            BivariateCumulativeNormalDistributionDr78 CBND2(rho2);

            long double delta = (1.0-(a>0.0 ? 1.0: -1.0)*(b>0.0 ? 1.0: -1.0))/4.0;

            result= bivCumNormalDist(a, 0.0) + CBND2(b, 0.0) - delta;
        } else {
            QL_FAIL("case not handled");
        }

        return result;
    }


    // West 2004

    namespace {

        class eqn3 { /* Relates to eqn3 Genz 2004 */
          public:
            eqn3(long double h, long double k, long double asr) {
                hk_ = h * k;
                hs_  = (h * h + k * k) / 2;
                asr_ = asr;
            }
            long double operator()(long double x) const {
                long double sn = std::sin(asr_ * (-x + 1) * 0.5);
                return std::exp((sn * hk_ - hs_) / (1.0 - sn * sn));
            }
          private:
            long double hk_, asr_, hs_;
        };

        class eqn6 { /* Relates to eqn6 Genz 2004 */
          public:
            eqn6(long double a, long double c, long double d, long double bs, long double hk)
            : a_(a), c_(c), d_(d), bs_(bs), hk_(hk) {}
            long double operator()(long double x) const {
                long double xs = a_ * (-x + 1);
                xs = std::fabs(xs*xs);
                long double rs = std::sqrt(1 - xs);
                long double asr = -(bs_ / xs + hk_) / 2;
                if (asr > -100.0) {
                    return (a_ * std::exp(asr) *
                            (std::exp(-hk_ * (1 - rs) / (2 * (1 + rs))) / rs -
                             (1 + c_ * xs * (1 + d_ * xs))));
                } else {
                    return 0.0;
                }
            }
          private:
            long double a_, c_, d_, bs_, hk_;
        };

    }

    BivariateCumulativeNormalDistributionWe04DP::
    BivariateCumulativeNormalDistributionWe04DP(long double rho)
    : correlation_(rho) {

        QL_REQUIRE(rho>=-1.0,
                   "rho must be >= -1.0 (" << rho << " not allowed)");
        QL_REQUIRE(rho<=1.0,
                   "rho must be <= 1.0 (" << rho << " not allowed)");
    }


    long double BivariateCumulativeNormalDistributionWe04DP::operator()(
                                                       long double x, long double y) const {

        /* The implementation is described at section 2.4 "Hybrid
           Numerical Integration Algorithms" of "Numerical Computation
           of Rectangular Bivariate an Trivariate Normal and t
           Probabilities", Genz (2004), Statistics and Computing 14,
           151-160. (available at
           www.sci.wsu.edu/math/faculty/henz/homepage)

           The Gauss-Legendre quadrature have been extracted to
           TabulatedGaussLegendre (x,w zero-based)

           Tthe functions ot be integrated numerically have been moved
           to classes eqn3 and eqn6

           Change some magic numbers to M_PI */

        TabulatedGaussLegendre gaussLegendreQuad(20);
        if (std::fabs(correlation_) < 0.3) {
            gaussLegendreQuad.order(6);
        } else if (std::fabs(correlation_) < 0.75) {
            gaussLegendreQuad.order(12);
        }

        long double h = -x;
        long double k = -y;
        long double hk = h * k;
        long double BVN = 0.0;

        if (std::fabs(correlation_) < 0.925)
        {
            if (std::fabs(correlation_) > 0)
            {
                long double asr = std::asin(correlation_);
                eqn3 f(h,k,asr);
                BVN = gaussLegendreQuad(f);
                BVN *= asr * (0.25 / M_PI);
            }
            BVN += cumnorm_(-h) * cumnorm_(-k);
        }
        else
        {
            if (correlation_ < 0)
            {
                k *= -1;
                hk *= -1;
            }
            if (std::fabs(correlation_) < 1)
            {
                long double Ass = (1 - correlation_) * (1 + correlation_);
                long double a = std::sqrt(Ass);
                long double bs = (h-k)*(h-k);
                long double c = (4 - hk) / 8;
                long double d = (12 - hk) / 16;
                long double asr = -(bs / Ass + hk) / 2;
                if (asr > -100)
                {
                    BVN = a * std::exp(asr) *
                        (1 - c * (bs - Ass) * (1 - d * bs / 5) / 3 +
                         c * d * Ass * Ass / 5);
                }
                if (-hk < 100)
                {
                    long double B = std::sqrt(bs);
                    BVN -= std::exp(-hk / 2) * 2.506628274631 *
                        cumnorm_(-B / a) * B *
                        (1 - c * bs * (1 - d * bs / 5) / 3);
                }
                a /= 2;
                eqn6 f(a,c,d,bs,hk);
                BVN += gaussLegendreQuad(f);
                BVN /= (-2.0 * M_PI);
            }

            if (correlation_ > 0) {
                BVN += cumnorm_(-std::max(h, k));
            } else {
                BVN *= -1;
                if (k > h) {
                    // evaluate cumnorm where it is most precise, that
                    // is in the lower tail because of double accuracy
                    // around 0.0 vs around 1.0
                    if (h >= 0) {
                        BVN += cumnorm_(-h) - cumnorm_(-k);
                    } else {
                        BVN += cumnorm_(k) - cumnorm_(h);
                    }
                }
            }
        }
        return BVN;
    }

}
