/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Klaus Spanderen

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

/*! \file modifiedbessel.cpp
    \brief modified Bessel functions of first and second kind
*/

#include <ql/math/modifiedbessel.hpp>
#include <ql/math/distributions/gammadistribution.hpp>

#include <cmath>
#include <limits>
#include <array>

namespace {
    static const std::array<long double, 7> c1 = {-1.142022680371168e0, 6.5165112670737e-3,
                                                  3.087090173086e-4, -3.4706269649e-6, 6.9437664e-9, 3.67795e-11,
                                                  -1.356e-13};
    static const std::array<long double, 8> c2 = {1.843740587300905e0, -7.68528408447867e-2,
                                                  1.2719271366546e-3, -4.9717367042e-6, -3.31261198e-8, 2.423096e-10,
                                                  -1.702e-13, -1.49e-15};

// Chebyshev series.
    template<size_t size>
    inline long double chebev(const std::array<long double, size> &c, long double x) {
        long double sv;
        long double d = 0.0, dd = 0.0;
        for (size_t j = size - 1; j > 0; --j) {
            sv = d;
            d = 2.0 * x * d - dd + c[j];
            dd = sv;
        }
        return x * d - dd + 0.5 * c[0];
    }
// Mostly adapted from `Numerical Recipes` 3rd edition by Press et al.
    long double besselik_impl(const long double nu, const long double x, bool compute_i) {
        QL_REQUIRE(x > 0.0 && nu >= 0.0, "bad arguments in modified Bessel function");

        const int MAXIT = 10000;
        const long double EPS = std::numeric_limits<long double>::epsilon();
        const long double FPMIN = std::numeric_limits<long double>::min() / EPS;
        const long double XMIN = 2.0, PI = 3.141592653589793238462643383279502884L;
        long double a, a1, delh, dels, e, f, fact, fact2, ff,
                gammi, gampl, p, pimu, q, q1, q2, qnew, ril, ril1, rimu, ripl,
                ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2, xx;
        int i;
        int nl = static_cast<int>(nu + 0.5);
        long double xmu = nu - nl;
        long double xmu2 = xmu * xmu;
        long double h = std::max(nu / x, FPMIN) ;
        long double b = 2.0 * nu / x;
        long double d = 0.0;
        long double c = h;
        for (i = 0; i < MAXIT; i++) {
            b += 2.0 / x;
            d = 1.0 / (b + d);
            c = b + 1.0 / c;
            h = c * d * h;
            if (std::abs(c * d - 1.0) <= EPS) break;
        }
        QL_REQUIRE(i < MAXIT, "x too large in modified Bessel function");
        ril = FPMIN;
        ripl = h * ril;
        ril1 = ril;
        fact = nu / x;
        for (int l = nl - 1; l >= 0; --l) {
            ritemp = fact * ril + ripl;
            fact -= 1.0 / x;
            ripl = fact * ritemp + ril;
            ril = ritemp;
        }
        f = ripl / ril;
        if (x < XMIN) {
            x2 = 0.5 * x;
            pimu = PI * xmu;
            fact = (std::abs(pimu) < EPS ? 1.0 : pimu / std::sin(pimu));
            d = -std::log(x2);
            e = xmu * d;
            fact2 = (std::abs(e) < EPS ? 1.0 : std::sinh(e) / e);
            xx = 8.0 * xmu * xmu - 1.0;
            long double gam1 = chebev(c1, xx);
            long double gam2 = chebev(c2, xx);
            gampl = gam2 - xmu * gam1;
            gammi = gam2 + xmu * gam1;
            ff = fact * (gam1 * std::cosh(e) + gam2 * fact2 * d);
            sum = ff;
            e = std::exp(e);
            p = 0.5 * e / gampl;
            q = 0.5 / (e * gammi);
            c = 1.0;
            d = x2 * x2;
            sum1 = p;
            for (i = 1; i <= MAXIT; i++) {
                ff = (i * ff + p + q) / (i * i - xmu2);
                c *= (d / i);
                p /= (i - xmu);
                q /= (i + xmu);
                sum += c * ff;
                sum1 += c * (p - i * ff);
                if (std::abs(c * ff) < std::abs(sum) * EPS) break;
            }
            QL_REQUIRE(i < MAXIT, "Modified Bessel function fails to converge");
            rkmu = sum;
            rk1 = sum1 * 2.0 / x;
        } else {
            b = 2.0 * (1.0 + x);
            d = 1.0 / b;
            h = delh = d;
            q1 = 0.0;
            q2 = 1.0;
            a1 = 0.25 - xmu2;
            q = c = a1;
            a = -a1;
            s = 1.0 + q * delh;
            for (i = 1; i < MAXIT; ++i) {
                a -= 2 * i;
                c = -a * c / (i + 1.0);
                qnew = (q1 - b * q2) / a;
                q1 = q2;
                q2 = qnew;
                q += c * qnew;
                b += 2.0;
                d = 1.0 / (b + a * d);
                delh = (b * d - 1.0) * delh;
                h += delh;
                dels = q * delh;
                s += dels;
                if (std::abs(dels / s) <= EPS) break;
            }
            QL_REQUIRE(i < MAXIT, "Modified Bessel function fails to converge in cf2");
            h = a1 * h;
            rkmu = std::sqrt(PI / (2.0 * x)) * std::exp(-x) / s;
            rk1 = rkmu * (xmu + x + 0.5 - h) / x;
        }
        rkmup = xmu * rkmu / x - rk1;
        rimu = 1.0 / ((f * rkmu - rkmup) * x);
        if (compute_i)
            return (rimu * ril1) / ril;
        else {
            for (i = 1; i <= nl; i++) {
                rktemp = (xmu + i) * 2.0 * rk1 / x + rkmu;
                rkmu = rk1;
                rk1 = rktemp;
            }
            return rkmu;
        }
    }

}

namespace QuantLib {

    namespace {

        template<class T>
        struct I {
        };

        template<>
        struct I<Real> {
            Real value() { return 0.0; }
        };

        template<>
        struct I<std::complex<Real> > {
            std::complex<Real> value() { return std::complex<Real>(0.0, 1.0); }
        };

        template<class T>
        struct Unweighted {
            T weightSmallX(const T &x) { return 1.0; }

            T weight1LargeX(const T &x) { return std::exp(x); }

            T weight2LargeX(const T &x) { return std::exp(-x); }
        };

        template<class T>
        struct ExponentiallyWeighted {
            T weightSmallX(const T &x) { return std::exp(-x); }

            T weight1LargeX(const T &x) { return 1.0; }

            T weight2LargeX(const T &x) { return std::exp(-2.0 * x); }
        };

        template<class T, template<class> class W>
        T modifiedBesselFunction_i_impl(Real nu, const T &x) {
            if (std::abs(x) < 13.0) {
                const T alpha = std::pow(0.5 * x, nu)
                                / std::tgamma(1.0 + nu);
                const T Y = 0.25 * x * x;
                Size k = 1;
                T sum = alpha, B_k = alpha;

                while (std::abs(B_k *= Y / (k * (k + nu))) > std::abs(sum) * QL_EPSILON) {
                    sum += B_k;
                    QL_REQUIRE(++k < 1000, "max iterations exceeded");
                }
                return sum * W<T>().weightSmallX(x);
            } else {
                Real na_k = 1.0, sign = 1.0;
                T da_k = T(1.0);

                T s1 = T(1.0), s2 = T(1.0);
                for (Size k = 1; k < 30; ++k) {
                    sign *= -1;
                    na_k *= (4.0 * nu * nu -
                             (2.0 * static_cast<Real>(k) - 1.0) *
                             (2.0 * static_cast<Real>(k) - 1.0));
                    da_k *= (8.0 * k) * x;
                    const T a_k = na_k / da_k;

                    s2 += a_k;
                    s1 += sign * a_k;
                }

                const T i = I<T>().value();
                return 1.0 / std::sqrt(2 * M_PI * x) *
                       (W<T>().weight1LargeX(x) * s1 +
                        i * std::exp(i * nu * M_PI) * W<T>().weight2LargeX(x) * s2);
            }
        }

        template<class T, template<class> class W>
        T modifiedBesselFunction_k_impl(Real nu, const T &x) {
            return M_PI_2 * (modifiedBesselFunction_i_impl<T, W>(-nu, x) -
                             modifiedBesselFunction_i_impl<T, W>(nu, x)) /
                   std::sin(M_PI * nu);
        }
    }

    Real modifiedBesselFunction_i(Real nu, Real x) {
        QL_REQUIRE(x >= 0.0, "negative argument requires complex version of "
                "modifiedBesselFunction");
        return modifiedBesselFunction_i_impl<Real, Unweighted>(nu, x);
    }

    std::complex<Real> modifiedBesselFunction_i(Real nu,
                                                const std::complex<Real> &z) {
        return modifiedBesselFunction_i_impl<
                std::complex<Real>, Unweighted>(nu, z);
    }

    Real modifiedBesselFunction_k(Real nu, Real x) {
        return modifiedBesselFunction_k_impl<Real, Unweighted>(nu, x);
    }

    std::complex<Real> modifiedBesselFunction_k(Real nu,
                                                const std::complex<Real> &z) {
        return modifiedBesselFunction_k_impl<
                std::complex<Real>, Unweighted>(nu, z);
    }

    Real modifiedBesselFunction_i_exponentiallyWeighted(Real nu, Real x) {
        QL_REQUIRE(x >= 0.0, "negative argument requires complex version of "
                "modifiedBesselFunction");
        return modifiedBesselFunction_i_impl<Real, ExponentiallyWeighted>(
                nu, x);
    }

    std::complex<Real> modifiedBesselFunction_i_exponentiallyWeighted(
            Real nu, const std::complex<Real> &z) {
        return modifiedBesselFunction_i_impl<
                std::complex<Real>, ExponentiallyWeighted>(nu, z);
    }

    Real modifiedBesselFunction_k_exponentiallyWeighted(Real nu, Real x) {
        return modifiedBesselFunction_k_impl<Real, ExponentiallyWeighted>(
                nu, x);
    }

    std::complex<Real> modifiedBesselFunction_k_exponentiallyWeighted(
            Real nu, const std::complex<Real> &z) {
        return modifiedBesselFunction_k_impl<
                std::complex<Real>, ExponentiallyWeighted>(nu, z);
    }

    long double modifiedBesselFunction_i_Press(long double nu, long double x) {
        return besselik_impl(nu, x, true);
    }

    long double modifiedBesselFunction_k_Press(long double nu, long double x) {
        return besselik_impl(nu, x, false);
    }

}
