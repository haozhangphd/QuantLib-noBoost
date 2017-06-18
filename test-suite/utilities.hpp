/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003, 2004, 2008 StatPro Italia srl

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

#ifndef quantlib_test_utilities_hpp
#define quantlib_test_utilities_hpp

#include "catch.hpp"
#include <ql/instruments/payoffs.hpp>
#include <ql/exercise.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/quote.hpp>
#include <ql/patterns/observable.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <functional>
#include <vector>
#include <string>
#include <numeric>
#include <iomanip>
#include <chrono>
#include <iostream>

// This makes it easier to use array literals (alas, no std::vector literals)
#define LENGTH(a) (sizeof(a)/sizeof(a[0]))

#define QL_FIXED std::fixed
#define QL_SCIENTIFIC std::scientific

/* the following displays the elapsed time for the test if
   QL_DISPLAY_TEST_TIME is defined. */
#if defined(QL_DISPLAY_TEST_TIME)
class progress_timer {
public:
    progress_timer( std::ostream & os = std::cout ) : startT{std::chrono::steady_clock::now()}, os_{os} {}

    progress_timer(const progress_timer&) = delete;

    ~progress_timer() {
        try {
            std::chrono::time_point<std::chrono::steady_clock> endT = std::chrono::steady_clock::now();
            double elapsed = static_cast<double>((endT - startT).count()) / 1.0e9;

            std::istream::fmtflags old_flags = os_.setf(std::istream::fixed,
                                                         std::istream::floatfield);
            std::streamsize old_prec = os_.precision(2);
            os_ << elapsed << " s\n" << std::endl;
            os_.flags(old_flags);
            os_.precision(old_prec);
        }
        catch (...) {}
    }

private:
    std::chrono::time_point<std::chrono::steady_clock> startT;
    std::ostream &os_;
};

#define QL_TEST_START_TIMING progress_timer t;
#else
#define QL_TEST_START_TIMING
#endif

namespace QuantLib {

    std::string payoffTypeToString(const std::shared_ptr<Payoff> &);

    std::string exerciseTypeToString(const std::shared_ptr<Exercise> &);


    std::shared_ptr<YieldTermStructure>
    flatRate(const Date &today,
             const std::shared_ptr<Quote> &forward,
             const DayCounter &dc);

    std::shared_ptr<YieldTermStructure>
    flatRate(const Date &today,
             Rate forward,
             const DayCounter &dc);

    std::shared_ptr<YieldTermStructure>
    flatRate(const std::shared_ptr<Quote> &forward,
             const DayCounter &dc);

    std::shared_ptr<YieldTermStructure>
    flatRate(Rate forward,
             const DayCounter &dc);


    std::shared_ptr<BlackVolTermStructure>
    flatVol(const Date &today,
            const std::shared_ptr<Quote> &volatility,
            const DayCounter &dc);

    std::shared_ptr<BlackVolTermStructure>
    flatVol(const Date &today,
            Volatility volatility,
            const DayCounter &dc);

    std::shared_ptr<BlackVolTermStructure>
    flatVol(const std::shared_ptr<Quote> &volatility,
            const DayCounter &dc);

    std::shared_ptr<BlackVolTermStructure>
    flatVol(Volatility volatility,
            const DayCounter &dc);


    Real relativeError(Real x1, Real x2, Real reference);

    //bool checkAbsError(Real x1, Real x2, Real tolerance){
    //    return std::fabs(x1 - x2) < tolerance;
    //};

    class Flag : public QuantLib::Observer {
    private:
        bool up_;
    public:
        Flag() : up_(false) {}

        void raise() { up_ = true; }

        void lower() { up_ = false; }

        bool isUp() const { return up_; }

        void update() { raise(); }
    };

    template<class Iterator>
    Real norm(const Iterator &begin, const Iterator &end, Real h) {
        // squared values
        std::vector<Real> f2(end - begin);
        std::transform(begin, end, begin, f2.begin(),
                       std::multiplies<Real>());
        // numeric integral of f^2
        Real I = h * (std::accumulate(f2.begin(), f2.end(), Real(0.0))
                      - 0.5 * f2.front() - 0.5 * f2.back());
        return std::sqrt(I);
    }


    // this cleans up index-fixing histories when destroyed
    class IndexHistoryCleaner {
    public:
        IndexHistoryCleaner();

        ~IndexHistoryCleaner();
    };


    // Allow streaming vectors to error messages.

    // The standard forbids defining new overloads in the std
    // namespace, so we have to use a wrapper instead of overloading
    // operator<< to send a vector to the stream directly.

    template<class T>
    struct vector_streamer {
        vector_streamer(const std::vector<T> &v) : v(v) {}

        std::vector<T> v;
    };

    template<class T>
    vector_streamer<T> to_stream(const std::vector<T> &v) {
        return vector_streamer<T>(v);
    }

    template<class T>
    std::ostream &operator<<(std::ostream &out, const vector_streamer<T> &s) {
        out << "{ ";
        if (!s.v.empty()) {
            for (size_t n = 0; n < s.v.size() - 1; ++n)
                out << s.v[n] << ", ";
            out << s.v.back();
        }
        out << " }";
        return out;
    }


}

#endif
