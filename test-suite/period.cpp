/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007, 2014 Ferdinando Ametrano

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
#include "ql/time/period.hpp"

using namespace QuantLib;


TEST_CASE("Period_YearsMonthsAlgebra", "[Period]") {

    INFO("Testing period algebra on years/months...");

    Period OneYear(1, Years);
    Period SixMonths(6, Months);
    Period ThreeMonths(3, Months);

    Integer n = 4;
    if (OneYear/n!=ThreeMonths)
        FAIL_CHECK("division error: " << OneYear << "/" << n <<
                    " not equal to " << ThreeMonths);
    n = 2;
    if (OneYear/n!=SixMonths)
        FAIL_CHECK("division error: " << OneYear << "/" << n <<
                    " not equal to " << SixMonths);

    Period sum=ThreeMonths;
    sum+=SixMonths;
    if (sum!=Period(9, Months))
        FAIL_CHECK("sum error: " << ThreeMonths <<
                    " + " << SixMonths <<
                    " != " << Period(9, Months));

    sum+=OneYear;
    if (sum!=Period(21, Months))
        FAIL_CHECK("sum error: " << ThreeMonths <<
                    " + " << SixMonths <<
                    " + " << OneYear <<
                    " != " << Period(21, Months));

    Period TwelveMonths(12, Months);
    if (TwelveMonths.length()!=12)
        FAIL_CHECK("normalization error: TwelveMonths.length()" <<
                    " is " << TwelveMonths.length() <<
                    " instead of 12");
    if (TwelveMonths.units()!=Months)
        FAIL_CHECK("normalization error: TwelveMonths.units()" <<
                    " is " << TwelveMonths.units() <<
                    " instead of " << Months);

    Period NormalizedTwelveMonths(12, Months);
    NormalizedTwelveMonths.normalize();
    if (NormalizedTwelveMonths.length()!=1)
        FAIL_CHECK("normalization error: NormalizedTwelveMonths.length()" <<
                    " is " << NormalizedTwelveMonths.length() <<
                    " instead of 1");
    if (NormalizedTwelveMonths.units()!=Years)
        FAIL_CHECK("normalization error: NormalizedTwelveMonths.units()" <<
                    " is " << NormalizedTwelveMonths.units() <<
                    " instead of " << Years);
}

TEST_CASE("Period_WeeksDaysAlgebra", "[Period]") {

    INFO("Testing period algebra on weeks/days...");

    Period TwoWeeks(2, Weeks);
    Period OneWeek(1, Weeks);
    Period ThreeDays(3, Days);
    Period OneDay(1, Days);

    Integer n = 2;
    if (TwoWeeks/n!=OneWeek)
        FAIL_CHECK("division error: " << TwoWeeks << "/" << n <<
                    " not equal to " << OneWeek);
    n = 7;
    if (OneWeek/n!=OneDay)
        FAIL_CHECK("division error: " << OneWeek << "/" << n <<
                    " not equal to " << OneDay);

    Period sum=ThreeDays;
    sum+=OneDay;
    if (sum!=Period(4, Days))
        FAIL_CHECK("sum error: " << ThreeDays <<
                    " + " << OneDay <<
                    " != " << Period(4, Days));

    sum+=OneWeek;
    if (sum!=Period(11, Days))
        FAIL_CHECK("sum error: " << ThreeDays <<
                    " + " << OneDay <<
                    " + " << OneWeek <<
                    " != " << Period(11, Days));

    Period SevenDays(7, Days);
    if (SevenDays.length()!=7)
        FAIL_CHECK("normalization error: SevenDays.length()" <<
                    " is " << SevenDays.length() <<
                    " instead of 7");
    if (SevenDays.units()!=Days)
        FAIL_CHECK("normalization error: SevenDays.units()" <<
                    " is " << SevenDays.units() <<
                    " instead of " << Days);
}
