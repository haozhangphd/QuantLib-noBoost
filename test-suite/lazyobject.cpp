/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016 StatPro Italia srl

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
#include <ql/instruments/stock.hpp>
#include <ql/quotes/simplequote.hpp>

using namespace QuantLib;

using std::shared_ptr;

TEST_CASE("LazyObject_DiscardingNotifications", "[LazyObject]") {

    INFO(
        "Testing that lazy objects discard notifications after the first...");

    std::shared_ptr<SimpleQuote> q = std::make_shared<SimpleQuote>(0.0);
    std::shared_ptr<Instrument> s = std::make_shared<Stock>(Handle<Quote>(q));

    Flag f;
    f.registerWith(s);
    
    s->NPV();
    q->setValue(1.0);
    if (!f.isUp())
        FAIL("Observer was not notified of change");
    
    f.lower();
    q->setValue(2.0);
    if (f.isUp())
        FAIL("Observer was notified of second change");

    f.lower();
    s->NPV();
    q->setValue(3.0);
    if (!f.isUp())
        FAIL("Observer was not notified of change after recalculation");
}


TEST_CASE("LazyObject_ForwardingNotifications", "[LazyObject]") {

    INFO(
        "Testing that lazy objects forward all notifications when told...");

    std::shared_ptr<SimpleQuote> q = std::make_shared<SimpleQuote>(0.0);
    std::shared_ptr<Instrument> s = std::make_shared<Stock>(Handle<Quote>(q));

    s->alwaysForwardNotifications();

    Flag f;
    f.registerWith(s);
    
    s->NPV();
    q->setValue(1.0);
    if (!f.isUp())
        FAIL("Observer was not notified of change");
    
    f.lower();
    q->setValue(2.0);
    if (!f.isUp())
        FAIL("Observer was not notified of second change");
}
