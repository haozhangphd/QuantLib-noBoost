/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2004 StatPro Italia srl

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
#include <ql/money.hpp>
#include <ql/currencies/asia.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/currencies/america.hpp>
#include <ql/currencies/exchangeratemanager.hpp>

using namespace QuantLib;


TEST_CASE("Money_None", "[Money]") {

    INFO("Testing money arithmetic without conversions...");

    Currency EUR = EURCurrency();

    Money m1 = 50000.0 * EUR;
    Money m2 = 100000.0 * EUR;
    Money m3 = 500000.0 * EUR;

    Money::conversionType = Money::NoConversion;

    Money calculated = m1 * 3.0 + 2.5 * m2 - m3 / 5.0;
    Decimal x = m1.value() * 3.0 + 2.5 * m2.value() - m3.value() / 5.0;
    Money expected(x, EUR);

    if (calculated != expected) {
        FAIL("Wrong result: \n"
                     << "    expected:   " << expected << "\n"
                     << "    calculated: " << calculated);
    }
}


TEST_CASE("Money_BaseCurrency", "[Money]") {

    INFO("Testing money arithmetic with conversion "
                 "to base currency...");

    Currency EUR = EURCurrency(), GBP = GBPCurrency(), USD = USDCurrency();

    Money m1 = 50000.0 * GBP;
    Money m2 = 100000.0 * EUR;
    Money m3 = 500000.0 * USD;

    ExchangeRateManager::instance().clear();
    ExchangeRate eur_usd = ExchangeRate(EUR, USD, 1.2042);
    ExchangeRate eur_gbp = ExchangeRate(EUR, GBP, 0.6612);
    ExchangeRateManager::instance().add(eur_usd);
    ExchangeRateManager::instance().add(eur_gbp);

    Money::conversionType = Money::BaseCurrencyConversion;
    Money::baseCurrency = EUR;

    Money calculated = m1 * 3.0 + 2.5 * m2 - m3 / 5.0;

    Rounding round = Money::baseCurrency.rounding();
    Decimal x = round(m1.value() * 3.0 / eur_gbp.rate()) + 2.5 * m2.value()
                - round(m3.value() / (5.0 * eur_usd.rate()));
    Money expected(x, EUR);

    Money::conversionType = Money::NoConversion;

    if (calculated != expected) {
        FAIL("Wrong result: \n"
                     << "    expected:   " << expected << "\n"
                     << "    calculated: " << calculated);
    }
}


TEST_CASE("Money_Automated", "[Money]") {

    INFO("Testing money arithmetic with automated conversion...");

    Currency EUR = EURCurrency(), GBP = GBPCurrency(), USD = USDCurrency();

    Money m1 = 50000.0 * GBP;
    Money m2 = 100000.0 * EUR;
    Money m3 = 500000.0 * USD;

    ExchangeRateManager::instance().clear();
    ExchangeRate eur_usd = ExchangeRate(EUR, USD, 1.2042);
    ExchangeRate eur_gbp = ExchangeRate(EUR, GBP, 0.6612);
    ExchangeRateManager::instance().add(eur_usd);
    ExchangeRateManager::instance().add(eur_gbp);

    Money::conversionType = Money::AutomatedConversion;

    Money calculated = (m1 * 3.0 + 2.5 * m2) - m3 / 5.0;

    Rounding round = m1.currency().rounding();
    Decimal x = m1.value() * 3.0 + round(2.5 * m2.value() * eur_gbp.rate())
                - round((m3.value() / 5.0) * eur_gbp.rate() / eur_usd.rate());
    Money expected(x, GBP);

    Money::conversionType = Money::NoConversion;

    if (calculated != expected) {
        FAIL("Wrong result: \n"
                     << "    expected:   " << expected << "\n"
                     << "    calculated: " << calculated);
    }
}

TEST_CASE("Money_Formatting", "[.]") {

    INFO("Testing formatting money with different currency...");

    Currency EUR = EURCurrency(), GBP = GBPCurrency(), USD = USDCurrency();
    Currency CNY = CNYCurrency(), VND = VNDCurrency(), KWD = KWDCurrency();

    Money gbp = 50000.0 * GBP;
    Money eur = 100000.0 * EUR;
    Money usd = 500000.0 * USD;
    Money cny = 1234.5678 * CNY;
    Money vnd = 9876.54321 * VND;
    Money kwd = 3.1415926 * KWD;

    std::ostringstream gbp_str, eur_str, usd_str, cny_str, vnd_str, kwd_str;

    gbp_str << gbp;
    CHECK(gbp_str.str() == "£ 50,000.00");

    eur_str << eur;
    CHECK(eur_str.str() == "EUR 100,000.00");

    usd_str << usd;
    CHECK(usd_str.str() == "$ 500,000.00");

    cny_str << cny;
    CHECK(cny_str.str() == "¥ 1,234.57");

    vnd_str << vnd;
    CHECK(vnd_str.str() == "9,877 ₫");

    kwd_str << kwd;
    CHECK(kwd_str.str() == "KD 3.142");


}
