/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2004, 2005 StatPro Italia srl

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
    Data from http://fx.sauder.ubc.ca/currency_table.html
    and http://www.thefinancials.com/vortex/CurrencyFormats.html
*/

#include <ql/currencies/asia.hpp>

namespace QuantLib {

    // Bangladesh taka
    /* The ISO three-letter code is BDT; the numeric code is 50.
       It is divided in 100 paisa.
    */
    BDTCurrency::BDTCurrency() {
        static std::shared_ptr<Data> bdtData =
                std::make_shared<Data>("Bangladesh taka", "BDT", 50,
                                       "Bt", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = bdtData;
    }

    // Chinese yuan
    /* The ISO three-letter code is CNY; the numeric code is 156.
       It is divided in 100 fen.
    */
    CNYCurrency::CNYCurrency() {
        static std::shared_ptr<Data> cnyData =
                std::make_shared<Data>("Chinese yuan", "CNY", 156,
                                       u8"¥", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = cnyData;
    }

    // Hong Kong dollar
    /* The ISO three-letter code is HKD; the numeric code is 344.
      It is divided in 100 cents.
    */
    HKDCurrency::HKDCurrency() {
        static std::shared_ptr<Data> hkdData =
                std::make_shared<Data>("Hong Kong dollar", "HKD", 344,
                                       "HK$", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = hkdData;
    }

    // Indonesian Rupiah
    /* The ISO three-letter code is IDR; the numeric code is 360.
       It is divided in 100 sen.
    */
    IDRCurrency::IDRCurrency() {
        static std::shared_ptr<Data> idrData =
                std::make_shared<Data>("Indonesian Rupiah", "IDR", 360,
                                       "Rp", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = idrData;
    }

    // Israeli shekel
    /* The ISO three-letter code is ILS; the numeric code is 376.
      It is divided in 100 agorot.
    */
    ILSCurrency::ILSCurrency() {
        static std::shared_ptr<Data> ilsData =
                std::make_shared<Data>("Israeli shekel", "ILS", 376,
                                       "NIS", "", 100,
                                       Rounding(),
                                       "%1$.2f %3%");
        data_ = ilsData;
    }

    // Indian rupee
    /* The ISO three-letter code is INR; the numeric code is 356.
       It is divided in 100 paise.
    */
    INRCurrency::INRCurrency() {
        static std::shared_ptr<Data> inrData =
                std::make_shared<Data>("Indian rupee", "INR", 356,
                                       "Rs", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = inrData;
    }

    // Iraqi dinar
    /* The ISO three-letter code is IQD; the numeric code is 368.
       It is divided in 1000 fils.
    */
    IQDCurrency::IQDCurrency() {
        static std::shared_ptr<Data> iqdData =
                std::make_shared<Data>("Iraqi dinar", "IQD", 368,
                                       "ID", "", 1000,
                                       Rounding(),
                                       "%2% %1$.3f");
        data_ = iqdData;
    }

    // Iranian rial
    /* The ISO three-letter code is IRR; the numeric code is 364.
       It has no subdivisions.
    */
    IRRCurrency::IRRCurrency() {
        static std::shared_ptr<Data> irrData =
                std::make_shared<Data>("Iranian rial", "IRR", 364,
                                       "Rls", "", 1,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = irrData;
    }

    // Japanese yen
    /* The ISO three-letter code is JPY; the numeric code is 392.
       It is divided into 100 sen.
    */
    JPYCurrency::JPYCurrency() {
        static std::shared_ptr<Data> jpyData =
                std::make_shared<Data>("Japanese yen", "JPY", 392,
                                       u8"¥", "", 100,
                                       Rounding(),
                                       "%3% %1$.0f");
        data_ = jpyData;
    }

    // South-Korean won
    /* The ISO three-letter code is KRW; the numeric code is 410.
       It is divided in 100 chon.
    */
    KRWCurrency::KRWCurrency() {
        static std::shared_ptr<Data> krwData =
                std::make_shared<Data>("South-Korean won", "KRW", 410,
                                       u8"₩", "", 100,
                                       Rounding(),
                                       "%3% %1$.0f");
        data_ = krwData;
    }

    // Kuwaiti dinar
    /* The ISO three-letter code is KWD; the numeric code is 414.
       It is divided in 1000 fils.
    */
    KWDCurrency::KWDCurrency() {
        static std::shared_ptr<Data> kwdData =
                std::make_shared<Data>("Kuwaiti dinar", "KWD", 414,
                                       "KD", "", 1000,
                                       Rounding(),
                                       "%3% %1$.3f");
        data_ = kwdData;
    }

    // Malaysian Ringgit
    /* The ISO three-letter code is MYR; the numeric code is 458.
       It is divided in 100 sen.
    */
    MYRCurrency::MYRCurrency() {
        static std::shared_ptr<Data> myrData =
                std::make_shared<Data>("Malaysian Ringgit",
                                       "MYR", 458,
                                       "RM", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = myrData;
    }

    // Nepal rupee
    /* The ISO three-letter code is NPR; the numeric code is 524.
       It is divided in 100 paise.
    */
    NPRCurrency::NPRCurrency() {
        static std::shared_ptr<Data> nprData =
                std::make_shared<Data>("Nepal rupee", "NPR", 524,
                                       "NRs", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = nprData;
    }

    // Pakistani rupee
    /* The ISO three-letter code is PKR; the numeric code is 586.
       It is divided in 100 paisa.
    */
    PKRCurrency::PKRCurrency() {
        static std::shared_ptr<Data> pkrData =
                std::make_shared<Data>("Pakistani rupee", "PKR", 586,
                                       "Rs", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = pkrData;
    }

    // Saudi riyal
    /* The ISO three-letter code is SAR; the numeric code is 682.
       It is divided in 100 halalat.
    */
    SARCurrency::SARCurrency() {
        static std::shared_ptr<Data> sarData =
                std::make_shared<Data>("Saudi riyal", "SAR", 682,
                                       "SRls", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = sarData;
    }

    // %Singapore dollar
    /* The ISO three-letter code is SGD; the numeric code is 702.
       It is divided in 100 cents.
    */
    SGDCurrency::SGDCurrency() {
        static std::shared_ptr<Data> sgdData =
                std::make_shared<Data>("Singapore dollar", "SGD", 702,
                                       "S$", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = sgdData;
    }

    // Thai baht
    /* The ISO three-letter code is THB; the numeric code is 764.
       It is divided in 100 stang.
    */
    THBCurrency::THBCurrency() {
        static std::shared_ptr<Data> thbData =
                std::make_shared<Data>("Thai baht", "THB", 764,
                                       "Bht", "", 100,
                                       Rounding(),
                                       "%1$.2f %3%");
        data_ = thbData;
    }

    // %Taiwan dollar
    /* The ISO three-letter code is TWD; the numeric code is 901.
       It is divided in 100 cents.
    */
    TWDCurrency::TWDCurrency() {
        static std::shared_ptr<Data> twdData =
                std::make_shared<Data>("Taiwan dollar", "TWD", 901,
                                       "NT$", "", 100,
                                       Rounding(),
                                       "%3% %1$.2f");
        data_ = twdData;
    }

    // Vietnamese Dong
    /* The ISO three-letter code is VND; the numeric code is 704.
       It was divided in 100 xu.
    */
    VNDCurrency::VNDCurrency() {
        static std::shared_ptr<Data> twdData =
                std::make_shared<Data>("Vietnamese Dong", "VND", 704,
                                       u8"₫", "", 100,
                                       Rounding(),
                                       "%1$.0f %3%");
        data_ = twdData;
    }

}

