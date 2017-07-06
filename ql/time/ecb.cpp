/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2009, 2011 Ferdinando Ametrano
 Copyright (C) 2015 Paolo Mazzocchi

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

#include <ql/time/ecb.hpp>
#include <ql/settings.hpp>
#include <ql/utilities/dataparsers.hpp>
#include <ql/utilities/stringutils.hpp>
#include <algorithm>

using std::string;

namespace QuantLib {

    static std::set<Date> knownDateSet;

    const std::set<Date>& ECB::knownDates() {

        // one-off inizialization
        static const Date::serial_type knownDatesArray[] = {
              38370, 38390, 38419, 38454, 38482, 38510, 38545, 38573, 38601, 38636, 38664, 38691 // 2005
            , 38734, 38755, 38783, 38818, 38846, 38882, 38909, 38937, 38965, 39000, 39028, 39063 // 2006
            , 39098, 39126, 39154, 39189, 39216, 39245, 39273, 39301, 39336, 39364, 39399, 39427 // 2007
            , 39462, 39490, 39518, 39553, 39581, 39609, 39637, 39672, 39700, 39728, 39763, 39791 // 2008
            , 39833, 39854, 39882, 39910, 39945, 39973, 40001, 40036, 40064, 40099, 40127, 40154 // 2009
            , 40197, 40218, 40246, 40281, 40309, 40344, 40372, 40400, 40428, 40463, 40491, 40519 // 2010
            , 40561, 40582, 40610, 40645, 40673, 40708, 40736, 40764, 40799, 40827, 40855, 40890 // 2011
            // http://www.ecb.europa.eu/press/pr/date/2011/html/pr110520.en.html
            , 40925, 40953, 40981, 41009, 41037, 41072, 41100, 41128, 41163, 41191, 41226, 41254 // 2012
            , 41289, 41317, 41345, 41373, 41401, 41436, 41464, 41492, 41527, 41555, 41590, 41618 // 2013
            // http://www.ecb.europa.eu/press/pr/date/2013/html/pr130610.en.html
            , 41653, 41681, 41709, 41737, 41772, 41800, 41828, 41863, 41891, 41919, 41954, 41982 // 2014
            // http://www.ecb.europa.eu/press/pr/date/2014/html/pr140717_1.en.html
            , 42031, 42073, 42115, 42164, 42206, 42255, 42304, 42346// 2015
            // https://www.ecb.europa.eu/press/pr/date/2015/html/pr150622.en.html
            , 42395, 42444, 42486, 42528, 42577, 42626, 42668, 42717 // 2016
            // https://www.ecb.europa.eu/press/calendars/reserve/html/index.en.html
            , 42759, 42808, 42857, 42899, 42941, 42990, 43039, 43088 //2017
        };
        if (knownDateSet.empty()) {
            Size n = sizeof(knownDatesArray)/sizeof(Date::serial_type);
            for (Size i=0; i<n; ++i)
                knownDateSet.insert(Date(knownDatesArray[i]));
        }

        return knownDateSet;
    }

    void ECB::addDate(const Date& d) {
        knownDates(); // just to ensure inizialization
        knownDateSet.insert(d);
    }

    void ECB::removeDate(const Date& d) {
        knownDates(); // just to ensure inizialization
        knownDateSet.erase(d);
    }

    Date ECB::date(const string& ecbCode,
                   const Date& refDate) {

        QL_REQUIRE(isECBcode(ecbCode),
                   ecbCode << " is not a valid ECB code");

        string code = to_upper_copy(ecbCode);
        string monthString = code.substr(0, 3);
        Month m;
        if (monthString=="JAN")      m = January;
        else if (monthString=="FEB") m = February;
        else if (monthString=="MAR") m = March;
        else if (monthString=="APR") m = April;
        else if (monthString=="MAY") m = May;
        else if (monthString=="JUN") m = June;
        else if (monthString=="JUL") m = July;
        else if (monthString=="AUG") m = August;
        else if (monthString=="SEP") m = September;
        else if (monthString=="OCT") m = October;
        else if (monthString=="NOV") m = November;
        else if (monthString=="DEC") m = December;
        else QL_FAIL("not an ECB month (and it should have been)");

        Year y = io::to_integer(code.substr(3, 2));
        Date referenceDate = (refDate != Date() ?
                              refDate :
                              Date(Settings::instance().evaluationDate()));
        Year referenceYear = (referenceDate.year() % 100);
        y += referenceDate.year() - referenceYear;
        if (y<Date::minDate().year())
            return ECB::nextDate(Date::minDate());

        return ECB::nextDate(Date(1, m, y) - 1);
    }

    string ECB::code(const Date& ecbDate) {

        QL_REQUIRE(isECBdate(ecbDate),
                   ecbDate << " is not a valid ECB date");

        std::ostringstream ECBcode;
        unsigned int y = ecbDate.year() % 100;
        string padding;
        if (y < 10)
            padding = "0";
        switch(ecbDate.month()) {
          case January:
            ECBcode << "JAN" << padding << y;
            break;
          case February:
            ECBcode << "FEB" << padding << y;
            break;
          case March:
            ECBcode << "MAR" << padding << y;
            break;
          case April:
            ECBcode << "APR" << padding << y;
            break;
          case May:
            ECBcode << "MAY" << padding << y;
            break;
          case June:
            ECBcode << "JUN" << padding << y;
            break;
          case July:
            ECBcode << "JUL" << padding << y;
            break;
          case August:
            ECBcode << "AUG" << padding << y;
            break;
          case September:
            ECBcode << "SEP" << padding << y;
            break;
          case October:
            ECBcode << "OCT" << padding << y;
            break;
          case November:
            ECBcode << "NOV" << padding << y;
            break;
          case December:
            ECBcode << "DEC" << padding << y;
            break;
          default:
            QL_FAIL("not an ECB month (and it should have been)");
        }

        #if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_ENSURE(isECBcode(ECBcode.str()),
                  "the result " << ECBcode.str() <<
                  " is an invalid ECB code");
        #endif
        return ECBcode.str();
    }



    Date ECB::nextDate(const Date& date) {
        Date d = (date == Date() ?
                  Settings::instance().evaluationDate() :
                  date);

        std::set<Date>::const_iterator i =
            std::upper_bound(knownDates().begin(), knownDates().end(), d);

        QL_REQUIRE(i!=knownDates().end(),
                   "ECB dates after " << *(--knownDates().end()) << " are unknown");
        return Date(*i);
    }

    std::vector<Date> ECB::nextDates(const Date& date) {
        Date d = (date == Date() ?
                  Settings::instance().evaluationDate() :
                  date);

        std::set<Date>::const_iterator i =
            std::upper_bound(knownDates().begin(), knownDates().end(), d);

        QL_REQUIRE(i!=knownDates().end(),
                   "ECB dates after " << *knownDates().end() << " are unknown");
        return std::vector<Date>(i, knownDates().end());
    }


    bool ECB::isECBcode(const std::string& ecbCode) {

        if (ecbCode.length() != 5)
            return false;

        string code = to_upper_copy(ecbCode);

        string str1("0123456789");
        string::size_type loc = str1.find(code.substr(3, 1), 0);
        if (loc == string::npos)
            return false;
        loc = str1.find(code.substr(4, 1), 0);
        if (loc == string::npos)
            return false;

        string monthString = code.substr(0, 3);
        if (monthString=="JAN")      return true;
        else if (monthString=="FEB") return true;
        else if (monthString=="MAR") return true;
        else if (monthString=="APR") return true;
        else if (monthString=="MAY") return true;
        else if (monthString=="JUN") return true;
        else if (monthString=="JUL") return true;
        else if (monthString=="AUG") return true;
        else if (monthString=="SEP") return true;
        else if (monthString=="OCT") return true;
        else if (monthString=="NOV") return true;
        else if (monthString=="DEC") return true;
        else return false;
    }

    string ECB::nextCode(const std::string& ecbCode) {
        QL_REQUIRE(isECBcode(ecbCode),
                   ecbCode << " is not a valid ECB code");

        string code = to_upper_copy(ecbCode);
        std::ostringstream result;

        string monthString = code.substr(0, 3);
        if (monthString=="JAN")      result << "FEB" << code.substr(3, 2);
        else if (monthString=="FEB") result << "MAR" << code.substr(3, 2);
        else if (monthString=="MAR") result << "APR" << code.substr(3, 2);
        else if (monthString=="APR") result << "MAY" << code.substr(3, 2);
        else if (monthString=="MAY") result << "JUN" << code.substr(3, 2);
        else if (monthString=="JUN") result << "JUL" << code.substr(3, 2);
        else if (monthString=="JUL") result << "AUG" << code.substr(3, 2);
        else if (monthString=="AUG") result << "SEP" << code.substr(3, 2);
        else if (monthString=="SEP") result << "OCT" << code.substr(3, 2);
        else if (monthString=="OCT") result << "NOV" << code.substr(3, 2);
        else if (monthString=="NOV") result << "DEC" << code.substr(3, 2);
        else if (monthString=="DEC") {
            unsigned int y = (io::to_integer(code.substr(3, 2)) + 1) % 100;
            string padding;
            if (y < 10)
                padding = "0";

            result << "JAN" << padding << y;
        } else QL_FAIL("not an ECB month (and it should have been)");


        #if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_ENSURE(isECBcode(result.str()),
                  "the result " << result.str() <<
                  " is an invalid ECB code");
        #endif
        return result.str();
    }

}
