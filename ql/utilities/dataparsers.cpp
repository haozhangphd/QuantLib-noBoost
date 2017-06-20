/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 Decillion Pty(Ltd)
 Copyright (C) 2006 Joseph Wang
 Copyright (2) 2009 Mark Joshi
 Copyright (2) 2009 StatPro Italia srl

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

#include <ql/utilities/dataparsers.hpp>
#include <ql/utilities/null.hpp>
#include <ql/time/period.hpp>
#include <ql/errors.hpp>

#if defined(__GNUC__) && (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif
#ifndef QL_PATCH_SOLARIS
#endif
#if defined(__GNUC__) && (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ > 4))
#pragma GCC diagnostic pop
#endif

#include <locale>
#include <cctype>
#include <iomanip>

namespace QuantLib {

    namespace io {

        Integer to_integer(const std::string &str) {
            return std::atoi(str.c_str());
        }

    }

    Period PeriodParser::parse(const std::string &str) {
        QL_REQUIRE(str.length() > 1, "period string length must be at least 2");

        std::vector<std::string> subStrings;
        std::string reducedString = str;

        Size iPos, reducedStringDim = 100000, max_iter = 0;
        while (reducedStringDim > 0) {
            iPos = reducedString.find_first_of("DdWwMmYy");
            Size subStringDim = iPos + 1;
            reducedStringDim = reducedString.length() - subStringDim;
            subStrings.emplace_back(reducedString.substr(0, subStringDim));
            reducedString = reducedString.substr(iPos + 1, reducedStringDim);
            ++max_iter;
            QL_REQUIRE(max_iter < str.length(), "unknown '" << str << "' unit");
        }

        Period result = parseOnePeriod(subStrings[0]);
        for (Size i = 1; i < subStrings.size(); ++i)
            result += parseOnePeriod(subStrings[i]);
        return result;
    }

    Period PeriodParser::parseOnePeriod(const std::string &str) {
        QL_REQUIRE(str.length() > 1, "single period require a string of at "
                "least 2 characters");

        Size iPos = str.find_first_of("DdWwMmYy");
        QL_REQUIRE(iPos == str.length() - 1, "unknown '" <<
                                                         str.substr(str.length() - 1, str.length()) << "' unit");
        TimeUnit units = Days;
        char abbr = static_cast<char>(std::toupper(str[iPos]));
        if (abbr == 'D') units = Days;
        else if (abbr == 'W') units = Weeks;
        else if (abbr == 'M') units = Months;
        else if (abbr == 'Y') units = Years;

        Size nPos = str.find_first_of("-+0123456789");
        QL_REQUIRE(nPos < iPos, "no numbers of " << units << " provided");
        Integer n;
        try {
            n = io::to_integer(str.substr(nPos, iPos));
        } catch (std::exception &e) {
            QL_FAIL("unable to parse the number of units of " << units <<
                                                              " in '" << str << "'. Error:" << e.what());
        }

        return Period(n, units);
    }


#ifdef QL_HIGH_RESOLUTION_DATE
    Date DateParser::parseFormatted(const std::string &str,
                                    const std::string &fmt) {
        std::istringstream is(str);
        is.imbue(std::locale());
        std::tm t = {};
        is >> std::get_time(&t, fmt.c_str());
        return Date(std::chrono::time_point_cast<std::chrono::microseconds>(
                std::chrono::system_clock::from_time_t(timegm(&t))));
    }
#else
    Date DateParser::parseFormatted(const std::string &str,
                                    const std::string &fmt) {
        std::istringstream is(str);
        is.imbue(std::locale());
        std::tm t = {};
        is >> std::get_time(&t, fmt.c_str());
        return Date(t.tm_mday, static_cast<Month>(t.tm_mon+1), t.tm_year+1900);
    }
#endif

    Date DateParser::parseISO(const std::string &str) {
        QL_REQUIRE(str.size() == 10 && str[4] == '-' && str[7] == '-',
                   "invalid format");
        Integer year = io::to_integer(str.substr(0, 4));
        Month month = static_cast<Month>(io::to_integer(str.substr(5, 2)));
        Integer day = io::to_integer(str.substr(8, 2));

        return Date(day, month, year);
    }

}
