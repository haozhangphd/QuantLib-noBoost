/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2004, 2007 StatPro Italia srl

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

#include <ql/currency.hpp>
#include <iomanip>

namespace QuantLib {

    std::string currency_format(const Currency& c, Decimal value) {

        std::ostringstream ret;
        ret.imbue(std::locale(""));

        std::string fmt_str{c.format()};

        switch (fmt_str[1]) {
            case '1':
                switch (fmt_str[8]) {
                    case '2':
                        ret << std::fixed << std::setprecision(int{fmt_str[4]} - '0') << value << " " << c.code();
                        break;
                    case '3':
                        ret << std::fixed << std::setprecision(int{fmt_str[4]} - '0') << value << " " << c.symbol();
                        break;
                    default:
                        QL_FAIL("\"Could not parse format string: \" << fmt_str");
                }
                break;
            case '2':
                ret << c.code() << " " << std::fixed << std::setprecision(int{fmt_str[8]} - '0') << value;
                break;
            case '3':
                ret << c.symbol() << " " << std::fixed << std::setprecision(int{fmt_str[8]} - '0') << value;
                break;
            default:
                QL_FAIL("\"Could not parse format string: \" << fmt_str");
        }
        return ret.str();
    }

    std::ostream &operator<<(std::ostream &out, const Currency &c) {
        if (!c.empty())
            return out << c.code();
        else
            return out << "null currency";
    }

    Currency::Data::Data(const std::string &name,
                         const std::string &code,
                         Integer numericCode,
                         const std::string &symbol,
                         const std::string &fractionSymbol,
                         Integer fractionsPerUnit,
                         const Rounding &rounding,
                         const std::string &formatString,
                         const Currency &triangulationCurrency)
            : name(name), code(code), numeric(numericCode),
              symbol(symbol), fractionSymbol(fractionSymbol),
              fractionsPerUnit(fractionsPerUnit), rounding(rounding),
              triangulated(triangulationCurrency),
              formatString(formatString) {}

}

