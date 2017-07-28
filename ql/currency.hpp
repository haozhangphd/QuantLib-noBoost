/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2004 Decillion Pty(Ltd)
 Copyright (C) 2004, 2005, 2006, 2007 StatPro Italia srl

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

/*! \file currency.hpp
    \brief Currency specification
*/

#ifndef quantlib_currency_hpp
#define quantlib_currency_hpp

#include <ql/math/rounding.hpp>
#include <ql/errors.hpp>
#include <iosfwd>

namespace QuantLib {

    //! %Currency specification
    class Currency {
    public:
        //! default constructor
        /*! Instances built via this constructor have undefined
            behavior. Such instances can only act as placeholders
            and must be reassigned to a valid currency before being
            used.
        */
        Currency();

        //! \name Inspectors
        //@{
        //! currency name, e.g, "U.S. Dollar"
        const std::string &name() const;

        //! ISO 4217 three-letter code, e.g, "USD"
        const std::string &code() const;

        //! ISO 4217 numeric code, e.g, "840"
        Integer numericCode() const;

        //! symbol, e.g, "$"
        const std::string &symbol() const;

        //! fraction symbol, e.g, "¢"
        const std::string &fractionSymbol() const;

        //! number of fractionary parts in a unit, e.g, 100
        Integer fractionsPerUnit() const;

        //! rounding convention
        const Rounding &rounding() const;
        //! output format
        /*! The format will be fed three positional parameters,
            namely, value, code, and symbol, in this order.
        */
        std::string format() const;

        //@}
        //! \name Other information
        //@{
        //! is this a usable instance?
        bool empty() const;

        //! currency used for triangulated exchange when required
        const Currency &triangulationCurrency() const;
        //@}
    protected:
        struct Data;
        std::shared_ptr<Data> data_;
    };

    struct Currency::Data {
        std::string name_, code_;
        Integer numeric_;
        std::string symbol_, fractionSymbol_;
        Integer fractionsPerUnit_;
        Rounding rounding_;
        Currency triangulated_;
        std::string formatString_;

        Data(const std::string &name,
             const std::string &code,
             Integer numericCode,
             const std::string &symbol,
             const std::string &fractionSymbol,
             Integer fractionsPerUnit,
             const Rounding &rounding,
             const std::string &formatString,
             const Currency &triangulationCurrency = Currency());
    };

    /*! \relates Currency */
    bool operator==(const Currency &,
                    const Currency &);

    /*! \relates Currency */
    bool operator!=(const Currency &,
                    const Currency &);

    std::string currency_format(const Currency & c, Decimal value);

    /*! \relates Currency */
    std::ostream &operator<<(std::ostream &,
                             const Currency &);


    // inline definitions

    inline Currency::Currency() {}

    inline const std::string &Currency::name() const {
        return data_->name_;
    }

    inline const std::string &Currency::code() const {
        return data_->code_;
    }

    inline Integer Currency::numericCode() const {
        return data_->numeric_;
    }

    inline const std::string &Currency::symbol() const {
        return data_->symbol_;
    }

    inline const std::string &Currency::fractionSymbol() const {
        return data_->fractionSymbol_;
    }

    inline Integer Currency::fractionsPerUnit() const {
        return data_->fractionsPerUnit_;
    }

    inline const Rounding &Currency::rounding() const {
        return data_->rounding_;
    }

    inline std::string Currency::format() const {
        return data_->formatString_;
    }

    inline bool Currency::empty() const {
        return !data_;
    }

    inline const Currency &Currency::triangulationCurrency() const {
        return data_->triangulated_;
    }

    inline bool operator==(const Currency &c1, const Currency &c2) {
        return c1.name() == c2.name();
    }

    inline bool operator!=(const Currency &c1, const Currency &c2) {
        return !(c1 == c2);
    }

}


#endif
