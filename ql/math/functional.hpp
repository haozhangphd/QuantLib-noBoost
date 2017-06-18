/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 RiskMap srl

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

/*! \file functional.hpp
    \brief functionals and combinators not included in the STL
*/

#ifndef quantlib_functional_hpp
#define quantlib_functional_hpp

#include <ql/types.hpp>
#include <cmath>
#include <functional>

namespace QuantLib {

    auto constant = [](auto c){return [c](auto x){return c;};};
    auto identity = [](auto x) { return x; };
    auto square = [](auto x) { return x * x; };
    auto cube = [](auto x) { return x * x * x; };
    auto fourth_power = [](auto x) { return x * x * x * x; };
    auto everywhere = [](auto x) { return true; };
    auto nowhere = [](auto x) { return false; };
    auto equal_whithin = [](auto eps){return [eps](auto a, auto b){return std::fabs(a-b) <= eps;};};



}


#endif
