/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2009 Klaus Spanderen

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

/*! \file ranluxuniformrng.hpp
    \brief "Luxury" random number generator.
*/

#ifndef quantlib_ranlux_uniform_rng_h
#define quantlib_ranlux_uniform_rng_h

#include <ql/methods/montecarlo/sample.hpp>
#include <tr1/random.h>


using ranlux64_3_01 = std::tr1::discard_block<std::tr1::ranlux64_base_01, 223, 24>;
using ranlux64_4_01 = std::tr1::discard_block<std::tr1::ranlux64_base_01, 389, 24>;

namespace QuantLib {

    //! Uniform random number generator
    /*! M. Luescher's "luxury" random number generator

        Implementation is a proxy for the corresponding tr1 random
        number generator. For more detail see the tr1 documentation and:
          M.Luescher, A portable high-quality random number generator for
          lattice field theory simulations, Comp. Phys. Comm. 79 (1994) 100
          
        Available luxury levels:
        Ranlux3: Any theoretically possible correlations have very small change
                 of being observed.
        Ranlux4: highest possible luxury.         
    */
    class Ranlux3UniformRng {
      public:
        typedef Sample<Real> sample_type;
 
        explicit Ranlux3UniformRng(Size seed = 19780503u)
        : ranlux3_(std::tr1::ranlux64_base_01(seed)){ }

        sample_type next() const {
            return sample_type(ranlux3_(), 1.0);
        }

      private:
        mutable ranlux64_3_01 ranlux3_;
    };

    class Ranlux4UniformRng {
      public:
        typedef Sample<Real> sample_type;
 
        explicit Ranlux4UniformRng(Size seed = 19780503u)
        : ranlux4_(std::tr1::ranlux64_base_01(seed)){ }

        sample_type next() const {
            return sample_type(ranlux4_(), 1.0);
        }

      private:
        mutable ranlux64_4_01 ranlux4_;
    };
}


#endif
