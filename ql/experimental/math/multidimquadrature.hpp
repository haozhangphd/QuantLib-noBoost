/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Jose Aparicio

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

#ifndef quantlib_math_multidimquadrature_hpp
#define quantlib_math_multidimquadrature_hpp

#include <ql/qldefines.hpp>

/* Currently, this doesn't compile under Sun C++ (see
   https://github.com/lballabio/QuantLib/issues/223).  Until that's
   fixed, we disable it so that the rest of the library can be built.
*/

#ifndef QL_PATCH_SOLARIS

#include <ql/math/integrals/gaussianquadratures.hpp>
#include <functional>

namespace QuantLib {

    namespace detail {
        typedef std::vector<Real> DispArray;
    }

    /*! \brief Integrates a vector or scalar function of vector domain. 
        
        A template recursion along dimensions avoids calling depth 
        test or virtual functions.

        \todo Add coherence test between the integrand function dimensions (the
        vector size) and the declared dimension in the constructor.

        \todo Split into integrator classes for functions returning scalar and 
            vector?
    */
    class GaussianQuadMultidimIntegrator {
    private:
        // Vector integration. Quadrature to functions returning a vector of 
        // real numbers, turns 1D quadratures into ND
        class VectorIntegrator : public GaussHermiteIntegration {
        public:
            VectorIntegrator(Size n, Real mu = 0.0) 
            : GaussHermiteIntegration(n, mu) {}

            template <class F> // todo: fix copies.
            detail::DispArray operator()(const F& f) const {
                //first one, we do not know the size of the vector returned by f
                Integer i = order()-1;
                std::vector<Real> term = f(x_[i]);// potential copy! @#$%^!!!
                std::for_each(term.begin(), term.end(),
                              [this, i](Real r){return this->w_[i]* r;});
                std::vector<Real> sum = term;
           
                for (i--; i >= 0; --i) {
                    term = f(x_[i]);// potential copy! @#$%^!!!
                    // sum[j] += term[j] * w_[i];
                    std::transform(term.begin(), term.end(), sum.begin(), 
                        sum.begin(),[this, i](Real r1, Real r2) {return this->w_[i] * r1 + r2;});
                }
                return sum;
            }
        };

    public:
        /*!
            @param dimension The number of dimensions of the argument of the 
            function we want to integrate.
            @param quadOrder Quadrature order.
            @param mu Parameter in the Gauss Hermite weight (i.e. points load).
        */
        GaussianQuadMultidimIntegrator(Size dimension, Size quadOrder, 
            Real mu = 0.);
        //! Integration quadrature order.
        Size order() const {return integralV_.order();}

        //! Integrates function f over \f$ R^{dim} \f$
        /* This function is just syntax since the only thing it does is calling 
        to integrate<RetType> which has to exist for the type returned by the 
        function. So theres one redundant call but there should not be any extra 
        cost... up to the compiler. It can not be templated all the way since
        the integration entries functions can not be templates.
        Most times integrands will return a scalar or vector but could be a 
        matrix too. Also vectors might be returned as vector.
         */
        template<class RetType_T>
        RetType_T operator()(std::function<RetType_T (
            const std::vector<Real>& arg)> f) const
        {
            return integrate<RetType_T>(f);
        }


        //---------------------------------------------------------
        // Declare, spezializations follow.
        template<class RetType_T>
        RetType_T integrate(std::function<RetType_T (
            const std::vector<Real>& v1)> f) const;

    private:
        /* The maximum number of dimensions of the integration variable domain
            A higher than this number of dimension would presumably be 
           impractical and another integration algorithm (MC) should be 
           considered.
           \to do Consider moving it to a library configuration variable.
        */
        static const Size maxDimensions_ = 15;

        //! \name Integration entry points generation
        //@{
        //! Recursive template methods to statically generate (at this 
        //    class construction time) handles to the integration entry points
        template<Size levelSpawn>
        void spawnFcts() const {

            integrationEntries_[levelSpawn-1] =
                    [this](auto i1, auto i2){return this->scalarIntegrator<levelSpawn>(i1, i2);};

            integrationEntriesVR_[levelSpawn-1] =
                    [this](auto i1, auto i2){return this->vectorIntegratorVR<levelSpawn>(i1, i2);};

            spawnFcts<levelSpawn-1>();
        }
        //@}

        //---------------------------------------------------------

        template <int intgDepth>
        Real scalarIntegrator(
            std::function<Real (const std::vector<Real>& arg1)> f,
            const Real mFctr) const
        {
            varBuffer_[intgDepth-1] = mFctr;
            return integral_(
            [this, &f](Real r){return this->scalarIntegrator<intgDepth-1>(f, r);});
        }

        template <int intgDepth>
        detail::DispArray vectorIntegratorVR(
            std::function<detail::DispArray(const std::vector<Real>& arg1)> f,
            const Real mFctr) const 
        {
            varBuffer_[intgDepth-1] = mFctr;
            return 
              integralV_([this,f](auto r){return this->vectorIntegratorVR<intgDepth-1>(f, r);});
        }
    private:
        // Same object for all dimensions poses problems when using the 
        //   parallelized integrals version.
        //! The actual integrators.
        GaussHermiteIntegration integral_;
        VectorIntegrator integralV_;

        //! Buffer to allow acces to integrations. We do not know at which 
        //    level/dimension we are going to start integration
        // \todo Declare typedefs for traits
        mutable std::vector<
        std::function<Real (std::function<Real (
            const std::vector<Real>&)>,
            const Real)> > integrationEntries_;
        mutable std::vector<
        std::function<detail::DispArray (std::function<detail::DispArray(
            const std::vector<Real>&)>,
            const Real)> > integrationEntriesVR_;

        Size dimension_;
        // integration veriable buffer
        mutable std::vector<Real> varBuffer_;
    };


    // Template specializations ---------------------------------------------

    template<>
    inline Real GaussianQuadMultidimIntegrator::operator()(
        std::function<Real (const std::vector<Real>& v1)> f) const
    {
        return integral_(
        [this, &f](auto r){return integrationEntries_[dimension_-1](f, r);});
    }

    // Scalar integrand version (merge with vector case?)
    template<>
    inline Real GaussianQuadMultidimIntegrator::integrate<Real>(
        std::function<Real (const std::vector<Real>& v1)> f) const
    {
        // integration variables
        // call vector quadrature integration with the function and start 
        // values, kicks in recursion over the dimensions of the integration
        // variable.
        return integral_(
       [this, &f](auto r){return this->integrationEntries_[dimension_-1](f, r);});
    }

    // Vector integrand version
    template<>
    inline detail::DispArray GaussianQuadMultidimIntegrator::integrate<detail::DispArray>(
        std::function<detail::DispArray(const std::vector<Real>&)> f) const
    {
        return integralV_([this, &f](auto v){return (this->integrationEntriesVR_[dimension_-1])(f,v);});
    }

    //! Terminal integrand; scalar function version
    template<> 
    inline Real GaussianQuadMultidimIntegrator::scalarIntegrator<1>(
        std::function<Real (const std::vector<Real>& arg1)> f,
        const Real mFctr) const
    {
        varBuffer_[0] = mFctr;
        return f(varBuffer_);
    }

    //! Terminal integrand; vector function version
    template<>
    inline detail::DispArray
        GaussianQuadMultidimIntegrator::vectorIntegratorVR<1>(
        std::function<detail::DispArray (const std::vector<Real>& arg1)> f,
        const Real mFctr) const 
    {
        varBuffer_[0] = mFctr;
        return f(varBuffer_);
    }

    //! Terminal level:
    template<>
    inline void GaussianQuadMultidimIntegrator::spawnFcts<1>() const {
        integrationEntries_[0] =
                [this](auto i1, auto i2){return this->scalarIntegrator<1>(i1, i2);};
        integrationEntriesVR_[0] =
                [this](auto i1, auto i2){return this->vectorIntegratorVR<1>(i1, i2);};
    }

}

#endif

#endif
