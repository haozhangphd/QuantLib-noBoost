/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2013, 2015 Peter Caspers

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

/*! \file gaussian1dmodel.hpp
    \brief basic interface for one factor interest rate models
*/


#ifndef quantlib_gaussian1dmodel_hpp
#define quantlib_gaussian1dmodel_hpp

#include <functional>
#include <ql/models/model.hpp>
#include <ql/models/parameter.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/indexes/swapindex.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/time/date.hpp>
#include <ql/time/period.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/utilities/null.hpp>
#include <ql/patterns/lazyobject.hpp>
#include <unordered_map>

namespace {
	template <class T>
	inline void hash_combine(std::size_t& seed, const T& v)
	{
		std::hash<T> hasher;
		seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
	}
}

namespace QuantLib {

/*! One factor interest rate model interface class
    The only methods that must be implemented by subclasses
    are the numeraire and zerobond methods for an input array
    of state variable values. The variable $y$ is understood
    to be the standardized (zero mean, unit variance) version
    of the model's original state variable $x$.

    \warning the variance of the state process conditional on
    $x(t)=x$ must be independent of the value of $x$

*/

class Gaussian1dModel : public TermStructureConsistentModel, public LazyObject {
  public:
    const std::shared_ptr<StochasticProcess1D> stateProcess() const;

    Real numeraire(const Time t, const Real y = 0.0,
                   const Handle<YieldTermStructure> &yts =
                             Handle<YieldTermStructure>()) const;

    Real zerobond(
        const Time T, const Time t = 0.0, const Real y = 0.0,
        const Handle<YieldTermStructure> &yts = Handle<YieldTermStructure>()) const;

    Real numeraire(const Date &referenceDate, const Real y = 0.0,
                   const Handle<YieldTermStructure> &yts =
                             Handle<YieldTermStructure>()) const;

    Real zerobond(
        const Date &maturity, const Date &referenceDate = Null<Date>(),
        const Real y = 0.0,
        const Handle<YieldTermStructure> &yts = Handle<YieldTermStructure>()) const;

    Real zerobondOption(
        const Option::Type &type, const Date &expiry, const Date &valueDate,
        const Date &maturity, const Rate strike,
        const Date &referenceDate = Null<Date>(), const Real y = 0.0,
        const Handle<YieldTermStructure> &yts = Handle<YieldTermStructure>(),
        const Real yStdDevs = 7.0, const Size yGridPoints = 64,
        const bool extrapolatePayoff = true,
        const bool flatPayoffExtrapolation = false) const;

    Real forwardRate(
        const Date &fixing, const Date &referenceDate = Null<Date>(),
        const Real y = 0.0,
        std::shared_ptr<IborIndex> iborIdx = std::shared_ptr<IborIndex>()) const;

    Real swapRate(
        const Date &fixing, const Period &tenor,
        const Date &referenceDate = Null<Date>(), const Real y = 0.0,
        std::shared_ptr<SwapIndex> swapIdx = std::shared_ptr<SwapIndex>()) const;

    Real swapAnnuity(
        const Date &fixing, const Period &tenor,
        const Date &referenceDate = Null<Date>(), const Real y = 0.0,
        std::shared_ptr<SwapIndex> swapIdx = std::shared_ptr<SwapIndex>()) const;

    /*! Computes the integral
    \f[ {2\pi}^{-0.5} \int_{a}^{b} p(x) \exp{-0.5*x*x} \mathrm{d}x \f]
    with
    \f[ p(x) = ax^4+bx^3+cx^2+dx+e \f].
    */
    static Real gaussianPolynomialIntegral(const Real a, const Real b,
                                           const Real c, const Real d,
                                           const Real e, const Real x0,
                                           const Real x1);

    /*! Computes the integral
    \f[ {2\pi}^{-0.5} \int_{a}^{b} p(x) \exp{-0.5*x*x} \mathrm{d}x \f]
    with
    \f[ p(x) = a(x-h)^4+b(x-h)^3+c(x-h)^2+d(x-h)+e \f].
    */
    static Real
    gaussianShiftedPolynomialIntegral(const Real a, const Real b, const Real c,
                                      const Real d, const Real e, const Real h,
                                      const Real x0, const Real x1);

    /*! Generates a grid of values for the standardized state variable $y$
       at time $T$
        conditional on $y(t)=y$, covering yStdDevs standard deviations
       consisting of
        2*gridPoints+1 points */

    Array yGrid(const Real yStdDevs, const int gridPoints,
                                  const Real T = 1.0, const Real t = 0,
                                  const Real y = 0) const;

  private:
    // It is of great importance for performance reasons to cache underlying
    // swaps generated from indexes. In addition the indexes may only be given
    // as templates for the conventions with the tenor replaced by the actual
    // one later on.

    struct CachedSwapKey {
        const std::shared_ptr<SwapIndex> index;
        const Date fixing;
        const Period tenor;
        bool operator==(const CachedSwapKey &o) const {
            return index->name() == o.index->name() && fixing == o.fixing &&
                   tenor == o.tenor;
        }
    };

    struct CachedSwapKeyHasher {

        std::size_t operator()(const CachedSwapKey &x) const {
            std::size_t seed = 0;
            hash_combine(seed, x.index->name());
            hash_combine(seed, x.fixing.serialNumber());
            hash_combine(seed, x.tenor.length());
            hash_combine(seed, x.tenor.units());
            return seed;
        }
    };

    typedef std::unordered_map<CachedSwapKey, std::shared_ptr<VanillaSwap>,
                                 CachedSwapKeyHasher> CacheType;

    mutable CacheType swapCache_;

  protected:
    // we let derived classes register with the termstructure
    Gaussian1dModel(const Handle<YieldTermStructure> &yieldTermStructure)
        : TermStructureConsistentModel(yieldTermStructure) {
        registerWith(Settings::instance().evaluationDate());
    }

    virtual ~Gaussian1dModel() {}

    virtual Real
    numeraireImpl(const Time t, const Real y,
                  const Handle<YieldTermStructure> &yts) const = 0;

    virtual Real zerobondImpl(const Time T, const Time t, const Real y,
                              const Handle<YieldTermStructure> &yts) const = 0;

    void performCalculations() const {
        evaluationDate_ = Settings::instance().evaluationDate();
        enforcesTodaysHistoricFixings_ =
            Settings::instance().enforcesTodaysHistoricFixings();
    }

    void generateArguments() {
        calculate();
        notifyObservers();
    }

    // retrieve underlying swap from cache if possible, otherwise
    // create it and store it in the cache
    std::shared_ptr<VanillaSwap>
    underlyingSwap(const std::shared_ptr<SwapIndex> &index,
                   const Date &expiry, const Period &tenor) const {

        CachedSwapKey k = {index, expiry, tenor};
        CacheType::iterator i = swapCache_.find(k);
        if (i == swapCache_.end()) {
            std::shared_ptr<VanillaSwap> underlying =
                index->clone(tenor)->underlyingSwap(expiry);
            swapCache_.insert(std::make_pair(k, underlying));
            return underlying;
        }
        return i->second;
    }

    std::shared_ptr<StochasticProcess1D> stateProcess_;
    mutable Date evaluationDate_;
    mutable bool enforcesTodaysHistoricFixings_;
};

inline const std::shared_ptr<StochasticProcess1D>
Gaussian1dModel::stateProcess() const {

    QL_REQUIRE(stateProcess_ != NULL, "state process not set");
    return stateProcess_;
}

inline Real
Gaussian1dModel::numeraire(const Time t, const Real y,
                           const Handle<YieldTermStructure> &yts) const {

    return numeraireImpl(t, y, yts);
}

inline Real
Gaussian1dModel::zerobond(const Time T, const Time t, const Real y,
                          const Handle<YieldTermStructure> &yts) const {
    return zerobondImpl(T, t, y, yts);
}

inline Real
Gaussian1dModel::numeraire(const Date &referenceDate, const Real y,
                           const Handle<YieldTermStructure> &yts) const {

    return numeraire(termStructure()->timeFromReference(referenceDate), y, yts);
}

inline Real
Gaussian1dModel::zerobond(const Date &maturity, const Date &referenceDate,
                          const Real y, const Handle<YieldTermStructure> &yts) const {

    return zerobond(termStructure()->timeFromReference(maturity),
                    referenceDate != Null<Date>()
                        ? termStructure()->timeFromReference(referenceDate)
                        : 0.0,
                    y, yts);
}
}

#endif
