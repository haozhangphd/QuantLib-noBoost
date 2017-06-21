/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Johannes Göttker-Schnetmann
 Copyright (C) 2015 Klaus Spanderen

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


#include <ql/timegrid.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/models/equity/hestonmodel.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/integrals/discreteintegrals.hpp>
#include <ql/math/interpolations/bilinearinterpolation.hpp>
#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/fixedlocalvolsurface.hpp>
#include <ql/methods/finitedifferences/schemes/all.hpp>
#include <ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/meshers/predefined1dmesher.hpp>
#include <ql/methods/finitedifferences/meshers/concentrating1dmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/utilities/fdmmesherintegral.hpp>
#include <ql/experimental/models/hestonslvfdmmodel.hpp>
#include <ql/experimental/finitedifferences/fdmhestonfwdop.hpp>
#include <ql/experimental/finitedifferences/localvolrndcalculator.hpp>
#include <ql/experimental/finitedifferences/squarerootprocessrndcalculator.hpp>

#include <memory>
#include <algorithm>
#include <memory>
#include <algorithm>

#include <functional>


namespace QuantLib {

    namespace {
        std::shared_ptr<Fdm1dMesher> varianceMesher(
            const SquareRootProcessRNDCalculator& rnd,
            Time t0, Time t1, Size vGrid,
            Real v0, const HestonSLVFokkerPlanckFdmParams& params) {

            std::vector<std::tuple<Real, Real, bool> > cPoints;

            const Real v0Density = params.v0Density;
            const Real upperBoundDensity = params.vUpperBoundDensity;
            const Real lowerBoundDensity = params.vLowerBoundDensity;

            Real lowerBound = Null<Real>(), upperBound = -Null<Real>();

            for (Size i=0; i <= 10; ++i) {
                const Time t = t0 + i/10.0*(t1-t0);
                lowerBound = std::min(
                    lowerBound, rnd.invcdf(params.vLowerEps, t));
                upperBound = std::max(
                    upperBound, rnd.invcdf(1.0-params.vUpperEps, t));
            }

            lowerBound = std::max(lowerBound, params.vMin);
            switch (params.trafoType) {
                case FdmSquareRootFwdOp::Log:
                  {
                    lowerBound = std::log(lowerBound);
                    upperBound = std::log(upperBound);

                    const Real v0Center = std::log(v0);

                    cPoints = {
                        std::make_tuple(lowerBound, lowerBoundDensity, false),
                        std::make_tuple(v0Center, v0Density, true),
                        std::make_tuple(upperBound, upperBoundDensity, false)};

                    return std::make_shared<Concentrating1dMesher>(
                        lowerBound, upperBound, vGrid, cPoints, 1e-8);
                  }
                break;
                case FdmSquareRootFwdOp::Plain:
                  {
                      const Real v0Center = v0;

                      cPoints = {
                          std::make_tuple(lowerBound, lowerBoundDensity, false),
                          std::make_tuple(v0Center, v0Density, true),
                          std::make_tuple(upperBound, upperBoundDensity, false)};

                      return std::make_shared<Concentrating1dMesher>(
                          lowerBound, upperBound, vGrid, cPoints, 1e-8);
                  }
                break;
                case FdmSquareRootFwdOp::Power:
                {
                    const Real v0Center = v0;

                    cPoints = {
                        std::make_tuple(lowerBound, lowerBoundDensity, false),
                        std::make_tuple(v0Center, v0Density, true),
                        std::make_tuple(upperBound, upperBoundDensity, false)};

                    return std::make_shared<Concentrating1dMesher>(
                        lowerBound, upperBound, vGrid, cPoints, 1e-8);
                }
                break;
                default:
                    QL_FAIL("transformation type is not implemented");
            }
        }

        Real integratePDF(const Array& p,
                          const std::shared_ptr<FdmMesherComposite>& mesher,
                          FdmSquareRootFwdOp::TransformationType trafoType,
                          Real alpha) {

            if (trafoType != FdmSquareRootFwdOp::Power) {
                return FdmMesherIntegral(
                        mesher, DiscreteSimpsonIntegral()).integrate(p);
            }
            else {
                Array tp(p.size());
                const FdmLinearOpIterator end = mesher->layout()->end();
                for (FdmLinearOpIterator iter = mesher->layout()->begin();
                    iter != end; ++iter) {
                    const Size idx = iter.index();
                    const Real nu = mesher->location(iter, 1);

                    tp[idx] = p[idx]*std::pow(nu, alpha-1);
                }

                return FdmMesherIntegral(
                        mesher, DiscreteSimpsonIntegral()).integrate(tp);
            }
        }


        Array rescalePDF(
            const Array& p,
            const std::shared_ptr<FdmMesherComposite>& mesher,
            FdmSquareRootFwdOp::TransformationType trafoType, Real alpha) {

            Array retVal = p/integratePDF(p, mesher, trafoType, alpha);

            return retVal;
        }


        template <class Interpolator>
        Array reshapePDF(
            const Array& p,
            const std::shared_ptr<FdmMesherComposite>& oldMesher,
            const std::shared_ptr<FdmMesherComposite>& newMesher,
            const Interpolator& interp = Interpolator()) {

            const std::shared_ptr<FdmLinearOpLayout> oldLayout
                = oldMesher->layout();
            const std::shared_ptr<FdmLinearOpLayout> newLayout
                = newMesher->layout();

            QL_REQUIRE(   oldLayout->size() == newLayout->size()
                       && oldLayout->size() == p.size(),
                       "inconsistent mesher or vector size given");

            Matrix m(oldLayout->dim()[1], oldLayout->dim()[0]);
            for (Size i=0; i < m.rows(); ++i) {
                std::copy(p.begin() + i*m.columns(),
                          p.begin() + (i+1)*m.columns(), m.row_begin(i));
            }
            const Interpolation2D interpol = interp.interpolate(
                oldMesher->getFdm1dMeshers()[0]->locations().begin(),
                oldMesher->getFdm1dMeshers()[0]->locations().end(),
                oldMesher->getFdm1dMeshers()[1]->locations().begin(),
                oldMesher->getFdm1dMeshers()[1]->locations().end(), m);

            Array pNew(p.size());
            const FdmLinearOpIterator endIter = newLayout->end();
            for (FdmLinearOpIterator iter = newLayout->begin();
                iter != endIter; ++iter) {
                const Real x = newMesher->location(iter, 0);
                const Real v = newMesher->location(iter, 1);

                if (   x > interpol.xMax() || x < interpol.xMin()
                    || v > interpol.yMax() || v < interpol.yMin() ) {
                    pNew[iter.index()] = 0;
                }
                else {
                    pNew[iter.index()] = interpol(x, v);
                }
            }

            return pNew;
        }

        class FdmScheme {
          public:
            virtual ~FdmScheme() {}
            virtual void step(Array& a, Time t) = 0;
            virtual void setStep(Time dt) = 0;
        };

        template <class T>
        class FdmSchemeWrapper : public FdmScheme {
          public:
            explicit FdmSchemeWrapper(T* scheme)
            : scheme_(scheme) { }

            void step(Array& a, Time t) {
                scheme_->step(a, t);
            }
            void setStep(Time dt) {
                scheme_->setStep(dt);
            }

          private:
            const std::unique_ptr<T> scheme_;
        };

        std::shared_ptr<FdmScheme> fdmSchemeFactory(
            const FdmSchemeDesc desc,
            const std::shared_ptr<FdmLinearOpComposite>& op) {

            switch (desc.type) {
              case FdmSchemeDesc::HundsdorferType:
                  return std::shared_ptr<FdmScheme>(
                      new FdmSchemeWrapper<HundsdorferScheme>(
                          new HundsdorferScheme(desc.theta, desc.mu, op)));
              case FdmSchemeDesc::DouglasType:
                  return std::shared_ptr<FdmScheme>(
                      new FdmSchemeWrapper<DouglasScheme>(
                          new DouglasScheme(desc.theta, op)));
              case FdmSchemeDesc::CraigSneydType:
                  return std::shared_ptr<FdmScheme>(
                      new FdmSchemeWrapper<CraigSneydScheme>(
                          new CraigSneydScheme(desc.theta, desc.mu, op)));
              case FdmSchemeDesc::ModifiedCraigSneydType:
                  return std::shared_ptr<FdmScheme>(
                     new FdmSchemeWrapper<ModifiedCraigSneydScheme>(
                          new ModifiedCraigSneydScheme(
                              desc.theta, desc.mu, op)));
              case FdmSchemeDesc::ImplicitEulerType:
                  return std::shared_ptr<FdmScheme>(
                      new FdmSchemeWrapper<ImplicitEulerScheme>(
                          new ImplicitEulerScheme(op)));
              case FdmSchemeDesc::ExplicitEulerType:
                  return std::shared_ptr<FdmScheme>(
                      new FdmSchemeWrapper<ExplicitEulerScheme>(
                          new ExplicitEulerScheme(op)));
              default:
                  QL_FAIL("Unknown scheme type");
            }
        }
    }

    HestonSLVFDMModel::HestonSLVFDMModel(
        const Handle<LocalVolTermStructure>& localVol,
        const Handle<HestonModel>& hestonModel,
        const Date& endDate,
        const HestonSLVFokkerPlanckFdmParams& params,
        const bool logging,
        const std::vector<Date>& mandatoryDates)
    : localVol_(localVol),
      hestonModel_(hestonModel),
      endDate_(endDate),
      params_(params),
      mandatoryDates_(mandatoryDates),
      logging_(logging) {

        registerWith(localVol_);
        registerWith(hestonModel_);
    }

    std::shared_ptr<HestonProcess> HestonSLVFDMModel::hestonProcess() const {
        return hestonModel_->process();
    }

    std::shared_ptr<LocalVolTermStructure> HestonSLVFDMModel::localVol() const {
        return localVol_.currentLink();
    }

    std::shared_ptr<LocalVolTermStructure>
    HestonSLVFDMModel::leverageFunction() const {
        calculate();

        return leverageFunction_;
    }

    void HestonSLVFDMModel::performCalculations() const {
        logEntries_.clear();

        const std::shared_ptr<HestonProcess> hestonProcess
            = hestonModel_->process();
        const std::shared_ptr<Quote> spot
            = hestonProcess->s0().currentLink();
        const std::shared_ptr<YieldTermStructure> rTS
            = hestonProcess->riskFreeRate().currentLink();
        const std::shared_ptr<YieldTermStructure> qTS
            = hestonProcess->dividendYield().currentLink();

        const Real v0    = hestonProcess->v0();
        const Real kappa = hestonProcess->kappa();
        const Real theta = hestonProcess->theta();
        const Real sigma = hestonProcess->sigma();
        const Real alpha = 2*kappa*theta/(sigma*sigma);

        const Size xGrid = params_.xGrid;
        const Size vGrid = params_.vGrid;

        const DayCounter dc = rTS->dayCounter();
        const Date referenceDate = rTS->referenceDate();

        const Time T = dc.yearFraction(referenceDate, endDate_);

        QL_REQUIRE(referenceDate < endDate_,
            "reference date must be smaller than final calibration date");

        QL_REQUIRE(localVol_->maxTime() >= T,
            "final calibration maturity exceeds local volatility surface");

        // set-up exponential time step scheme
        const Time maxDt = 1.0/params_.tMaxStepsPerYear;
        const Time minDt = 1.0/params_.tMinStepsPerYear;

        Time tIdx=0.0;
        std::vector<Time> times(1, tIdx);
        times.reserve(Size(T*params_.tMinStepsPerYear));
        while (tIdx < T) {
            const Real decayFactor = std::exp(-params_.tStepNumberDecay*tIdx);
            const Time dt = maxDt*decayFactor + minDt*(1.0-decayFactor);

            times.emplace_back(std::min(T, tIdx+=dt));
        }

        for (Size i=0; i < mandatoryDates_.size(); ++i) {
            times.emplace_back(
                dc.yearFraction(referenceDate, mandatoryDates_[i]));
        }

        const std::shared_ptr<TimeGrid> timeGrid(
            new TimeGrid(times.begin(), times.end()));

        // build 1d meshers
        const LocalVolRNDCalculator localVolRND(
            spot, rTS, qTS, localVol_.currentLink(),
            timeGrid, xGrid,
            params_.x0Density,
            params_.localVolEpsProb,
            params_.maxIntegrationIterations);

        const std::vector<Size> rescaleSteps
            = localVolRND.rescaleTimeSteps();

        const SquareRootProcessRNDCalculator squareRootRnd(
            v0, kappa, theta, sigma);

        const FdmSquareRootFwdOp::TransformationType trafoType
          = params_.trafoType;

        std::vector<std::shared_ptr<Fdm1dMesher> > xMesher, vMesher;
        xMesher.reserve(timeGrid->size());
        vMesher.reserve(timeGrid->size());

        xMesher.emplace_back(localVolRND.mesher(0.0));
        vMesher.emplace_back(std::make_shared<Predefined1dMesher>(
            std::vector<Real>(vGrid, v0)));

        Size rescaleIdx = 0;
        for (Size i=1; i < timeGrid->size(); ++i) {
            xMesher.emplace_back(localVolRND.mesher(timeGrid->at(i)));

            if (i == rescaleSteps[rescaleIdx]) {
                ++rescaleIdx;
                vMesher.emplace_back(varianceMesher(squareRootRnd,
                    timeGrid->at(rescaleSteps[rescaleIdx-1]),
                    (rescaleIdx < rescaleSteps.size())
                        ? timeGrid->at(rescaleSteps[rescaleIdx])
                        : timeGrid->back(),
                    vGrid, v0, params_));
            }
            else
                vMesher.emplace_back(vMesher.back());
        }

        // start probability distribution
        std::shared_ptr<FdmMesherComposite> mesher
            = std::make_shared<FdmMesherComposite>(
                xMesher.at(1), vMesher.at(1));

        const Volatility lv0
            = localVol_->localVol(0.0, spot->value())/std::sqrt(v0);

        std::shared_ptr<Matrix> L(new Matrix(xGrid, timeGrid->size()));

        const Real l0 = lv0;
        std::fill(L->column_begin(0),L->column_end(0), l0);
        std::fill(L->column_begin(1),L->column_end(1), l0);

        // create strikes from meshers
        std::vector<std::shared_ptr<std::vector<Real> > > vStrikes(
            timeGrid->size());

        for (Size i=0; i < timeGrid->size(); ++i) {
            vStrikes[i] = std::make_shared<std::vector<Real> >(xGrid);
            std::transform(xMesher[i]->locations().begin(),
                           xMesher[i]->locations().end(),
                           vStrikes[i]->begin(),
                           std::ptr_fun<Real, Real>(std::exp));
        }

        const std::shared_ptr<FixedLocalVolSurface> leverageFct(
            new FixedLocalVolSurface(referenceDate, times, vStrikes, L, dc));

        std::shared_ptr<FdmLinearOpComposite> hestonFwdOp(
            new FdmHestonFwdOp(mesher, hestonProcess, trafoType, leverageFct));

        Array p = FdmHestonGreensFct(mesher, hestonProcess, trafoType, lv0)
            .get(timeGrid->at(1), params_.greensAlgorithm);

        if (logging_) {
            const LogEntry entry = { timeGrid->at(1),
                std::shared_ptr<Array>(new Array(p)), mesher };
            logEntries_.emplace_back(entry);
        }

        for (Size i=2; i < times.size(); ++i) {
            const Time t = timeGrid->at(i);
            const Time dt = t - timeGrid->at(i-1);

            if (   mesher->getFdm1dMeshers()[0] != xMesher[i]
                || mesher->getFdm1dMeshers()[1] != vMesher[i]) {
                const std::shared_ptr<FdmMesherComposite> newMesher(
                    new FdmMesherComposite(xMesher[i], vMesher[i]));

                p = reshapePDF<Bilinear>(p, mesher, newMesher);
                mesher = newMesher;

                p = rescalePDF(p, mesher, trafoType, alpha);

                hestonFwdOp = std::shared_ptr<FdmLinearOpComposite>(
                                new FdmHestonFwdOp(mesher, hestonProcess,
                                               trafoType, leverageFct));
            }

            Array pn = p;
            const Array x(Exp(
                Array(mesher->getFdm1dMeshers()[0]->locations().begin(),
                      mesher->getFdm1dMeshers()[0]->locations().end())));
            const Array v(
                    mesher->getFdm1dMeshers()[1]->locations().begin(),
                    mesher->getFdm1dMeshers()[1]->locations().end());

            // predictor corrector steps
            for (Size r=0; r < params_.predictionCorretionSteps; ++r) {
                const std::shared_ptr<FdmScheme> fdmScheme(
                    fdmSchemeFactory(params_.schemeDesc, hestonFwdOp));

                for (Size j=0; j < x.size(); ++j) {
                    Array pSlice(vGrid);
                    for (Size k=0; k < vGrid; ++k)
                        pSlice[k] = pn[j + k*xGrid];

                    const Real pInt = (trafoType == FdmSquareRootFwdOp::Power)
                       ? DiscreteSimpsonIntegral()(v, Pow(v, alpha-1)*pSlice)
                       : DiscreteSimpsonIntegral()(v, pSlice);

                    const Real vpInt = (trafoType == FdmSquareRootFwdOp::Log)
                      ? DiscreteSimpsonIntegral()(v, Exp(v)*pSlice)
                      : (trafoType == FdmSquareRootFwdOp::Power)
                      ? DiscreteSimpsonIntegral()(v, Pow(v, alpha)*pSlice)
                      : DiscreteSimpsonIntegral()(v, v*pSlice);

                    const Real scale = pInt/vpInt;
                    const Volatility localVol = localVol_->localVol(t, x[j]);

                    const Real l = (scale >= 0.0)
                      ? localVol*std::sqrt(scale) : 1.0;

                    (*L)[j][i] = std::min(50.0, std::max(0.001, l));

                    leverageFct->setInterpolation(Linear());
                }

                const Real sLowerBound = std::max(x.front(),
                    std::exp(localVolRND.invcdf(
                        params_.leverageFctPropEps, t)));
                const Real sUpperBound = std::min(x.back(),
                    std::exp(localVolRND.invcdf(
                        1.0-params_.leverageFctPropEps, t)));

                const Real lowerL = leverageFct->localVol(t, sLowerBound);
                const Real upperL = leverageFct->localVol(t, sUpperBound);

                for (Size j=0; j < x.size(); ++j) {
                    if (x[j] < sLowerBound)
                        std::fill(L->row_begin(j)+i,
                          std::min(L->row_begin(j)+i+1, L->row_end(j)),
                          lowerL);
                    else if (x[j] > sUpperBound)
                        std::fill(L->row_begin(j)+i,
                          std::min(L->row_begin(j)+i+1, L->row_end(j)),
                          upperL);
                    else if ((*L)[j][i] == Null<Real>())
                        QL_FAIL("internal error");
                }
                leverageFct->setInterpolation(Linear());

                pn = p;

                fdmScheme->setStep(dt);
                fdmScheme->step(pn, t);
            }
            p = pn;
            p = rescalePDF(p, mesher, trafoType, alpha);

            if (logging_) {
                const LogEntry entry
                    = { t, std::shared_ptr<Array>(new Array(p)), mesher };
                logEntries_.emplace_back(entry);
            }
        }

        leverageFunction_ = leverageFct;
    }

    const std::list<HestonSLVFDMModel::LogEntry>& HestonSLVFDMModel::logEntries()
    const {
        performCalculations();
        return logEntries_;
    }
}

