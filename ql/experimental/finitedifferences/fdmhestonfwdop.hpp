/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2012, 2013, 2015 Klaus Spanderen
 Copyright (C) 2014, 2015 Johannes Göttker-Schnetmann

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

/*! \file fdmhestonfwdop.hpp
    \brief Heston Fokker-Planck forward operator
*/

#ifndef quantlib_fdm_heston_fwd_op_hpp
#define quantlib_fdm_heston_fwd_op_hpp

#include <ql/types.hpp>
#include <ql/experimental/finitedifferences/fdmsquarerootfwdop.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp>
#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>

namespace QuantLib {

    class YieldTermStructure;
    class FdmMesher;
    class HestonProcess;
    class FirstDerivativeOp;
    class TripleBandLinearOp;
    class NinePointLinearOp;

    class FdmHestonFwdOp : public FdmLinearOpComposite {
      public:
        FdmHestonFwdOp(
            const std::shared_ptr<FdmMesher>& mesher,
            const std::shared_ptr<HestonProcess>& process,
            FdmSquareRootFwdOp::TransformationType type 
                = FdmSquareRootFwdOp::Plain,
            const std::shared_ptr<LocalVolTermStructure> & leverageFct
                = std::shared_ptr<LocalVolTermStructure>());

        Size size() const;
        void setTime(Time t1, Time t2);

        Array apply(const Array& r) const;
        Array apply_mixed(const Array& r) const;

        Array apply_direction(Size direction,
                                          const Array& r) const;
        Array solve_splitting(Size direction,
                                          const Array& r, Real s) const;
        Array preconditioner(const Array& r, Real s) const;

#if !defined(QL_NO_UBLAS_SUPPORT)
        std::vector<SparseMatrix>  toMatrixDecomp() const;
#endif
      private:
        Array getLeverageFctSlice(Time t1, Time t2) const;
        const FdmSquareRootFwdOp::TransformationType type_;
        const Real kappa_, theta_, sigma_, rho_, v0_;

        const std::shared_ptr<YieldTermStructure> rTS_;
        const std::shared_ptr<YieldTermStructure> qTS_;

        const Array varianceValues_;
        const std::shared_ptr<FirstDerivativeOp> dxMap_;
        const std::shared_ptr<ModTripleBandLinearOp> dxxMap_;
        const std::shared_ptr<ModTripleBandLinearOp> boundary_;
        Array L_;

        const std::shared_ptr<TripleBandLinearOp> mapX_;
        const std::shared_ptr<FdmSquareRootFwdOp> mapY_;

        const std::shared_ptr<NinePointLinearOp> correlation_;

        const std::shared_ptr<LocalVolTermStructure> leverageFct_;
        const std::shared_ptr<FdmMesher> mesher_;
        const Array x_;
    };
}

#endif
