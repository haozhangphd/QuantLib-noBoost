/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Andreas Gaida
 Copyright (C) 2008 Ralph Schreyer
 Copyright (C) 2008 Klaus Spanderen

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

#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/operators/ninepointlinearop.hpp>

namespace QuantLib {

    NinePointLinearOp::NinePointLinearOp(
        Size d0, Size d1,
        const std::shared_ptr<FdmMesher>& mesher)
    : d0_(d0), d1_(d1),
      i00_(mesher->layout()->size()),
      i10_(mesher->layout()->size()),
      i20_(mesher->layout()->size()),
      i01_(mesher->layout()->size()),
      i21_(mesher->layout()->size()),
      i02_(mesher->layout()->size()),
      i12_(mesher->layout()->size()),
      i22_(mesher->layout()->size()),
      a00_(mesher->layout()->size()),
      a10_(mesher->layout()->size()),
      a20_(mesher->layout()->size()),
      a01_(mesher->layout()->size()),
      a11_(mesher->layout()->size()),
      a21_(mesher->layout()->size()),
      a02_(mesher->layout()->size()),
      a12_(mesher->layout()->size()),
      a22_(mesher->layout()->size()),
      mesher_(mesher) {

        QL_REQUIRE(   d0_ != d1_
            && d0_ < mesher->layout()->dim().size()
            && d1_ < mesher->layout()->dim().size(),
            "inconsistent derivative directions");

        const std::shared_ptr<FdmLinearOpLayout> layout = mesher->layout();
        const FdmLinearOpIterator endIter = layout->end();

        for (FdmLinearOpIterator iter = layout->begin(); iter!=endIter; ++iter) {
            const Size i = iter.index();

            i10_[i] = layout->neighbourhood(iter, d1_, -1);
            i01_[i] = layout->neighbourhood(iter, d0_, -1);
            i21_[i] = layout->neighbourhood(iter, d0_,  1);
            i12_[i] = layout->neighbourhood(iter, d1_,  1);
            i00_[i] = layout->neighbourhood(iter, d0_, -1, d1_, -1);
            i20_[i] = layout->neighbourhood(iter, d0_,  1, d1_, -1);
            i02_[i] = layout->neighbourhood(iter, d0_, -1, d1_,  1);
            i22_[i] = layout->neighbourhood(iter, d0_,  1, d1_,  1);
        }
    }

    NinePointLinearOp::NinePointLinearOp(const NinePointLinearOp& m)
    : i00_(m.mesher_->layout()->size()),
      i10_(m.mesher_->layout()->size()),
      i20_(m.mesher_->layout()->size()),
      i01_(m.mesher_->layout()->size()),
      i21_(m.mesher_->layout()->size()),
      i02_(m.mesher_->layout()->size()),
      i12_(m.mesher_->layout()->size()),
      i22_(m.mesher_->layout()->size()),
      a00_(m.mesher_->layout()->size()),
      a10_(m.mesher_->layout()->size()),
      a20_(m.mesher_->layout()->size()),
      a01_(m.mesher_->layout()->size()),
      a11_(m.mesher_->layout()->size()),
      a21_(m.mesher_->layout()->size()),
      a02_(m.mesher_->layout()->size()),
      a12_(m.mesher_->layout()->size()),
      a22_(m.mesher_->layout()->size()),
      mesher_(m.mesher_) {

        const Size size = mesher_->layout()->size();
        std::copy(m.i00_.begin(), m.i00_.end(), i00_.begin());
        std::copy(m.i10_.begin(), m.i10_.end(), i10_.begin());
        std::copy(m.i20_.begin(), m.i20_.end(), i20_.begin());
        std::copy(m.i01_.begin(), m.i01_.end(), i01_.begin());
        std::copy(m.i21_.begin(), m.i21_.end(), i21_.begin());
        std::copy(m.i02_.begin(), m.i02_.end(), i02_.begin());
        std::copy(m.i12_.begin(), m.i12_.end(), i12_.begin());
        std::copy(m.i22_.begin(), m.i22_.end(), i22_.begin());
        std::copy(m.a00_.begin(), m.a00_.end(), a00_.begin());
        std::copy(m.a10_.begin(), m.a10_.end(), a10_.begin());
        std::copy(m.a20_.begin(), m.a20_.end(), a20_.begin());
        std::copy(m.a01_.begin(), m.a01_.end(), a01_.begin());
        std::copy(m.a11_.begin(), m.a11_.end(), a11_.begin());
        std::copy(m.a21_.begin(), m.a21_.end(), a21_.begin());
        std::copy(m.a02_.begin(), m.a02_.end(), a02_.begin());
        std::copy(m.a12_.begin(), m.a12_.end(), a12_.begin());
        std::copy(m.a22_.begin(), m.a22_.end(), a22_.begin());
    }

    NinePointLinearOp& NinePointLinearOp::operator=(
        const NinePointLinearOp& m) {
        NinePointLinearOp temp(m);
        swap(temp);
        return *this;
    }

    NinePointLinearOp& NinePointLinearOp::operator=(
        const Disposable<NinePointLinearOp>& m) {
        swap(const_cast<Disposable<NinePointLinearOp>&>(m));
        return *this;
    }

    NinePointLinearOp::NinePointLinearOp(
        const Disposable<NinePointLinearOp>& from) {
        swap(const_cast<Disposable<NinePointLinearOp>&>(from));
    }

    Disposable<Array> NinePointLinearOp::apply(const Array& u)
        const {

        const std::shared_ptr<FdmLinearOpLayout> index=mesher_->layout();
        QL_REQUIRE(u.size() == index->size(),"inconsistent length of r "
                    << u.size() << " vs " << index->size());

        Array retVal(u.size());
        // #pragma omp parallel for
        for (Size i=0; i < retVal.size(); ++i) {
            retVal[i] =   a00_[i]*u[i00_[i]]
                        + a01_[i]*u[i01_[i]]
                        + a02_[i]*u[i02_[i]]
                        + a10_[i]*u[i10_[i]]
                        + a11_[i]*u[i]
                        + a12_[i]*u[i12_[i]]
                        + a20_[i]*u[i20_[i]]
                        + a21_[i]*u[i21_[i]]
                        + a22_[i]*u[i22_[i]];
        }
        return retVal;
    }

#if !defined(QL_NO_UBLAS_SUPPORT)
    Disposable<SparseMatrix> NinePointLinearOp::toMatrix() const {
        const std::shared_ptr<FdmLinearOpLayout> index = mesher_->layout();
        const Size n = index->size();

        SparseMatrix retVal(n, n, 9*n);
        for (Size i=0; i < index->size(); ++i) {
            retVal(i, i00_[i]) += a00_[i];
            retVal(i, i01_[i]) += a01_[i];
            retVal(i, i02_[i]) += a02_[i];
            retVal(i, i10_[i]) += a10_[i];
            retVal(i, i      ) += a11_[i];
            retVal(i, i12_[i]) += a12_[i];
            retVal(i, i20_[i]) += a20_[i];
            retVal(i, i21_[i]) += a21_[i];
            retVal(i, i22_[i]) += a22_[i];
        }

        return retVal;
    }
#endif


    Disposable<NinePointLinearOp>
        NinePointLinearOp::mult(const Array & u) const {

        NinePointLinearOp retVal(d0_, d1_, mesher_);
        const Size size = mesher_->layout()->size();

        // #pragma omp parallel for
        for (Size i=0; i < size; ++i) {
            const Real s = u[i];
            retVal.a11_[i]=a11_[i]*s; retVal.a00_[i]=a00_[i]*s;
            retVal.a01_[i]=a01_[i]*s; retVal.a02_[i]=a02_[i]*s;
            retVal.a10_[i]=a10_[i]*s; retVal.a20_[i]=a20_[i]*s;
            retVal.a21_[i]=a21_[i]*s; retVal.a12_[i]=a12_[i]*s;
            retVal.a22_[i]=a22_[i]*s;
        }

        return retVal;
    }

    void NinePointLinearOp::swap(NinePointLinearOp& m) {
        std::swap(d0_, m.d0_);
        std::swap(d1_, m.d1_);

        i00_.swap(m.i00_); i10_.swap(m.i10_); i20_.swap(m.i20_);
        i01_.swap(m.i01_); i21_.swap(m.i21_); i02_.swap(m.i02_);
        i12_.swap(m.i12_); i22_.swap(m.i22_);
        a00_.swap(m.a00_); a10_.swap(m.a10_); a20_.swap(m.a20_);
        a01_.swap(m.a01_); a21_.swap(m.a21_); a02_.swap(m.a02_);
        a12_.swap(m.a12_); a22_.swap(m.a22_); a11_.swap(m.a11_);

        std::swap(mesher_, m.mesher_);
    }
}
