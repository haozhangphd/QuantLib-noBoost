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

    Array NinePointLinearOp::apply(const Array& u)
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

    SparseMatrix NinePointLinearOp::toMatrix() const {
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


    NinePointLinearOp
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
