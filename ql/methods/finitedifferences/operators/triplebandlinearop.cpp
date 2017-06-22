/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Andreas Gaida
 Copyright (C) 2008 Ralph Schreyer
 Copyright (C) 2008 Klaus Spanderen
 Copyright (C) 2014 Johannes Göttker-Schnetmann

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
#include <ql/methods/finitedifferences/tridiagonaloperator.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/operators/triplebandlinearop.hpp>

namespace QuantLib {

    TripleBandLinearOp::TripleBandLinearOp(
            Size direction,
            const std::shared_ptr<FdmMesher> &mesher)
            : direction_(direction),
              i0_(mesher->layout()->size()),
              i2_(mesher->layout()->size()),
              reverseIndex_(mesher->layout()->size()),
              lower_(mesher->layout()->size()),
              diag_(mesher->layout()->size()),
              upper_(mesher->layout()->size()),
              mesher_(mesher) {

        const std::shared_ptr<FdmLinearOpLayout> layout = mesher->layout();
        const FdmLinearOpIterator endIter = layout->end();

        std::vector<Size> newDim(layout->dim());
        std::iter_swap(newDim.begin(), newDim.begin() + direction_);
        std::vector<Size> newSpacing = FdmLinearOpLayout(newDim).spacing();
        std::iter_swap(newSpacing.begin(), newSpacing.begin() + direction_);

        for (FdmLinearOpIterator iter = layout->begin(); iter != endIter; ++iter) {
            const Size i = iter.index();

            i0_[i] = layout->neighbourhood(iter, direction, -1);
            i2_[i] = layout->neighbourhood(iter, direction, 1);

            const std::vector<Size> &coordinates = iter.coordinates();
            const Size newIndex =
                    std::inner_product(coordinates.begin(), coordinates.end(),
                                       newSpacing.begin(), Size(0));
            reverseIndex_[newIndex] = i;
        }
    }

    void TripleBandLinearOp::swap(TripleBandLinearOp &m) {
        std::swap(mesher_, m.mesher_);
        std::swap(direction_, m.direction_);

        i0_.swap(m.i0_);
        i2_.swap(m.i2_);
        reverseIndex_.swap(m.reverseIndex_);
        lower_.swap(m.lower_);
        diag_.swap(m.diag_);
        upper_.swap(m.upper_);
    }

    void TripleBandLinearOp::axpyb(const Array &a,
                                   const TripleBandLinearOp &x,
                                   const TripleBandLinearOp &y,
                                   const Array &b) {
        const Size size = mesher_->layout()->size();


        if (a.empty()) {
            if (b.empty()) {
                // #pragma omp parallel for
                for (Size i = 0; i < size; ++i) {
                    diag_[i] = y.diag_[i];
                    lower_[i] = y.lower_[i];
                    upper_[i] = y.upper_[i];
                }
            } else {
                if (b.size() > 1) {
                    // #pragma omp parallel for
                    for (Size i = 0; i < size; ++i) {
                        diag_[i] = y.diag_[i] + b[i];
                        lower_[i] = y.lower_[i];
                        upper_[i] = y.upper_[i];
                    }

                } else {
                    // #pragma omp parallel for
                    const Real s = b[0];
                    for (Size i = 0; i < size; ++i) {
                        diag_[i] = y.diag_[i] + s;
                        lower_[i] = y.lower_[i];
                        upper_[i] = y.upper_[i];
                    }
                }
            }
        } else if (b.empty()) {
            if (a.size() > 1) {
                // #pragma omp parallel for
                for (Size i = 0; i < size; ++i) {
                    const Real s = a[i];
                    diag_[i] = y.diag_[i] + s * x.diag_[i];
                    lower_[i] = y.lower_[i] + s * x.lower_[i];
                    upper_[i] = y.upper_[i] + s * x.upper_[i];
                }
            } else {
                const Real s = a[0];
                // #pragma omp parallel for
                for (Size i = 0; i < size; ++i) {
                    diag_[i] = y.diag_[i] + s * x.diag_[i];
                    lower_[i] = y.lower_[i] + s * x.lower_[i];
                    upper_[i] = y.upper_[i] + s * x.upper_[i];
                }
            }

        } else {
            if (b.size() > 1 && a.size() > 1) {
                // #pragma omp parallel for
                for (Size i = 0; i < size; ++i) {
                    const Real s = a[i];
                    diag_[i] = y.diag_[i] + s * x.diag_[i] + b[i];
                    lower_[i] = y.lower_[i] + s * x.lower_[i];
                    upper_[i] = y.upper_[i] + s * x.upper_[i];
                }

            } else if (b.size() > 1) {

                const Real s = a[0];
                // #pragma omp parallel for
                for (Size i = 0; i < size; ++i) {
                    diag_[i] = y.diag_[i] + s * x.diag_[i] + b[i];
                    lower_[i] = y.lower_[i] + s * x.lower_[i];
                    upper_[i] = y.upper_[i] + s * x.upper_[i];
                }

            } else if (a.size() > 1) {
                const Real sb = b[0];
                // #pragma omp parallel for
                for (Size i = 0; i < size; ++i) {
                    const Real sa = a[i];
                    diag_[i] = y.diag_[i] + sa * x.diag_[i] + sb;
                    lower_[i] = y.lower_[i] + sa * x.lower_[i];
                    upper_[i] = y.upper_[i] + sa * x.upper_[i];
                }

            } else {
                const Real sa = a[0];
                const Real sb = b[0];
                // #pragma omp parallel for
                for (Size i = 0; i < size; ++i) {
                    diag_[i] = y.diag_[i] + sa * x.diag_[i] + sb;
                    lower_[i] = y.lower_[i] + sa * x.lower_[i];
                    upper_[i] = y.upper_[i] + sa * x.upper_[i];
                }

            }
        }
    }

    TripleBandLinearOp
    TripleBandLinearOp::add(const TripleBandLinearOp &m) const {

        TripleBandLinearOp retVal(direction_, mesher_);
        const Size size = mesher_->layout()->size();
        // #pragma omp parallel for
        for (Size i = 0; i < size; ++i) {
            retVal.lower_[i] = lower_[i] + m.lower_[i];
            retVal.diag_[i] = diag_[i] + m.diag_[i];
            retVal.upper_[i] = upper_[i] + m.upper_[i];
        }

        return retVal;
    }


    TripleBandLinearOp TripleBandLinearOp::mult(const Array &u) const {

        TripleBandLinearOp retVal(direction_, mesher_);

        const Size size = mesher_->layout()->size();
        // #pragma omp parallel for
        for (Size i = 0; i < size; ++i) {
            const Real s = u[i];
            retVal.lower_[i] = lower_[i] * s;
            retVal.diag_[i] = diag_[i] * s;
            retVal.upper_[i] = upper_[i] * s;
        }

        return retVal;
    }

    TripleBandLinearOp TripleBandLinearOp::multR(const Array &u) const {
        const std::shared_ptr<FdmLinearOpLayout> layout = mesher_->layout();
        const Size size = layout->size();
        QL_REQUIRE(u.size() == size, "inconsistent size of rhs");
        TripleBandLinearOp retVal(direction_, mesher_);

        // #pragma omp parallel for
        for (Size i = 0; i < size; ++i) {
            const Real sm1 = i > 0 ? u[i - 1] : 1.0;
            const Real s0 = u[i];
            const Real sp1 = i < size - 1 ? u[i + 1] : 1.0;
            retVal.lower_[i] = lower_[i] * sm1;
            retVal.diag_[i] = diag_[i] * s0;
            retVal.upper_[i] = upper_[i] * sp1;
        }

        return retVal;
    }

    TripleBandLinearOp TripleBandLinearOp::add(const Array &u) const {

        TripleBandLinearOp retVal(direction_, mesher_);

        const Size size = mesher_->layout()->size();
        // #pragma omp parallel for
        for (Size i = 0; i < size; ++i) {
            retVal.lower_[i] = lower_[i];
            retVal.upper_[i] = upper_[i];
            retVal.diag_[i] = diag_[i] + u[i];
        }

        return retVal;
    }

    Array TripleBandLinearOp::apply(const Array &r) const {
        const std::shared_ptr<FdmLinearOpLayout> index = mesher_->layout();

        QL_REQUIRE(r.size() == index->size(), "inconsistent length of r");


        array_type retVal(r.size());
        // #pragma omp parallel for
        for (Size i = 0; i < index->size(); ++i) {
            retVal[i] = r[i0_[i]] * lower_[i] + r[i] * diag_[i] + r[i2_[i]] * upper_[i];
        }

        return retVal;
    }

    SparseMatrix TripleBandLinearOp::toMatrix() const {
        const std::shared_ptr<FdmLinearOpLayout> index = mesher_->layout();
        const Size n = index->size();

        SparseMatrix retVal(n, n, 3 * n);
        for (Size i = 0; i < n; ++i) {
            retVal(i, i0_[i]) += lower_[i];
            retVal(i, i) += diag_[i];
            retVal(i, i2_[i]) += upper_[i];
        }

        return retVal;
    }


    Array
    TripleBandLinearOp::solve_splitting(const Array &r, Real a, Real b) const {
        const std::shared_ptr<FdmLinearOpLayout> layout = mesher_->layout();
        QL_REQUIRE(r.size() == layout->size(), "inconsistent size of rhs");

#ifdef QL_EXTRA_SAFETY_CHECKS
        for (FdmLinearOpIterator iter = layout->begin();
             iter!=layout->end(); ++iter) {
            const std::vector<Size>& coordinates = iter.coordinates();
            QL_REQUIRE(   coordinates[direction_] != 0
                       || lower_[iter.index()] == 0,"removing non zero entry!");
            QL_REQUIRE(   coordinates[direction_] != layout->dim()[direction_]-1
                       || upper_[iter.index()] == 0,"removing non zero entry!");
        }
#endif

        Array retVal(r.size()), tmp(r.size());


        // Thomson algorithm to solve a tridiagonal system.
        // Example code taken from Tridiagonalopertor and
        // changed to fit for the triple band operator.
        Size rim1 = reverseIndex_[0];
        Real bet = 1.0 / (a * diag_[rim1] + b);
        QL_REQUIRE(bet != 0.0, "division by zero");
        retVal[reverseIndex_[0]] = r[rim1] * bet;

        for (Size j = 1; j <= layout->size() - 1; j++) {
            const Size ri = reverseIndex_[j];
            tmp[j] = a * upper_[rim1] * bet;

            bet = b + a * (diag_[ri] - tmp[j] * lower_[ri]);
            QL_ENSURE(bet != 0.0, "division by zero");
            bet = 1.0 / bet;

            retVal[ri] = (r[ri] - a * lower_[ri] * retVal[rim1]) * bet;
            rim1 = ri;
        }
        // cannot be j>=0 with Size j
        for (Size j = layout->size() - 2; j > 0; --j)
            retVal[reverseIndex_[j]] -= tmp[j + 1] * retVal[reverseIndex_[j + 1]];
        retVal[reverseIndex_[0]] -= tmp[1] * retVal[reverseIndex_[1]];

        return retVal;
    }
}
