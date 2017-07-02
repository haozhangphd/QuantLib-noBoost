/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Hao Zhang
 Copyright (C) 2007, 2008 Klaus Spanderen

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

/*! \file matrix.hpp
    \brief matrix used in linear algebra.
*/

#include <ql/math/matrix.hpp>
#ifdef QL_USE_MKL
#include <mkl.h>
#endif

#ifndef QL_USE_MKL
namespace {
    using QuantLib::Matrix;
    using QuantLib::Real;

    // inspired by Numerical Recipes 3rd edition by Press et al.
    std::tuple<Matrix, std::vector<int>, Real> luDecomposition(const Matrix &m) {
        Matrix ret(m);
        int size = ret.row_size();
        std::vector<Real> vv(size, 1);
        std::vector<int> index(size);
        Real perm = 1;
        for (int i; i < size; ++i) {
            Real largest_element = *std::max_element(ret.row_begin(i), ret.row_end(i),
                                                     [](Real x, Real y) { return std::abs(x) < std::abs(y); });
            if (largest_element == 0.0)
                return std::make_tuple(Matrix(size, size), std::vector<int>(size), 0);
            vv[i] /= largest_element;
        }

        for (int k = 0; k < size; ++k) {
            std::vector<Real> temp(size-k);
            std::transform(vv.begin()+k, vv.end(), ret.column_begin(k)+k,temp.begin(),[](Real x, Real y){return x * std::abs(y);});
            int largest_element_index = std::distance(temp.begin(), std::max_element(temp.begin(), temp.end())) + k;
            Real largest_element = temp[largest_element - k];
            if (largest_element_index != k) {
                std::swap_ranges(ret.row_begin(largest_element_index), ret.row_end(largest_element_index), ret.row_begin(k));
                perm = -perm;
                vv[largest_element_index] = vv[k];
            }

            index[k] = largest_element_index;
            temp.pop_back();
            std::for_each(ret.column_begin(k) + k + 1, ret.column_end(k), [&ret, &k](Real& x){x/=ret.diagonal()[k];});
            std::copy(ret.column_begin(k)+k + 1, ret.column_end(k), temp.begin());
            for (int i = k + 1; i < size; i++) {
                for (int j = k + 1; j < size; j++)
                    ret[i][j] -= temp[i-k-1] * ret[k][j];
            }
        }
        return std::make_tuple(ret, index, perm);
    }
}
#endif
namespace QuantLib {
    Matrix inverse(const Matrix &m) {
        const size_t size = m.rows();
        QL_REQUIRE(size == m.columns(), "matrix is not square");

#ifdef QL_USE_MKL
        std::vector<Real> lu(size * size);
        std::copy(m.begin(), m.end(), lu.begin());

        std::vector<int> ipiv(size);
        int info;

        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, size, size, lu.data(), size, ipiv.data());

        QL_REQUIRE(info == 0, "Failed to obtain LU decomposition of the matrix.");

        LAPACKE_dgetri(LAPACK_ROW_MAJOR, size, lu.data(), size, ipiv.data());

        QL_REQUIRE(info == 0, "Could not invert the matrix.");


        Matrix retVal(size, size);
        std::copy(lu.begin(), lu.end(), retVal.begin());

        return retVal;
#else
        auto [lu, index, _] = luDecomposition(m);
        Matrix ret(size, size);

        for (int k = 0; k < size; ++k) {
            std::vector<Real> ret_column(size, 0);
            ret_column[k] = 1;
            int ii=0;
            Real sum;
            for (int i = 0; i < size; ++i) {
                sum = ret_column[index[i]];
                ret_column[index[i]] = ret_column[i];
                if (ii != 0)
                    for (int j = ii - 1; j < i; ++j)
                        sum -= lu[i][j] * ret_column[j];
                else if (sum != 0.0)
                    ii = i + 1;
                ret_column[i] = sum;
            }
            for (int i = size - 1; i >= 0; --i) {
                sum = ret_column[i];
                for (int j = i + 1; j < size; j++)
                    sum -= lu[i][j] * ret_column[j];
                ret_column[i] = sum / lu[i][i];
            }
            std::copy(ret_column.begin(), ret_column.end(), ret.column_begin(k));
        }
        return ret;
#endif
    }

    Real determinant(const Matrix &m) {
        const size_t size = m.rows();
        QL_REQUIRE(size == m.columns(), "matrix is not square");
#ifdef QL_USE_MKL
        std::vector<Real> lu(size * size);

        std::copy(m.begin(), m.end(), lu.begin());

        std::vector<int> ipiv(size);
        int info;

        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, size, size, lu.data(), size, ipiv.data());

        if (info < 0)
            QL_FAIL("Failed to obtain LU decomposition of the matrix.");
        else if (info > 0)
            return 0;

        Real retVal = 1.0;

        for (Size i = 0; i < size; ++i) {
            if (ipiv[i] != i + 1)
                retVal *= -lu[i * (size + 1)];
            else
                retVal *= lu[i * (size + 1)];
        }
        return retVal;
#else
        auto [lu, _, ret] = luDecomposition(m);
        Array diag = lu.diagonal();
        return std::accumulate(diag.begin(), diag.end(), ret, std::multiplies<Real>());
#endif
    }

}
