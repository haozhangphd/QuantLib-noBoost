/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2017 Hao Zhang
Copyright (C) 2012 Klaus Spanderen

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

/*! \file sparsematrix.hpp
*/

#ifndef quantlib_sparse_matrix_hpp
#define quantlib_sparse_matrix_hpp

#include <ql/math/array.hpp>

#ifdef QL_USE_MKL
#include <mkl.h>
#endif

namespace QuantLib {

    template<typename T>
    class element_proxy;

    // A sparse matrix class implemented in terms of compressed sparse row (CSR) scheme (3-vector variation)
    // one-based indexing is used INTERNALLY, for interfacing with MKL (mkl_dcsradd and mkl_dcsrmultcsr)
    // externally still zero-based indexing
    template<typename T>
    class SparseMatrixGeneral {
    public:
        //typedefs
        using size_type = int;

        //constructors
        SparseMatrixGeneral() = default;

        SparseMatrixGeneral(int row_size, int column_size, int non_zeroes = 0)
                : values_(), columns_(), rowIndex_(row_size + 1, 1),
                  column_size_(column_size), row_size_(row_size), filled_row_until_(1)
        {
            values_.reserve(non_zeroes);
            columns_.reserve(non_zeroes);
        }

        //takes the stored data directly
        SparseMatrixGeneral(int row_size, int column_size, std::vector<T> &&values, std::vector<int> &&columns,
                            std::vector<int> &&rowIndex, int filled_row_until)
                : values_(std::move(values)), columns_(std::move(columns)), rowIndex_(std::move(rowIndex)),
                  column_size_(column_size), row_size_(row_size), filled_row_until_(filled_row_until) {}

        //copying these vectors are too expensive
        SparseMatrixGeneral(int row_size, int column_size, const std::vector<T>& values, const std::vector<int>& columns,
                            const std::vector<int>& rowIndex, int filled_row_until) = delete;

        //matrix info
	int filled_size() const {return values_.size();}
        int row_size() const {return row_size_;}
        int column_size() const {return column_size_;}

        //element access
        element_proxy<T> operator() (int m, int n){
            return element_proxy<T>(m, n, this);
        }
        const T operator() (int m, int n) const{
            std::pair<int, int> index = findIndex(m, n);
            return index.first == 0? 0: values_[index.second];
        }

        //Matrix-scaler oprations
        SparseMatrixGeneral<T> &operator*=(const T &x);

        SparseMatrixGeneral<T> operator*(const T &x) const;

        SparseMatrixGeneral<T> operator-() const;

        template<typename TT>
        friend SparseMatrixGeneral<TT> operator*(const TT &x, const SparseMatrixGeneral<TT> &y);

        //Matrix-array operations
        Array operator*(const Array &x) const;

        //Matrix-matrix operations
        SparseMatrixGeneral<T> &operator+=(const SparseMatrixGeneral<T> &x);

        SparseMatrixGeneral<T> operator+(const SparseMatrixGeneral<T> &x) const;

        SparseMatrixGeneral<T> &operator-=(const SparseMatrixGeneral<T> &x);

        SparseMatrixGeneral<T> operator-(const SparseMatrixGeneral<T> &x) const;

    private:
        std::vector<T> values_;
        std::vector<int> columns_;
        std::vector<int> rowIndex_;
        int column_size_;
        int row_size_;
        // this variable stores the index beyond the last filled rowIndex_
        // rowIndex_.begin() + filled_row_until to rowIndex_.end() is zero until absolutely necessary
        // lazy evaluation. incrementing all elements in a large rowIndex_ is extremely expensive
        int filled_row_until_;


        //the first number gives 1 if element found, 0 if not found
        //the second number gives the index of the element, or the index
        //of the element would be inserted
        const std::pair<int, int> findIndex(int row, int column) const;

        void insert(int index, int row, int column, const T&x);
        friend class element_proxy<T>;
    };

    //Matrix-scaler oprations
    template<typename T>
    SparseMatrixGeneral<T> &SparseMatrixGeneral<T>::operator*=(const T &x) {
        std::for_each(values_.begin(), values_.end(), [x](T& y) {  y *= x; });
    }

    template<typename T>
    SparseMatrixGeneral<T> SparseMatrixGeneral<T>::operator*(const T &x) const {
        SparseMatrixGeneral<T> ret(*this);
        ret *= x;
        return ret;
    }

    template<typename T>
    SparseMatrixGeneral<T> operator*(const T &x, const SparseMatrixGeneral<T> &y) {
        return y * x;
    }

    template<typename T>
    SparseMatrixGeneral<T> SparseMatrixGeneral<T>::operator-() const {
        SparseMatrixGeneral<T> ret(*this);
        std::for_each(ret.values_.begin(), ret.values_.end(), [](T& x) {  x = -x; });
        return ret;
    }

    //Matrix-array operations
    template<typename T>
    Array SparseMatrixGeneral<T>::operator*(const Array &x) const {
        QL_REQUIRE(x.size() == column_size_,
                   "Array size of " << x.size() << " mismatched with the numbers of columns " << column_size_
                                    << " of the matrix.");
#ifdef QL_USE_MKL
        const char transa = 'N';

        Array ret(column_size_);
        mkl_dcsrgemv(&transa, &row_size_, values_.data(), rowIndex_.data(),
                   columns_.data(), x.data(), ret.data());
        return ret;
#else
        Array ret(row_size_);
        for (int i = 0; i < filled_row_until_ - 1; ++i)
        {
            Real sum = 0.0;
            for (int j = rowIndex_[i] - 1; j < rowIndex_[i+1] - 1; ++j)
            {
                sum += values_[j] * x[columns_[j] - 1];
            }
            ret[i] = sum;
        }
        return ret;
#endif
    }

//Matrix-matrix operations
    template<typename T>
    SparseMatrixGeneral<T> &SparseMatrixGeneral<T>::operator+=(const SparseMatrixGeneral<T> &x) {
        QL_REQUIRE(x.column_size_ == column_size_ && x.row_size_ == row_size_,
                   "Matrices of different dimensions cannot be added.");
#ifdef QL_USE_MKL
        const char transa = 'N';
        const int request = 0;
        const int sort = 0;
        const double beta = 1;
        const int nzmax = column_size_ * row_size_;

        std::vector<int> rowIndex_copy(rowIndex_);
        std::vector<int> rowIndex_copy_other(x.rowIndex_);

        int temp1 = rowIndex_[filled_row_until_ - 1];
        std::for_each(rowIndex_copy.begin() + filled_row_until_, rowIndex_copy.end(), [&temp1](int &a) { a = temp1; });
        int temp2 = x.rowIndex_[x.filled_row_until_ - 1];
        std::for_each(rowIndex_copy_other.begin() + x.filled_row_until_, rowIndex_copy_other.end(), [&temp2](int &a) { a = temp2; });

        //return values
        std::vector<T> c(nzmax);
        std::vector<int> jc(nzmax);
        std::vector<int> ic(row_size_+1);
        int info;

        mkl_dcsradd(&transa, &request, &sort, &row_size_, &column_size_, values_.data(), columns_.data(),
                    rowIndex_copy.data(), &beta, x.values_.data(), x.columns_.data(), rowIndex_copy.data(), c.data(), jc.data(),
                    ic.data(), &nzmax,  &info);
        QL_REQUIRE(info == 0, "Matrix addition cannot be performed.");
        jc.erase(std::find(jc.begin(), jc.end(), 0), jc.end());
        c.erase(c.begin() + jc.size(), c.end());
        std::vector<int>::iterator last = std::find(ic.begin(), ic.end(), jc.size() + 1);
        filled_row_until_ = std::distance(ic.begin(), last) + 1;

        values_ = std::move(c);
        columns_ = std::move(jc);
        rowIndex_ = std::move(ic);

        return *this;
#else
        // inplace addition is not performed, because doing so may require multiple insertion in the middle of
        // values_ and columns_, which is slower than constructing a new matrix, if the number of columns is large
        std::vector<T> values_ret;
        values_ret.reserve(values_.size());
        std::vector<int> columns_ret;
        columns_ret.reserve(columns_.size());
        int row_end = std::min(filled_row_until_, x.filled_row_until_) - 1;
        int original_elements;
        if (filled_row_until_ > x.filled_row_until_)
            original_elements = rowIndex_[x.filled_row_until_ - 1];
        else if (filled_row_until_ < x.filled_row_until_)
            original_elements = x.rowIndex_[filled_row_until_ - 1];
        for (int row = 0; row < row_end; ++row) {
            int column_left = rowIndex_[row] - 1;
            int column_right = x.rowIndex_[row] - 1;
            int column_left_end = rowIndex_[row + 1] - 1;
            int column_right_end = x.rowIndex_[row + 1] - 1;
            while (column_left < column_left_end && column_right < column_right_end) {
                if (columns_[column_left] < x.columns_[column_right]) {
                    values_ret.emplace_back(values_[column_left]);
                    columns_ret.emplace_back(columns_[column_left]);
                    ++column_left;
                } else if (columns_[column_left] > x.columns_[column_right]) {
                    values_ret.emplace_back(x.values_[column_right]);
                    columns_ret.emplace_back(x.columns_[column_right]);
                    ++column_right;
                } else {
                    values_ret.emplace_back(values_[column_left] + x.values_[column_right]);
                    columns_ret.emplace_back(columns_[column_left]);
                    ++column_left;
                    ++column_right;
                }

            }
            if (column_left < column_left_end) {
                int size = values_ret.size() + column_left_end - column_left;
                values_ret.reserve(size);
                columns_ret.reserve(size);

                values_ret.insert(values_ret.end(), values_.begin() + column_left, values_.begin() + column_left_end);
                columns_ret.insert(columns_ret.end(), columns_.begin() + column_left,
                                   columns_.begin() + column_left_end);
            } else if (column_right < column_right_end) {
                int size = values_ret.size() + column_right_end - column_right;
                values_ret.reserve(size);
                columns_ret.reserve(size);

                values_ret.insert(values_ret.end(), x.values_.begin() + column_right,
                                  x.values_.begin() + column_right_end);
                columns_ret.insert(columns_ret.end(), x.columns_.begin() + column_right,
                                   x.columns_.begin() + column_right_end);
            }
            rowIndex_[row + 1] = values_ret.size() + 1;
        }
        if (filled_row_until_ > x.filled_row_until_) {
            values_ret.insert(values_ret.end(), values_.begin() + rowIndex_[x.filled_row_until_] - 1, values_.end());
            columns_ret.insert(columns_ret.end(), columns_.begin() + rowIndex_[x.filled_row_until_] - 1, columns_.end());
            int new_elements = rowIndex_[x.filled_row_until_ - 1] - original_elements;
            std::for_each(rowIndex_.begin()+x.filled_row_until_, rowIndex_.begin() + filled_row_until_, [&new_elements](int & x){x += new_elements;} );
        }
        else if (filled_row_until_ < x.filled_row_until_) {
            values_ret.insert(values_ret.end(), x.values_.begin() + x.rowIndex_[filled_row_until_] - 1, x.values_.end());
            columns_ret.insert(columns_ret.end(), x.columns_.begin() + x.rowIndex_[filled_row_until_] - 1, x.columns_.end());
            int new_elements = rowIndex_[filled_row_until_ - 1] - original_elements;
            std::transform(x.rowIndex_.begin() + filled_row_until_, x.rowIndex_.begin() + x.filled_row_until_,
                           rowIndex_.begin() + filled_row_until_, [&new_elements](int x) { return x + new_elements; });
            filled_row_until_ = x.filled_row_until_;
        }
        values_ = std::move(values_ret);
        columns_ = std::move(columns_ret);
#endif
    }

    template<typename T>
    SparseMatrixGeneral<T> SparseMatrixGeneral<T>::operator+(const SparseMatrixGeneral<T> &x) const {
        // operator+= is not used here to avoid one extra copy
        QL_REQUIRE(x.column_size_ == column_size_ && x.row_size_ == row_size_,
                   "Matrices of different dimensions cannot be added.");
#ifdef QL_USE_MKL
        const char transa = 'N';
        const int request = 0;
        const int sort = 0;
        const double beta = 1;
        const int nzmax = column_size_ * row_size_;

        std::vector<int> rowIndex_copy(rowIndex_);
        std::vector<int> rowIndex_copy_other(x.rowIndex_);

        int temp1 = rowIndex_[filled_row_until_ - 1];
        std::for_each(rowIndex_copy.begin() + filled_row_until_, rowIndex_copy.end(), [&temp1](int &a) { a = temp1; });
        int temp2 = x.rowIndex_[x.filled_row_until_ - 1];
        std::for_each(rowIndex_copy_other.begin() + x.filled_row_until_, rowIndex_copy_other.end(), [&temp2](int &a) { a = temp2; });

        //return values
        std::vector<T> c(nzmax);
        std::vector<int> jc(nzmax);
        std::vector<int> ic(row_size_+1);
        int info;

        //since sort is disabled, the const variables are not modified
        mkl_dcsradd(&transa, &request, &sort, &row_size_, &column_size_, const_cast<T *>(values_.data()),
                    const_cast<int *>(columns_.data()), rowIndex_copy.data(), &beta, const_cast<T *>(x.values_.data()),
                    const_cast<int *>(x.columns_.data()), rowIndex_copy_other.data(), c.data(), jc.data(), ic.data(), &nzmax, &info);
        QL_REQUIRE(info == 0, "Matrix addition cannot be performed.");

        jc.erase(std::find(jc.begin(), jc.end(), 0), jc.end());
        c.erase(c.begin() + jc.size(), c.end());
        std::vector<int>::iterator last = std::find(ic.begin(), ic.end(), jc.size() + 1);
        int filled_row_until = std::distance(ic.begin(), last) + 1;

        return SparseMatrixGeneral<T>(row_size_,column_size_, std::move(c), std::move(jc), std::move(ic), filled_row_until);
#else
        // inplace addition is not performed, because doing so may require multiple insertion in the middle of
        // values_ and columns_, which is slower than constructing a new matrix, if the number of columns is large
        std::vector<T> values_ret;
        values_ret.reserve(values_.size());
        std::vector<int> columns_ret;
        columns_ret.reserve(columns_.size());
        std::vector<int> rowIndex_ret;
        rowIndex_ret.reserve(row_size_ + 1);
        rowIndex_ret.emplace_back(1);
        int row_end = std::min(filled_row_until_, x.filled_row_until_) - 1;
        int original_elements;
        if (filled_row_until_ > x.filled_row_until_)
            original_elements = rowIndex_[x.filled_row_until_ - 1];
        else if (filled_row_until_ < x.filled_row_until_)
            original_elements = x.rowIndex_[filled_row_until_ - 1];
        for (int row = 0; row < row_end; ++row) {
            int column_left = rowIndex_[row] - 1;
            int column_right = x.rowIndex_[row] - 1;
            int column_left_end = rowIndex_[row + 1] - 1;
            int column_right_end = x.rowIndex_[row + 1] - 1;
            while (column_left < column_left_end && column_right < column_right_end) {
                if (columns_[column_left] < x.columns_[column_right]) {
                    values_ret.emplace_back(values_[column_left]);
                    columns_ret.emplace_back(columns_[column_left]);
                    ++column_left;
                } else if (columns_[column_left] > x.columns_[column_right]) {
                    values_ret.emplace_back(x.values_[column_right]);
                    columns_ret.emplace_back(x.columns_[column_right]);
                    ++column_right;
                } else {
                    values_ret.emplace_back(values_[column_left] + x.values_[column_right]);
                    columns_ret.emplace_back(columns_[column_left]);
                    ++column_left;
                    ++column_right;
                }

            }
            if (column_left < column_left_end) {
                int size = values_ret.size() + column_left_end - column_left;
                values_ret.reserve(size);
                columns_ret.reserve(size);

                values_ret.insert(values_ret.end(), values_.begin() + column_left, values_.begin() + column_left_end);
                columns_ret.insert(columns_ret.end(), columns_.begin() + column_left,
                                   columns_.begin() + column_left_end);
            } else if (column_right < column_right_end) {
                int size = values_ret.size() + column_right_end - column_right;
                values_ret.reserve(size);
                columns_ret.reserve(size);

                values_ret.insert(values_ret.end(), x.values_.begin() + column_right,
                                  x.values_.begin() + column_right_end);
                columns_ret.insert(columns_ret.end(), x.columns_.begin() + column_right,
                                   x.columns_.begin() + column_right_end);
            }
            rowIndex_ret.emplace_back(values_ret.size() + 1);
        }
        int filled_row_until_ret;
        if (filled_row_until_ > x.filled_row_until_) {
            values_ret.insert(values_ret.end(), values_.begin() + rowIndex_[x.filled_row_until_] - 1, values_.end());
            columns_ret.insert(columns_ret.end(), columns_.begin() + rowIndex_[x.filled_row_until_] - 1, columns_.end());
            int new_elements = rowIndex_ret[x.filled_row_until_ - 1] - original_elements;
            std::transform(rowIndex_.begin()+x.filled_row_until_, rowIndex_.begin() + filled_row_until_,
                           std::back_inserter(rowIndex_ret), [&new_elements](int x) { return x + new_elements; });
            filled_row_until_ret = filled_row_until_;
        }
        else if (filled_row_until_ < x.filled_row_until_) {
            values_ret.insert(values_ret.end(), x.values_.begin() + x.rowIndex_[filled_row_until_] - 1, x.values_.end());
            columns_ret.insert(columns_ret.end(), x.columns_.begin() + x.rowIndex_[filled_row_until_] - 1, x.columns_.end());
            int new_elements = rowIndex_[filled_row_until_ - 1] - original_elements;
            std::transform(x.rowIndex_.begin() + filled_row_until_, x.rowIndex_.begin() + x.filled_row_until_,
                           std::back_inserter(rowIndex_ret), [&new_elements](int x) { return x + new_elements; });
            filled_row_until_ret = x.filled_row_until_;
        }
        else
            filled_row_until_ret = filled_row_until_;
	return SparseMatrixGeneral<T>(row_size_, column_size_, std::move(values_ret), std::move(columns_ret),
                                      std::move(rowIndex_ret), filled_row_until_ret);
#endif
    }

    template<typename T>
    SparseMatrixGeneral<T> &SparseMatrixGeneral<T>::operator-=(const SparseMatrixGeneral<T> &x) {
        QL_REQUIRE(x.column_size_ == column_size_ && x.row_size_ == row_size_,
                   "Matrices of different dimensions cannot be subtracted.");
#ifdef QL_USE_MKL
        const char transa = 'N';
        const int request = 0;
        const int sort = 0;
        const double beta = 1;
        const int nzmax = column_size_ * row_size_;

        std::vector<int> rowIndex_copy(rowIndex_);
        std::vector<int> rowIndex_copy_other(x.rowIndex_);
        std::vector<T> values_copy_other(x.values_);

        int temp1 = rowIndex_[filled_row_until_ - 1];
        std::for_each(rowIndex_copy.begin() + filled_row_until_, rowIndex_copy.end(), [&temp1](int &a) { a = temp1; });
        int temp2 = x.rowIndex_[x.filled_row_until_ - 1];
        std::for_each(rowIndex_copy_other.begin() + x.filled_row_until_, rowIndex_copy_other.end(), [&temp2](int &a) { a = temp2; });
        std::for_each(values_copy_other.begin(), values_copy_other.end(), [](T& x) {  x = -x; });

        //return values
        std::vector<T> c(nzmax);
        std::vector<int> jc(nzmax);
        std::vector<int> ic(row_size_+1);
        int info;

        mkl_dcsradd(&transa, &request, &sort, &row_size_, &column_size_, values_.data(), columns_.data(),
                    rowIndex_copy.data(), &beta, values_copy_other.data(), x.columns_.data(), rowIndex_copy_other.data(), c.data(), jc.data(),
                    ic.data(), &nzmax,  &info);
        QL_REQUIRE(info == 0, "Matrix addition cannot be performed.");
        jc.erase(std::find(jc.begin(), jc.end(), 0), jc.end());
        c.erase(c.begin() + jc.size(), c.end());
        std::vector<int>::iterator last = std::find(ic.begin(), ic.end(), jc.size() + 1);
        filled_row_until_ = std::distance(ic.begin(), last) + 1;

        values_ = std::move(c);
        columns_ = std::move(jc);
        rowIndex_ = std::move(ic);

        return *this;
#else
        // inplace subtraction is not performed, because doing so may require multiple insertion in the middle of
        // values_ and columns_, which is slower than constructing a new matrix, if the number of columns is large
        std::vector<T> values_ret;
        values_ret.reserve(values_.size());
        std::vector<int> columns_ret;
        columns_ret.reserve(columns_.size());
        int row_end = std::min(filled_row_until_, x.filled_row_until_) - 1;
        int original_elements;
        if (filled_row_until_ > x.filled_row_until_)
            original_elements = rowIndex_[x.filled_row_until_ - 1];
        else if (filled_row_until_ < x.filled_row_until_)
            original_elements = x.rowIndex_[filled_row_until_ - 1];
        for (int row = 0; row < row_end; ++row) {
            int column_left = rowIndex_[row] - 1;
            int column_right = x.rowIndex_[row] - 1;
            int column_left_end = rowIndex_[row + 1] - 1;
            int column_right_end = x.rowIndex_[row + 1] - 1;
            while (column_left < column_left_end && column_right < column_right_end) {
                if (columns_[column_left] < x.columns_[column_right]) {
                    values_ret.emplace_back(values_[column_left]);
                    columns_ret.emplace_back(columns_[column_left]);
                    ++column_left;
                } else if (columns_[column_left] > x.columns_[column_right]) {
                    values_ret.emplace_back(-x.values_[column_right]);
                    columns_ret.emplace_back(x.columns_[column_right]);
                    ++column_right;
                } else {
                    values_ret.emplace_back(values_[column_left] - x.values_[column_right]);
                    columns_ret.emplace_back(columns_[column_left]);
                    ++column_left;
                    ++column_right;
                }

            }
            if (column_left < column_left_end) {
                int size = values_ret.size() + column_left_end - column_left;
                values_ret.reserve(size);
                columns_ret.reserve(size);

                values_ret.insert(values_ret.end(), values_.begin() + column_left, values_.begin() + column_left_end);
                columns_ret.insert(columns_ret.end(), columns_.begin() + column_left,
                                   columns_.begin() + column_left_end);
            } else if (column_right < column_right_end) {
                int size = values_ret.size() + column_right_end - column_right;
                values_ret.reserve(size);
                columns_ret.reserve(size);

                std::transform(x.values_.begin() + column_right,
                               x.values_.begin() + column_right_end, std::back_inserter(values_ret),
                               [](T x) { return -x; });
                columns_ret.insert(columns_ret.end(), x.columns_.begin() + column_right,
                                   x.columns_.begin() + column_right_end);
            }
            rowIndex_[row + 1] = values_ret.size() + 1;
        }
        if (filled_row_until_ > x.filled_row_until_) {
            values_ret.insert(values_ret.end(), values_.begin() + rowIndex_[x.filled_row_until_] - 1, values_.end());
            columns_ret.insert(columns_ret.end(), columns_.begin() + rowIndex_[x.filled_row_until_] - 1, columns_.end());
            int new_elements = rowIndex_[x.filled_row_until_ - 1] - original_elements;
            std::for_each(rowIndex_.begin()+x.filled_row_until_, rowIndex_.begin() + filled_row_until_, [&new_elements](int & x){x += new_elements;} );
        }
        else if (filled_row_until_ < x.filled_row_until_) {
            std::transform(x.values_.begin() + x.rowIndex_[filled_row_until_] - 1, x.values_.end(),
                           std::back_inserter(values_ret), [](Real &x) { return -x; });
            columns_ret.insert(columns_ret.end(), x.columns_.begin() + x.rowIndex_[filled_row_until_] - 1, x.columns_.end());
            int new_elements = rowIndex_[filled_row_until_ - 1] - original_elements;
            std::transform(x.rowIndex_.begin() + filled_row_until_, x.rowIndex_.begin() + x.filled_row_until_,
                           rowIndex_.begin() + filled_row_until_, [&new_elements](int x) { return x + new_elements; });
            filled_row_until_ = x.filled_row_until_;
        }
        values_ = std::move(values_ret);
        columns_ = std::move(columns_ret);
#endif
    }

    template<typename T>
    SparseMatrixGeneral<T> SparseMatrixGeneral<T>::operator-(const SparseMatrixGeneral<T> &x) const {
        QL_REQUIRE(x.column_size_ == column_size_ && x.row_size_ == row_size_,
                   "Matrices of different dimensions cannot be subtracted.");
#ifdef QL_USE_MKL
        const char transa = 'N';
        const int request = 0;
        const int sort = 0;
        const double beta = 1;
        const int nzmax = column_size_ * row_size_;

        std::vector<int> rowIndex_copy(rowIndex_);
        std::vector<int> rowIndex_copy_other(x.rowIndex_);
        std::vector<T> values_copy_other(x.values_);

        int temp1 = rowIndex_[filled_row_until_ - 1];
        std::for_each(rowIndex_copy.begin() + filled_row_until_, rowIndex_copy.end(), [&temp1](int &a) { a = temp1; });
        int temp2 = x.rowIndex_[x.filled_row_until_ - 1];
        std::for_each(rowIndex_copy_other.begin() + x.filled_row_until_, rowIndex_copy_other.end(), [&temp2](int &a) { a = temp2; });
        std::for_each(values_copy_other.begin(), values_copy_other.end(), [](T& x) {  x = -x; });

        //return values
        std::vector<T> c(nzmax);
        std::vector<int> jc(nzmax);
        std::vector<int> ic(row_size_+1);
        int info;

        //since sort is disabled, the const variables are not modified
        mkl_dcsradd(&transa, &request, &sort, &row_size_, &column_size_, const_cast<T *>(values_.data()),
                    const_cast<int *>(columns_.data()), rowIndex_copy.data(), &beta, values_copy_other.data()),
                    const_cast<int *>(x.columns_.data(), rowIndex_copy_other.data(), c.data(), jc.data(), ic.data(), &nzmax, &info);
        QL_REQUIRE(info == 0, "Matrix addition cannot be performed.");

        jc.erase(std::find(jc.begin(), jc.end(), 0), jc.end());
        c.erase(c.begin() + jc.size(), c.end());
        std::vector<int>::iterator last = std::find(ic.begin(), ic.end(), jc.size() + 1);
        int filled_row_until = std::distance(ic.begin(), last) + 1;

        return SparseMatrixGeneral<T>(row_size_,column_size_, std::move(c), std::move(jc), std::move(ic), filled_row_until);
#else
        // inplace subtraction is not performed, because doing so may require multiple insertion in the middle of
        // values_ and columns_, which is slower than constructing a new matrix, if the number of columns is large
        std::vector<T> values_ret;
        values_ret.reserve(values_.size());
        std::vector<int> columns_ret;
        columns_ret.reserve(columns_.size());
        std::vector<int> rowIndex_ret;
        rowIndex_ret.reserve(row_size_ + 1);
        rowIndex_ret.emplace_back(1);
        int row_end = std::min(filled_row_until_, x.filled_row_until_) - 1;
        int original_elements;
        if (filled_row_until_ > x.filled_row_until_)
            original_elements = rowIndex_[x.filled_row_until_ - 1];
        else if (filled_row_until_ < x.filled_row_until_)
            original_elements = x.rowIndex_[filled_row_until_ - 1];
        for (int row = 0; row < row_end; ++row) {
            int column_left = rowIndex_[row] - 1;
            int column_right = x.rowIndex_[row] - 1;
            int column_left_end = rowIndex_[row + 1] - 1;
            int column_right_end = x.rowIndex_[row + 1] - 1;
            while (column_left < column_left_end && column_right < column_right_end) {
                if (columns_[column_left] < x.columns_[column_right]) {
                    values_ret.emplace_back(values_[column_left]);
                    columns_ret.emplace_back(columns_[column_left]);
                    ++column_left;
                } else if (columns_[column_left] > x.columns_[column_right]) {
                    values_ret.emplace_back(-x.values_[column_right]);
                    columns_ret.emplace_back(x.columns_[column_right]);
                    ++column_right;
                } else {
                    values_ret.emplace_back(values_[column_left] - x.values_[column_right]);
                    columns_ret.emplace_back(columns_[column_left]);
                    ++column_left;
                    ++column_right;
                }

            }
            if (column_left < column_left_end) {
                int size = values_ret.size() + column_left_end - column_left;
                values_ret.reserve(size);
                columns_ret.reserve(size);

                values_ret.insert(values_ret.end(), values_.begin() + column_left, values_.begin() + column_left_end);
                columns_ret.insert(columns_ret.end(), columns_.begin() + column_left,
                                   columns_.begin() + column_left_end);
            } else if (column_right < column_right_end) { 
                int size = values_ret.size() + column_right_end - column_right;
                values_ret.reserve(size);
                columns_ret.reserve(size);

                std::transform(x.values_.begin() + column_right,
                               x.values_.begin() + column_right_end, std::back_inserter(values_ret),
                               [](T x) { return -x; });
                columns_ret.insert(columns_ret.end(), x.columns_.begin() + column_right,
                                   x.columns_.begin() + column_right_end);
            }
            rowIndex_ret.emplace_back(values_ret.size() + 1);
        }
        int filled_row_until_ret;
        if (filled_row_until_ > x.filled_row_until_) {
            values_ret.insert(values_ret.end(), values_.begin() + rowIndex_[x.filled_row_until_] - 1, values_.end());
            columns_ret.insert(columns_ret.end(), columns_.begin() + rowIndex_[x.filled_row_until_] - 1, columns_.end());
            int new_elements = rowIndex_ret[x.filled_row_until_ - 1] - original_elements;
            std::transform(rowIndex_.begin()+x.filled_row_until_, rowIndex_.begin() + filled_row_until_,
                           std::back_inserter(rowIndex_ret), [&new_elements](int x) { return x + new_elements; });
            filled_row_until_ret = filled_row_until_;
        }
        else if (filled_row_until_ < x.filled_row_until_) {
            std::transform(x.values_.begin() + x.rowIndex_[filled_row_until_] - 1, x.values_.end(),
                           std::back_inserter(values_ret), [](Real &x) { return -x; });
            columns_ret.insert(columns_ret.end(), x.columns_.begin() + x.rowIndex_[filled_row_until_] - 1, x.columns_.end());
            int new_elements = rowIndex_ret[filled_row_until_ - 1] - original_elements;
            std::transform(x.rowIndex_.begin() + filled_row_until_, x.rowIndex_.begin() + x.filled_row_until_,
                           std::back_inserter(rowIndex_ret), [&new_elements](int x) { return x + new_elements; });
            filled_row_until_ret = x.filled_row_until_;
        }
        else
            filled_row_until_ret = filled_row_until_;
        return SparseMatrixGeneral<T>(row_size_, column_size_, std::move(values_ret), std::move(columns_ret),
                                      std::move(rowIndex_ret), filled_row_until_ret);
#endif
    }

    template<typename T>
    const std::pair<int, int> SparseMatrixGeneral<T>::findIndex(int row, int column) const {
        // in a non-filled row
        if (row + 1 >= filled_row_until_)
            return std::make_pair(0, rowIndex_[filled_row_until_ - 1] - 1);
        auto index = std::equal_range(columns_.begin() + rowIndex_[row] - 1, columns_.begin() + rowIndex_[row+1] - 1, column + 1);
        int i = std::distance(columns_.begin(), index.first);
        if (index.first == index.second)
            return std::make_pair(0, i);
        else
            return std::make_pair(1, i);
    }

    template<typename T>
    void SparseMatrixGeneral<T>::insert(int index, int row, int column, const T &x) {
        columns_.insert(columns_.begin() + index, column + 1);
        values_.insert(values_.begin() + index, x);

        //if brand new row
        if (row + 1 >= filled_row_until_) {
            int temp = rowIndex_[filled_row_until_ - 1];
            std::for_each(rowIndex_.begin() + filled_row_until_, rowIndex_.begin() + row + 1, [&temp](int &a) { a = temp; });
            filled_row_until_ = row + 2;
            rowIndex_[row + 1] = columns_.size() + 1;
        }
        //insert into existing row
        else
            std::for_each(rowIndex_.begin() + row + 1, rowIndex_.begin() + filled_row_until_, [](int &a) { ++a; });
    }

    // a proxy class used to get and set values in the SparseMatrix
    // makes m(a, b) = c possible
    template<typename T>
    class element_proxy {
    public:
        element_proxy(int m, int n, SparseMatrixGeneral<T> *parent) : m_(m), n_(n), index_(parent->findIndex(m, n)), parent_(parent) {}

        const element_proxy &operator=(const T &x) {
            //if value not already present
            if (index_.first == 0) {
                parent_->insert(index_.second, m_, n_, x);
            }
                // if value present
            else {
                parent_->values_[index_.second] = x;
            }
            return *this;
        }

        const element_proxy &operator+=(const T &x) {
            //if value not already present
            if (index_.first == 0) {
                parent_->insert(index_.second, m_, n_, x);
            }
                // if value present
            else {
                parent_->values_[index_.second] += x;
            }
            return *this;
        }

        const element_proxy &operator-=(const T &x) {
            //if value not already present
            if (index_.first == 0) {
                parent_->insert(index_.second, m_, n_, -x);
            }
                // if value present
            else {
                parent_->values_[index_.second] -= x;
            }
            return *this;
        }

        operator const T() const {
            return index_.first == 0? 0: parent_->values_[index_.second];
        }

    private:
        const int m_;
        const int n_;
        const std::pair<int, int> index_;
        SparseMatrixGeneral <T> *parent_;
    };

    using SparseMatrix = SparseMatrixGeneral<Real>;
    using IntegerSparseMatrix = SparseMatrixGeneral<int>;

    // This is not used ANYWHERE in the QuantLib library
    class SparseMatrixReference {
    public:
        explicit SparseMatrixReference(SparseMatrix &m) : data_(m) {}

        operator const SparseMatrix &() const { return data_; }

        const Real operator()(int i, int j) const {
            return Real(data_(i, j));
        }

        element_proxy<Real> operator()(int i, int j) {
            return data_(i, j);
        }

        SparseMatrix operator+(const SparseMatrix &right) {
            return data_ + right;
        }

        friend SparseMatrix operator+(const SparseMatrix &left, const SparseMatrixReference &right) {
            return left + SparseMatrix(right);
        }

        Array operator*(const Array &A) {
            return data_ * A;
        }
    private:
        SparseMatrix &data_;

    };

    template<typename T>
    SparseMatrixGeneral<T> identity_matrix(int s) {
        std::vector<T> values(s, 1);
        std::vector<int> columns(s);
        std::iota(columns.begin(), columns.end(), 1);
        std::vector<int> rowIndex(columns);
        rowIndex.push_back(rowIndex.size());
        return SparseMatrixGeneral<T>(s, s, std::move(values), std::move(columns), std::move(rowIndex), s + 1);
    }

}

#endif
