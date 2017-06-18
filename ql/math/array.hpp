/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2009 StatPro Italia srl
 Copyright (C) 2004 Ferdinando Ametrano

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

/*! \file array.hpp
	\brief 1-D array used in linear algebra.
*/

#ifndef quantlib_array_hpp
#define quantlib_array_hpp

#include <ql/types.hpp>
#include <ql/errors.hpp>
#include <ql/utilities/disposable.hpp>
#include <ql/utilities/null.hpp>
#include <vector>
#include <type_traits>
#include <functional>
#include <numeric>
#include <iomanip>

namespace QuantLib {

    //! 1-D array used in linear algebra.
    /*! This class implements the concept of vector as used in linear
        algebra.
        As such, it is <b>not</b> meant to be used as a container -
        <tt>std::vector</tt> should be used instead.

        \test construction of arrays is checked in a number of cases
    */
    class Array {
    public:
        //! \name Constructors, destructor, and assignment
        //@{
        //! creates the array with the given dimension
        explicit Array(Size size = 0);

        //! creates the array and fills it with <tt>value</tt>
        Array(Size size, Real value);

        /*! \brief creates the array and fills it according to
            \f$ a_{0} = value, a_{i}=a_{i-1}+increment \f$
        */
        Array(Size size, Real value, Real increment);

        Array(const std::vector<Real> &);

        Array(const Array &);

        Array(const Disposable<Array> &);

        //! creates the array from an iterable sequence
        template<class ForwardIterator>
        Array(ForwardIterator begin, ForwardIterator end);

        Array &operator=(const Array &);

        Array &operator=(const Disposable<Array> &);

        bool operator==(const Array &) const;

        bool operator!=(const Array &) const;
        //@}
        /*! \name Vector algebra

            <tt>v += x</tt> and similar operation involving a scalar value
            are shortcuts for \f$ \forall i : v_i = v_i + x \f$

            <tt>v *= w</tt> and similar operation involving two vectors are
            shortcuts for \f$ \forall i : v_i = v_i \times w_i \f$

            \pre all arrays involved in an algebraic expression must have
            the same size.
        */
        //@{
        const Array &operator+=(const Array &);

        const Array &operator+=(Real);

        const Array &operator-=(const Array &);

        const Array &operator-=(Real);

        const Array &operator*=(const Array &);

        const Array &operator*=(Real);

        const Array &operator/=(const Array &);

        const Array &operator/=(Real);

        //@}
        //! \name Element access
        //@{
        //! read-only
        Real operator[](Size) const;

        Real at(Size) const;

        Real front() const;

        Real back() const;

        //! read-write
        Real &operator[](Size);

        Real &at(Size);

        Real &front();

        Real &back();

        //@}
        //! \name Inspectors
        //@{
        //! dimension of the array
        Size size() const;

        //! whether the array is empty
        bool empty() const;

        //@}
        typedef Size size_type;
        typedef Real value_type;
        typedef std::vector<Real>::iterator iterator;
        typedef std::vector<Real>::const_iterator const_iterator;
        typedef std::vector<Real>::reverse_iterator reverse_iterator;
        typedef std::vector<Real>::const_reverse_iterator const_reverse_iterator;

        //! \name Iterator access
        //@{
        const_iterator begin() const;

        iterator begin();

        const_iterator end() const;

        iterator end();

        const_reverse_iterator rbegin() const;

        reverse_iterator rbegin();

        const_reverse_iterator rend() const;

        reverse_iterator rend();

        //@}
        //! \name Utilities
        //@{
        void swap(Array &) noexcept;  // never throws
        //@}

    private:
        std::vector<Real> data_;
        Size n_;
    };

    //! specialization of null template for this class
    template<>
    class Null<Array> {
    public:
        Null() {}

        operator Array() const { return Array(); }
    };


    /*! \relates Array */
    Real DotProduct(const Array &, const Array &);

    // unary operators
    /*! \relates Array */
    const Disposable<Array> operator+(const Array &v);

    /*! \relates Array */
    const Disposable<Array> operator-(const Array &v);

    // binary operators
    /*! \relates Array */
    const Disposable<Array> operator+(const Array &, const Array &);

    /*! \relates Array */
    const Disposable<Array> operator+(const Array &, Real);

    /*! \relates Array */
    const Disposable<Array> operator+(Real, const Array &);

    /*! \relates Array */
    const Disposable<Array> operator-(const Array &, const Array &);

    /*! \relates Array */
    const Disposable<Array> operator-(const Array &, Real);

    /*! \relates Array */
    const Disposable<Array> operator-(Real, const Array &);

    /*! \relates Array */
    const Disposable<Array> operator*(const Array &, const Array &);

    /*! \relates Array */
    const Disposable<Array> operator*(const Array &, Real);

    /*! \relates Array */
    const Disposable<Array> operator*(Real, const Array &);

    /*! \relates Array */
    const Disposable<Array> operator/(const Array &, const Array &);

    /*! \relates Array */
    const Disposable<Array> operator/(const Array &, Real);

    /*! \relates Array */
    const Disposable<Array> operator/(Real, const Array &);

    // math functions
    /*! \relates Array */
    const Disposable<Array> Abs(const Array &);

    /*! \relates Array */
    const Disposable<Array> Sqrt(const Array &);

    /*! \relates Array */
    const Disposable<Array> Log(const Array &);

    /*! \relates Array */
    const Disposable<Array> Exp(const Array &);

    /*! \relates Array */
    const Disposable<Array> Pow(const Array &, Real);

    // utilities
    /*! \relates Array */
    void swap(Array &, Array &) noexcept;

    // format
    /*! \relates Array */
    std::ostream &operator<<(std::ostream &, const Array &);


    // inline definitions

    inline Array::Array(Size size)
            : data_(size), n_(size) {}

    inline Array::Array(Size size, Real value)
            : data_(size, value), n_(size) {
    }

    inline Array::Array(Size size, Real value, Real increment)
            : data_(size), n_(size) {
        for (iterator i = begin(); i != end(); ++i, value += increment)
            *i = value;
    }

    inline Array::Array(const std::vector<Real> &vec)
            : data_(vec), n_(vec.size()) {}

    inline Array::Array(const Array &from)
            : data_(from.n_), n_(from.n_) {
#if defined(QL_PATCH_MSVC) && defined(QL_DEBUG)
        if (n_)
#endif
        std::copy(from.begin(), from.end(), begin());
    }

    inline Array::Array(const Disposable<Array> &from)
            : data_(std::vector<Real>()), n_(0) {
        swap(const_cast<Disposable<Array> &>(from));
    }

    namespace detail {

        template<class I>
        inline void _fill_array_(Array &a,
                                 std::vector<Real> &data_,
                                 Size &n_,
                                 I begin, I end,
                                 const std::true_type &) {
            // we got redirected here from a call like Array(3, 4)
            // because it matched the constructor below exactly with
            // ForwardIterator = int.  What we wanted was fill an
            // Array with a given value, which we do here.
            Size n = begin;
            Real value = end;
            n_ = n;
            if (n_)
                data_.assign(n, value);
            else
                data_.clear();
        }

        template<class I>
        inline void _fill_array_(Array &a,
                                 std::vector<Real> &data_,
                                 Size &n_,
                                 I begin, I end,
                                 const std::false_type &) {
            // true iterators
            Size n = std::distance(begin, end);
            n_ = n;
            if (n_) {
                data_ = std::vector<Real>(n_);
                std::copy(begin, end, data_.begin());
            }
            else data_.clear();
        }

    }

    template<class ForwardIterator>
    inline Array::Array(ForwardIterator begin, ForwardIterator end) {
        // Unfortunately, calls such as Array(3, 4) match this constructor.
        // We have to detect integral types and dispatch.
        detail::_fill_array_(*this, data_, n_, begin, end,
                             std::is_integral<ForwardIterator>());
    }

    inline Array &Array::operator=(const Array &from) {
        // strong guarantee
        Array temp(from);
        swap(temp);
        return *this;
    }

    inline bool Array::operator==(const Array &to) const {
        return (n_ == to.n_) && std::equal(begin(), end(), to.begin());
    }

    inline bool Array::operator!=(const Array &to) const {
        return !(this->operator==(to));
    }

    inline Array &Array::operator=(const Disposable<Array> &from) {
        swap(const_cast<Disposable<Array> &>(from));
        return *this;
    }

    inline const Array &Array::operator+=(const Array &v) {
        QL_REQUIRE(n_ == v.n_,
                   "arrays with different sizes (" << n_ << ", "
                                                   << v.n_ << ") cannot be added");
        std::transform(begin(), end(), v.begin(), begin(),
                       [](Real i, Real j) { return i + j; });
        return *this;
    }


    inline const Array &Array::operator+=(Real x) {
        std::transform(begin(), end(), begin(),
                       [x](Real i) { return i + x; });
        return *this;
    }

    inline const Array &Array::operator-=(const Array &v) {
        QL_REQUIRE(n_ == v.n_,
                   "arrays with different sizes (" << n_ << ", "
                                                   << v.n_ << ") cannot be subtracted");
        std::transform(begin(), end(), v.begin(), begin(),
                       [](Real i, Real j) { return i + j; });
        return *this;
    }

    inline const Array &Array::operator-=(Real x) {
        std::transform(begin(), end(), begin(),
                       [x](Real i) { return i - x; });
        return *this;
    }

    inline const Array &Array::operator*=(const Array &v) {
        QL_REQUIRE(n_ == v.n_,
                   "arrays with different sizes (" << n_ << ", "
                                                   << v.n_ << ") cannot be multiplied");
        std::transform(begin(), end(), v.begin(), begin(),
                       [](Real i, Real j) { return i * j; });
        return *this;
    }

    inline const Array &Array::operator*=(Real x) {
        std::transform(begin(), end(), begin(),
                       [x](Real i) { return i * x; });
        return *this;
    }

    inline const Array &Array::operator/=(const Array &v) {
        QL_REQUIRE(n_ == v.n_,
                   "arrays with different sizes (" << n_ << ", "
                                                   << v.n_ << ") cannot be divided");
        std::transform(begin(), end(), v.begin(), begin(),
                       [](Real i, Real j) { return i / j; });
        return *this;
    }

    inline const Array &Array::operator/=(Real x) {
        std::transform(begin(), end(), begin(),
                       [x](Real i) { return i / x; });
        return *this;
    }

    inline Real Array::operator[](Size i) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(i < n_,
            "index (" << i << ") must be less than " << n_ <<
            ": array access out of range");
#endif
        return data_[i];
    }

    inline Real Array::at(Size i) const {
        QL_REQUIRE(i < n_,
                   "index (" << i << ") must be less than " << n_ <<
                             ": array access out of range");
        return data_[i];
    }

    inline Real Array::front() const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(n_ > 0, "null Array: array access out of range");
#endif
        return data_.front();
    }

    inline Real Array::back() const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(n_ > 0, "null Array: array access out of range");
#endif
        return data_.back();
    }

    inline Real &Array::operator[](Size i) {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(i < n_,
            "index (" << i << ") must be less than " << n_ <<
            ": array access out of range");
#endif
        return data_[i];
    }

    inline Real &Array::at(Size i) {
        QL_REQUIRE(i < n_,
                   "index (" << i << ") must be less than " << n_ <<
                             ": array access out of range");
        return data_[i];
    }

    inline Real &Array::front() {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(n_ > 0, "null Array: array access out of range");
#endif
        return data_.front();
    }

    inline Real &Array::back() {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(n_ > 0, "null Array: array access out of range");
#endif
        return data_.back();
    }

    inline Size Array::size() const {
        return n_;
    }

    inline bool Array::empty() const {
        return n_ == 0;
    }

    inline Array::const_iterator Array::begin() const {
        return data_.begin();
    }

    inline Array::iterator Array::begin() {
        return data_.begin();
    }

    inline Array::const_iterator Array::end() const {
        return data_.end();
    }

    inline Array::iterator Array::end() {
        return data_.end();
    }

    inline Array::const_reverse_iterator Array::rbegin() const {
        return data_.rbegin();
    }

    inline Array::reverse_iterator Array::rbegin() {
        return data_.rbegin();
    }

    inline Array::const_reverse_iterator Array::rend() const {
        return data_.rend();
    }

    inline Array::reverse_iterator Array::rend() {
        return data_.rend();
    }

    inline void Array::swap(Array &from) noexcept {
        using std::swap;
        data_.swap(from.data_);
        swap(n_, from.n_);
    }


    // dot product

    inline Real DotProduct(const Array &v1, const Array &v2) {
        QL_REQUIRE(v1.size() == v2.size(),
                   "arrays with different sizes (" << v1.size() << ", "
                                                   << v2.size() << ") cannot be multiplied");
        return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
    }

    // overloaded operators

    // unary

    inline const Disposable<Array> operator+(const Array &v) {
        Array result = v;
        return result;
    }

    inline const Disposable<Array> operator-(const Array &v) {
        Array result(v.size());
        std::transform(v.begin(), v.end(), result.begin(),
                       [](Real i) { return -i; });
        return result;
    }


    // binary operators

    inline const Disposable<Array> operator+(const Array &v1,
                                             const Array &v2) {
        QL_REQUIRE(v1.size() == v2.size(),
                   "arrays with different sizes (" << v1.size() << ", "
                                                   << v2.size() << ") cannot be added");
        Array result(v1.size());
        std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(),
                       [](Real i, Real j) { return i + j; });
        return result;
    }

    inline const Disposable<Array> operator+(const Array &v1, Real a) {
        Array result(v1.size());
        std::transform(v1.begin(), v1.end(), result.begin(),
                       [a](Real i) { return i + a; });
        return result;
    }

    inline const Disposable<Array> operator+(Real a, const Array &v2) {
        Array result(v2.size());
        std::transform(v2.begin(), v2.end(), result.begin(),
                       [a](Real i) { return i + a; });
        return result;
    }

    inline const Disposable<Array> operator-(const Array &v1,
                                             const Array &v2) {
        QL_REQUIRE(v1.size() == v2.size(),
                   "arrays with different sizes (" << v1.size() << ", "
                                                   << v2.size() << ") cannot be subtracted");
        Array result(v1.size());
        std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(),
                       [](Real i, Real j) { return i - j; });
        return result;
    }

    inline const Disposable<Array> operator-(const Array &v1, Real a) {
        Array result(v1.size());
        std::transform(v1.begin(), v1.end(), result.begin(),
                       [a](Real i) { return i - a; });
        return result;
    }

    inline const Disposable<Array> operator-(Real a, const Array &v2) {
        Array result(v2.size());
        std::transform(v2.begin(), v2.end(), result.begin(),
                       [a](Real i) { return a - i; });
        return result;
    }

    inline const Disposable<Array> operator*(const Array &v1,
                                             const Array &v2) {
        QL_REQUIRE(v1.size() == v2.size(),
                   "arrays with different sizes (" << v1.size() << ", "
                                                   << v2.size() << ") cannot be multiplied");
        Array result(v1.size());
        std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(),
                       [](Real i, Real j) { return i * j; });
        return result;
    }

    inline const Disposable<Array> operator*(const Array &v1, Real a) {
        Array result(v1.size());
        std::transform(v1.begin(), v1.end(), result.begin(),
                       [a](Real i) { return i * a; });
        return result;
    }

    inline const Disposable<Array> operator*(Real a, const Array &v2) {
        Array result(v2.size());
        std::transform(v2.begin(), v2.end(), result.begin(),
                       [a](Real i) { return i * a; });
        return result;
    }

    inline const Disposable<Array> operator/(const Array &v1,
                                             const Array &v2) {
        QL_REQUIRE(v1.size() == v2.size(),
                   "arrays with different sizes (" << v1.size() << ", "
                                                   << v2.size() << ") cannot be divided");
        Array result(v1.size());
        std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(),
                       [](Real i, Real j) { return i / j; });
        return result;
    }

    inline const Disposable<Array> operator/(const Array &v1, Real a) {
        Array result(v1.size());
        std::transform(v1.begin(), v1.end(), result.begin(),
                       [a](Real i) { return i / a; });
        return result;
    }

    inline const Disposable<Array> operator/(Real a, const Array &v2) {
        Array result(v2.size());
        std::transform(v2.begin(), v2.end(), result.begin(),
                       [a](Real i) { return a / i; });
        return result;
    }

    // functions

    inline const Disposable<Array> Abs(const Array &v) {
        Array result(v.size());
        std::transform(v.begin(), v.end(), result.begin(),
                       [](Real i) { return std::fabs(i); });
        return result;
    }

    inline const Disposable<Array> Sqrt(const Array &v) {
        Array result(v.size());
        std::transform(v.begin(), v.end(), result.begin(),
                       [](Real i) { return std::sqrt(i); });
        return result;
    }

    inline const Disposable<Array> Log(const Array &v) {
        Array result(v.size());
        std::transform(v.begin(), v.end(), result.begin(),
                       [](Real i) { return std::log(i); });
        return result;
    }

    inline const Disposable<Array> Exp(const Array &v) {
        Array result(v.size());
        std::transform(v.begin(), v.end(), result.begin(),
                       [](Real i) { return std::exp(i); });
        return result;
    }

    inline const Disposable<Array> Pow(const Array &v, Real alpha) {
        Array result(v.size());
        std::transform(v.begin(), v.end(), result.begin(),
                       [alpha](Real i) { return pow(i, alpha); });

        return result;
    }


    inline void swap(Array &v, Array &w) noexcept {
        v.swap(w);
    }

    inline std::ostream &operator<<(std::ostream &out, const Array &a) {
        std::streamsize width = out.width();
        out << "[ ";
        if (!a.empty()) {
            for (Size n = 0; n < a.size() - 1; ++n)
                out << std::setw(int(width)) << a[n] << "; ";
            out << std::setw(int(width)) << a.back();
        }
        out << " ]";
        return out;
    }

}


#endif
