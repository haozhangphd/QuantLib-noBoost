/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Hao Zhang

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

#ifndef quantlib_transformiterator_hpp
#define quantlib_transformiterator_hpp

#include <ql/types.hpp>

namespace {
    using QuantLib::Real;
    std::function<Real *(Real &)> address_of = [](Real &x) { return &x; };
    std::function<const Real *(const Real &)>
            const_address_of = [](const Real &x) -> const Real * { return &x; };

    template<typename Func, typename Iter>
    class transformIterator {
    public:
        using difference_type = typename Iter::difference_type;
        using value_type = typename Iter::value_type;
        using reference = typename Iter::reference;
        using pointer = typename Iter::pointer;
        using iterator_category = typename Iter::iterator_category;

        transformIterator() = default;

        explicit transformIterator(const Iter &other, const Func &func) : it_(other),
                                                                          func_(func) {}

        inline typename Func::result_type operator*() const { return func_(*it_); }

        inline typename Func::result_type operator[](difference_type i) const {
            return func_(it_[i]);
        }

        void increment() {
            it_.increment();
        }

        void decrement() {
            it_.decrement();
        }

        void advance(difference_type n) {
            it_.advance(n);
        }

        difference_type
        distance_to(const transformIterator &i) const {
            return i.it_ - it_;
        }

        bool operator==(const transformIterator &other) const { return it_ == other.it_; };

        bool operator!=(const transformIterator &other) const { return it_ != other.it_; }

        bool operator<(const transformIterator &other) const { return it_ < other.it_; }

        bool operator>(const transformIterator &other) const { return it_ > other.it_; }

        bool operator<=(const transformIterator &other) const { return it_ <= other.it_; }

        bool operator>=(const transformIterator &other) const { return it_ >= other.it_; }

        transformIterator &operator++() {
            ++it_;
            return *this;
        }

        transformIterator operator++(int) {
            transformIterator<Func, Iter> new_iter(*this);
            ++(this->it_);
            return new_iter;
        }

        transformIterator &operator--() {
            --it_;
            return *this;
        }

        transformIterator operator--(int) {
            transformIterator<Func, Iter> new_iter(*this);
            --(this->it_);
            return new_iter;
        }

        transformIterator &operator+=(difference_type d) {
            it_ += d;
            return *this;
        }

        transformIterator &operator-=(difference_type d) {
            it_ -= d;
            return *this;
        }


    protected:
        Iter it_;
        Func func_;

    };

// transformIterator binary operators
    template<typename Func, typename Iter>
    transformIterator<Func, Iter> operator+(const transformIterator<Func, Iter> &iter,
                                            typename Iter::difference_type d) {
        transformIterator<Func, Iter> new_iter(iter);
        new_iter += d;
        return new_iter;
    }

    template<typename Func, typename Iter>
    transformIterator<Func, Iter> operator+(typename Iter::difference_type d,
                                            const transformIterator<Func, Iter> &iter) {
        transformIterator<Func, Iter> new_iter(iter);
        new_iter.advance(d);
        return new_iter;
    }

    template<typename Func, typename Iter>
    transformIterator<Func, Iter> operator-(const transformIterator<Func, Iter> &iter,
                                            typename Iter::difference_type d) {
        transformIterator<Func, Iter> new_iter(iter);
        new_iter -= d;
        return new_iter;
    }

    template<typename Func, typename Iter>
    typename Iter::difference_type operator-(const transformIterator<Func, Iter> &iter,
                                             const transformIterator<Func, Iter> &other) {
        return -iter.distance_to(other);
    }


    template<typename BaseIter, typename IndexIter>
    class permutationIterator :
            public transformIterator<std::function<typename BaseIter::reference(int)>, IndexIter> {
    public:
        explicit permutationIterator(const BaseIter &base, const IndexIter &index) :
                baseIter_(base), transformIterator<std::function<typename BaseIter::reference(int)>, IndexIter>() {
            this->it_ = index;
            this->func_ = taking_index;
        }

    private:
        BaseIter baseIter_;
        std::function<typename BaseIter::reference(int)> taking_index =
                [this](int i) -> typename BaseIter::reference { return baseIter_[i]; };
    };

    //permutation iterator binary operators
    template<typename BaseIter, typename IndexIter>
    permutationIterator<BaseIter, IndexIter> operator+(const permutationIterator<BaseIter, IndexIter> &iter,
                                                       typename IndexIter::difference_type d) {
        permutationIterator<BaseIter, IndexIter> new_iter(iter);
        new_iter += d;
        return new_iter;
    }

    template<typename BaseIter, typename IndexIter>
    permutationIterator<BaseIter, IndexIter> operator+(typename IndexIter::difference_type d,
                                                       const permutationIterator<BaseIter, IndexIter> &iter) {
        permutationIterator<BaseIter, IndexIter> new_iter(iter);
        new_iter.advance(d);
        return new_iter;
    }

    template<typename BaseIter, typename IndexIter>
    permutationIterator<BaseIter, IndexIter> operator-(const permutationIterator<BaseIter, IndexIter> &iter,
                                                       typename IndexIter::difference_type d) {
        permutationIterator<BaseIter, IndexIter> new_iter(iter);
        new_iter -= d;
        return new_iter;
    }

    template<typename BaseIter, typename IndexIter>
    typename IndexIter::difference_type
    operator-(const permutationIterator<BaseIter, IndexIter> &iter,
              const permutationIterator<BaseIter, IndexIter> &other) {
        return -iter.distance_to(other);
    }

}
#endif
