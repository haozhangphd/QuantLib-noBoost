/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifndef quantlib_test_circular_buffer_hpp
#define quantlib_test_circular_buffer_hpp

#include <deque>
#include <iostream>
#include <ql/types.hpp>

// a VERY minimum circular buffer class, not used in QuantLib
// only used in one of the test cases
using QuantLib::Real;

namespace {
    class circular_buffer {
    public:

        using value_type = Real;
        using reference = Real &;
        using const_reference = const Real &;
        using pointer = Real *;
        using const_pointer = const Real *;
        using size_type = size_t;
        using iterator = std::deque<Real>::iterator;
        using const_iterator = std::deque<Real>::const_iterator;


        circular_buffer() : data_(), size_(0) {}

        explicit circular_buffer(
                size_t capacity)
                : data_(capacity), size_(capacity) {}

        template<class InputIterator>
        circular_buffer(InputIterator first, InputIterator last): data_(first, last),
                                                                  size_(std::distance(first, last)) {}

        circular_buffer(std::initializer_list<Real> init) : data_(init), size_(init.size()) {}

        reference operator[](size_t n) { return data_[n]; }

        const_reference operator[](size_type n) const { return data_[n]; }

        iterator begin() noexcept { return data_.begin(); }

        const_iterator begin() const noexcept { return data_.begin(); }

        iterator end() noexcept { return data_.end(); }

        const_iterator end() const noexcept { return data_.end(); }

        void push_back(Real v) {
            if (full()) {
                data_.pop_front();
                data_.push_back(v);
            } else
                data_.push_back(v);

        }


        void push_front(Real v) {
            if (full()) {
                data_.pop_back();
                data_.push_front(v);
            } else
                data_.push_front(v);

        }

        size_t size() const {
            return data_.size();
        }

        bool full() {
            return data_.size() == size_;
        }

        bool operator==(const circular_buffer &other) const { return data_ == other.data_; }

        friend std::ostream &operator<<(std::ostream &os, const circular_buffer &c) {
            os << "{ ";
            for (auto i = c.data_.begin(); i != std::prev(c.data_.end()); ++i)
                os << *i << ", ";
            os << c.data_.back() << " }";
            return os;
        }


    private:
        std::deque<Real> data_;
        size_t size_;
    };


}
#endif

