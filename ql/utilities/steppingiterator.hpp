#ifndef quantlib_stepping_iterator_hpp
#define quantlib_stepping_iterator_hpp

#include <ql/types.hpp>


namespace {
	using QuantLib::BigInteger, QuantLib::Size;

	template <class Iter>
	class step_iterator {
	public:
        
		using difference_type = typename Iter::difference_type;
		using value_type = typename Iter::value_type;
		using reference = typename Iter::reference;
		using pointer = typename Iter::pointer;
		using iterator_category = typename Iter::iterator_category;

		explicit step_iterator(const Iter& other, BigInteger step) :
			step_(step), it_(other) {}

		Size step() const { return static_cast<Size>(this->step_); }
		// iterator adapter interface

		void increment() {
			std::advance(it_, step_);
		}
		void decrement() {
			std::advance(it_, -step_);
		}
		void advance(difference_type n) {
			it_ += n*step_;
		}
		difference_type
			distance_to(const step_iterator& i) const {
			return (i.it_ - it_) / step_;
		}



		bool operator==(const step_iterator & other) const { return it_ == other.it_; };

        bool operator!=(const step_iterator & other) const { return it_ != other.it_;}

        bool operator<(const step_iterator & other) const { return it_ < other.it_; }
        bool operator>(const step_iterator & other) const { return it_ > other.it_; }
        bool operator<=(const step_iterator & other) const { return it_ <= other.it_; }
        bool operator>=(const step_iterator & other) const { return it_ >= other.it_; }

        step_iterator &operator++() { it_ += step_;  return *this; }
        step_iterator operator++(int) { step_iterator<Iter>new_iter(*this); this->it_ += step_; return new_iter;}
        step_iterator &operator--() { it_ -= step_;  return *this; }
        step_iterator operator--(int) { step_iterator<Iter>new_iter(*this); this->it_ -= step_; return new_iter;}
        step_iterator &operator+=(difference_type d) { it_ += d * step_; return *this; }
        step_iterator &operator-=(difference_type d) { it_ -= d * step_; return *this; }
        step_iterator operator+(difference_type d) const {
            step_iterator<Iter>new_iter(*this); new_iter += d; return new_iter;
        }
        friend step_iterator operator+(difference_type d, const step_iterator & iter) { step_iterator<Iter>new_iter(iter); new_iter.advance(d); return new_iter; }
		step_iterator operator-(difference_type d) const {
		step_iterator<Iter>new_iter(*this); new_iter -= d; return new_iter;
		}
        difference_type operator-(const step_iterator& other) const { return -distance_to(other); }

		reference operator*() const { return *it_; }
		pointer operator->() const { return &(*it_); }

        reference operator[](difference_type d) const { return it_[d * step_]; }

	private:
		BigInteger step_;
		Iter it_;

	};
  }
#endif
