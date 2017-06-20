#ifndef quantlib_stepping_iterator_hpp
#define quantlib_stepping_iterator_hpp

#include <ql/types.hpp>


namespace QuantLib {

	template <class Iter>
	class step_iterator {
	public:
        
		using difference_type = typename Iter::difference_type;
		using value_type = typename Iter::value_type;
		using reference = typename Iter::reference;
		using pointer = typename Iter::pointer;
		using iterator_category = std::random_access_iterator_tag;

		explicit step_iterator(const Iter& other, BigInteger step) :
			step_(step), super_t(other) {}

		template <class IterOther, typename = std::enable_if_t
		<std::is_convertible<IterOther, Iter>::value>>
		explicit step_iterator(const step_iterator<IterOther>& i)
			:  step_(static_cast<BigInteger>(i.step_)), super_t(i.super_t) {}

		Size step() const { return static_cast<Size>(this->step_); }
		// iterator adapter interface

		void increment() {
			std::advance(super_t, step_);
		}
		void decrement() {
			std::advance(super_t, -step_);
		}
		void advance(difference_type n) {
			super_t += n*step_;
		}
		difference_type
			distance_to(const step_iterator& i) const {
			return (i.super_t - super_t) / step_;
		}



		bool operator==(const step_iterator & other) const { return super_t == other.super_t; };

        bool operator!=(const step_iterator & other) const { return super_t != other.super_t;}

        bool operator<(const step_iterator & other) const { return super_t < other.super_t; }
        bool operator>(const step_iterator & other) const { return super_t > other.super_t; }
        bool operator<=(const step_iterator & other) const { return super_t <= other.super_t; }
        bool operator>=(const step_iterator & other) const { return super_t >= other.super_t; }

        step_iterator &operator++() { super_t += step_;  return *this; }
        step_iterator operator++(int) { step_iterator<Iter>new_iter(*this); this->super_t += step_; return new_iter;}
        step_iterator &operator--() { super_t -= step_;  return *this; }
        step_iterator operator--(int) { step_iterator<Iter>new_iter(*this); this->super_t -= step_; return new_iter;}
        step_iterator &operator+=(difference_type d) { super_t += d * step_; return *this; }
        step_iterator &operator-=(difference_type d) { super_t -= d * step_; return *this; }
        step_iterator operator+(difference_type d) const {
            step_iterator<Iter>new_iter(*this); new_iter += d; return new_iter;
        }
        friend step_iterator operator+(difference_type d, const step_iterator & iter) { step_iterator<Iter>new_iter(iter); new_iter.advance(d); return new_iter; }
		step_iterator operator-(difference_type d) const {
		step_iterator<Iter>new_iter(*this); new_iter -= d; return new_iter;
		}
        difference_type operator-(const step_iterator& other) const { return -distance_to(other); }

		reference operator*() const { return *super_t; }
		pointer operator->() const { return &(*super_t); }

        reference operator[](difference_type d) const { return super_t[d * step_]; }

	private:
		BigInteger step_;
		Iter super_t;

	};
  }
#endif
