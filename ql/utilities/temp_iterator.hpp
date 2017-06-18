
#ifndef quantlib_temp_iterator_hpp
#define quantlib_temp_iterator_hpp

#include <iterator>
#include <vector>
#include <ql/types.hpp>
#include <functional>


namespace QuantLib {

	template <class valueT>
	class temp_iterator {
	public:

		typedef typename std::vector<valueT>::difference_type difference_type;
		typedef typename std::vector<valueT>::value_type value_type;
		typedef typename std::vector<valueT>::reference reference;
		typedef valueT* pointer;
		typedef std::random_access_iterator_tag iterator_category;

		temp_iterator() = default;
		explicit temp_iterator(const typename std::vector<valueT>::iterator& other, BigInteger step) :
			step_(step), super_t(other) {}

		template <class valueTO, typename = std::enable_if_t
		<std::is_convertible<valueTO, valueT>::value>>
		explicit temp_iterator(const temp_iterator<valueT>& i)
			:  step_(static_cast<BigInteger>(i.step_)), super_t(i.super_t) {}

		Size step() const { return static_cast<Size>(this->step_); }
		// iterator adapter interface

		void increment() {
			std::advance(super_t, step_);
		}
		void decrement() {
			std::advance(super_t, -step_);
		}
		void advance(typename std::vector<valueT>::iterator::difference_type n) {
			super_t += n*step_;
		}

		typename std::vector<valueT>::iterator::difference_type
			distance_to(const temp_iterator& i) const {
			return (i.super_t - super_t) / step_;
		}

		temp_iterator(const temp_iterator &) = default;

		~temp_iterator() = default;

		temp_iterator &operator=(const temp_iterator &) = default;

		bool operator==(const temp_iterator & other) const { return super_t == other.super_t && step_ == other.step_; };

		bool operator!=(const temp_iterator & other) const { return super_t != other.super_t; }

		//        bool operator<(const temp_iterator &) const; //optional
		//        bool operator>(const temp_iterator &) const; //optional
		//        bool operator<=(const temp_iterator &) const; //optional
		//        bool operator>=(const temp_iterator &) const; //optional

		temp_iterator &operator++() { return super_t++; }

		temp_iterator &operator--() { return super_t--; }
		temp_iterator &operator+=(size_type); //optional
		temp_iterator operator+(size_type) const; //optional
		friend temp_iterator operator+(size_type, const temp_iterator &); //optional
		temp_iterator &operator-=(size_type); //optional
		temp_iterator operator-(size_type) const; //optional
		difference_type operator-(temp_iterator) const; //optional

		reference operator*() const { return *super_t; }
		pointer operator->() const { return &super_t; }

	private:
		BigInteger step_;
		std::vector<valueT>::iterator super_t;

	};
	
	/*   class const_temp_iterator {
	   public:
		   typedef std::vector<valueT>::difference_type difference_type;
		   typedef std::vector<valueT>::value_type value_type;
		   typedef std::vector<valueT>::reference const_reference;
		   typedef std::vector<valueT>::pointer const_pointer;
		   typedef std::random_access_iterator_tag iterator_category; //or another tag

		   const_iterator();

		   const_iterator(const const_iterator &);

		   const_iterator(const const_temp_iterator &);

		   ~const_iterator();

		   const_iterator &operator=(const const_iterator &);

		   bool operator==(const const_iterator &) const;

		   bool operator!=(const const_iterator &) const;

		   bool operator<(const const_iterator &) const; //optional
		   bool operator>(const const_iterator &) const; //optional
		   bool operator<=(const const_iterator &) const; //optional
		   bool operator>=(const const_iterator &) const; //optional

		   const_iterator &operator++();

		   const_iterator operator++(int); //optional
		   const_iterator &operator--(); //optional
		   const_iterator operator--(int); //optional
		   const_iterator &operator+=(size_type); //optional
		   const_iterator operator+(size_type) const; //optional
		   friend const_iterator operator+(size_type, const const_iterator &); //optional
		   const_iterator &operator-=(size_type); //optional
		   const_iterator operator-(size_type) const; //optional
		   difference_type operator-(const_iterator) const; //optional

		   const_reference operator*() const;

		   const_pointer operator->() const;

		   const_reference operator[](size_type) const; //optional
	   };

	   typedef std::reverse_iterator<temp_iterator> reverse_temp_iterator; //optional
	   typedef std::reverse_iterator<const_temp_iterator> const_reverse_temp_iterator; //optional */
   }
#endif
