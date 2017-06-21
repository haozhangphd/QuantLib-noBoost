/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*! \file threedimensionalarray.hpp
    \brief simple three dimensional array with minimum functionality.
*/

#ifndef quantlib_threedimensionalarray_hpp
#define quantlib_threedimensionalarray_hpp

#include <ql/utilities/steppingiterator.hpp>

namespace QuantLib {

    // An iterator for 2d matrix, requires 2 indices, It[i] gives a *pointer* to the first element of ith column

    template<typename Iter>
    class table_iterator_adaptor {
    public:
        using difference_type = typename Iter::difference_type;
        using value_type = typename Iter::value_type;
        using reference = typename Iter::reference;
        using pointer = typename Iter::pointer;
        using iterator_category = std::random_access_iterator_tag;

        explicit table_iterator_adaptor(const Iter &other, BigInteger step) : it_(other, step) {}

        typename step_iterator<Iter>::pointer
        operator[](typename step_iterator<Iter>::difference_type i) { return &it_[i]; }

    private:
        step_iterator <Iter> it_;

    };

    class threeDimensionalArray {

    private:
        using vecit = typename std::vector<Real>::iterator;
        using convecit = typename std::vector<Real>::const_iterator;

        std::vector<Real> data_;
        Size tables_, rows_, columns_;

    public:
        // Constructors, destructor, and assignment
        threeDimensionalArray();

        threeDimensionalArray(Size tables, Size rows, Size columns);

        threeDimensionalArray(Size tables, Size rows, Size columns, Real value);

        template<class Iterator>
        threeDimensionalArray(Size tables, Size rows, Size columns, Iterator begin, Iterator end);

        // Algebraic operators
        /* all matrices involved in an algebraic expression must have
                 the same size.
        */
        threeDimensionalArray &operator+=(const threeDimensionalArray &);

        threeDimensionalArray &operator-=(const threeDimensionalArray &);

        threeDimensionalArray &operator*=(Real);

        threeDimensionalArray &operator/=(Real);

        using iterator = std::vector<Real>::iterator;
        using const_iterator = std::vector<Real>::const_iterator;
        using reverse_iterator = std::vector<Real>::reverse_iterator;
        using const_reverse_iterator = std::vector<Real>::const_reverse_iterator;

        using table_iterator = table_iterator_adaptor<vecit>;
        using const_table_iterator = table_iterator_adaptor<convecit>;
        using reverse_table_iterator = std::reverse_iterator<table_iterator>;
        using const_reverse_table_iterator = std::reverse_iterator<const_table_iterator>;

        // Iterator access
        const_iterator begin() const;
        iterator begin();
        const_iterator end() const;
        iterator end();
        const_reverse_iterator rbegin() const;
        reverse_iterator rbegin();
        const_reverse_iterator rend() const;
        reverse_iterator rend();

        const_table_iterator table_begin(Size table) const;
        table_iterator table_begin(Size table);
        const_table_iterator table_end(Size table) const;
        table_iterator table_end(Size table);
        const_reverse_table_iterator table_rbegin(Size table) const;
        reverse_table_iterator table_rbegin(Size table);
        const_reverse_table_iterator table_rend(Size table) const;
        reverse_table_iterator table_rend(Size table);

        // Element access
        const_table_iterator operator[](Size) const;
        const_table_iterator at(Size) const;
        table_iterator operator[](Size);
        table_iterator at(Size);
        Real &operator()(Size i, Size j, Size k);

        // Inspectors
        Size tables() const;
        Size rows() const;
        Size columns() const;
        bool empty() const;
        Size size1() const;
        Size size2() const;
        Size size3() const;

        // Utilities
        void swap(threeDimensionalArray &) noexcept;

    };

    // algebraic operators

    threeDimensionalArray operator+(const threeDimensionalArray &, const threeDimensionalArray &);

    threeDimensionalArray operator-(const threeDimensionalArray &, const threeDimensionalArray &);

    threeDimensionalArray operator*(const threeDimensionalArray &, Real);

    threeDimensionalArray operator*(Real, const threeDimensionalArray &);

    threeDimensionalArray operator/(const threeDimensionalArray &, Real);

    // misc. operations

    void swap(threeDimensionalArray &, threeDimensionalArray &) noexcept;

    // inline definitions

    inline threeDimensionalArray::threeDimensionalArray()
            : data_(std::vector<Real>()), tables_(0), rows_(0), columns_(0) {}

    inline threeDimensionalArray::threeDimensionalArray(Size tables, Size rows, Size columns)
            : data_(tables * rows * columns),
              tables_(tables), rows_(rows), columns_(columns) {}

    inline threeDimensionalArray::threeDimensionalArray(Size tables, Size rows, Size columns, Real value)
            : data_(tables * rows * columns, value),
              tables_(tables), rows_(rows), columns_(columns) {}

    template<class Iterator>
    inline threeDimensionalArray::threeDimensionalArray(Size tables, Size rows, Size columns,
                                                        Iterator begin, Iterator end)
            : data_(tables * rows * columns),
              tables_(tables), rows_(rows), columns_(columns) {
        std::copy(begin, end, this->begin());
    }


    inline void threeDimensionalArray::swap(threeDimensionalArray &from) noexcept {
        using std::swap;
        data_.swap(from.data_);
        swap(tables_, from.tables_);
        swap(rows_, from.rows_);
        swap(columns_, from.columns_);
    }

    inline threeDimensionalArray &threeDimensionalArray::operator+=(const threeDimensionalArray &m) {
        QL_REQUIRE(tables_ == m.tables_ && rows_ == m.rows_ && columns_ == m.columns_,
                   "matrices with different sizes (" <<
                                                     m.tables_ << "x" << m.rows_ << "x" << m.columns_ << ", " <<
                                                     tables_ << "x" << rows_ << "x" << columns_ << ") cannot be "
                                                             "added");
        std::transform(begin(), end(), m.begin(),
                       begin(), std::plus<Real>());
        return *this;
    }

    inline threeDimensionalArray &threeDimensionalArray::operator-=(const threeDimensionalArray &m) {
        QL_REQUIRE(tables_ == m.tables_ && rows_ == m.rows_ && columns_ == m.columns_,
                   "matrices with different sizes (" <<
                                                     m.tables_ << "x" << m.rows_ << "x" << m.columns_ << ", " <<
                                                     tables_ << "x" << rows_ << "x" << columns_ << ") cannot be "
                                                             "subtracted");
        std::transform(begin(), end(), m.begin(), begin(),
                       std::minus<Real>());
        return *this;
    }

    inline threeDimensionalArray &threeDimensionalArray::operator*=(Real x) {
        std::transform(begin(), end(), begin(),
                       [x](Real r) { return x * r; });
        return *this;
    }

    inline threeDimensionalArray &threeDimensionalArray::operator/=(Real x) {
        std::transform(begin(), end(), begin(),
                       [x](Real r) { return r / x; });
        return *this;
    }

    inline threeDimensionalArray::const_iterator threeDimensionalArray::begin() const {
        return data_.begin();
    }

    inline threeDimensionalArray::iterator threeDimensionalArray::begin() {
        return data_.begin();
    }

    inline threeDimensionalArray::const_iterator threeDimensionalArray::end() const {
        return data_.end();
    }

    inline threeDimensionalArray::iterator threeDimensionalArray::end() {
        return data_.end();
    }

    inline threeDimensionalArray::const_reverse_iterator threeDimensionalArray::rbegin() const {
        return data_.rbegin();
    }

    inline threeDimensionalArray::reverse_iterator threeDimensionalArray::rbegin() {
        return data_.rbegin();
    }

    inline threeDimensionalArray::const_reverse_iterator threeDimensionalArray::rend() const {
        return data_.rend();
    }

    inline threeDimensionalArray::reverse_iterator threeDimensionalArray::rend() {
        return data_.rend();
    }

    inline threeDimensionalArray::const_table_iterator
    threeDimensionalArray::table_begin(Size table) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(table < tables_,
            "table index (" << table << ") must be less than " << tables_ <<
            ": 3d array cannot be accessed out of range");
#endif
        return const_table_iterator(data_.begin() + (rows_ * columns_ * table), columns_);
    }

    inline threeDimensionalArray::table_iterator threeDimensionalArray::table_begin(Size table) {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(table < tables_,
            "table index (" << table << ") must be less than " << tables_ <<
            ": 3d array cannot be accessed out of range");
#endif
        return table_iterator(data_.begin() + (rows_ * columns_ * table), columns_);
    }

    inline threeDimensionalArray::const_table_iterator threeDimensionalArray::table_end(Size table) const {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(table < tables_,
            "table index (" << table << ") must be less than " << tables_ <<
            ": 3d array cannot be accessed out of range");
#endif
        return const_table_iterator(data_.begin() + (rows_ * columns_ * (table + 1)), columns_);
    }

    inline threeDimensionalArray::table_iterator threeDimensionalArray::table_end(Size table) {
#if defined(QL_EXTRA_SAFETY_CHECKS)
        QL_REQUIRE(table < tables_,
            "table index (" << table << ") must be less than " << tables_ <<
            ": 3d array cannot be accessed out of range");
#endif
        return table_iterator(data_.begin() + rows_ * columns_ * (table + 1), columns_);
    }

    inline threeDimensionalArray::const_reverse_table_iterator
    threeDimensionalArray::table_rbegin(Size table) const {
        return const_reverse_table_iterator(table_end(table));
    }

    inline threeDimensionalArray::reverse_table_iterator threeDimensionalArray::table_rbegin(Size table) {
        return reverse_table_iterator(table_end(table));
    }

    inline threeDimensionalArray::const_reverse_table_iterator
    threeDimensionalArray::table_rend(Size table) const {
        return const_reverse_table_iterator(table_begin(table));
    }

    inline threeDimensionalArray::reverse_table_iterator threeDimensionalArray::table_rend(Size table) {
        return reverse_table_iterator(table_begin(table));
    }

    inline threeDimensionalArray::const_table_iterator
    threeDimensionalArray::operator[](Size table) const {
        return table_begin(table);
    }

    inline threeDimensionalArray::const_table_iterator
    threeDimensionalArray::at(Size table) const {
        QL_REQUIRE(table < tables_, "3d array access out of range");
        return table_begin(table);
    }

    inline threeDimensionalArray::table_iterator threeDimensionalArray::operator[](Size table) {
        return table_begin(table);
    }

    inline threeDimensionalArray::table_iterator threeDimensionalArray::at(Size table) {
        QL_REQUIRE(table < tables_, "3d array access out of range");
        return table_begin(table);
    }

    inline Real &threeDimensionalArray::operator()(Size i, Size j, Size k) {
        return data_[i * columns_ * rows_ + j * columns_ + k];
    }

    inline Size threeDimensionalArray::tables() const {
        return tables_;
    }

    inline Size threeDimensionalArray::rows() const {
        return rows_;
    }

    inline Size threeDimensionalArray::columns() const {
        return columns_;
    }

    inline Size threeDimensionalArray::size1() const {
        return tables_;
    }

    inline Size threeDimensionalArray::size2() const {
        return rows_;
    }

    inline Size threeDimensionalArray::size3() const {
        return columns_;
    }

    inline bool threeDimensionalArray::empty() const {
        return tables_ == 0 || rows_ == 0 || columns_ == 0;
    }

    inline threeDimensionalArray operator+(const threeDimensionalArray &m1,
                                           const threeDimensionalArray &m2) {
        QL_REQUIRE(m1.tables() == m2.tables() &&
                   m1.rows() == m2.rows() &&
                   m1.columns() == m2.columns(),
                   "matrices with different sizes (" <<
                                                     m1.tables() << "x" << m1.rows() << "x" << m1.columns() << ", " <<
                                                     m2.tables() << "x" << m2.rows() << "x" << m2.columns()
                                                     << ") cannot be "
                                                             "added");
        threeDimensionalArray temp(m1.tables(), m1.rows(), m1.columns());
        std::transform(m1.begin(), m1.end(), m2.begin(), temp.begin(),
                       std::plus<Real>());
        return temp;
    }

    inline threeDimensionalArray operator-(const threeDimensionalArray &m1,
                                           const threeDimensionalArray &m2) {
        QL_REQUIRE(m1.tables() == m2.tables() &&
                   m1.rows() == m2.rows() &&
                   m1.columns() == m2.columns(),
                   "matrices with different sizes (" <<
                                                     m1.tables() << "x" << m1.rows() << "x" << m1.columns() << ", " <<
                                                     m2.tables() << "x" << m2.rows() << "x" << m2.columns()
                                                     << ") cannot be "
                                                             "subtracted");
        threeDimensionalArray temp(m1.tables(), m1.rows(), m1.columns());
        std::transform(m1.begin(), m1.end(), m2.begin(), temp.begin(),
                       std::minus<Real>());
        return temp;
    }

    inline threeDimensionalArray operator*(const threeDimensionalArray &m, Real x) {
        threeDimensionalArray temp(m.tables(), m.rows(), m.columns());
        std::transform(m.begin(), m.end(), temp.begin(),
                       [x](Real r) { return x * r; });
        return temp;
    }

    inline threeDimensionalArray operator*(Real x, const threeDimensionalArray &m) {
        threeDimensionalArray temp(m.tables(), m.rows(), m.columns());
        std::transform(m.begin(), m.end(), temp.begin(),
                       [x](Real r) { return x * r; });
        return temp;
    }

    inline threeDimensionalArray operator/(const threeDimensionalArray &m, Real x) {
        threeDimensionalArray temp(m.tables(), m.rows(), m.columns());
        std::transform(m.begin(), m.end(), temp.begin(),
                       [x](Real r) { return r / x; });
        return temp;
    }


    inline void swap(threeDimensionalArray &m1, threeDimensionalArray &m2) noexcept {
        m1.swap(m2);
    }

}

#endif
