/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007 StatPro Italia srl
 Copyright (C) 2004, 2005, 2006 Ferdinando Ametrano
 Copyright (C) 2006 Katiuscia Manzoni
 Copyright (C) 2006 Toyin Akin
 Copyright (C) 2015 Klaus Spanderen

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

#include <ql/time/date.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/errors.hpp>

#include <iomanip>
#include <ctime>
#include <chrono>
#include <cmath>
#include <ratio>

namespace std {
    using ::time; using ::time_t; using ::tm;
}

#ifdef QL_HIGH_RESOLUTION_DATE
using std::chrono::time_point;
using std::chrono::system_clock;
using std::chrono::duration;
using std::chrono::floor;
using std::chrono::time_point_cast;
#ifdef _MSC_VER
#define timegm _mkgmtime
#endif
#endif


namespace QuantLib {
	// Year 1900 IS NOT a leap year
    bool Date::isLeap(Year y) {
        return (!(y % 4)) && ((y % 100) || (!(y % 400)));
    }

    Day Date::monthLength(Month m, bool leapYear) {
        if (m == 4 || m == 6 || m == 9 || m == 11)
            return 30;
        else if (m == 2)
            return leapYear ? 29 : 28;
        else return 31;
    }

#ifndef QL_HIGH_RESOLUTION_DATE

    // constructors
    Date::Date()
            : serialNumber_(Date::serial_type(0)) {}

    Date::Date(Date::serial_type serialNumber)
            : serialNumber_(serialNumber) {
        checkSerialNumber(serialNumber);
    }

    Date::Date(Day d, Month m, Year y) {
        QL_REQUIRE(y > 1900 && y < 2200,
                   "year " << y << " out of bound. It must be in [1901,2199]");
        QL_REQUIRE(Integer(m) > 0 && Integer(m) < 13,
                   "month " << Integer(m)
                            << " outside January-December range [1,12]");

        bool leap = isLeap(y);
        Day len = monthLength(m, leap), offset = monthOffset(m, leap);
        QL_REQUIRE(d <= len && d > 0,
                   "day outside month (" << Integer(m) << ") day-range "
                                         << "[1," << len << "]");

        serialNumber_ = d + offset + yearOffset(y);
    }

    Month Date::month() const {
        Day d = dayOfYear(); // dayOfYear is 1 based
        Integer m = d / 30 + 1;
        bool leap = isLeap(year());
        while (d <= monthOffset(Month(m), leap))
            --m;
        while (d > monthOffset(Month(m + 1), leap))
            ++m;
        return Month(m);
    }

    Year Date::year() const {
        Year y = (serialNumber_ / 365) + 1900;
        // yearOffset(y) is December 31st of the preceding year
        if (serialNumber_ <= yearOffset(y))
            --y;
        return y;
    }

    Date &Date::operator+=(Date::serial_type days) {
        Date::serial_type serial = serialNumber_ + days;
        checkSerialNumber(serial);
        serialNumber_ = serial;
        return *this;
    }

    Date &Date::operator+=(const Period &p) {
        advance(p.length(), p.units());
        return *this;
    }

    Date &Date::operator-=(Date::serial_type days) {
        Date::serial_type serial = serialNumber_ - days;
        checkSerialNumber(serial);
        serialNumber_ = serial;
        return *this;
    }

    Date &Date::operator-=(const Period &p) {
        advance(-p.length(), p.units());
        return *this;
    }

    Date &Date::operator++() {
        Date::serial_type serial = serialNumber_ + 1;
        checkSerialNumber(serial);
        serialNumber_ = serial;
        return *this;
    }

    Date &Date::operator--() {
        Date::serial_type serial = serialNumber_ - 1;
        checkSerialNumber(serial);
        serialNumber_ = serial;
        return *this;
    }

    void Date::advance(Integer n, TimeUnit units) {
        switch (units) {
            case Days:
                serialNumber_ += n;
                break;
            case Weeks:
                serialNumber_ += 7 * n;
                break;
            case Months: {
                Day d = dayOfMonth();
                Integer m = (Integer(month()) + n) % 12;
                if (m <= 0)
                    m = 12 + m;
                Year y = year() + floor((Integer(month()) + n - 1) / 12.0);

                QL_ENSURE(y >= 1900 && y <= 2199,
                          "year " << y << " out of bounds. "
                                  << "It must be in [1901,2199]");

                Integer length = monthLength(Month(m), isLeap(y));
                d = std::min(length, d);

                serialNumber_ = d + monthOffset(static_cast<Month>(m), isLeap(y)) + yearOffset(y);
                break;
            }
            case Years: {
                Day d = dayOfMonth();
                Month m = month();
                Year y = year() + n;

                QL_ENSURE(y >= 1900 && y <= 2199,
                          "year " << y << " out of bounds. "
                                  << "It must be in [1901,2199]");

                if (d == 29 && m == February && !isLeap(y))
                    d = 28;

                serialNumber_ = d + monthOffset(m, isLeap(y)) + yearOffset(y);
                break;
            }
            default:
                QL_FAIL("undefined time units");
        }
    }

    Integer Date::monthOffset(Month m, bool leapYear) {
        static const Integer MonthOffset[] = {
                0, 31, 59, 90, 120, 151,   // Jan - Jun
                181, 212, 243, 273, 304, 334,   // Jun - Dec
                365     // used in dayOfMonth to bracket day
        };
        static const Integer MonthLeapOffset[] = {
                0, 31, 60, 91, 121, 152,   // Jan - Jun
                182, 213, 244, 274, 305, 335,   // Jun - Dec
                366     // used in dayOfMonth to bracket day
        };
        return (leapYear ? MonthLeapOffset[m - 1] : MonthOffset[m - 1]);
    }


    Date::serial_type Date::yearOffset(Year y) {
        // the list of all December 31st in the preceding year
        // e.g. for 1901 yearOffset[1] is 365, that is, December 31 1900
        static const Date::serial_type YearOffset[] = {
                // 1900-1909
                0, 365, 730, 1095, 1460, 1826, 2191, 2556, 2921, 3287,
                // 1910-1919
                3652, 4017, 4382, 4748, 5113, 5478, 5843, 6209, 6574, 6939,
                // 1920-1929
                7304, 7670, 8035, 8400, 8765, 9131, 9496, 9861, 10226, 10592,
                // 1930-1939
                10957, 11322, 11687, 12053, 12418, 12783, 13148, 13514, 13879, 14244,
                // 1940-1949
                14609, 14975, 15340, 15705, 16070, 16436, 16801, 17166, 17531, 17897,
                // 1950-1959
                18262, 18627, 18992, 19358, 19723, 20088, 20453, 20819, 21184, 21549,
                // 1960-1969
                21914, 22280, 22645, 23010, 23375, 23741, 24106, 24471, 24836, 25202,
                // 1970-1979
                25567, 25932, 26297, 26663, 27028, 27393, 27758, 28124, 28489, 28854,
                // 1980-1989
                29219, 29585, 29950, 30315, 30680, 31046, 31411, 31776, 32141, 32507,
                // 1990-1999
                32872, 33237, 33602, 33968, 34333, 34698, 35063, 35429, 35794, 36159,
                // 2000-2009
                36524, 36890, 37255, 37620, 37985, 38351, 38716, 39081, 39446, 39812,
                // 2010-2019
                40177, 40542, 40907, 41273, 41638, 42003, 42368, 42734, 43099, 43464,
                // 2020-2029
                43829, 44195, 44560, 44925, 45290, 45656, 46021, 46386, 46751, 47117,
                // 2030-2039
                47482, 47847, 48212, 48578, 48943, 49308, 49673, 50039, 50404, 50769,
                // 2040-2049
                51134, 51500, 51865, 52230, 52595, 52961, 53326, 53691, 54056, 54422,
                // 2050-2059
                54787, 55152, 55517, 55883, 56248, 56613, 56978, 57344, 57709, 58074,
                // 2060-2069
                58439, 58805, 59170, 59535, 59900, 60266, 60631, 60996, 61361, 61727,
                // 2070-2079
                62092, 62457, 62822, 63188, 63553, 63918, 64283, 64649, 65014, 65379,
                // 2080-2089
                65744, 66110, 66475, 66840, 67205, 67571, 67936, 68301, 68666, 69032,
                // 2090-2099
                69397, 69762, 70127, 70493, 70858, 71223, 71588, 71954, 72319, 72684,
                // 2100-2109
                73049, 73414, 73779, 74144, 74509, 74875, 75240, 75605, 75970, 76336,
                // 2110-2119
                76701, 77066, 77431, 77797, 78162, 78527, 78892, 79258, 79623, 79988,
                // 2120-2129
                80353, 80719, 81084, 81449, 81814, 82180, 82545, 82910, 83275, 83641,
                // 2130-2139
                84006, 84371, 84736, 85102, 85467, 85832, 86197, 86563, 86928, 87293,
                // 2140-2149
                87658, 88024, 88389, 88754, 89119, 89485, 89850, 90215, 90580, 90946,
                // 2150-2159
                91311, 91676, 92041, 92407, 92772, 93137, 93502, 93868, 94233, 94598,
                // 2160-2169
                94963, 95329, 95694, 96059, 96424, 96790, 97155, 97520, 97885, 98251,
                // 2170-2179
                98616, 98981, 99346, 99712, 100077, 100442, 100807, 101173, 101538, 101903,
                // 2180-2189
                102268, 102634, 102999, 103364, 103729, 104095, 104460, 104825, 105190, 105556,
                // 2190-2199
                105921, 106286, 106651, 107017, 107382, 107747, 108112, 108478, 108843, 109208,
                // 2200
                109573
        };
        return YearOffset[y - 1900];
    }

#else

    namespace {
        time_point<system_clock, std::chrono::microseconds> dateToTimePoint(Day d, Month m, Year y) {
            tm timeinfo = {};
            timeinfo.tm_year = y - 1900;
            timeinfo.tm_mon = m - 1;
            timeinfo.tm_mday = d;
            return time_point_cast<std::chrono::microseconds>(
                    system_clock::from_time_t(timegm(&timeinfo)));

        }

        std::chrono::time_point<std::chrono::system_clock, std::chrono::microseconds> ref_date =
                dateToTimePoint(31, December, 1899);

    }

    void Date::advance(Integer n, TimeUnit units) {

        switch (units) {
            case Days:
                dateTime_ += duration<long int, std::ratio<86400>>{n};
                break;
            case Weeks:
                dateTime_ += duration<long int, std::ratio<604800>>{n};
                break;
            case Months: {
                time_t t{system_clock::to_time_t(dateTime_)};
                tm *tt = gmtime(&t);
                Integer month_bare = (tt->tm_mon + n) % 12;
                if (month_bare < 0)
                    month_bare = 12 + month_bare;
                Integer year = tt->tm_year + static_cast<Integer>(floor((tt->tm_mon + n) / 12.0));

                tt->tm_mon = month_bare;
                tt->tm_year = year;

                tt->tm_mday = std::min(tt->tm_mday,
                                       Date::monthLength(static_cast<Month>(tt->tm_mon + 1),
                                                         Date::isLeap(tt->tm_year)));
                dateTime_ = time_point_cast<std::chrono::microseconds>(
                        system_clock::from_time_t(timegm(tt)));

                break;
            }
            case Years : {
                time_t t{system_clock::to_time_t(dateTime_)};
                tm *tt = gmtime(&t);
                tt->tm_year += n;
                tt->tm_mday = std::min(tt->tm_mday, Date::monthLength(static_cast<Month>(tt->tm_mon + 1),
                                                                      Date::isLeap(tt->tm_year)));
                dateTime_ = time_point_cast<std::chrono::microseconds>(
                        system_clock::from_time_t(timegm(tt)));

                break;
            }

            default:
                QL_FAIL("undefined time units");
        }
    }


    Date::Date()
            : dateTime_(ref_date) {}

    Date::Date(const std::chrono::time_point<std::chrono::system_clock, std::chrono::microseconds> &time)
            : dateTime_(time) {}

    Date::Date(Day d, Month m, Year y)
            : dateTime_(dateToTimePoint(d, m, y)) {}

    Date::Date(Day d, Month m, Year y,
               Hour hours, Minute minutes, Second seconds,
               Millisecond millisec, Microsecond microsec)
            : dateTime_(
            dateToTimePoint(d, m, y) + std::chrono::hours{hours} + std::chrono::minutes{minutes} +
            std::chrono::seconds{seconds} + std::chrono::milliseconds{millisec} +
            +std::chrono::microseconds{microsec}) {}

    Date::Date(Date::serial_type serialNumber)
            : dateTime_(ref_date +
                        std::chrono::hours(24 * serialNumber)) {
        checkSerialNumber(serialNumber);
    }

    Weekday Date::weekday() const {
        time_t t{system_clock::to_time_t(dateTime_)};
        tm *tt = gmtime(&t);
        return static_cast<Weekday>(tt->tm_wday + 1);

    }

    Day Date::dayOfMonth() const {
        time_t t{system_clock::to_time_t(dateTime_)};
        tm *tt = gmtime(&t);
        return tt->tm_mday;
    }

    Day Date::dayOfYear() const {
        time_t t{system_clock::to_time_t(dateTime_)};
        tm *tt = gmtime(&t);
        return tt->tm_yday + 1;
    }

    Month Date::month() const {
        time_t t{system_clock::to_time_t(dateTime_)};
        tm *tt = gmtime(&t);
        return static_cast<Month>(tt->tm_mon + 1);
    }

    Year Date::year() const {
        time_t t{system_clock::to_time_t(dateTime_)};
        tm *tt = gmtime(&t);
        return tt->tm_year + 1900;
    }

    Hour Date::hours() const {
        time_point<system_clock, duration<long int, std::ratio<86400>>> time_round =
                floor<duration<long int, std::ratio<86400>>>(dateTime_);
        return floor<std::chrono::hours>(dateTime_ - time_round).count();
    }

    Minute Date::minutes() const {
        time_point<system_clock, std::chrono::microseconds> time_round = floor<std::chrono::hours>(dateTime_);
        return floor<std::chrono::minutes>(dateTime_ - time_round).count();
    }

    Second Date::seconds() const {
        time_point<system_clock, std::chrono::microseconds> time_round = floor<std::chrono::minutes>(dateTime_);
        return floor<std::chrono::seconds>(dateTime_ - time_round).count();
    }

    Time Date::fractionOfDay() const {
        time_point<system_clock, std::chrono::microseconds> time_round = floor<duration<long int, std::ratio<86400>>>(
                dateTime_);
        return (std::chrono::duration_cast<duration<Time, std::ratio<86400>>>(dateTime_ - time_round)).count();
    }

    Time Date::fractionOfSecond() const {
        time_point<system_clock, std::chrono::microseconds> time_round = floor<std::chrono::seconds>(
                dateTime_);
        return (std::chrono::duration_cast<duration<Time>>(dateTime_ - time_round)).count();
    }

    Millisecond Date::milliseconds() const {
        time_point<system_clock, std::chrono::microseconds> time_round = floor<std::chrono::seconds>(dateTime_);
        return floor<std::chrono::milliseconds>(dateTime_ - time_round).count();
    }

    Microsecond Date::microseconds() const {
        time_point<system_clock, std::chrono::microseconds> time_round = floor<std::chrono::milliseconds>(dateTime_);
        return floor<std::chrono::microseconds>(dateTime_ - time_round).count();
    }

    long Date::ticksPerSecond() {
        return system_clock::period::den;
    }

    Date::serial_type Date::serialNumber() const {
        const Date::serial_type n = floor<duration<long int, std::ratio<86400>>>(dateTime_
                                                                                 - ref_date).count();
        checkSerialNumber(n);

        return n;
    }

    const time_point<system_clock, std::chrono::microseconds> &Date::dateTime() const { return dateTime_; }

    Date &Date::operator+=(Date::serial_type d) {
        dateTime_ += duration<long int, std::ratio<86400>>{d};
        return *this;
    }

    Date &Date::operator+=(const Period &p) {
        advance(p.length(), p.units());
        return *this;
    }

    Date &Date::operator-=(Date::serial_type d) {
        dateTime_ -= duration<long int, std::ratio<86400>>{d};
        return *this;
    }

    Date &Date::operator-=(const Period &p) {
        advance(-p.length(), p.units());
        return *this;
    }

    Date &Date::operator++() {
        dateTime_ += duration<long int, std::ratio<86400>>{1};
        return *this;
    }

    Date &Date::operator--() {
        dateTime_ -= duration<long int, std::ratio<86400>>{1};
        return *this;
    }

    Date Date::operator+(Date::serial_type days) const {
        Date retVal(*this);
        retVal += days;

        return retVal;
    }

    Date Date::operator-(Date::serial_type days) const {
        Date retVal(*this);
        retVal -= days;

        return retVal;
    }

    Date Date::operator+(const Period &p) const {
        Date retVal(*this);
        retVal += p;

        return retVal;
    }

    Date Date::operator-(const Period &p) const {
        Date retVal(*this);
        retVal -= p;

        return retVal;
    }

    Date Date::endOfMonth(const Date &d) {
        const Month m = d.month();
        const Year y = d.year();
        const Day eoM = Date::monthLength(m, Date::isLeap(y));

        return Date(eoM, m, y);
    }

    bool Date::isEndOfMonth(const Date &d) {
        return d.dayOfMonth() == Date::monthLength(d.month(), Date::isLeap(d.year()));
    }


    Date::serial_type operator-(const Date &d1, const Date &d2) {
        return floor<duration<long int, std::ratio<86400>>>(d1.dateTime() - d2.dateTime()).count();
    }

    Time daysBetween(const Date &d1, const Date &d2) {
        return (std::chrono::duration_cast<duration<Time, std::ratio<86400>>>(d2.dateTime() - d1.dateTime())).count();
    }

    bool operator<(const Date &d1, const Date &d2) {
        return (d1.dateTime() < d2.dateTime());
    }

    bool operator<=(const Date &d1, const Date &d2) {
        return (d1.dateTime() <= d2.dateTime());
    }

    bool operator>(const Date &d1, const Date &d2) {
        return (d1.dateTime() > d2.dateTime());
    }

    bool operator>=(const Date &d1, const Date &d2) {
        return (d1.dateTime() >= d2.dateTime());
    }

    bool operator==(const Date &d1, const Date &d2) {
        return (d1.dateTime() == d2.dateTime());
    }

    bool operator!=(const Date &d1, const Date &d2) {
        return (d1.dateTime() != d2.dateTime());
    }

#endif

    Date::serial_type Date::minimumSerialNumber() {
        return 366;       // Jan 1st, 1901
    }

    Date::serial_type Date::maximumSerialNumber() {
        return 109573;    // Dec 31st, 2199
    }

    void Date::checkSerialNumber(Date::serial_type serialNumber) {
        QL_REQUIRE(serialNumber >= minimumSerialNumber() &&
                   serialNumber <= maximumSerialNumber(),
                   "Date's serial number (" << serialNumber << ") outside "
                           "allowed range [" << minimumSerialNumber() <<
                                            "-" << maximumSerialNumber() << "], i.e. [" <<
                                            minDate() << "-" << maxDate() << "]");
    }

    Date Date::minDate() {
        static const Date minimumDate(minimumSerialNumber());
        return minimumDate;
    }

    Date Date::maxDate() {
        static const Date maximumDate(maximumSerialNumber());
        return maximumDate;
    }

    Date Date::operator++(int) {
        Date old(*this);
        ++*this; // use the pre-increment
        return old;
    }

    Date Date::operator--(int) {
        Date old(*this);
        --*this; // use the pre-decrement
        return old;
    }

    Date Date::todaysDate() {
        std::time_t t;

        if (std::time(&t) == std::time_t(-1)) // -1 means time() didn't work
            return Date();
        std::tm *lt = std::localtime(&t);
        return Date(Day(lt->tm_mday),
                    Month(lt->tm_mon + 1),
                    Year(lt->tm_year + 1900));
    }

    Date Date::nextWeekday(const Date &d, Weekday dayOfWeek) {
        Weekday wd = d.weekday();
        return d + ((wd > dayOfWeek ? 7 : 0) - wd + dayOfWeek);
    }

    Date Date::nthWeekday(Size nth, Weekday dayOfWeek,
                          Month m, Year y) {
        QL_REQUIRE(nth > 0,
                   "zeroth day of week in a given (month, year) is undefined");
        QL_REQUIRE(nth < 6,
                   "no more than 5 weekday in a given (month, year)");
        Weekday first = Date(1, m, y).weekday();
        Size skip = nth - (dayOfWeek >= first ? 1 : 0);
        return Date((1 + dayOfWeek + skip * 7) - first, m, y);
    }



// month formatting

    std::ostream &operator<<(std::ostream &out, Month m) {
        switch (m) {
            case January:
                return out << "January";
            case February:
                return out << "February";
            case March:
                return out << "March";
            case April:
                return out << "April";
            case May:
                return out << "May";
            case June:
                return out << "June";
            case July:
                return out << "July";
            case August:
                return out << "August";
            case September:
                return out << "September";
            case October:
                return out << "October";
            case November:
                return out << "November";
            case December:
                return out << "December";
            default:
                QL_FAIL("unknown month (" << Integer(m) << ")");
        }
    }


// date formatting

    std::ostream &operator<<(std::ostream &out, const Date &d) {
        return out << io::long_date(d);
    }

    namespace detail {

        struct FormatResetter {
            // An instance of this object will have undefined behaviour
            // if the object out passed in the constructor is destroyed
            // before this instance
            struct nopunct : std::numpunct<char> {
                std::string do_grouping() const { return ""; }
            };

            explicit FormatResetter(std::ostream &out)
                    : out_(&out), flags_(out.flags()), filler_(out.fill()),
                      loc_(out.getloc()) {
                std::locale loc(out.getloc(), new nopunct);
                out.imbue(loc);
                out << std::resetiosflags(
                        std::ios_base::adjustfield | std::ios_base::basefield |
                        std::ios_base::floatfield | std::ios_base::showbase |
                        std::ios_base::showpos | std::ios_base::uppercase);
                out << std::right;
            }

            ~FormatResetter() {
                out_->flags(flags_);
                out_->fill(filler_);
                out_->imbue(loc_);
            }

            std::ostream *out_;
            std::ios_base::fmtflags flags_;
            char filler_;
            std::locale loc_;
        };

        std::ostream &operator<<(std::ostream &out,
                                 const short_date_holder &holder) {
            const Date &d = holder.d;
            if (d == Date()) {
                out << "null date";
            } else {
                FormatResetter resetter(out);
                Integer dd = d.dayOfMonth(), mm = Integer(d.month()),
                        yyyy = d.year();
                char filler = out.fill();
                out << std::setw(2) << std::setfill('0') << mm << "/";
                out << std::setw(2) << std::setfill('0') << dd << "/";
                out << yyyy;
                out.fill(filler);
            }
            return out;
        }

        std::ostream &operator<<(std::ostream &out,
                                 const long_date_holder &holder) {
            const Date &d = holder.d;
            if (d == Date()) {
                out << "null date";
            } else {
                FormatResetter resetter(out);
                out << d.month() << " ";
                out << io::ordinal(d.dayOfMonth()) << ", ";
                out << d.year();
            }
            return out;
        }

        std::ostream &operator<<(std::ostream &out,
                                 const iso_date_holder &holder) {
            const Date &d = holder.d;
            if (d == Date()) {
                out << "null date";
            } else {
                FormatResetter resetter(out);
                Integer dd = d.dayOfMonth(), mm = Integer(d.month()),
                        yyyy = d.year();
                out << yyyy << "-";
                out << std::setw(2) << std::setfill('0') << mm << "-";
                out << std::setw(2) << std::setfill('0') << dd;
            }
            return out;
        }


#ifdef QL_HIGH_RESOLUTION_DATE

        std::ostream &operator<<(std::ostream &out,
                                 const iso_datetime_holder &holder) {
            const Date &d = holder.d;

            out << io::iso_date(d) << "T";
            FormatResetter resetter(out);
            const Hour hh = d.hours();
            const Minute mm = d.minutes();
            const Second s = d.seconds();
            const Millisecond millis = d.milliseconds();
            const Microsecond micros = d.microseconds();

            out << std::setw(2) << std::setfill('0') << hh << ":"
                << std::setw(2) << std::setfill('0') << mm << ":"
                << std::setw(2) << std::setfill('0') << s << ","
                << std::setw(3) << std::setfill('0') << millis
                << std::setw(3) << std::setfill('0') << micros;

            return out;
        }

        std::ostream &operator<<(std::ostream &out,
                                 const formatted_date_holder &holder) {
            const Date &d = holder.d;
            if (d == Date()) {
                out << "null date";
            } else {
                FormatResetter resetter(out);
                out.imbue(std::locale());
                std::time_t t{system_clock::to_time_t(d.dateTime())};
                std::tm *tm = std::gmtime(&t);
                out << std::put_time(tm, holder.f.c_str());
            }
            return out;
        }

#else

        std::ostream &operator<<(std::ostream &out,
                                 const formatted_date_holder &holder) {
            const Date &d = holder.d;
            if (d == Date()) {
                out << "null date";
            } else {
                FormatResetter resetter(out);
                out.imbue(std::locale());
                std::tm tm = {};
                tm.tm_year = d.year() - 1900;
                tm.tm_mon = static_cast<Integer>(d.month()) - 1;
                tm.tm_mday = d.dayOfMonth();
                out << std::put_time(&tm, holder.f.c_str());
            }
            return out;
        }

#endif
    }

    namespace io {
        detail::short_date_holder short_date(const Date &d) {
            return detail::short_date_holder(d);
        }

        detail::long_date_holder long_date(const Date &d) {
            return detail::long_date_holder(d);
        }

        detail::iso_date_holder iso_date(const Date &d) {
            return detail::iso_date_holder(d);
        }

        detail::formatted_date_holder formatted_date(const Date &d,
                                                     const std::string &f) {
            return detail::formatted_date_holder(d, f);
        }

#ifdef QL_HIGH_RESOLUTION_DATE

        detail::iso_datetime_holder iso_datetime(const Date &d) {
            return detail::iso_datetime_holder(d);
        }

#endif
    }
}
