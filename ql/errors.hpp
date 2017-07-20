/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005 StatPro Italia srl

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

/*! \file errors.hpp
    \brief Classes and functions for error handling.
*/

#ifndef quantlib_errors_hpp
#define quantlib_errors_hpp

#include <ql/qldefines.hpp>
#include <memory>
#include <algorithm>
#include <exception>
#include <sstream>
#include <string>

namespace QuantLib {

    //! Base error class
    class Error : public std::exception {
      public:
        /*! The explicit use of this constructor is not advised.
            Use the QL_FAIL macro instead.
        */
        Error(const std::string& file,
              long line,
              const std::string& functionName,
              const std::string& message = "");
        //! returns the error message.
        const char* what() const noexcept;
      private:
        std::shared_ptr<std::string> message_;
    };

}

#define MULTILINE_MACRO_BEGIN do {
    #define MULTILINE_MACRO_END } while(false)
#endif


/*! \def QL_FAIL
    \brief throw an error (possibly with file and line information)
*/
#ifndef _MSC_VER
#define QL_FAIL(message) \
MULTILINE_MACRO_BEGIN \
    std::ostringstream _ql_msg_stream; \
    _ql_msg_stream << message; \
    throw QuantLib::Error(__FILE__,__LINE__, \
                          __PRETTY_FUNCTION__,_ql_msg_stream.str()); \
MULTILINE_MACRO_END
#else
#define QL_FAIL(message) \
MULTILINE_MACRO_BEGIN \
    std::ostringstream _ql_msg_stream; \
    _ql_msg_stream << message; \
    throw QuantLib::Error(__FILE__,__LINE__, \
                          __FUNCSIG__,_ql_msg_stream.str()); \
MULTILINE_MACRO_END
#endif

/*! \def QL_ASSERT
    \brief throw an error if the given condition is not verified
*/
#ifndef _MSC_VER 
#define QL_ASSERT(condition,message) \
if (!(condition)) { \
    std::ostringstream _ql_msg_stream; \
    _ql_msg_stream << message; \
    throw QuantLib::Error(__FILE__,__LINE__, \
                          __PRETTY_FUNCTION__,_ql_msg_stream.str()); \
 } else
#else
#define QL_ASSERT(condition,message) \
if (!(condition)) { \
    std::ostringstream _ql_msg_stream; \
    _ql_msg_stream << message; \
    throw QuantLib::Error(__FILE__,__LINE__, \
                          __FUNCSIG__,_ql_msg_stream.str()); \
 } else 
#endif



/*! \def QL_REQUIRE
    \brief throw an error if the given pre-condition is not verified
*/
#ifndef _MSC_VER
#define QL_REQUIRE(condition,message) \
if (!(condition)) { \
    std::ostringstream _ql_msg_stream; \
    _ql_msg_stream << message; \
    throw QuantLib::Error(__FILE__,__LINE__, \
                          __PRETTY_FUNCTION__,_ql_msg_stream.str()); \
 } else
#else
#define QL_REQUIRE(condition,message) \
if (!(condition)) { \
    std::ostringstream _ql_msg_stream; \
    _ql_msg_stream << message; \
    throw QuantLib::Error(__FILE__,__LINE__, \
                          __FUNCSIG__,_ql_msg_stream.str()); \
 } else
#endif


/*! \def QL_ENSURE
    \brief throw an error if the given post-condition is not verified
*/
#ifndef _MSC_VER
#define QL_ENSURE(condition,message) \
if (!(condition)) { \
    std::ostringstream _ql_msg_stream; \
    _ql_msg_stream << message; \
    throw QuantLib::Error(__FILE__,__LINE__, \
                          __PRETTY_FUNCTION__,_ql_msg_stream.str()); \
 } else
#else
#define QL_ENSURE(condition,message) \
if (!(condition)) { \
    std::ostringstream _ql_msg_stream; \
    _ql_msg_stream << message; \
    throw QuantLib::Error(__FILE__,__LINE__, \
                          __FUNCSIG__,_ql_msg_stream.str()); \
 } else
#endif

