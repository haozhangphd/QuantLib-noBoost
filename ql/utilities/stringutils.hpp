#ifndef quantlib_stringutils_hpp
#define quantlib_stringutils_hpp

#include <string>
#include <algorithm>

namespace {

    inline std::string to_upper_copy(const std::string &input) {
        std::string ret;
        std::for_each(input.cbegin(), input.cend(), [&ret](const char c) { ret.push_back(std::toupper(c)); });
        return ret;
    }

    inline std::string to_lower_copy(const std::string &input) {
        std::string ret;
        std::for_each(input.cbegin(), input.cend(), [&ret](const char c) { ret.push_back(std::tolower(c)); });
        return ret;
    }

}

#endif
