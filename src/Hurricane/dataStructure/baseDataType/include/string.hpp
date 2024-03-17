/*!
 * \file string.hpp
 * \brief Header of string
 * \author Rao Sihang
 * \version V2.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of Hurricane project (Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation)
 * \copyright Copyright (C) 2019-2023, Prof. Xu Xu's group at Beihang University.
 *
 * License
 *		This file is part of OpenHurricane
 *
 *		OpenHurricane is free software: you can redistribute it and/or modify it
 *		under the terms of the GNU General Public License as published by
 *		the Free Software Foundation, either version 3 of the License, or
 *		(at your option) any later version.
 *
 *		OpenHurricane is distributed in the hope that it will be useful, but WITHOUT
 *		ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *		FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *		for more details.
 *
 *		You should have received a copy of the GNU General Public License
 *		along with OpenHurricane.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */
#pragma once
#include "chars.hpp"
#include <regex>

namespace OpenHurricane {
    template <class Cmpt> hur_nodiscard inline std::string toString(const Cmpt &) {
        return "Undefined toString()";
    }

    inline void stringToUpperCase(std::string &str) {
        std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    }

    inline std::string stringToUpperCase(const std::string &str) {
        std::string uppStr(str);
        std::transform(uppStr.begin(), uppStr.end(), uppStr.begin(), ::toupper);
        return uppStr;
    }

    inline void stringToLowerCase(std::string &str) {
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    }

    inline std::string stringToLowerCase(const std::string &str) {
        std::string lowStr(str);
        std::transform(lowStr.begin(), lowStr.end(), lowStr.begin(), ::tolower);
        return lowStr;
    }

    /*!\brief Find comment flag in a string line
        The default comment flags are:
            (1) single comment mark: "#" and "!"
        (2) multi comment mark: "//"
        user can change comment marks
        cS for single comment mark
        cM for multi comment mark*/
    std::string::size_type findComment(std::string &str, const std::string cS = "#!",
                                       const std::string cM = "//");

    inline void replaceAllMarks(std::string &str, const std::string &marks, const std::string &re) {
        size_t mSize = marks.size();
        size_t pos0 = 0;
        while (1) {
            size_t pos = str.find(marks, pos0);
            if (pos == std::string::npos) {
                return;
            }
            str.replace(pos, mSize, re);
            pos0 = pos + re.size();
        }
    }

    inline int countMarks(const std::string &str, const std::string &marks) {
        int count = 0;
        size_t index = 0;
        while (1) {
            index = str.find(marks, index);
            if (index == std::string::npos) {
                break;
            }
            count++;
            index++;
        }
        return count;
    }

    /**\brief Trim string.*/
    void trim(std::string &);

    /**\brief Trim and copy a string.*/
    hur_nodiscard std::string trimCopy(const std::string &);

    /**\brief Erase the blank line in a string.*/
    void eraseNullBlank(std::string &str);

    /**\brief get the part of string in front of "_".*/
    hur_nodiscard std::string getPrefix(const std::string &str);

    static const std::string nullString;

    /**
     * \brief The class of string.
     */
    class string : public std::string {
    public:
        using Base = std::string;

    public:
        static const string null;

        // Constructors

        inline string() : Base(){};
        inline string(const string &str) : Base(str) {}
        inline string &operator=(const string &str) {
            if (this != std::addressof(str)) {
                Base::operator=(str);
            }
            return *this;
        }

        inline string(string &&str) noexcept : Base() { Base::swap(str); }
        inline string &operator=(string &&str) noexcept {
            Base::swap(str);
            str.clear();
            return *this;
        }

        inline string(const char c) : Base(&c) {}

        inline string(const Base &str) : Base(str) {}
        inline string &operator=(const std::string &str) {
            Base::operator=(str);
            return *this;
        }

        inline string(Base &&str) noexcept : Base() { Base::swap(str); }
        inline string &operator=(Base &&str) noexcept {
            Base::swap(str);
            str.clear();
            return *this;
        }

        inline string(const char *c) : Base(c) {}

        inline string(const char *c, const size_type n) : Base(c, n) {}
        inline string &operator=(const char *c) {
            Base::operator=(c);
            return *this;
        }
    };

    template <> hur_nodiscard inline std::string toString(const string &w) {
        return static_cast<std::string>(w);
    }

} // namespace OpenHurricane