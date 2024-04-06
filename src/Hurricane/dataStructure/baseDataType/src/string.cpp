/*!
 * \file string.cpp
 * \brief The subroutines of string
 * \author Rao Sihang
 * \version V2.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of Hurricane project (Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation)
 * \copyright Copyright (C) 2019-2024, Prof. Xu Xu's group at Beihang University.
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
#include "string.hpp"

#include "Lists.hpp"
#include <cstring>

namespace OpenHurricane {

    std::string::size_type findComment(std::string &str, const std::string cS,
                                       const std::string cM) {
        std::string::size_type posS, posM;
        posS = str.find_first_of(cS);
        posM = str.find(cM);
        if ((posS == std::string::npos) && (posM == std::string::npos)) {
            return std::string::npos;
        } else if (posS <= posM) {
            return posS;
        } else {
            return posM;
        }
    }

    void trim(std::string &str) {
        if (str.empty()) {
            return;
        }
        str.erase(0, str.find_first_not_of(" "));
        str.erase(str.find_last_not_of(" ") + 1);
    }

    hur_nodiscard std::string trimCopy(const std::string &str) {
        std::string s = str;
        trim(s);
        return s;
    }

    void eraseNullBlank(std::string &str) {
        auto num = std::strlen(str.c_str());
        integer count = Zero;

        for (size_t i = 0; i < num; ++i) {
            if (str[i] == '\n') {
                count++;
            }
        }

        count++;

        size_t lastPos = 0;
        size_t currentPos = 0;
        std::string temStr;
        for (integer i = 0; i < count; ++i) {
            currentPos = str.substr(lastPos).find_first_of('\n');
            currentPos += lastPos;

            std::string lineStr = trimCopy(str.substr(lastPos, currentPos - lastPos));
            if (!lineStr.empty()) {
                temStr += lineStr;
                temStr += '\n';
            }

            lastPos = ++currentPos;
        }

        str = temStr;
    }

    hur_nodiscard std::string getPrefix(const std::string &str) {
        stdStringList rll;
        split(str, rll, "_");
        return rll[0];
    }

    const string string::null;

} // namespace OpenHurricane