/*!
 * \file Lists.cpp
 * \brief The subroutines of lists
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

#include "Lists.hpp"
#include <cstring>
#include <iostream>

hur_nodiscard OpenHurricane::stdStringList OpenHurricane::stringToLineList(const std::string &str) {
    auto num = strlen(str.c_str());
    integer count = Zero;

    for (size_t i = 0; i < num; ++i) {
        if (str[i] =='\n') {
            count++;
        }
    }
    if (count == 0) {
        stdStringList stl(1, str);
        return stl;
    } else if (str[num - 1] !='\n') {
        count++;
    }
    stdStringList stl(count);

    size_t lastPos = 0;
    size_t currentPos = 0;

    for (integer i = 0; i < count; ++i) {
        currentPos = str.substr(lastPos).find_first_of('\n');
        currentPos += lastPos;
        stl[i] = str.substr(lastPos, currentPos - lastPos);

        // Replace "\r" with a space.
        replaceAllMarks(stl[i], "\r", " ");
        lastPos = ++currentPos;
    }

    return stl;
}

hur_nodiscard OpenHurricane::stringList OpenHurricane::stringToLineList(const string &str) {
    auto num = strlen(str.c_str());
    integer count = Zero;

    for (size_t i = 0; i < num; ++i) {
        if (str[i] =='\n') {
            count++;
        }
    }
    if (count == 0) {
        stringList stl(1, str);
        return stl;
    } else if (str[num - 1] !='\n') {
        count++;
    }
    stringList stl(count);

    size_t lastPos = 0;
    size_t currentPos = 0;

    for (integer i = 0; i < count; ++i) {
        currentPos = str.substr(lastPos).find_first_of('\n');
        currentPos += lastPos;
        stl[i] = str.substr(lastPos, currentPos - lastPos);

        // Replace "\r" with a space.
        replaceAllMarks(stl[i], "\r", " ");
        lastPos = ++currentPos;
    }

    return stl;
}

void OpenHurricane::split(const std::string &str, stdStringList &taken, const char *c) {
    std::string::size_type lastPos = str.find_first_not_of(c, 0);
    std::string::size_type pos = str.find_first_of(c, lastPos);
    while (pos != std::string::npos || lastPos != std::string::npos) {
        taken.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(c, pos);
        pos = str.find_first_of(c, lastPos);
    }
}

void OpenHurricane::split(const std::string &str, stdStringList &taken, const std::string &c) {
    split(str, taken, c.c_str());
}

void OpenHurricane::split(const std::string &str, stringList &taken, const char *c) {
    std::string::size_type lastPos = str.find_first_not_of(c, 0);
    std::string::size_type pos = str.find_first_of(c, lastPos);
    while (pos != std::string::npos || lastPos != std::string::npos) {
        taken.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(c, pos);
        pos = str.find_first_of(c, lastPos);
    }
}

void OpenHurricane::split(const std::string &str, stringList &taken, const std::string &c) {
    split(str, taken, c.c_str());
}
