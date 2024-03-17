/*!
 * \file version.cpp
 * \brief The subroutines and functions of version.
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

#include "version.hpp"
#include "Lists.hpp"
#include "errorAbort.hpp"
#include "integer.hpp"
#include <iomanip>
#include <iostream>

OpenHurricane::versionData::versionData(const string &version)
    : majorVersion(HURRICANE_VER_MAJOR), minorVersion(HURRICANE_VER_MINOR),
      subminorVersion(HURRICANE_VER_SUBMINOR) {
    stringList taken;
    split(version, taken, ".");

    if (taken.size() == 3) {
        char *endptr;
        majorVersion = static_cast<int>(strtol(taken[0].c_str(), &endptr, 10));
        minorVersion = static_cast<int>(strtol(taken[1].c_str(), &endptr, 10));
        subminorVersion = static_cast<int>(strtol(taken[2].c_str(), &endptr, 10));
    } else {
        LFatal("Invalid version data");
    }
}

hur_nodiscard bool OpenHurricane::versionData::operator==(const versionData &v) const {
    return (majorVersion == v.majorVersion) && (minorVersion == v.minorVersion) &&
           (subminorVersion == v.subminorVersion);
}

hur_nodiscard bool OpenHurricane::versionData::operator!=(const versionData &v) const {
    return !operator==(v);
}

hur_nodiscard bool OpenHurricane::versionData::operator==(const string &v) const {
    return operator==(versionData(v));
}

hur_nodiscard bool OpenHurricane::versionData::operator!=(const string &v) const {
    return !operator==(v);
}

const OpenHurricane::string OpenHurricane::programName::name = "OpenHurricane";
const OpenHurricane::versionData OpenHurricane::programName::version;

hur_nodiscard OpenHurricane::string OpenHurricane::programName::getVersion() {
    string vv;
    vv = "V" + OpenHurricane::toString(integer(HURRICANE_VER_MAJOR)) + "." +
         OpenHurricane::toString(integer(HURRICANE_VER_MINOR)) + "." +
         OpenHurricane::toString(integer(HURRICANE_VER_SUBMINOR));
    return vv;
}

void OpenHurricane::programName::printVersion() {
    std::string butTime = "     Built on " + std::string(__TIME__) + " " + std::string(__DATE__);
    string name = programName::getVersion();

#ifdef HURRICANE_DP
    std::string pName = "     Built in double precision,";
#endif // HURRICANE_DP

#ifdef HURRICANE_SP
    std::string pName = "     Built in single precision";
#endif // HURRICANE_SP

#ifdef MPI_PARALLEL
    pName += " parallel";
#else
    pName += " serial";
#endif // MPI_PARALLEL

#ifdef HUR_DEBUG
    std::string verStr = "     CFD computation tools: " + programName::name + " debug, " + name;
#else
    std::string verStr = "     CFD computation tools: " + programName::name + " release, " + name;
#endif // HUR_DEBUG

    std::string welInfo = "  Welcome to use " + programName::name + " software";
    std::string welInfo2 = "     Highly Universal Rocket & Ramjet sImulation "
                           "Codes for ANalysis and Evaluation";

    std::cout << "    *" << std::setw(85) << std::setfill('*') << "*" << std::endl;
    std::cout << std::left << "    *" << std::setfill(' ') << std::setw(84) << welInfo.c_str()
              << "*" << std::endl;
    std::cout << std::left << "    *" << std::setfill(' ') << std::setw(84) << welInfo2.c_str()
              << "*" << std::endl;
    std::cout << "    *" << std::setw(84) << "     Copyright(c) 2019-2022 BUAA XX Research Group"
              << "*" << std::endl;
    std::cout << "    *" << std::setw(84) << butTime.c_str() << "*" << std::endl;
    std::cout << "    *" << std::setw(84) << pName.c_str() << "*" << std::endl;
    std::cout << "    *" << std::setw(84) << verStr.c_str() << "*" << std::endl;
    std::cout << "    *" << std::setw(85) << std::setfill('*') << "*" << std::endl;
    std::cout << std::right << std::setfill(' ');
}
