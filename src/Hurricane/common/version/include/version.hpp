/*!
 * \file version.hpp
 * \brief Header of version
 *       The subroutines and functions are in the <i>version.cpp</i> file.
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

#include "OpenHurricaneConfig.hpp"
#include "preset.hpp"

/* OpenHurricane's dimension */
#define DIMENSIONSET 3

namespace OpenHurricane {
    class string;

    class versionData {
    public:
        int majorVersion;
        int minorVersion;
        int subminorVersion;

    public:
        versionData()
            : majorVersion(HURRICANE_VER_MAJOR), minorVersion(HURRICANE_VER_MINOR),
              subminorVersion(HURRICANE_VER_SUBMINOR) {}

        versionData(const int verMajor, const int verMinor, const int verSubminor)
            : majorVersion(verMajor), minorVersion(verMinor), subminorVersion(verSubminor) {}

        versionData(const string &version);

        inline ~versionData() noexcept {}

        hur_nodiscard bool operator==(const versionData &v) const;
        hur_nodiscard bool operator!=(const versionData &v) const;

        hur_nodiscard bool operator==(const string &v) const;
        hur_nodiscard bool operator!=(const string &v) const;
    };

    /*!\brief Program name and version.*/
    class programName {
    public:
        static const string name;
        static const versionData version;

        static hur_nodiscard string getVersion();

        static void printVersion();
    };
} // namespace OpenHurricane