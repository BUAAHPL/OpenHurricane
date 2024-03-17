/*!
 * \file hdf5I.cpp
 * \brief Main subroutines of the <i>hdf5I.hpp</i> file.
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
#include "hdf5I.hpp"

void OpenHurricane::hdf5I::open(const unsigned int flg) {
    std::ifstream myin(filename_, std::ios::binary);
    if (!myin.good()) {
        LFatal("File: \"%s\" is not exist", filename_.c_str());
    }
    myin.close();
    if (flg != H5F_ACC_RDONLY) {
        LFatal("Must be read-only");
    }
    if (openOption_ == ONLY_MASTER) {
        if (HurMPIBase::master()) {
            hdf5IO::open(flg);
        }
    } else {
        hdf5IO::open(flg);
    }
}

void OpenHurricane::hdf5I::readString(std::string &str, const string &dataName) const {
    if (openOption_ == ONLY_MASTER) {
        if (HurMPIBase::master()) {
            hdf5IO::readString(str, dataName);
        }
    } else {
        hdf5IO::readString(str, dataName);
    }
}

void OpenHurricane::hdf5I::readString(std::string &str, const string &groupName,
                                  const string &dataName) const {
    if (openOption_ == ONLY_MASTER) {
        if (HurMPIBase::master()) {
            hdf5IO::readString(str, groupName, dataName);
        }
    } else {
        hdf5IO::readString(str, groupName, dataName);
    }
}

hur_nodiscard bool OpenHurricane::hdf5I::exist(const string &name) const {
    if (openOption_ == ONLY_MASTER && !HurMPIBase::master()) {
        errorAbortStr(("Called on process: " + toString(HurMPIBase::getProcRank()) +
                       ", while it is opened in \"ONLY_MASTER\""));
    }
    return hdf5IO::exist(name);
}

hur_nodiscard bool OpenHurricane::hdf5I::exist(const string &groupName, const string &name) const {
    if (openOption_ == ONLY_MASTER && !HurMPIBase::master()) {
        errorAbortStr(("Called on process: " + toString(HurMPIBase::getProcRank()) +
                       ", while it is opened in \"ONLY_MASTER\""));
    }
    return hdf5IO::exist(groupName, name);
}
