/*!
 * \file hdf5O.cpp
 * \brief Main subroutines of the <i>hdf5O.hpp</i> file.
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
#include "hdf5O.hpp"

void OpenHurricane::hdf5O::open(const unsigned int flg) {
    if (flg == H5F_ACC_RDONLY) {
        LFatal("Cannot be read-only in write data section");
    }
    if (openOption_ == ONLY_MASTER) {
        if (HurMPIBase::master()) {
            hdf5IO::open(flg);
        }
    } else if (openOption_ == CREATE_ALL) {
        auto pathOut = filename_.parentPath();
        auto fname = filename_.name(true);
        string fext;
        fext += filename_.ext();
        fname += "-";
        fname += toString(HurMPIBase::getProcRank());
        filename_ = fname + fext;
        filename_ = pathOut / filename_;

        hdf5IO::open(flg);
    }
}

void OpenHurricane::hdf5O::writeString(const std::string &str, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeString(str, dataName);
        }
    } else {
        hdf5IO::writeString(str, dataName);
    }
}

void OpenHurricane::hdf5O::writeString(const std::string &str, const string &groupName,
                                   const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeString(str, groupName, dataName);
        }
    } else {
        hdf5IO::writeString(str, groupName, dataName);
    }
}

template <> void OpenHurricane::hdf5O::write(const List<std::string> &str, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write(str, dataName);
        }
    } else {
        hdf5IO::write(str, dataName);
    }
}

template <>
void OpenHurricane::hdf5O::write(const List<std::string> &str, const string &groupName,
                             const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write(str, groupName, dataName);
        }
    } else {
        hdf5IO::write(str, groupName, dataName);
    }
}

template <> void OpenHurricane::hdf5O::write(const List<string> &str, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write(str, dataName);
        }
    } else {
        hdf5IO::write(str, dataName);
    }
}

template <>
void OpenHurricane::hdf5O::write(const List<string> &str, const string &groupName,
                             const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write(str, groupName, dataName);
        }
    } else {
        hdf5IO::write(str, groupName, dataName);
    }
}