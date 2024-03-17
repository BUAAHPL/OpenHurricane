#include "hdf5I.hpp"
/*!
 * \file hdf5I.inl
 * \brief In-Line subroutines of the <i>hdf5I.hpp</i> file.
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

inline OpenHurricane::hdf5I::hdf5I() : hdf5IO(), openOption_(ONLY_MASTER) {
    flag_ = H5F_ACC_RDONLY;
}

inline OpenHurricane::hdf5I::hdf5I(const fileName &fN) : hdf5IO(fN), openOption_(ONLY_MASTER) {
    flag_ = H5F_ACC_RDONLY;
}

inline OpenHurricane::hdf5I::~hdf5I() noexcept {}

inline void OpenHurricane::hdf5I::open(const fileName &fN) {
    filename_ = fN;
    open(flag_);
}

inline void OpenHurricane::hdf5I::open() {
    open(flag_);
}

inline bool OpenHurricane::hdf5I::onlyMaster() const noexcept {
    return openOption_ == ONLY_MASTER;
}

inline void OpenHurricane::hdf5I::readTypeFromAttribute(const H5::DataSet &dataset,
                                                    int &dataTypr) const {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::readTypeFromAttribute(dataset, dataTypr);
        }
    } else {
        hdf5IO::readTypeFromAttribute(dataset, dataTypr);
    }
}

template <class Type>
inline void OpenHurricane::hdf5I::readSingleCmpt(List<Type> &l, const string &dataName) const {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::readSingleCmpt<Type>(l, dataName);
        }
    } else {
        hdf5IO::readSingleCmpt<Type>(l, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5I::readSingleCmpt(List<Type> &l, const string &groupName,
                                             const string &dataName) const {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::readSingleCmpt<Type>(l, groupName, dataName);
        }
    } else {
        hdf5IO::readSingleCmpt<Type>(l, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5I::readMultipleCmpt(List<Type> &l, const string &dataName) const {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::readMultipleCmpt<Type>(l, dataName);
        }
    } else {
        hdf5IO::readMultipleCmpt<Type>(l, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5I::readMultipleCmpt(List<Type> &l, const string &groupName,
                                               const string &dataName) const {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::readMultipleCmpt<Type>(l, groupName, dataName);
        }
    } else {
        hdf5IO::readMultipleCmpt<Type>(l, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5I::read(List<Type> &l, const string &dataName) const {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::read<Type>(l, dataName);
        }
    } else {
        hdf5IO::read<Type>(l, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5I::read(List<Type> &l, const string &groupName,
                                   const string &dataName) const {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::read<Type>(l, groupName, dataName);
        }
    } else {
        hdf5IO::read<Type>(l, groupName, dataName);
    }
}

template <template <typename T> class Form, typename Type>
inline void OpenHurricane::hdf5I::readArrayArray(Form<Type> &l, const string &dataName) const {

    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::readArrayArray<Form, Type>(l, dataName);
        }
    } else {
        hdf5IO::readArrayArray<Form, Type>(l, dataName);
    }
}

template <template <typename T> class Form, typename Type>
inline void OpenHurricane::hdf5I::readArrayArray(Form<Type> &l, const string &groupName,
                                             const string &dataName) const {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::readArrayArray<Form, Type>(l, groupName, dataName);
        }
    } else {
        hdf5IO::readArrayArray<Form, Type>(l, groupName, dataName);
    }
}
