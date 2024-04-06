#include "hdf5O.hpp"
/*!
 * \file hdf5O.inl
 * \brief In-Line subroutines of the <i>hdf5O.hpp</i> file.
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
#pragma once

inline OpenHurricane::hdf5O::hdf5O() : hdf5IO(), openOption_(ONLY_MASTER) {
    flag_ = H5F_ACC_TRUNC;
}

inline OpenHurricane::hdf5O::hdf5O(const fileName &fN) : hdf5IO(fN), openOption_(ONLY_MASTER) {
    flag_ = H5F_ACC_TRUNC;
}

inline OpenHurricane::hdf5O::~hdf5O() noexcept {}

inline bool OpenHurricane::hdf5O::onlyMaster() const noexcept {
    return openOption_ == ONLY_MASTER;
}

inline void OpenHurricane::hdf5O::open(const fileName &fN) {
    filename_ = fN;
    open(flag_);
}

inline void OpenHurricane::hdf5O::open() {
    open(flag_);
}

template <class Type>
inline void OpenHurricane::hdf5O::writeSingleCmpt(const List<Type> &l, const Type factor,
                                              const integer size, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeSingleCmpt<Type>(l, factor, size, dataName);
        }
    } else {
        hdf5IO::writeSingleCmpt<Type>(l, factor, size, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeSingleCmpt(const List<Type> &l, const Type factor,
                                              const integer size, const string &groupName,
                                              const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeSingleCmpt<Type>(l, factor, size, groupName, dataName);
        }
    } else {
        hdf5IO::writeSingleCmpt<Type>(l, factor, size, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeSingleCmpt(const List<Type> &l, const integer size,
                                              const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeSingleCmpt<Type>(l, size, dataName);
        }
    } else {
        hdf5IO::writeSingleCmpt<Type>(l, size, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeSingleCmpt(const List<Type> &l, const integer size,
                                              const string &groupName, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeSingleCmpt<Type>(l, size, groupName, dataName);
        }
    } else {
        hdf5IO::writeSingleCmpt<Type>(l, size, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeSingleCmpt(const List<Type> &l, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeSingleCmpt<Type>(l, dataName);
        }
    } else {
        hdf5IO::writeSingleCmpt<Type>(l, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeSingleCmpt(const List<Type> &l, const string &groupName,
                                              const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeSingleCmpt<Type>(l, groupName, dataName);
        }
    } else {
        hdf5IO::writeSingleCmpt<Type>(l, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeMultipleCmpt(const List<Type> &l,
                                                const typename feature<Type>::elementType factor,
                                                const integer size, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeMultipleCmpt<Type>(l, factor, size, dataName);
        }
    } else {
        hdf5IO::writeMultipleCmpt<Type>(l, factor, size, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeMultipleCmpt(const List<Type> &l,
                                                const typename feature<Type>::elementType factor,
                                                const integer size, const string &groupName,
                                                const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeMultipleCmpt<Type>(l, factor, size, groupName, dataName);
        }
    } else {
        hdf5IO::writeMultipleCmpt<Type>(l, factor, size, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeMultipleCmpt(const List<Type> &l, const integer size,
                                                const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeMultipleCmpt<Type>(l, size, dataName);
        }
    } else {
        hdf5IO::writeMultipleCmpt<Type>(l, size, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeMultipleCmpt(const List<Type> &l, const integer size,
                                                const string &groupName, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeMultipleCmpt<Type>(l, size, groupName, dataName);
        }
    } else {
        hdf5IO::writeMultipleCmpt<Type>(l, size, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeMultipleCmpt(const List<Type> &l, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeMultipleCmpt<Type>(l, dataName);
        }
    } else {
        hdf5IO::writeMultipleCmpt<Type>(l, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::writeMultipleCmpt(const List<Type> &l, const string &groupName,
                                                const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeMultipleCmpt<Type>(l, groupName, dataName);
        }
    } else {
        hdf5IO::writeMultipleCmpt<Type>(l, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::write(const List<Type> &l,
                                    const typename feature<Type>::elementType factor,
                                    const integer size, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write<Type>(l, factor, size, dataName);
        }
    } else {
        hdf5IO::write<Type>(l, factor, size, dataName);
    }
}

template <class Type>
inline void
OpenHurricane::hdf5O::write(const List<Type> &l, const typename feature<Type>::elementType factor,
                        const integer size, const string &groupName, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write<Type>(l, factor, size, groupName, dataName);
        }
    } else {
        hdf5IO::write<Type>(l, factor, size, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::write(const List<Type> &l, const integer size,
                                    const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write<Type>(l, size, dataName);
        }
    } else {
        hdf5IO::write<Type>(l, size, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::write(const List<Type> &l, const integer size,
                                    const string &groupName, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write<Type>(l, size, groupName, dataName);
        }
    } else {
        hdf5IO::write<Type>(l, size, groupName, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::write(const List<Type> &l, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write<Type>(l, dataName);
        }
    } else {
        hdf5IO::write<Type>(l, dataName);
    }
}

template <class Type>
inline void OpenHurricane::hdf5O::write(const List<Type> &l, const string &groupName,
                                    const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::write<Type>(l, groupName, dataName);
        }
    } else {
        hdf5IO::write<Type>(l, groupName, dataName);
    }
}

template <template <typename Type> class Form, typename Type>
inline void OpenHurricane::hdf5O::writeArrayArray(const Form<Type> &l, const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeArrayArray<Form, Type>(l, dataName);
        }
    } else {
        hdf5IO::writeArrayArray<Form, Type>(l, dataName);
    }
}

template <template <typename Type> class Form, typename Type>
inline void OpenHurricane::hdf5O::writeArrayArray(const Form<Type> &l, const string &groupName,
                                              const string &dataName) {
    if (onlyMaster()) {
        if (HurMPIBase::master()) {
            hdf5IO::writeArrayArray<Form, Type>(l, groupName, dataName);
        }
    } else {
        hdf5IO::writeArrayArray<Form, Type>(l, groupName, dataName);
    }
}
