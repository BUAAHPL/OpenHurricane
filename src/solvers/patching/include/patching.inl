#include "patching.hpp"
/*!
 * \file patching.inl
 * \brief The In-Line functions of the <i>patching.hpp</i> file.
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

inline const OpenHurricane::runtimeMesh &OpenHurricane::patching::mesh() const {
    return mesh_;
}

inline const OpenHurricane::controller &OpenHurricane::patching::cont() const {
    return cont_;
}

inline const std::string &OpenHurricane::patching::regionName() const noexcept {
    return regionName_;
}

inline const std::string &OpenHurricane::patching::typeName() const noexcept {
    return typeName_;
}

inline const OpenHurricane::stdStringList &OpenHurricane::patching::varName() const {
    return varName_;
}

inline OpenHurricane::integer OpenHurricane::patching::startStep() const noexcept {
    return startStep_;
}

inline OpenHurricane::integer OpenHurricane::patching::stayStep() const noexcept {
    return stayStep_;
}

inline void OpenHurricane::patching::patchRegion(const integer istep) const {
    if (typeName_ == "uniform") {
        uniformPatch(istep);
    } else if (typeName_ == "distributed") {
        distributePatch(istep);
    } else {
        LFatal("Invalid type: %s", typeName_.c_str());
    }
}

inline void OpenHurricane::patching::patchRegion() const {
    if (typeName_ == "uniform") {
        uniformPatch();
    } else if (typeName_ == "distributed") {
        distributePatch();
    } else {
        LFatal("Invalid type: %s", typeName_.c_str());
    }
}

inline void OpenHurricane::patching::uniformPatch(const integer istep) const {
    if (istep >= startStep_ && istep <= (startStep_ + stayStep_)) {
        uniformPatch();
    }
}