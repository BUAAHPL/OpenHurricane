#include "searchProcedures.hpp"
/*!
 * \file searchProcedures.inl
 * \brief In-Line subroutines of the <i>searchProcedures.hpp</i> file.
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
hur_nodiscard inline OpenHurricane::vectorArray &OpenHurricane::faceZonePackage::faceCentre() noexcept {
    return faceCentre_;
}

hur_nodiscard inline const OpenHurricane::vectorArray &
OpenHurricane::faceZonePackage::faceCentre() const noexcept {
    return faceCentre_;
}

hur_nodiscard inline OpenHurricane::vectorArray &OpenHurricane::faceZonePackage::faceArea() noexcept {
    return faceArea_;
}

hur_nodiscard inline const OpenHurricane::vectorArray &
OpenHurricane::faceZonePackage::faceArea() const noexcept {
    return faceArea_;
}

hur_nodiscard inline OpenHurricane::vectorArrayArray &OpenHurricane::faceZonePackage::fp() noexcept {
    return fp_;
}

hur_nodiscard inline const OpenHurricane::vectorArrayArray &
OpenHurricane::faceZonePackage::fp() const noexcept {
    return fp_;
}

inline void OpenHurricane::faceZonePackage::clear() noexcept {
    faceArea_.clear();
    faceCentre_.clear();
    fp_.clear();
}

inline bool OpenHurricane::searchProcedures::getDistance(cellRealArray &y) {
    return getDistance(y, const_cast<cellVectorArray &>(cellVectorArray::nullObject()));
}