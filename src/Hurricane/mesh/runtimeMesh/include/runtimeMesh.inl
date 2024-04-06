#include "runtimeMesh.hpp"
/*!
 * \file runtimeMesh.inl
 * \brief In-Line subroutines of the <i>runtimeMesh.hpp</i> file.
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

inline OpenHurricane::runtimeMesh::runtimeMesh(const object &ob, const controller &cont)
    : geometryMesh(ob, cont), regions_(), regionMap_() {
    checkAndRportMeshQuality();
    createRegions(cont);
}

inline OpenHurricane::runtimeMesh::runtimeMesh(object &&ob, const controller &cont) noexcept
    : geometryMesh(std::move(ob), cont), regions_(), regionMap_() {
    checkAndRportMeshQuality();
    createRegions(cont);
}

inline OpenHurricane::runtimeMesh::runtimeMesh(object &&ob, const controller &cac,
                                           const std::string &meshStr) noexcept
    : geometryMesh(std::move(ob), cac, meshStr), regions_(), regionMap_() {
    checkAndRportMeshQuality();
    createRegions(cac);
}

hur_nodiscard inline const OpenHurricane::realArray &OpenHurricane::runtimeMesh::cellVol() const {
    return cellVolume();
}

hur_nodiscard inline const OpenHurricane::vectorArray &OpenHurricane::runtimeMesh::cellCntr() const {
    return cellCentre();
}

hur_nodiscard inline const OpenHurricane::vectorArray &OpenHurricane::runtimeMesh::fA() const {
    return faceArea();
}

hur_nodiscard inline const OpenHurricane::vectorArray &OpenHurricane::runtimeMesh::faceCntr() const {
    return faceCentre();
}

hur_nodiscard inline const OpenHurricane::realArray &OpenHurricane::runtimeMesh::faceWgt() const {
    return faceWeight();
}

hur_nodiscard inline const OpenHurricane::sharedPtrList<OpenHurricane::markRegion> &
OpenHurricane::runtimeMesh::regions() const noexcept{
    return regions_;
}

hur_nodiscard inline const OpenHurricane::markRegion &
OpenHurricane::runtimeMesh::region(const std::string &regionName) const {
    auto iter = regionMap_.find(regionName);
    if (iter == regionMap_.end()) {
        LFatal("Cannot find region: %s", regionName.c_str());
    }
    return *(iter->second);
}
