#include "globalMesh.hpp"
/*!
 * \file globalMesh.inl
 * \brief In-Line subroutines of the <i>globalMesh.hpp</i> file.
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

inline OpenHurricane::globalMesh::globalMesh(const geometryMesh &gMesh)
    : mesh_(gMesh), globalCellIndecesPtr_(nullptr), globalFaceIndecesPtr_(nullptr),
      globalPointIndecesPtr_(nullptr) {
    globalCellIndecesPtr_ = new globalCellIndex(gMesh);
}

hur_nodiscard inline const OpenHurricane::geometryMesh &OpenHurricane::globalMesh::mesh() const noexcept {
    return mesh_;
}

inline OpenHurricane::globalMesh::globalMesh(const geometryMesh &gMesh, OpenHurricane::globalCellIndex &gCI,
                                         OpenHurricane::globalFaceIndex &gFI,
                                         OpenHurricane::globalPointIndex &gPI)
    : mesh_(gMesh), globalCellIndecesPtr_(new globalCellIndex(gCI)),
      globalFaceIndecesPtr_(new globalFaceIndex(gFI)),
      globalPointIndecesPtr_(new globalPointIndex(gPI)) {
    globalFaceIndecesPtr_->updatePairMap(*globalCellIndecesPtr_);
}

inline OpenHurricane::globalMesh::globalMesh(const geometryMesh &gMesh,
                                         const decomposingMeshByMetis &dcM)
    : mesh_(gMesh), globalCellIndecesPtr_(new globalCellIndex(gMesh)),
      globalFaceIndecesPtr_(new globalFaceIndex(gMesh, dcM.facePairMaps(), dcM.faceOffSet())),
      globalPointIndecesPtr_(new globalPointIndex(gMesh, dcM.pointPairMaps(), dcM.pointOffSet())) {
    globalFaceIndecesPtr_->updatePairMap(*globalCellIndecesPtr_);
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalMesh::toGlobalCellIndex(const integer zoneid, const integer localIndex) const {
    return globalCellIndecesPtr_->toGlobalIndex(zoneid, localIndex);
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalMesh::toGlobalCellIndex(const integer localIndex) const {
    return globalCellIndecesPtr_->toGlobalIndex(localIndex);
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalMesh::toGlobalFaceIndex(const integer zoneId, const integer localId) const {
    return globalFaceIndecesPtr_->toGlobalIndex(zoneId, localId);
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalMesh::toGlobalPointIndex(const integer zoneId, const integer localId) const {
    return globalPointIndecesPtr_->toGlobalIndex(zoneId, localId);
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalMesh::toGlobalPointIndex(const integer localId) const {
    return globalPointIndecesPtr_->toGlobalIndex(localId);
}

hur_nodiscard inline const OpenHurricane::globalCellIndex &
OpenHurricane::globalMesh::globalCellIndeces() const noexcept {
    return *globalCellIndecesPtr_;
}

hur_nodiscard inline const OpenHurricane::globalFaceIndex &
OpenHurricane::globalMesh::globalFaceIndeces() const noexcept {
    return *globalFaceIndecesPtr_;
}

hur_nodiscard inline const OpenHurricane::globalPointIndex &
OpenHurricane::globalMesh::globalPointIndeces() const noexcept {
    return *globalPointIndecesPtr_;
}

inline void OpenHurricane::globalMesh::gatherFaceZone(const integer zonei, faceZone &gfz,
                                                  integerArrayArray &faceConnect,
                                                  integer &_offset) const {
    globalFaceIndecesPtr_->gatherFaceZone(zonei, gfz, faceConnect, *globalCellIndecesPtr_,
                                          *globalPointIndecesPtr_, mesh_.faces(), _offset);
}
