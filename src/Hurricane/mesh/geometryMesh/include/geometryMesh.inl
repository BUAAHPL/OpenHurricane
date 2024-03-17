/*!
 * \file geometryMesh.inl
 * \brief In-Line subroutines of the <i>geometryMesh.hpp</i> file.
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
#include "geometryMesh.hpp"

hur_nodiscard inline const OpenHurricane::geometryMesh &OpenHurricane::geometryMesh::nullObject() {
    return NullRefObj::nullRef<geometryMesh>();
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::geometryMesh::totalCutFaces() const {
    if (!totalCutFacesPtr_) {
        calcNumberOfCutFaces();
    }
    if (totalCutFacesPtr_) {
        return *totalCutFacesPtr_;
    }
    return 0;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::geometryMesh::totalPerFaces() const {
    if (!totalPerFacesPtr_) {
        calcNumberOfPerFaces();
    }
    if (totalPerFacesPtr_) {
        return *totalPerFacesPtr_;
    }
    return 0;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::geometryMesh::internalArraySize() const {
    return baseMesh::nCells();
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::geometryMesh::size() const noexcept {
    return baseMesh::nTotalCells();
}

hur_nodiscard inline bool OpenHurricane::geometryMesh::isInteriorFace(const integer i) const noexcept {
    return i < baseMesh::nInteriorFaces();
}

hur_nodiscard inline bool OpenHurricane::geometryMesh::isBoundaryFace(const integer i) const noexcept {
    return !isInteriorFace(i);
}

hur_nodiscard inline bool OpenHurricane::geometryMesh::isInternalCell(const integer i) const noexcept {
    return i < baseMesh::nCells();
}

hur_nodiscard inline bool OpenHurricane::geometryMesh::isDummyCell(const integer i) const noexcept {
    return !isInternalCell(i);
}

hur_nodiscard inline const OpenHurricane::integerArrayArray &
OpenHurricane::geometryMesh::tarOfProc() const noexcept {
    return *tarOfProc_;
}

hur_nodiscard inline const OpenHurricane::vectorArray &
OpenHurricane::geometryMesh::sorCellCentre() const noexcept {
    return *sorCellCentrePtr_;
}

hur_nodiscard inline const OpenHurricane::integerArrayArray &
OpenHurricane::geometryMesh::sorKNN() const noexcept {
    return *sorKNN_;
}

hur_nodiscard inline const OpenHurricane::integerArray &
OpenHurricane::geometryMesh::cellLoadWeights() const noexcept {
    if (!cellLoadWeightsPtr_) {
        cellLoadWeightsPtr_.reset(new integerArray(nCells(), Zero));
    }
    return *cellLoadWeightsPtr_;
}

hur_nodiscard inline OpenHurricane::integerArray &OpenHurricane::geometryMesh::cellLoadWeights() noexcept {
    if (!cellLoadWeightsPtr_) {
        cellLoadWeightsPtr_.reset(new integerArray(nCells(), Zero));
    }
    return *cellLoadWeightsPtr_;
}

hur_nodiscard inline bool OpenHurricane::geometryMesh::hasCellLoadWeights() const noexcept {
    return !(cellLoadWeightsPtr_.isNull());
}
