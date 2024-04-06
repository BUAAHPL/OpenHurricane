#include "baseMesh.hpp"
/*!
 * \file baseMesh.inl
 * \brief In-Line subroutines of the <i>baseMesh.hpp</i> file.
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

void OpenHurricane::baseMesh::setGhostCellsType(const short t) {
    std::map<short, Options::ghostCellLayers::ghostCellLayer>::iterator iter;
    iter = Options::ghostCellLayers::ghostCellLayerMap.find(t);
    if (iter == Options::ghostCellLayers::ghostCellLayerMap.end()) {
        std::string errMsg = "ghost cell layer type not found! type_ = ";
        errMsg = errMsg + std::to_string(t);
        errorAbortStr(errMsg);
    }
    ghostCellType_ = t;
}

inline void OpenHurricane::baseMesh::setOriginAndAxis(const point &_origin,
                                                  const point &_axis) noexcept {
    origin_ = _origin;
    axis_ = _axis;
}

hur_nodiscard inline const OpenHurricane::controller &OpenHurricane::baseMesh::cont() const noexcept {
    return cont_;
}

/*!\brief Return true if the ghost cells have been created.*/
hur_nodiscard inline bool OpenHurricane::baseMesh::isGhostCellCreated() const noexcept {
    return isGhostCellCreated_;
}

inline void OpenHurricane::baseMesh::setIsGhostCellCreated(const bool isSet) noexcept {
    isGhostCellCreated_ = isSet;
}

hur_nodiscard inline short OpenHurricane::baseMesh::ghostCellsType() const noexcept {
    return ghostCellType_;
}

/*!\brief Are all cells hexahedral in 3D or quadrilateral in 2D.*/
hur_nodiscard inline bool OpenHurricane::baseMesh::areAllHexOrQuadri() const noexcept {
    return areAllHexOrQuadri_;
}

hur_nodiscard inline short OpenHurricane::baseMesh::secondNeighbourCellsSize() const noexcept {
    return secondNeighbourCellsSize_;
}

inline void OpenHurricane::baseMesh::setSecondNeighbourCellsSize(const short s) {
    if (s < 0) {
        LFatal("Second neighbour cells size must large than 0");
    }

    if (s % 2 != 0) {
        LFatal("Second neighbour cells size must be even number");
    }

    secondNeighbourCellsSize_ = s;
}

