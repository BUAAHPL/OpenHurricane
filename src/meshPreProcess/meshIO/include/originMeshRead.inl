/*!
 * \file originMeshRead.inl
 * \brief In-Line subroutines of the <i>originMeshRead.hpp</i> file.
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

inline const OpenHurricane::vector &OpenHurricane::originMeshRead::origin() const noexcept {
    return origin_;
}

inline const OpenHurricane::vector &OpenHurricane::originMeshRead::axis() const noexcept {
    return axis_;
}

inline const OpenHurricane::fileName &OpenHurricane::originMeshRead::meshFileName() const {
    return fileName_;
}

inline int OpenHurricane::originMeshRead::numberOfProcessor() const {
    return nProcs_;
}

inline short OpenHurricane::originMeshRead::meshDimensionSet() const {
    return ND_;
}

inline OpenHurricane::integer OpenHurricane::originMeshRead::globalNodes() const {
    return globalNodes_;
}

inline OpenHurricane::integer OpenHurricane::originMeshRead::globalFaces() const {
    return globalFaces_;
}

inline OpenHurricane::integer OpenHurricane::originMeshRead::globalInterior() const {
    return globalInterior_;
}

inline OpenHurricane::integer OpenHurricane::originMeshRead::globalCells() const {
    return globalCells_;
}

inline const OpenHurricane::pointField &OpenHurricane::originMeshRead::points() const {
    return points_;
}

inline const OpenHurricane::faceList &OpenHurricane::originMeshRead::faces() const {
    return faces_;
}

inline const OpenHurricane::cellList &OpenHurricane::originMeshRead::cells() const {
    return cells_;
}

inline const OpenHurricane::pointZoneList &OpenHurricane::originMeshRead::pointZones() const {
    return pointZones_;
}

inline const OpenHurricane::faceZoneList &OpenHurricane::originMeshRead::faceZones() const {
    return faceZones_;
}

inline const OpenHurricane::cellZoneList &OpenHurricane::originMeshRead::cellZones() const {
    return cellZones_;
}

inline const OpenHurricane::periodicPairList &OpenHurricane::originMeshRead::periodicPairZones() const {
    return periodicPairZone_;
}

inline OpenHurricane::integerArrayArray &OpenHurricane::originMeshRead::originCellIndex() {
    return originCellIndex_;
}

inline bool OpenHurricane::originMeshRead::hasBeenRead() const {
    return hasBeenRead_;
}

inline bool OpenHurricane::originMeshRead::hasPeriodic() const {
    return hasPeriodic_;
}

inline bool OpenHurricane::originMeshRead::hasHangingNodes() const {
    return hasHangingNodes_;
}

inline bool OpenHurricane::originMeshRead::hasInterface() const {
    return hasInterface_;
}

inline std::map<OpenHurricane::integer, OpenHurricane::integer> &
OpenHurricane::originMeshRead::interiorWallFaceMap() {
    return interiorWallFaceMap_;
}

inline const std::map<OpenHurricane::integer, OpenHurricane::integer> &
OpenHurricane::originMeshRead::interiorWallFaceMap() const {
    return interiorWallFaceMap_;
}

hur_nodiscard inline bool OpenHurricane::originMeshRead::checkCoupledWall() const noexcept {
    return false;
}

hur_nodiscard inline bool OpenHurricane::originMeshRead::isHurricaneMesh() const noexcept {
    return false;
}

inline const OpenHurricane::integerArrayArray &OpenHurricane::originMeshRead::decomposeList() const {
    return decomposeList_;
}

inline OpenHurricane::integer OpenHurricane::originMeshRead::originMeshDecompSize() const noexcept {
    return originMeshDecompSize_;
}

hur_nodiscard inline const OpenHurricane::integerArray &
OpenHurricane::originMeshRead::cellLoadWeights() const noexcept {
    return cellLoadWeights_;
}
