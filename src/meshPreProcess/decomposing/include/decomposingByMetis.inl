/*!
 * \file decomposingByMetis.inl
 * \brief In-Line subroutines of the <i>decomposingByMetis.hpp</i> file.
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

#include "decomposingByMetis.hpp"
#pragma once

inline const OpenHurricane::integerArray &
OpenHurricane::decomposingMeshByMetis::partOfProcs() const noexcept {
    return partOfProcs_;
}

inline OpenHurricane::integerArray &OpenHurricane::decomposingMeshByMetis::partOfProcs() {
    return partOfProcs_;
}

inline void OpenHurricane::decomposingMeshByMetis::clear() noexcept {
    adjncy_.clear();
    xadj_.clear();
    HurDelete(loadWgtPtr_);
}

inline OpenHurricane::decomposingMeshByMetis::decomposingMeshByMetis(
    const originMeshRead &originMeshes, const integer nP, const controller &cont)
    : originMesh_(originMeshes), origin_(), axis_(), partOfProcs_(), xadj_(), adjncy_(),
      nparts_(nP), cells_(originMeshes.cells()), faces_(originMeshes.faces()),
      nodes_(originMeshes.points()), cellZones_(originMeshes.cellZones()),
      faceZones_(originMeshes.faceZones()), pointZones_(originMeshes.pointZones()), cutZones_(),
      faceShareZones_(), perZones_(), periodicPairSize_(0), facePairMaps_(), pointPairMaps_(),
      perPairMaps_(), cellOffSet_(nP, OpenHurricane::Zero), faceOffSet_(nP, OpenHurricane::Zero),
      pointOffSet_(nP, OpenHurricane::Zero), color_(), unStructed_(false), hasStructed_(false),
      twoGhostLayer_(true), decomposeList_(),
      originMeshDecompSize_(originMeshes.originMeshDecompSize()), loadWgtPtr_(nullptr) {
    if (HurMPI::master()) {
        origin_ = originMeshes.origin();
        axis_ = originMeshes.axis();
    }
    HurMPI::bcast(&originMeshDecompSize_, 1, feature<integer>::MPIType);
    HurMPI::bcast(&origin_[0], 1, feature<real>::MPIType);
    HurMPI::bcast(&origin_[1], 1, feature<real>::MPIType);
    HurMPI::bcast(&origin_[2], 1, feature<real>::MPIType);

    HurMPI::bcast(&axis_[0], 1, feature<real>::MPIType);
    HurMPI::bcast(&axis_[1], 1, feature<real>::MPIType);
    HurMPI::bcast(&axis_[2], 1, feature<real>::MPIType);

    if (originMeshDecompSize_ > 1) {
        if (HurMPI::master()) {
            decomposeList_.resize(originMeshes.decomposeList().size());
            for (integer i = 0; i < decomposeList_.size(); ++i) {
                decomposeList_[i].resize(originMeshes.decomposeList()[i].size());
                decomposeList_[i] = originMeshes.decomposeList()[i];
            }
        }
    }

    if (HurMPI::master()) {
        if (cont.found("loadBalancing")) {
            const auto &lbCont = cont.subController("loadBalancing");

            loadWgtPtr_ = new loadBalancingWeights(originMesh_, lbCont);
        }
    }
}

inline OpenHurricane::decomposingMeshByMetis::~decomposingMeshByMetis() noexcept {
    clear();
}

inline const OpenHurricane::originMeshRead &
OpenHurricane::decomposingMeshByMetis::originMesh() const {
    return originMesh_;
}

inline OpenHurricane::integer OpenHurricane::decomposingMeshByMetis::nProcessors() const {
    return nparts_;
}

inline bool OpenHurricane::decomposingMeshByMetis::isUnStructed() const {
    return unStructed_;
}

inline bool OpenHurricane::decomposingMeshByMetis::hasStructed() const {
    return hasStructed_;
}

inline bool OpenHurricane::decomposingMeshByMetis::twoGhostLayer() const {
    return twoGhostLayer_;
}

inline const OpenHurricane::vector &OpenHurricane::decomposingMeshByMetis::origin() const noexcept {
    return origin_;
}

inline const OpenHurricane::vector &OpenHurricane::decomposingMeshByMetis::axis() const noexcept {
    return axis_;
}

inline OpenHurricane::cellList &OpenHurricane::decomposingMeshByMetis::cells() {
    return cells_;
}

inline OpenHurricane::faceList &OpenHurricane::decomposingMeshByMetis::faces() {
    return faces_;
}

inline OpenHurricane::pointField &OpenHurricane::decomposingMeshByMetis::points() {
    return nodes_;
}

inline OpenHurricane::cellZoneList &OpenHurricane::decomposingMeshByMetis::cellZones() {
    return cellZones_;
}

inline OpenHurricane::faceZoneList &OpenHurricane::decomposingMeshByMetis::faceZones() {
    return faceZones_;
}

inline OpenHurricane::pointZoneList &OpenHurricane::decomposingMeshByMetis::pointZones() {
    return pointZones_;
}

inline OpenHurricane::perZoneList &OpenHurricane::decomposingMeshByMetis::perZones() {
    return perZones_;
}

inline OpenHurricane::cutZoneList &OpenHurricane::decomposingMeshByMetis::cutZones() {
    return cutZones_;
}

inline OpenHurricane::integer
OpenHurricane::decomposingMeshByMetis::periodicPairSize() const noexcept {
    return periodicPairSize_;
}

inline const std::map<OpenHurricane::integer, OpenHurricane::integerArray> &
OpenHurricane::decomposingMeshByMetis::facePairMaps() const {
    return facePairMaps_;
}

inline std::map<OpenHurricane::integer, OpenHurricane::integerArray> &
OpenHurricane::decomposingMeshByMetis::facePairMaps() {
    return facePairMaps_;
}

inline std::map<OpenHurricane::integer, OpenHurricane::integerArray> &
OpenHurricane::decomposingMeshByMetis::perPairMaps() {
    return perPairMaps_;
}

inline const std::map<OpenHurricane::integer, OpenHurricane::integerVectorList> &
OpenHurricane::decomposingMeshByMetis::pointPairMaps() const {
    return pointPairMaps_;
}

inline const OpenHurricane::integerArray &
OpenHurricane::decomposingMeshByMetis::cellOffSet() const {
    return cellOffSet_;
}

inline const OpenHurricane::integerArray &
OpenHurricane::decomposingMeshByMetis::faceOffSet() const {
    return faceOffSet_;
}

inline const OpenHurricane::integerArray &
OpenHurricane::decomposingMeshByMetis::pointOffSet() const {
    return pointOffSet_;
}