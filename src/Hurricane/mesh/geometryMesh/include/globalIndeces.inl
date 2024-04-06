/*!
 * \file globalIndeces.inl
 * \brief In-Line subroutines of the <i>globalIndeces.hpp</i> file.
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

inline void OpenHurricane::globalPointIndex::clear() noexcept {
    pointPairMap_.clear();
}

hur_nodiscard inline const OpenHurricane::geometryMesh &
OpenHurricane::globalPointIndex::mesh() const noexcept {
    return mesh_;
}

hur_nodiscard inline bool OpenHurricane::globalPointIndex::isSharedPoint(const integer i) const {
    pointPairMapType::const_iterator iter;
    iter = pointPairMap_.find(i);
    if (iter == pointPairMap_.end()) {
        return false;
    } else {
        return true;
    }
}

hur_nodiscard inline bool OpenHurricane::globalPointIndex::isSharedPoint(const integer i,
                                                                     bool &isMinId) const {
    pointPairMapType::const_iterator iter;
    iter = pointPairMap_.find(i);
    int pid = HurMPI::getProcRank();
    isMinId = true;
    if (iter == pointPairMap_.end()) {
        return false;
    } else {
        for (integer k = 0; k < iter->second.size(); ++k) {
            // If the process id is not the minimum of all processes that share the same point, then break.
            if (iter->second[k][0] < pid) {
                isMinId = false;
                break;
            }
        }
        return true;
    }
}

hur_nodiscard inline bool OpenHurricane::globalPointIndex::isSharedPoint(const integer i, bool &isMinId,
                                                                     int &minId) const {
    pointPairMapType::const_iterator iter;
    iter = pointPairMap_.find(i);
    minId = HurMPI::getProcRank();
    isMinId = true;
    if (iter == pointPairMap_.end()) {
        return false;
    } else {
        for (integer k = 0; k < iter->second.size(); ++k) {
            // If the process id is not the minimum of all processes that share the same point, then break.
            if (iter->second[k][0] < minId) {
                minId = iter->second[k][0];
                isMinId = false;
            }
        }
        return true;
    }
}

inline bool OpenHurricane::globalPointIndex::isSharedPoint(const integer i, bool &isMinId, int &minId,
                                                       integer &minLocalI) const {
    pointPairMapType::const_iterator iter;
    iter = pointPairMap_.find(i);
    minId = HurMPI::getProcRank();
    isMinId = true;
    if (iter == pointPairMap_.end()) {
        return false;
    } else {
        for (integer k = 0; k < iter->second.size(); ++k) {
            // If the process id is not the minimum of all processes that share the same point, then break.
            if (iter->second[k][0] < minId) {
                minId = iter->second[k][0];
                minLocalI = iter->second[k][2];
                isMinId = false;
            }
        }
        return true;
    }
}

hur_nodiscard inline bool
OpenHurricane::globalPointIndex::isSharedPoint(const integer i, integerVectorList &pairsList) const {
    pointPairMapType::const_iterator iter;
    iter = pointPairMap_.find(i);
    if (iter == pointPairMap_.end()) {
        return false;
    } else {
        pairsList.resize(iter->second.size());
        pairsList = iter->second;
        return true;
    }
}

hur_nodiscard inline bool
OpenHurricane::globalPointIndex::isSharedPoint(const integer i, bool &isMinId,
                                           integerVectorList &pairsList) const {
    pointPairMapType::const_iterator iter;
    iter = pointPairMap_.find(i);
    int pid = HurMPI::getProcRank();
    isMinId = true;
    if (iter == pointPairMap_.end()) {
        return false;
    } else {
        pairsList.resize(iter->second.size());
        pairsList = iter->second;
        for (integer k = 0; k < iter->second.size(); ++k) {
            // If the process id is not the minimum of all processes that share the same point, then break.
            if (iter->second[k][0] < pid) {
                isMinId = false;
                break;
            }
        }
        return true;
    }
}

hur_nodiscard inline bool
OpenHurricane::globalPointIndex::isSharedPoint(const integer i, bool &isMinId, int &minId,
                                           integerVectorList &pairsList) const {
    pointPairMapType::const_iterator iter;
    iter = pointPairMap_.find(i);
    minId = HurMPI::getProcRank();
    isMinId = true;
    if (iter == pointPairMap_.end()) {
        return false;
    } else {
        pairsList.resize(iter->second.size());
        pairsList = iter->second;
        for (integer k = 0; k < iter->second.size(); ++k) {
            // If the process id is not the minimum of all processes that share the same point, then break.
            if (iter->second[k][0] < minId) {
                minId = iter->second[k][0];
                isMinId = false;
            }
        }
        return true;
    }
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalPointIndex::toGlobalIndex(const integer zoneId, const integer localI) const {
    bool isMinId;
    int minId;
    integer minLocalI;
    isSharedPoint(localI, isMinId, minId, minLocalI);

    if (isMinId) {
        const auto counti = pointCountMap_.at(localI);

        return zoneOffset_[zoneId][minId] + processOffset_[zoneId][minId] + counti;

    } else {

        return zoneOffset_[zoneId][minId] + processOffset_[zoneId][minId] + minLocalI;
    }
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalPointIndex::toGlobalIndex(const integer localI) const {
    return toGlobalIndex(whichZone(localI), localI);
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalPointIndex::sharedPointSize() const noexcept {
    return static_cast<integer>(pointPairMap_.size());
}

inline void OpenHurricane::globalFaceIndex::clear() noexcept {
    facePairMap_.clear();
    zoneOffset_.clear();
    processOffset_.clear();
}

hur_nodiscard inline const OpenHurricane::geometryMesh &
OpenHurricane::globalFaceIndex::mesh() const noexcept {
    return mesh_;
}

/*!\brief Return true if the local face is shared with other process.*/
hur_nodiscard inline bool OpenHurricane::globalFaceIndex::isSharedFace(const integer localI) const {
    facePairMapType::const_iterator iter;
    iter = facePairMap_.find(localI);
    if (iter == facePairMap_.end()) {
        return false;
    } else {
        return true;
    }
}
/*!\brief Return true if the local face is shared with other process.*/
/*!\brief and return face pair info if true.*/
hur_nodiscard inline bool OpenHurricane::globalFaceIndex::isSharedFace(const integer localI,
                                                                   integerArray &pairs) const {
    facePairMapType::const_iterator iter;
    iter = facePairMap_.find(localI);
    if (iter == facePairMap_.end()) {
        return false;
    } else {
        pairs = iter->second;
        return true;
    }
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalFaceIndex::sharedFaceSize() const noexcept {
    return static_cast<integer>(facePairMap_.size());
}

hur_nodiscard inline const typename OpenHurricane::globalFaceIndex::facePairMapType &
OpenHurricane::globalFaceIndex::facePairMap() const noexcept {
    return facePairMap_;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::globalFaceZoneIndex::id() const noexcept {
    return id_;
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalFaceZoneIndex::totalNodes() const noexcept {
    return totalNodes_;
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::globalFaceZoneIndex::totalFaces() const noexcept {
    return totalFaces_;
}

hur_nodiscard inline const OpenHurricane::integerList &
OpenHurricane::globalFaceZoneIndex::writeNodeList() const noexcept {
    return writeNodeList_;
}

hur_nodiscard inline const std::unordered_map<OpenHurricane::integer, OpenHurricane::integer> &
OpenHurricane::globalFaceZoneIndex::antiOrderMap() const noexcept {
    return antiOrderMap_;
}

hur_nodiscard inline const OpenHurricane::pointField &
OpenHurricane::globalFaceZoneIndex::facePoints() const noexcept {
    return facePoints_;
}

inline void OpenHurricane::globalFaceZoneIndex::bcastFacePoints() const {
    integer ssz;
    if (HurMPIBase::master()) {
        ssz = facePoints_.size();
    }
    HurMPIBase::bcast(&ssz, 1, feature<integer>::MPIType, HurMPIBase::masterNo());
    if (!HurMPIBase::master()) {
        facePoints_.resize(ssz);
    }
    HurMPI::bcastVectorList(facePoints_, HurMPIBase::masterNo(), HurMPIBase::getComm());
}

hur_nodiscard inline const OpenHurricane::integerList &
OpenHurricane::globalFaceZoneIndex::faceDispls() const noexcept {
    return faceDispls_;
}

hur_nodiscard inline const OpenHurricane::integerList &
OpenHurricane::globalFaceZoneIndex::faceRecvcnt() const noexcept {
    return faceRecvcnt_;
}
