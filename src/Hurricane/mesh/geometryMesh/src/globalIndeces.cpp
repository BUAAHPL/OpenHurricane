/*!
 * \file globalIndeces.cpp
 * \brief Main subroutines for global indeces.
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

#include "globalIndeces.hpp"
#include "geometryMesh.hpp"
#include "globalMesh.hpp"

OpenHurricane::globalPointIndex::globalPointIndex(
    const geometryMesh &gMesh, const std::map<integer, integerVectorList> &pointPairMaps)
    : mesh_(gMesh), pointPairMap_(pointPairMaps) {
    setLocalSize(mesh_.nPoints());
    setOffset();
}

OpenHurricane::globalPointIndex::globalPointIndex(
    const geometryMesh &gMesh, const std::map<integer, integerVectorList> &pointPairMaps,
    const integerList &pointOffSet)
    : mesh_(gMesh), pointPairMap_(pointPairMaps) {
    check(mesh_.nPoints());
    setOffset();
}

OpenHurricane::globalPointIndex::globalPointIndex(const geometryMesh &gMesh) : mesh_(gMesh) {
    setOffset();
}

void OpenHurricane::globalPointIndex::setLocalSize(const integer nLocalPoints) {
    integer localNonSharedPointSize = check(nLocalPoints);
}

void OpenHurricane::globalPointIndex::setPointPairMap(const std::map<integer, integerVectorList> &pM) {
    pointPairMap_ = pM;
    setLocalSize(mesh_.nPoints());
}

OpenHurricane::integer OpenHurricane::globalPointIndex::check(const integer nLocalPoints) const {
    if (!HurMPI::parRun()) {
        return nLocalPoints;
    }
    integer localNonSharedPointSize = 0;
    integer sharedPointSizeInMap = integer(pointPairMap_.size());
    integer sharedPointSize = 0;

    pointPairMapType::const_iterator iter;

    for (integer nodeI = 0; nodeI < nLocalPoints; nodeI++) {
        iter = pointPairMap_.find(nodeI);
        if (iter == pointPairMap_.end()) {
            localNonSharedPointSize++;
        } else {
            sharedPointSize++;
        }
    }

    if (sharedPointSize != sharedPointSizeInMap) {
        LFatal(
            "The size of the shared points: %d is not equal to the size of the point pair map: %d",
            sharedPointSize, sharedPointSizeInMap);
    }

    if ((localNonSharedPointSize + sharedPointSize) != mesh_.nPoints()) {
        LFatal("The sum of the shared points and non-shared points is not equal to the total "
               "local faces");
    }

    return localNonSharedPointSize;
}

void OpenHurricane::globalPointIndex::setOffset() {
    const auto &pointZones = mesh_.pointZones();

    zoneOffset_.resize(pointZones.size() + 1, Zero);

    processOffset_.resize(pointZones.size());
    pSize_.resize(pointZones.size());
    for (integer i = 0; i < pointZones.size(); ++i) {
        processOffset_[i].resize(HurMPI::getProcSize() + 1, Zero);
        pSize_[i].resize(HurMPI::getProcSize(), Zero);
    }
    for (integer i = 0; i <= pointZones.size(); ++i) {
        zoneOffset_[i].resize(HurMPI::getProcSize(), Zero);
    }

    for (integer i = 0; i < pointZones.size(); ++i) {
        integerList sendSizeL(HurMPI::getProcSize(), Zero);
        integerList recvSizeL(HurMPI::getProcSize(), Zero);
        for (integer fi = pointZones[i].firstIndex(); fi <= pointZones[i].lastIndex(); ++fi) {
            bool minId;
            integerVectorList pairList;
            integer minPid;
            bool isShared = isSharedPoint(fi, minId, minPid, pairList);
            if (!isShared) {
            } else if (isShared && minId) {
                for (integer k = 0; k < pairList.size(); ++k) {
                    if (pairList[k][0] > HurMPI::getProcRank()) {
                        sendSizeL[pairList[k][0]]++;
                    }
                }
            } else {
                auto &pplist = pointPairMap_.at(fi);
                for (integer k = 0; k < pplist.size(); k++) {
                    if (pplist[k][0] == minPid) {
                        recvSizeL[pplist[k][0]]++;
                    }
                }
            }
        }
        integerListList sendL(HurMPI::getProcSize());
        integerListList recvL(HurMPI::getProcSize());
        for (integer ii = 0; ii < sendL.size(); ++ii) {
            if (sendSizeL[ii] != 0) {
                sendL[ii].resize(2 * sendSizeL[ii]);
            }
            if (recvSizeL[ii] != 0) {
                recvL[ii].resize(2 * recvSizeL[ii]);
            }
        }
        integerList countIP(HurMPI::getProcSize(), Zero);
        integer fs = 0;
        for (integer fi = pointZones[i].firstIndex(); fi <= pointZones[i].lastIndex(); ++fi) {
            bool minId;
            integerVectorList pairList;
            integer minPid;
            bool isShared = isSharedPoint(fi, minId, minPid, pairList);
            if (!isShared) {
                pointCountMap_.emplace(fi, fs);
                fs++;
            } else if (isShared && minId) {
                pointCountMap_.emplace(fi, fs);
                for (integer k = 0; k < pairList.size(); ++k) {
                    if (pairList[k][0] > HurMPI::getProcRank()) {
                        // send the output index of the shared point of the minimum process to the remote process
                        // Note the data structure of pairsList[k] is (remote process id, shared point's local index of remote process)
                        // Namely: pairsList[k][0] = remote process id
                        //		   pairsList[k][1] = shared point's local index of remote process

                        sendL[pairList[k][0]][countIP[pairList[k][0]]++] = pairList[k][1];
                        sendL[pairList[k][0]][countIP[pairList[k][0]]++] = fs;
                    }
                }
                fs++;
            } else {
            }
        }
        HurMPI::Status status;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
            if (sendSizeL[ip] != 0) {
                HurMPI::sendList(sendL[ip], ip, HurMPI::getProcRank(), HurMPI::getComm());
            }

            if (recvSizeL[ip] != 0) {
                HurMPI::recvList(recvL[ip], ip, ip, HurMPI::getComm(), &status);
            }
        }
        for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
            if (recvSizeL[ip] != 0) {
                for (integer ii = 0; ii < recvSizeL[ip]; ++ii) {
                    const auto fi = recvL[ip][2 * ii];
                    const auto id = recvL[ip][2 * ii + 1];
                    auto &pplist = pointPairMap_.at(fi);
                    for (integer k = 0; k < pplist.size(); k++) {
                        if (pplist[k][0] == ip) {
                            pplist[k][2] = id;
                        }
                    }
                }
            }
        }
        pSize_[i][HurMPI::getProcRank()] = fs;
        HurMPI::allGatherList(pSize_[i], HurMPI::getComm());

        //zoneOffset_[i][HurMPI::getProcRank()] -= pointZones[i].firstIndex();
        HurMPI::allGatherList(zoneOffset_[i], HurMPI::getComm());

        integer offset = 0;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
            integer oldOffset = offset;
            offset += pSize_[i][ip];
            if (offset < oldOffset) {
                LFatal("Overflow: sumof sizes: %d exceeds capability of int (%d). Please recompile "
                       "with larger datatype for int.",
                       offset, integerMax);
            }
            processOffset_[i][ip + 1] = offset;
        }
        zoneOffset_[i + 1][HurMPI::getProcRank()] = offset;
    }

    HurMPI::allGatherList(zoneOffset_[pointZones.size()], HurMPI::getComm());
}

void OpenHurricane::globalPointIndex::getNodeList(const integer zoneId, pointField &pp) const {
    if (!(HurMPI::parRun())) {
        return;
    }

    const auto &f = mesh_.faces();
    const auto &p = mesh_.points();
    const auto &pointZones = mesh_.pointZones();
    /*integer nP = 0;
    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        nP += pSize_[zoneId][ip];
    }
    realArray ppsf;
    if (HurMPI::master()) {
        pp.resize(nP, Zero);
        ppsf.resize(nP * vector::nElements_, Zero);
    }
    realArray psf(pSize_[zoneId][HurMPI::getProcRank()] * vector::nElements_, Zero);

    integer fs = 0;

    for (integer pi = pointZones[zoneId].firstIndex(); pi <= pointZones[zoneId].lastIndex(); ++pi) {
        bool minId;
        bool isShared = isSharedPoint(pi, minId);
        if (!isShared || (isShared && minId)) {
            const integer xi = 3 * fs + 0;
            const integer yi = 3 * fs + 1;
            const integer zi = 3 * fs + 2;
            psf[xi] = p[pi].x();
            psf[yi] = p[pi].y();
            psf[zi] = p[pi].z();
            fs++;
        }
    }

    integerList nSizeL(HurMPI::getProcSize(), Zero);
    integerList displs;
    if (HurMPI::master()) {
        displs.resize(HurMPI::getProcSize(), Zero);
    }
    nSizeL[HurMPI::getProcRank()] = psf.size();
    HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());
    integer allSize = 0;
    if (HurMPI::master()) {
        for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
            displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
        }
        for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
            allSize += nSizeL[ip];
        }
        Pout << " all size: " << allSize << std::endl;
    }

    HurMPI::Request request;
    auto ierr = HurMPI::igatherv(psf.data(), psf.size(), feature<point>::MPIType, ppsf.data(),
                                 nSizeL.data(), displs.data(), feature<point>::MPIType,
                                 HurMPI::masterNo(), HurMPI::getComm(), &request);
    if (ierr != MPI_SUCCESS) {
        std::string errMsg;
        HurMPI::errorString(ierr, errMsg);
        LFatal("gather data failed: %s", errMsg.c_str());
    }

    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

    if (HurMPI::master()) {
        for (integer j = 0; j < ppsf.size() / 3; ++j) {
            const integer xi = 3 * j + 0;
            const integer yi = 3 * j + 1;
            const integer zi = 3 * j + 2;
            pp[j][0] = ppsf[xi];
            pp[j][1] = ppsf[yi];
            pp[j][2] = ppsf[zi];
        }
    }*/

    integer nP = 0;
    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        nP += pSize_[zoneId][ip];
    }
    realArray ppsf;
    if (HurMPI::master()) {
        pp.resize(nP, Zero);
        ppsf.resize(nP, Zero);
    }
    realArray psf(pSize_[zoneId][HurMPI::getProcRank()], Zero);

    for (int ie = 0; ie < point::nElements_; ++ie) {
        integer fs = 0;

        for (integer pi = pointZones[zoneId].firstIndex(); pi <= pointZones[zoneId].lastIndex();
             ++pi) {
            bool minId;
            bool isShared = isSharedPoint(pi, minId);
            if (!isShared || (isShared && minId)) {
                psf[fs] = p[pi][ie];
                fs++;
            }
        }

        integerList nSizeL(HurMPI::getProcSize(), Zero);
        integerList displs;
        if (HurMPI::master()) {
            displs.resize(HurMPI::getProcSize(), Zero);
        }
        nSizeL[HurMPI::getProcRank()] = psf.size();
        HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());
        integer allSize = 0;
        if (HurMPI::master()) {
            for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
            }
            for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                allSize += nSizeL[ip];
            }
        }

        HurMPI::Request request;
        HurMPI::igatherv(psf.data(), psf.size(), feature<feature<point>::elementType>::MPIType,
                         ppsf.data(), nSizeL.data(), displs.data(),
                         feature<feature<point>::elementType>::MPIType, HurMPI::masterNo(),
                         HurMPI::getComm(), &request);

        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::master()) {
            for (integer j = 0; j < ppsf.size(); ++j) {
                pp[j][ie] = ppsf[j];
            }
        }
    }
}

hur_nodiscard OpenHurricane::integer OpenHurricane::globalPointIndex::whichZone(const integer pi) const {
    for (integer i = 0; i < mesh_.pointZones().size(); ++i) {
        if (pi >= mesh_.pointZones()[i].firstIndex() && pi <= mesh_.pointZones()[i].lastIndex()) {
            return i;
        }
    }
    LFatal("Current point does not belong to any pointZone.");
    return -1;
}

void OpenHurricane::globalCellIndex::setOffset() {
    const auto &cellZones = mesh_.cellZones();

    zoneOffset_.resize(cellZones.size() + 1, Zero);

    processOffset_.resize(cellZones.size());

    for (integer i = 0; i < cellZones.size(); ++i) {
        processOffset_[i].resize(HurMPI::getProcSize() + 1, Zero);
    }
    for (integer i = 0; i <= cellZones.size(); ++i) {
        zoneOffset_[i].resize(HurMPI::getProcSize(), Zero);
    }
    integerList pSize(HurMPI::getProcSize(), Zero);
    for (integer i = 0; i < cellZones.size(); ++i) {
        pSize[HurMPI::getProcRank()] = cellZones[i].size();
        HurMPI::allGatherList(pSize, HurMPI::getComm());
        zoneOffset_[i][HurMPI::getProcRank()] -= cellZones[i].firstIndex();
        HurMPI::allGatherList(zoneOffset_[i], HurMPI::getComm());
        integer offset = 0;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
            integer oldOffset = offset;
            offset += pSize[ip];
            if (offset < oldOffset) {
                LFatal("Overflow: sumof sizes: %d exceeds capability of int (%d). Please recompile "
                       "with larger datatype for int.",
                       offset, integerMax);
            }
            processOffset_[i][ip + 1] = offset;
        }
        zoneOffset_[i + 1] = offset;
    }

    HurMPI::allGatherList(zoneOffset_[cellZones.size()], HurMPI::getComm());
}

void OpenHurricane::globalCellIndex::checkIndex(const integer zoneId, const integer index) const {
    const auto &cellZones = mesh_.cellZones();
    if (index < cellZones[zoneId].firstIndex() || index > cellZones[zoneId].lastIndex()) {
        LFatal("The cell index: %d is not in cell zone: %d", index, zoneId);
    }
}

OpenHurricane::globalCellIndex::globalCellIndex(const geometryMesh &gMesh) : mesh_(gMesh) {
    setOffset();
}

hur_nodiscard OpenHurricane::integer OpenHurricane::globalCellIndex::whichZone(const integer ci) const {
    for (integer i = 0; i < mesh_.cellZones().size(); ++i) {
        if (ci >= mesh_.cellZones()[i].firstIndex() && ci <= mesh_.cellZones()[i].lastIndex()) {
            return i;
        }
    }
    return -1;
}

OpenHurricane::globalFaceIndex::globalFaceIndex(const geometryMesh &gMesh,
                                            const std::map<integer, integerArray> &facePairMaps)
    : mesh_(gMesh), facePairMap_(facePairMaps) {
    getCutZonesIdArrays();
    setLocalSize(mesh_.nFaces());
    setOffset();
}

OpenHurricane::globalFaceIndex::globalFaceIndex(const geometryMesh &gMesh,
                                            const std::map<integer, integerArray> &facePairMaps,
                                            const integerArray &faceOffSet)
    : mesh_(gMesh), facePairMap_(facePairMaps) {
    getCutZonesIdArrays();
    check(mesh_.nFaces());
    setOffset();
}

OpenHurricane::globalFaceIndex::globalFaceIndex(const geometryMesh &gMesh) : mesh_(gMesh) {
    getCutZonesIdArrays();
    setOffset();
}

OpenHurricane::globalFaceIndex::~globalFaceIndex() noexcept {}

inline void OpenHurricane::globalFaceIndex::setLocalSize(const integer nLocalFace) {
    integer localNonSharedFaceSize = check(nLocalFace);
}

void OpenHurricane::globalFaceIndex::setfacePairMap(const std::map<integer, integerArray> &fM) {
    facePairMap_ = fM;
    setLocalSize(mesh_.nFaces());
}

void OpenHurricane::globalFaceIndex::updatePairMap(const globalCellIndex &gci) {
    integerArray sendSizeL(HurMPI::getProcSize(), Zero);
    integerArray recvSizeL(HurMPI::getProcSize(), Zero);
    for (auto iter = facePairMap_.begin(); iter != facePairMap_.end(); ++iter) {
        if (iter->second[0] < HurMPI::getProcRank()) {
            sendSizeL[iter->second[0]]++;
        } else {
            recvSizeL[iter->second[0]]++;
        }
    }
    integerArrayArray sendL(HurMPI::getProcSize());
    integerArrayArray recvL(HurMPI::getProcSize());
    for (integer i = 0; i < sendL.size(); ++i) {
        if (sendSizeL[i] != 0) {
            sendL[i].resize(2 * sendSizeL[i]);
        }
        if (recvSizeL[i] != 0) {
            recvL[i].resize(2 * recvSizeL[i]);
        }
    }

    integerArray countIP(HurMPI::getProcSize(), Zero);
    HurMPI::Status status;
    for (auto iter = facePairMap_.begin(); iter != facePairMap_.end(); ++iter) {
        if (iter->second[0] < HurMPI::getProcRank()) {
            const auto cl = mesh_.faces()[iter->first].leftCell();
            auto gcl = gci.toGlobalIndex(cl);
            sendL[iter->second[0]][countIP[iter->second[0]]++] = iter->second[2];
            sendL[iter->second[0]][countIP[iter->second[0]]++] = gcl;
        }
    }

    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (sendSizeL[ip] != 0) {
            HurMPI::sendList(sendL[ip], ip, HurMPI::getProcRank(), HurMPI::getComm());
        }

        if (recvSizeL[ip] != 0) {
            HurMPI::recvList(recvL[ip], ip, ip, HurMPI::getComm(), &status);
        }
    }

    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (recvSizeL[ip] != 0) {
            for (integer i = 0; i < recvSizeL[ip]; ++i) {
                const auto fi = recvL[ip][2 * i];
                const auto gcr = recvL[ip][2 * i + 1];
                auto iter = facePairMap_.find(fi);
                if (iter == facePairMap_.end()) {
                    LFatal("Cannot find: face index: %d in shared face pair map", fi);
                }
                facePairMap_.at(fi).append(gcr);
            }
        }
    }
}

OpenHurricane::integer OpenHurricane::globalFaceIndex::check(const integer nLocalFace) const {
    if (!HurMPI::parRun()) {
        return nLocalFace;
    }
    integer localNonSharedFaceSize = 0;
    integer sharedFaceSizeInMap = integer(facePairMap_.size());
    integer sharedFaceSize = 0;

    if (sharedFaceSizeInMap == 0) {
        LFatal("The face pair map is not specified. Please check!");
    }

    facePairMapType::const_iterator iter;

    for (integer faceI = 0; faceI < nLocalFace; faceI++) {
        iter = facePairMap_.find(faceI);
        if (iter == facePairMap_.end()) {
            localNonSharedFaceSize++;
        } else {
            sharedFaceSize++;
        }
    }

    if (sharedFaceSize != sharedFaceSizeInMap) {
        LFatal("The size of the shared faces: %d is not equal to the size of the face pair map: %d",
               sharedFaceSize, sharedFaceSizeInMap);
    }

    if ((localNonSharedFaceSize + sharedFaceSize) != nLocalFace) {
        LFatal("The sum of the shared faces and non-shared faces: %d is not equal to the total "
               "local faces: %d",
               localNonSharedFaceSize + sharedFaceSize, nLocalFace);
    }

    return localNonSharedFaceSize;
}

void OpenHurricane::globalFaceIndex::setOffset() {
    // Face zone
    const auto &faceZones = mesh_.faceZones();

    zoneOffset_.resize(faceZones.size() + 1, Zero);

    processOffset_.resize(faceZones.size());

    for (integer i = 0; i < faceZones.size(); ++i) {
        processOffset_[i].resize(HurMPI::getProcSize() + 1, Zero);
    }
    for (integer i = 0; i <= faceZones.size(); ++i) {
        zoneOffset_[i].resize(HurMPI::getProcSize(), Zero);
    }
    integerArray pSize(HurMPI::getProcSize(), Zero);

    // Offset for all face
    integer offset0 = 0;
    for (integer i = 0; i < faceZones.size(); ++i) {
        if (faceZones[i].isCutFace()) {
            continue;
        }
        integer fs = 0;
        for (integer fi = faceZones[i].firstIndex(); fi <= faceZones[i].lastIndex(); ++fi) {
            bool minId;
            bool isShared = isSharedFace(fi, minId);
            if (!isShared || (isShared && minId)) {
                fs++;
            }
        }
        if (cutId_[i].size() != 0) {
            for (integer j = 0; j < cutId_[i].size(); ++j) {
                integer fzi = cutId_[i][j];
                for (integer fi = faceZones[fzi].firstIndex(); fi <= faceZones[fzi].lastIndex();
                     ++fi) {
                    bool minId;
                    bool isShared = isSharedFace(fi, minId);
                    if (!isShared || (isShared && minId)) {
                        fs++;
                    }
                }
            }
        }
        pSize[HurMPI::getProcRank()] = fs;
        HurMPI::allGatherList(pSize, HurMPI::getComm());

        zoneOffset_[i][HurMPI::getProcRank()] -= faceZones[i].firstIndex();
        HurMPI::allGatherList(zoneOffset_[i], HurMPI::getComm());

        // Offset for all face of this face zone in different process
        integer offset = 0;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
            integer oldOffset = offset;
            offset += pSize[ip];
            if (offset < oldOffset) {
                std::string errMsg;

                LFatal("Overflow: sumof sizes: %d exceeds capability of int (%d). Please recompile "
                       "with larger datatype for int.",
                       offset, integerMax);
            }
            processOffset_[i][ip + 1] = offset;
        }

        zoneOffset_[i + 1] = offset0 + offset;
        offset0 += offset;
    }
    HurMPI::allGatherList(zoneOffset_[faceZones.size()], HurMPI::getComm());
}

void OpenHurricane::globalFaceIndex::getCutZonesIdArrays() {
    cutId_.resize(mesh_.faceZones().size());
    cutLinkId_.resize(mesh_.faceZones().size());
    cutLinkId_ = -1;
    for (integer fzid = 0; fzid < mesh_.faceZones().size(); ++fzid) {
        if (mesh_.faceZones()[fzid].isCutFace()) {
            for (integer i = 0; i < mesh_.faceZones().size(); ++i) {
                if (mesh_.faceZones()[i].index() == -mesh_.faceZones()[fzid].index()) {
                    cutId_[i].append(fzid);
                    cutLinkId_[fzid] = i;
                }
            }
        }
    }
}

hur_nodiscard bool OpenHurricane::globalFaceIndex::isSharedFace(const integer i, bool &minId) const {
    facePairMapType::const_iterator iter;
    iter = facePairMap_.find(i);
    bool isShared = false;
    int pid = HurMPI::getProcRank();
    minId = true;
    if (iter == facePairMap_.end()) {
        isShared = false;
    } else {
        // If the process id is not the minimum of all processes that share the same point, then break.
        if (iter->second[0] < pid) {
            minId = false;
        }
        isShared = true;
    }
    return isShared;
}

bool OpenHurricane::globalFaceIndex::isSharedFace(const integer i, bool &isMinId, int &minId) const {
    facePairMapType::const_iterator iter;
    iter = facePairMap_.find(i);
    bool isShared = false;
    minId = HurMPI::getProcRank();
    isMinId = true;
    if (iter == facePairMap_.end()) {
        isShared = false;
    } else {
        // If the process id is not the minimum of all processes that share the same point, then break.
        if (iter->second[0] < minId) {
            isMinId = false;
            minId = iter->second[0];
        }
        isShared = true;
    }
    return isShared;
}

hur_nodiscard OpenHurricane::integer OpenHurricane::globalFaceIndex::toGlobalIndex(const integer zoneId,
                                                                           const integer i) const {
    bool isMinId;
    int minId;
    isSharedFace(i, isMinId, minId);
    integer zid = zoneId;
    if (mesh_.faceZones()[zoneId].isCutFace()) {
        zid = cutLinkId_[zoneId];
        return zoneOffset_[zid][minId] + processOffset_[zid][minId] +
               mesh_.faceZones()[zid].lastIndex() + 1 + i - mesh_.faceZones()[zoneId].firstIndex();
    }
    return zoneOffset_[zid][minId] + processOffset_[zid][minId] + i;
}

hur_nodiscard OpenHurricane::integer
OpenHurricane::globalFaceIndex::toGlobalIndex(const integer zoneId, const integer i,
                                          const integer processId) const {
    bool isMinId;
    int minId;
    integer zid = zoneId;
    if (isSharedFace(i, isMinId, minId)) {
        if (mesh_.faceZones()[zoneId].isCutFace()) {
            zid = cutLinkId_[zoneId];
        }
        return zoneOffset_[zid][minId] + processOffset_[zid][minId] +
               mesh_.faceZones()[zid].lastIndex() + 1 + i - mesh_.faceZones()[zoneId].firstIndex();
    }
    return zoneOffset_[zid][processId] + processOffset_[zid][processId] + i;
}

void OpenHurricane::globalFaceIndex::gatherFaceZone(const integer zonei, faceZone &gfz,
                                                integerArrayArray &faceConnect,
                                                const globalCellIndex &gci,
                                                const globalPointIndex &gpi, const faceList &fs,
                                                integer &_offset) const {

    if (mesh_.faceZones()[zonei].isCutFace()) {
        return;
    }

    integer zsize = mesh_.faceZones()[zonei].size();
    if (HurMPI::master()) {
        gfz.resetName(mesh_.faceZones()[zonei].name());
        gfz.setIndex(mesh_.faceZones()[zonei].index());
        gfz.setBcType(mesh_.faceZones()[zonei].bcType());
        gfz.setFaceType(mesh_.faceZones()[zonei].faceType());
    }
    if (cutId_[zonei].size() != 0) {
        integer fs = 0;
        for (integer j = 0; j < cutId_[zonei].size(); ++j) {
            integer fzi = cutId_[zonei][j];
            for (integer fi = mesh_.faceZones()[fzi].firstIndex();
                 fi <= mesh_.faceZones()[fzi].lastIndex(); ++fi) {
                bool minId;
                bool isShared = isSharedFace(fi, minId);
                if (!isShared || (isShared && minId)) {
                    fs++;
                }
            }
        }
        zsize += fs;
    }
    integerArrayArray tmpFaceConnect(zsize);
    HurMPI::reduce(zsize, MPI_SUM);

    if (HurMPI::master()) {
        gfz.setFirstIndex(_offset);
        gfz.setLastIndex(_offset + zsize - 1);
        _offset += zsize;
    }

    bool isMix = false;
    integer offset = 0;
    if (mesh_.faceZones()[zonei].isMixed()) {
        isMix = true;
        offset = 1;
    }

    integer count = 0;
    integer maxSize = 0;
    for (integer i = mesh_.faceZones()[zonei].firstIndex();
         i <= mesh_.faceZones()[zonei].lastIndex(); ++i) {
        tmpFaceConnect[count].resize(15);
        if (isMix) {
            tmpFaceConnect[count][0] = fs[i].size();
            maxSize = max(maxSize, fs[i].size() + 1 + 2);
        } else {
            maxSize = max(maxSize, fs[i].size() + 2);
        }

        if (maxSize > 15) {
            LFatal("Exceed the maximum limits");
        }

        for (integer j = 0; j < fs[i].size(); ++j) {
            tmpFaceConnect[count][j + offset] = gpi.toGlobalIndex(fs[i][j]);
        }
        tmpFaceConnect[count][fs[i].size() + offset] = gci.toGlobalIndex(fs[i].leftCell());

        if (fs[i].rightCell() >= mesh_.nCells()) {
            tmpFaceConnect[count][fs[i].size() + offset + 1] = -1;
        } else {
            tmpFaceConnect[count][fs[i].size() + offset + 1] = gci.toGlobalIndex(fs[i].rightCell());
        }
        count++;
    }

    if (cutId_[zonei].size() != 0) {
        for (integer j = 0; j < cutId_[zonei].size(); ++j) {
            integer fzi = cutId_[zonei][j];
            for (integer fi = mesh_.faceZones()[fzi].firstIndex();
                 fi <= mesh_.faceZones()[fzi].lastIndex(); ++fi) {

                bool minId = true;
                const auto pairL = facePairMap_.at(fi);

                if (pairL[0] < HurMPI::getProcRank()) {
                    minId = false;
                }
                if (minId) {
                    tmpFaceConnect[count].resize(15);
                    if (isMix) {
                        tmpFaceConnect[count][0] = fs[fi].size();
                        maxSize = max(maxSize, fs[fi].size() + 1 + 2);
                    } else {
                        maxSize = max(maxSize, fs[fi].size() + 2);
                    }
                    if (maxSize > 15) {
                        LFatal("Exceed the maximum limits");
                    }
                    for (integer k = 0; k < fs[fi].size(); ++k) {
                        tmpFaceConnect[count][k + offset] = gpi.toGlobalIndex(fs[fi][k]);
                    }
                    tmpFaceConnect[count][fs[fi].size() + offset] =
                        gci.toGlobalIndex(fs[fi].leftCell());
                    tmpFaceConnect[count][fs[fi].size() + offset + 1] = pairL[4];

                    count++;
                }
            }
        }
    }
    HurMPI::allReduce(maxSize, MPI_MAX);

    integerArray nSizeL(HurMPI::getProcSize(), Zero);
    integerArray displs;
    if (HurMPI::master()) {
        displs.resize(HurMPI::getProcSize(), Zero);
    }
    nSizeL[HurMPI::getProcRank()] = tmpFaceConnect.size();
    HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());
    integerArray fC(tmpFaceConnect.size());
    integerArray allFC;
    if (HurMPI::master()) {
        faceConnect.resize(zsize);
        for (integer i = 0; i < zsize; i++) {
            faceConnect[i].resize(maxSize);
        }
        allFC.resize(zsize);
        for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
            displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
        }
    }
    for (integer i = 0; i < maxSize; ++i) {
        for (integer j = 0; j < tmpFaceConnect.size(); ++j) {
            fC[j] = tmpFaceConnect[j][i];
        }
        HurMPI::Request request;
        HurMPI::igatherv(fC.data(), fC.size(), feature<integer>::MPIType, allFC.data(),
                         nSizeL.data(), displs.data(), feature<integer>::MPIType,
                         HurMPI::masterNo(), HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::master()) {
            for (integer j = 0; j < allFC.size(); ++j) {

                faceConnect[j][i] = allFC[j];
            }
        }
    }
}

void OpenHurricane::globalFaceZoneIndex::getNodeList() {
    // If it is not parallel task.
    if (!(HurMPI::parRun())) {
        getNodeListForSerialRun();
        return;
    }
    // The list of faces in this process.
    const auto &f = mesh_.faces();

    // The list of points in this process.
    const auto &p = mesh_.points();

    const auto &gm = mesh_.globalMeshInfo();

    const auto &gp = gm.globalPointIndeces();

    // To get the number of faces of this face zone.
    totalFaces_ = fZ().size();
    HurMPI::allReduce(totalFaces_, MPI_SUM);

    integerListList sendIPIP(HurMPI::getProcSize());
    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        sendIPIP[ip].resize(HurMPI::getProcSize(), Zero);
    }

    // flag = -1 -- point is not set yet
    // flag = 0  -- point is set and is not shared point
    // flag = 1  -- point is set and is shared point
    integerList flagIp(p.size(), -1);
    for (integer fi = fZ().firstIndex(); fi <= fZ().lastIndex(); ++fi) {
        for (integer pi = 0; pi < f[fi].size(); pi++) {
            // The index of pi_th point of face fi
            const auto j = f[fi][pi];
            if (flagIp[j] == -1) {
                integerVectorList pairsList;
                // To check if this point is shared by other processes.
                // And make sure that it is only counted
                // by the process having the min procRank among these shared processes.
                if (gp.isSharedPoint(j, pairsList)) {
                    auto pid = HurMPI::getProcRank();
                    for (integer k = 0; k < pairsList.size(); ++k) {
                        // Note the data structure of pairsList[k] is (remote process id, shared point's local index of remote process)
                        // Namely: pairsList[k][0] = remote process id (procRank)
                        //		   pairsList[k][1] = shared point's local index of remote process
                        if (pairsList[k][0] != pid) {
                            sendIPIP[pairsList[k][0]][pid]++;
                        }
                    }
                    flagIp[j] = 1;
                } else {
                    flagIp[j] = 0;
                }
            }
        }
    }

    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        HurMPI::allGatherList(sendIPIP[ip]);
    }
    integerListList sendIPIPList(HurMPI::getProcSize());
    integerListList recvIPIPList(HurMPI::getProcSize());
    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (!HurMPI::isThisProc(ip)) {
            if (sendIPIP[ip][HurMPI::getProcRank()] != 0) {
                sendIPIPList[ip].resize(2 * sendIPIP[ip][HurMPI::getProcRank()], -1);
            }
            if (sendIPIP[HurMPI::getProcRank()][ip] != 0) {
                recvIPIPList[ip].resize(2 * sendIPIP[HurMPI::getProcRank()][ip], -1);
            }
        }
    }
    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        sendIPIP[ip] = Zero;
    }
    flagIp = -1;
    for (integer fi = fZ().firstIndex(); fi <= fZ().lastIndex(); ++fi) {
        for (integer pi = 0; pi < f[fi].size(); pi++) {
            // The index of pi_th point of face fi
            const auto j = f[fi][pi];
            if (flagIp[j] == -1) {
                integerVectorList pairsList;
                // To check if this point is shared by other processes.
                // And make sure that it is only counted
                // by the process having the min procRank among these shared processes.
                if (gp.isSharedPoint(j, pairsList)) {
                    auto pid = HurMPI::getProcRank();
                    for (integer k = 0; k < pairsList.size(); ++k) {
                        // Note the data structure of pairsList[k] is (remote process id, shared point's local index of remote process)
                        // Namely: pairsList[k][0] = remote process id (procRank)
                        //		   pairsList[k][1] = shared point's local index of remote process
                        if (pairsList[k][0] != pid) {
                            auto &iid = sendIPIP[pairsList[k][0]][pid];
                            sendIPIPList[pairsList[k][0]][2 * iid] = pairsList[k][1];
                            sendIPIPList[pairsList[k][0]][2 * iid + 1] = j;
                            iid++;
                        }
                    }
                    flagIp[j] = 1;
                } else {
                    flagIp[j] = 0;
                }
            }
        }
    }

    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (sendIPIPList[ip].size() != 0) {
            HurMPI::sendList(sendIPIPList[ip], ip, HurMPI::getProcRank(), HurMPI::getComm());
        }

        if (recvIPIPList[ip].size() != 0) {
            HurMPI::Status status;
            HurMPI::recvList(recvIPIPList[ip], ip, ip, HurMPI::getComm(), &status);
        }
    }

    List<std::map<integer, integer>> faceZoneSharedPMap(HurMPI::getProcSize());
    for (integer ipp = 0; ipp < HurMPI::getProcSize(); ++ipp) {
        if (recvIPIPList[ipp].size() != 0) {
            for (integer ii = 0; ii < recvIPIPList[ipp].size() / 2; ++ii) {
                faceZoneSharedPMap[ipp].emplace(recvIPIPList[ipp][2 * ii],
                                                recvIPIPList[ipp][2 * ii + 1]);
            }
        }
    }

    // flag = -1 -- point is not set yet
    // flag = 0  -- point is set and is not shared point
    // flag = 1  -- point is set and is shared point and the process id
    //              of the current process is the minimum id in the shared list.
    // flag = 2  -- point is set and is shared point and the process id
    //              of the current process is not the minimum id in the shared list.
    integerList flag(p.size(), -1);
    // To get the numbe of points of this face zone.
    integer psize = 0;
    for (integer fi = fZ().firstIndex(); fi <= fZ().lastIndex(); ++fi) {
        for (integer pi = 0; pi < f[fi].size(); pi++) {
            // The index of pi_th point of face fi
            const auto j = f[fi][pi];
            if (flag[j] == -1) {
                integerVectorList pairsList;
                // To check if this point is shared by other processes.
                // And make sure that it is only counted
                // by the process having the min procRank among these shared processes.
                if (gp.isSharedPoint(j, pairsList)) {
                    auto pid = HurMPI::getProcRank();
                    bool minId = true;
                    for (integer k = 0; k < pairsList.size(); ++k) {
                        // Note the data structure of pairsList[k] is (remote process id, shared point's local index of remote process)
                        // Namely: pairsList[k][0] = remote process id (procRank)
                        //		   pairsList[k][1] = shared point's local index of remote process
                        if (pairsList[k][0] < pid) {
                            auto iter = faceZoneSharedPMap[pairsList[k][0]].find(j);
                            if (iter != faceZoneSharedPMap[pairsList[k][0]].end()) {
                                minId = false;
                                break;
                            }
                        }
                    }
                    if (minId) {
                        flag[j] = 1;
                        psize++;
                    } else {
                        flag[j] = 2;
                    }
                } else {
                    psize++;
                    flag[j] = 0;
                }
            }
        }
    }

    integerList offSet(HurMPI::getProcSize());
    offSet[HurMPI::getProcRank()] = psize;
    HurMPI::allGatherList(offSet, HurMPI::getComm());

    totalNodes_ = sumArray(offSet);
    writeNodeList_.resize(psize);
    integer tsize = 0;
    for (integer ip = 0; ip < HurMPI::getProcRank(); ++ip) {
        tsize += offSet[ip];
    }

    flag = -1;
    psize = 0;

    integerList sendSizeL(HurMPI::getProcSize(), Zero);
    integerList recvSizeL(HurMPI::getProcSize(), Zero);
    for (integer fi = fZ().firstIndex(); fi <= fZ().lastIndex(); ++fi) {
        for (integer pi = 0; pi < f[fi].size(); pi++) {
            const auto j = f[fi][pi];
            if (flag[j] == -1) {
                integerVectorList pairsList;
                if (gp.isSharedPoint(j, pairsList)) // For shared points
                {
                    // The index of the current process
                    auto pid = HurMPI::getProcRank();
                    bool minId = true;
                    for (integer k = 0; k < pairsList.size(); ++k) {
                        // If the process id is not the minimum of all processes that share the same point, then break.
                        if (pairsList[k][0] < pid) {
                            auto iter = faceZoneSharedPMap[pairsList[k][0]].find(j);
                            if (iter != faceZoneSharedPMap[pairsList[k][0]].end()) {
                                minId = false;
                                break;
                            }
                        }
                    }
                    if (minId) {
                        flag[j] = 1;
                        for (integer k = 0; k < pairsList.size(); ++k) {
                            if (pairsList[k][0] > pid) {
                                auto iter = faceZoneSharedPMap[pairsList[k][0]].find(j);
                                if (iter != faceZoneSharedPMap[pairsList[k][0]].end()) {
                                    sendSizeL[pairsList[k][0]]++;
                                }
                            }
                        }
                    } else {
                        flag[j] = 2;
                        integer minIp = HurMPI::getProcSize();
                        for (integer k = 0; k < pairsList.size(); ++k) {
                            auto iter = faceZoneSharedPMap[pairsList[k][0]].find(j);
                            if (iter != faceZoneSharedPMap[pairsList[k][0]].end()) {
                                minIp = min(minIp, pairsList[k][0]);
                            }
                        }
                        recvSizeL[minIp]++;
                    }
                } else // For not-shared points
                {
                    flag[j] = 0;
                }
            }
        }
    }

    integerListList sendL(HurMPI::getProcSize());
    integerListList recvL(HurMPI::getProcSize());
    for (integer ii = 0; ii < sendL.size(); ++ii) {
        if (sendSizeL[ii] != 0) {
            sendL[ii].resize(2 * sendSizeL[ii]);
        }
        if (recvSizeL[ii] != 0) {
            recvL[ii].resize(2 * recvSizeL[ii]);
        }
    }
    integerList countIP(HurMPI::getProcSize(), Zero);
    HurMPI::Status status;
    flag = -1;
    psize = 0;
    for (integer fi = fZ().firstIndex(); fi <= fZ().lastIndex(); ++fi) {
        for (integer pi = 0; pi < f[fi].size(); pi++) {
            const auto j = f[fi][pi];
            if (flag[j] == -1) {
                integerVectorList pairsList;
                if (gp.isSharedPoint(j, pairsList)) // For shared points
                {
                    // The index of the current process
                    auto pid = HurMPI::getProcRank();
                    bool minId = true;
                    for (integer k = 0; k < pairsList.size(); ++k) {
                        // If the process id is not the minimum of all processes that share the same point, then break.
                        if (pairsList[k][0] < pid) {
                            auto iter = faceZoneSharedPMap[pairsList[k][0]].find(j);
                            if (iter != faceZoneSharedPMap[pairsList[k][0]].end()) {
                                minId = false;
                                break;
                            }
                        }
                    }
                    if (minId) {
                        flag[j] = 1;
                        writeNodeList_[psize] = j;
                        integer id = tsize + psize;
                        antiOrderMap_.emplace(j, id);
                        for (integer k = 0; k < pairsList.size(); ++k) {
                            if (pairsList[k][0] > pid) {
                                // send the output index of the shared point of the minimum process to the remote process
                                // Note the data structure of pairsList[k] is (remote process id, shared point's local index of remote process)
                                // Namely: pairsList[k][0] = remote process id
                                //		   pairsList[k][1] = shared point's local index of remote process
                                auto iter = faceZoneSharedPMap[pairsList[k][0]].find(j);
                                if (iter != faceZoneSharedPMap[pairsList[k][0]].end()) {
                                    sendL[pairsList[k][0]][countIP[pairsList[k][0]]++] =
                                        pairsList[k][1];
                                    sendL[pairsList[k][0]][countIP[pairsList[k][0]]++] = id;
                                }
                            }
                        }

                        psize++;
                    } else {
                        flag[j] = 2;
                    }
                } else // For not-shared points
                {
                    writeNodeList_[psize] = j;
                    antiOrderMap_.emplace(j, tsize + psize);
                    psize++;
                    flag[j] = 0;
                }
            }
        }
    }
    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (sendSizeL[ip] != 0) {
            HurMPI::sendList(sendL[ip], ip, HurMPI::getProcRank(), HurMPI::getComm());
        }

        if (recvSizeL[ip] != 0) {
            HurMPI::recvList(recvL[ip], ip, ip, HurMPI::getComm(), &status);
        }
    }
    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (recvSizeL[ip] != 0) {
            for (integer ii = 0; ii < recvSizeL[ip]; ++ii) {
                const auto j = recvL[ip][2 * ii];
                const auto id = recvL[ip][2 * ii + 1];
                antiOrderMap_.emplace(j, id);
            }
        }
    }
}

void OpenHurricane::globalFaceZoneIndex::getNodeListForSerialRun() {
    const auto &fz = fZ();
    totalFaces_ = fz.size();
    const auto &f = mesh_.faces();
    const auto &p = mesh_.points();
    integerList flag(p.size(), -1);
    // The size of points of this face zone.
    integer psize = 0;

    for (integer fi = fz.firstIndex(); fi <= fz.lastIndex(); ++fi) {
        for (integer i = 0; i < f[fi].size(); ++i) {
            const integer k = f[fi][i];
            if (flag[k] == -1) {
                flag[k] = 0;
                psize++;
            }
        }
    }

    writeNodeList_.resize(psize, -1);
    totalNodes_ = psize;
    flag = -1;
    psize = 0;
    for (integer fi = fz.firstIndex(); fi <= fz.lastIndex(); ++fi) {
        for (integer i = 0; i < f[fi].size(); ++i) {
            const integer k = f[fi][i];
            if (flag[k] == -1) {
                flag[k] = 0;
                writeNodeList_[psize] = k;
                antiOrderMap_.emplace(k, psize);
                psize++;
            }
        }
    }
}

void OpenHurricane::globalFaceZoneIndex::getPoints() {
    const auto &p = mesh_.points();
    if (!HurMPI::parRun()) {
        facePoints_.resize(totalNodes_);
        for (integer j = 0; j < writeNodeList_.size(); ++j) {
            const integer k = writeNodeList_[j];
            facePoints_[j] = p[k];
        }
        return;
    }

    const auto &f = mesh_.faces();
    integerList nSizeL(HurMPI::getProcSize(), Zero);
    integerList displs;
    if (HurMPI::master()) {
        displs.resize(HurMPI::getProcSize(), Zero);
    }
    nSizeL[HurMPI::getProcRank()] = writeNodeList_.size();

    HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());

    if (HurMPI::master()) {
        facePoints_.resize(totalNodes_);
        for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
            displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
        }
    }

    realList pp(writeNodeList_.size(), Zero);
    realList rootPp;
    if (HurMPI::master()) {
        rootPp.resize(totalNodes_);
    }
    HurMPI::barrier(HurMPI::getComm());
    for (int i = 0; i < feature<vector>::nElements_; ++i) {
        for (integer j = 0; j < writeNodeList_.size(); ++j) {
            const integer k = writeNodeList_[j];
            pp[j] = p[k][i];
        }
        HurMPI::Request request;
        HurMPI::igatherv(pp.data(), pp.size(), feature<point>::MPIType, rootPp.data(),
                         nSizeL.data(), displs.data(), feature<point>::MPIType, HurMPI::masterNo(),
                         HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::master()) {
            for (integer j = 0; j < rootPp.size(); ++j) {
                facePoints_[j][i] = rootPp[j];
            }
        }
    }
}

OpenHurricane::globalFaceZoneIndex::globalFaceZoneIndex(const geometryMesh &mesh, const integer id)
    : mesh_(mesh), id_(id), totalNodes_(0), totalFaces_(0) {
    getNodeList();
    getPoints();

    faceRecvcnt_.resize(HurMPI::getProcSize(), Zero);
    faceRecvcnt_[HurMPI::getProcRank()] = fZ().size();
    HurMPI::gatherList(faceRecvcnt_, HurMPI::masterNo(), HurMPI::getComm());
    if (HurMPI::master()) {
        faceDispls_.resize(HurMPI::getProcSize(), Zero);
    }
    if (HurMPI::master()) {
        for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
            faceDispls_[ip] = faceDispls_[ip - 1] + faceRecvcnt_[ip - 1];
        }
    }
    if (!HurMPI::master()) {
        faceRecvcnt_.clear();
    }
}

OpenHurricane::globalFaceZoneIndex::~globalFaceZoneIndex() noexcept {}

hur_nodiscard const OpenHurricane::faceZone &OpenHurricane::globalFaceZoneIndex::fZ() const {
    const auto &fzl = mesh_.faceZones();
    for (integer i = 0; i < fzl.size(); ++i) {
        if (fzl[i].index() == id_) {
            return fzl[i];
        }
    }
    checkWarningStr(("Unknown face zone index: " + toString(id_)));
    return NullRefObj::nullRef<faceZone>();
}