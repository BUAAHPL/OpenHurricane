/*!
 * \file decomposingByMetis.cpp
 * \brief Main subroutines for decomposing mesh.
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

#include "decomposingByMetis.hpp"
#include "metis.h"

void OpenHurricane::decomposingMeshByMetis::computeArray() {
    //const integer ncell = cells_.size();
    const integer ncell = originMesh_.cells().size();
    //const integer mm = nEdges_;
    const integer mm = originMesh_.globalInterior();

    xadj_.resize(ncell + 1);
    adjncy_.resize(2 * mm);

    const cellList &cells = originMesh_.cells();
    const faceList &faces = originMesh_.faces();

    integer k = 0;
    xadj_[0] = 0;
    for (integer celli = 0; celli < ncell; celli++) {
        xadj_[celli + 1] = xadj_[celli];

        for (integer fci = 0; fci < cells[celli].faceSize(); fci++) {
            integer facei = cells[celli].facei(fci);
            const face &fc = faces[facei];

            if (fc.isInterior() && (!fc.isInterface())) {
                adjncy_[k] = fc.leftCell() + fc.rightCell() - celli;
                xadj_[celli + 1] += 1;
                k += 1;
            }
        }
    }
}

int OpenHurricane::decomposingMeshByMetis::decomposing() {
    if (!HurMPI::master()) {
        return 1;
    }
    if (loadWgtPtr_ != nullptr) {
        if (!loadWgtPtr_->isUnWeight()) {
            return decomposing(loadWgtPtr_->weights());
        }
    }
    const cellList &cells = originMesh_.cells();
    const faceList &faces = originMesh_.faces();
    const integer ncell = cells.size();
    partOfProcs_.resize(ncell);

    if (nparts_ == 1) {
        for (integer i = 0; i < partOfProcs_.size(); i++) {
            partOfProcs_[i] = 0;
        }
        return 1;
    } else if (nparts_ == originMeshDecompSize_) {
        Pout << "    Info: decomposing mesh..." << std::endl;
        decomposeFromOriginMesh();
        return 1;
    }

    Pout << "    Info: decomposing mesh by MetiS..." << std::endl;
    computeArray();

    idx_t *vwgt, *vsize;
    vwgt = NULL;
    vsize = NULL;

    idx_t *adjwgt;
    adjwgt = NULL;

    integer nvtxs, ncon;
    real_t *tpwgts, *ubvec;

    tpwgts = NULL;
    ubvec = NULL;

    nvtxs = ncell;
    ncon = 1;

    idx_t options0[METIS_NOPTIONS];

    idx_t objval;

    METIS_SetDefaultOptions(options0);

    // options[METIS OPTION NUMBERING]
    //	 Used to indicate which numbering scheme is used for the adjacency structure of a graph or the elementnode
    //	 structure of a mesh.The possible values are :
    //    0 C - style numbering is assumed that starts from 0.
    //	  1 Fortran - style numbering is assumed that starts from 1.
    options0[METIS_OPTION_NUMBERING] = 0;

    // METIS_OBJTYPE_CUT: Edge-cut minimization
    // METIS_OBJTYPE_VOL: Total communication volume minimization
    options0[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options0[METIS_OPTION_NITER] = 15;

    int returnFlag;
    if (nparts_ <= 8 && ncell <= 5000000) {
        options0[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
        returnFlag = METIS_PartGraphRecursive(&nvtxs, &ncon, xadj_.data(), adjncy_.data(), vwgt,
                                              vsize, adjwgt, &nparts_, tpwgts, ubvec, options0,
                                              &objval, partOfProcs_.data());
    } else {
        options0[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
        returnFlag =
            METIS_PartGraphKway(&nvtxs, &ncon, xadj_.data(), adjncy_.data(), vwgt, vsize, adjwgt,
                                &nparts_, tpwgts, ubvec, options0, &objval, partOfProcs_.data());
    }

    clear();

    /*!\brief Check the return flag of the metis function.*/
    if (returnFlag == METIS_OK) {
        return 1;
    } else if (returnFlag == METIS_ERROR_INPUT) {
        Pout << "An input error occur in metis function!" << std::endl;
        return 0;
    } else if (returnFlag == METIS_ERROR_MEMORY) {
        Pout << "Error: could not allocate the required memory in metis "
                "function!"
             << std::endl;
        return 0;
    } else if (returnFlag == METIS_ERROR) {
        Pout << "Error: Unknown error in metis functions!" << std::endl;
        return 0;
    }

    return 0;
}

void OpenHurricane::decomposingMeshByMetis::decomposeFromOriginMesh() {
    for (integer czi = 0; czi < cellZones_.size(); ++czi) {
        integer count = 0;
        for (integer pi = 0; pi < decomposeList_[czi].size(); ++pi) {
            for (integer i = 0; i < decomposeList_[czi][pi]; ++i) {
                integer cellIndex = cellZones_[czi].firstIndex() + i + count;
                partOfProcs_[cellIndex] = pi;
            }
            count += decomposeList_[czi][pi];
        }
    }
}

int OpenHurricane::decomposingMeshByMetis::decomposing(const integerArray &cellWeight) {
    if (!HurMPI::master()) {
        return 1;
    }
    const cellList &cells = originMesh_.cells();
    const faceList &faces = originMesh_.faces();
    const integer ncell = cells.size();
    partOfProcs_.resize(ncell);

    if (nparts_ == 1) {
        for (integer i = 0; i < partOfProcs_.size(); i++) {
            partOfProcs_[i] = 0;
        }
        return 1;
    }

    Pout << "    Info: decomposing mesh by MetiS..." << std::endl;
    computeArray();

    idx_t *vwgt, *vsize;
    vwgt = NULL;

    if (cellWeight.size() == ncell) {
        vwgt = new idx_t[ncell];

        for (integer i = 0; i < ncell; ++i) {
            vwgt[i] = cellWeight[i];
        }
    }

    vsize = NULL;

    idx_t *adjwgt;
    adjwgt = NULL;

    integer nvtxs, ncon;
    real_t *tpwgts, *ubvec;

    tpwgts = NULL;
    ubvec = NULL;

    nvtxs = ncell;
    ncon = 1;

    idx_t options0[METIS_NOPTIONS];

    idx_t objval;

    METIS_SetDefaultOptions(options0);

    // options[METIS OPTION NUMBERING]
    //	 Used to indicate which numbering scheme is used for the adjacency structure of a graph or the elementnode
    //	 structure of a mesh.The possible values are :
    //    0 C - style numbering is assumed that starts from 0.
    //	  1 Fortran - style numbering is assumed that starts from 1.
    options0[METIS_OPTION_NUMBERING] = 0;

    // METIS_OBJTYPE_CUT: Edge-cut minimization
    // METIS_OBJTYPE_VOL: Total communication volume minimization
    options0[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options0[METIS_OPTION_NITER] = 15;

    int returnFlag;
    if (nparts_ <= 8 && ncell <= 5000000) {
        options0[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
        returnFlag = METIS_PartGraphRecursive(&nvtxs, &ncon, xadj_.data(), adjncy_.data(), vwgt,
                                              vsize, adjwgt, &nparts_, tpwgts, ubvec, options0,
                                              &objval, partOfProcs_.data());
    } else {
        options0[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
        returnFlag =
            METIS_PartGraphKway(&nvtxs, &ncon, xadj_.data(), adjncy_.data(), vwgt, vsize, adjwgt,
                                &nparts_, tpwgts, ubvec, options0, &objval, partOfProcs_.data());
    }

    clear();

    if (vwgt != NULL) {
        delete[] vwgt;
    }

    /*!\brief Check the return flag of the metis function.*/
    if (returnFlag == METIS_OK) {
        return 1;
    } else if (returnFlag == METIS_ERROR_INPUT) {
        Pout << "An input error occur in metis function!" << std::endl;
        return 0;
    } else if (returnFlag == METIS_ERROR_MEMORY) {
        Pout << "Error: could not allocate the required memory in metis "
                "function!"
             << std::endl;
        return 0;
    } else if (returnFlag == METIS_ERROR) {
        Pout << "Error: Unknown error in metis functions!" << std::endl;
        return 0;
    }

    return 0;
}

void OpenHurricane::decomposingMeshByMetis::meshToProcessor() {
    if (HurMPI::getProcSize() == 1) {
        cellToProcessor();
        //isUnStructedMesh();
        integerArrayArray faceofproc_(nparts_);
        faceToProcessor(faceofproc_);
        //PointToProcessor(faceofproc_);
        PointToProcessorFix(faceofproc_);
        cellAssign();
        //isTwoLayer();
        //formCutPerZone(cutZones_);
        //formCutPerZone(perZones_);
        color_.clear();
    } else {
        Pout << "    Info: distributing mesh data..." << std::endl;
        cellToProcessor();
        //isUnStructedMesh();
        integerArrayArray faceofproc_(nparts_);
        faceToProcessor(faceofproc_);
        //PointToProcessor(faceofproc_);
        PointToProcessorFix(faceofproc_);
        cellAssign();
        //isTwoLayer();
        //formCutPerZone(cutZones_);
        //formCutPerZone(perZones_);
        color_.clear();
        Pout << "    Info: distributing mesh complete" << std::endl;
    }
}

void OpenHurricane::decomposingMeshByMetis::cellToProcessor() {
    if (HurMPI::master()) {
        color_.resize((originMesh_.cells().size()));
        cellZones_.resize(originMesh().cellZones().size());
        cellZones_ = originMesh().cellZones();
    }
    bcastZone(cellZones_);

    integerArrayArray cellOfProc_(nparts_);

    cellAlloc(cellOfProc_);

    cellScatter(cellOfProc_ /*, color*/);
}

void OpenHurricane::decomposingMeshByMetis::cellAlloc(integerArrayArray &cellOfProc) {
    integerArrayArray null;
    //integerArray nullColor;
    //integerArray size(nparts_, OpenHurricane::Zero);
    integer *size = new integer[nparts_]();

    if (HurMPI::master()) {
        for (integer i = 0; i < originMesh_.cells().size(); i++) {
            cellDistribute(null /*, nullColor*/, size, true, i);
        }

        for (integer ip = 0; ip < nparts_; ip++) {
            cellOfProc[ip].resize(size[ip]);
        }
    }

    integer localSize = size[HurMPI::getProcRank()];
    HurMPI::barrierDefaultComm();

    HurMPI::scatter(size, 1, feature<integer>::MPIType, &localSize, 1, feature<integer>::MPIType,
                    HurMPI::masterNo(), HurMPI::getComm());

    cells_.resize(localSize);

    if (HurMPI::barrierDefaultComm()) {
        delete[] size;
    }
}

void OpenHurricane::decomposingMeshByMetis::cellDistribute(integerArrayArray &cellOfProc,
                                                           //integerArray& color,
                                                           integer *size, const bool count,
                                                           const integer id) {
    integer ip = partOfProcs_[id];
    if (!count) {
        cellOfProc[ip][size[ip]] = id;
        color_[id] = size[ip];
    }
    size[ip]++;
}

void OpenHurricane::decomposingMeshByMetis::cellScatter(
    integerArrayArray &cellOfProc /*, integerArray& color*/) {
    integer *size = new integer[nparts_]();
    for (integer zoneI = 0; zoneI < cellZones_.size(); ++zoneI) {
        integer *oldSize = new integer[nparts_];
        for (integer ip = 0; ip < nparts_; ip++) {
            oldSize[ip] = size[ip];
        }

        integer *scatterInfo = new integer[1 + nparts_];
        integerArray sendBuf;
        integerArray displs;

        if (HurMPI::master()) {
            integer firstIndex = originMesh_.cellZones()[zoneI].firstIndex();
            integer lastIndex = originMesh_.cellZones()[zoneI].lastIndex();

            for (integer i = firstIndex; i < lastIndex + 1; i++) {
                cellDistribute(cellOfProc /*, color*/, size, false, i);
            }
            cellMember(cellOfProc, sendBuf, displs, oldSize, size, scatterInfo);
        }
        HurMPI::bcast(scatterInfo, 1 + nparts_, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());

        integerArray zFL(2 * nparts_);
        integerArray fl(2);
        integerArray dispfl(nparts_);
        integerArray sCountfl(nparts_, integer(2));
        for (integer ip = 0; ip < nparts_; ip++) {
            dispfl[ip] = 2 * ip;
            zFL[2 * ip] = oldSize[ip];
            zFL[2 * ip + 1] = size[ip];
        }
        HurMPI::scatterv(zFL.data(), sCountfl.data(), dispfl.data(), feature<integer>::MPIType,
                         fl.data(), sCountfl[HurMPI::getProcRank()], feature<integer>::MPIType,
                         HurMPI::masterNo(), HurMPI::getComm());

        integerArray sendCount(nparts_);
        integerArray recvBuf(scatterInfo[nparts_]);
        //short* recvBuf = new short[scatterInfo[nparts_]];

        for (integer ip = 0; ip < nparts_; ip++) {
            sendCount[ip] = scatterInfo[ip];
        }

        HurMPI::scatterv(sendBuf.data(), sendCount.data(), displs.data(), feature<integer>::MPIType,
                         recvBuf.data(), sendCount[HurMPI::getProcRank()],
                         feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

        cellReconstr(fl, recvBuf, sendCount, zoneI);

        delete[] oldSize;
        delete[] scatterInfo;
    }
    if (HurMPI::master()) {
        for (integer ip = 1; ip < nparts_; ip++) {
            cellOffSet_[ip] = cellOffSet_[ip - 1] + size[ip - 1];
        }
    }
    HurMPI::bcastList(cellOffSet_, HurMPI::masterNo(), HurMPI::getComm());
    delete[] size;
}

void OpenHurricane::decomposingMeshByMetis::cellMember(const integerArrayArray &cellOfProc,
                                                       integerArray &sendBuf, integerArray &displs,
                                                       integer *oldSize, integer *size,
                                                       integer *info) {
    integer total = 0;
    for (integer ip = 0; ip < nparts_; ip++) {
        integer i = size[ip] - oldSize[ip];
        info[ip] = 2 * i;
        total += 2 * i;
        //info[ip] = i;
        //total += i;
    }
    sendBuf.resize(total);

    displs.resize(nparts_);
    displs[0] = 0;
    for (integer ip = 1; ip < nparts_; ip++) {
        displs[ip] = info[ip - 1] + displs[ip - 1];
    }
    integer recvSize = *std::max_element(info, info + nparts_);
    info[nparts_] = recvSize;

    for (integer ip = 0; ip < nparts_; ip++) {
        integer start = displs[ip];
        for (integer i = oldSize[ip]; i < size[ip]; i++) {
            integer j = i - oldSize[ip];
            integer id = cellOfProc[ip][i];
            //sendBuf[2 * j + start] = originMesh_.cells()[id].wallFlag();
            //sendBuf[2 * j + 1 + start] = originMesh_.cells()[id].types();
            //sendBuf[j + start] = originMesh_.cells()[id].types();
            sendBuf[2 * j + start] = originMesh_.cells()[id].shapeType();
            sendBuf[2 * j + 1 + start] = originMesh_.cells()[id].orignIndex();
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::cellReconstr(const integerArray &fl,
                                                         const integerArray &recvBuf,
                                                         const integerArray &sendCount,
                                                         const integer zoneI) {
    cellZones_[zoneI].setFirstIndex(fl[0]);
    cellZones_[zoneI].setLastIndex(fl[1] - 1);
    integer j = 0;
    for (integer i = fl[0]; i < fl[1]; i++) {
        //cells_[i].setWallFlag(recvBuf[j++]);
        //cells_[i].setType(recvBuf[j++]);
        cells_[i].setType(short(recvBuf[j++]));
        cells_[i].setOrignIndex(recvBuf[j++]);
    }
}

void OpenHurricane::decomposingMeshByMetis::cellAssign() {
    originMeshRead::formingFaces(cells_.size(), faces_, cells_);
    originMeshRead::formingNodes(cells_.size(), faces_, cells_);
    originMeshRead::computingCellWall(nodes_.size(), cells_.size(), faceZones_, faces_, cells_);
}

void OpenHurricane::decomposingMeshByMetis::getPeriodicTransformVector(const periodicPair &pp,
                                                                       point &tranfV) const {
    if (HurMPI::master()) {
        const integer pf = pp[0].x();
        const integer sf = pp[0].y();

        const auto &pface = originMesh_.faces()[pf];
        const auto &sface = originMesh_.faces()[sf];

        pointField pPoints(pface.size());
        for (integer i = 0; i < pface.size(); ++i) {
            pPoints[i] = originMesh_.points()[pface[i]];
        }
        point pfA, pfC;
        face::areaAndCentre(pPoints, pfA, pfC);
        pointField sPoints(sface.size());
        for (integer i = 0; i < sface.size(); ++i) {
            sPoints[i] = originMesh_.points()[sface[i]];
        }
        point sfA, sfC;
        face::areaAndCentre(sPoints, sfA, sfC);

        tranfV = sfC - pfC;
    }
}

void OpenHurricane::decomposingMeshByMetis::getPeriodicRotationalMatrix(const periodicPair &pp,
                                                                        tensor &RM) const {
    if (HurMPI::master()) {
        const integer pf = pp[0].x();
        const integer sf = pp[0].y();

        const auto &pface = originMesh_.faces()[pf];
        const auto &sface = originMesh_.faces()[sf];

        pointField pPoints(pface.size());
        for (integer i = 0; i < pface.size(); ++i) {
            pPoints[i] = originMesh_.points()[pface[i]];
        }
        point pfA, pfC;
        face::areaAndCentre(pPoints, pfA, pfC);
        pointField sPoints(sface.size());
        for (integer i = 0; i < sface.size(); ++i) {
            sPoints[i] = originMesh_.points()[sface[i]];
        }
        point sfA, sfC;
        face::areaAndCentre(sPoints, sfA, sfC);
        real theta = fabs(cos(pfA, sfA));
        theta = acos(theta);

        vector direction = pfA ^ sfA;
        vector axisv = axis_ - origin_;

        real flag = direction * axisv;
        if (flag > 0) {
            theta = -theta;
        }

        calcRotationMatrix(axisv, theta, RM);
    }
}

void OpenHurricane::decomposingMeshByMetis::calcRotationMatrix(const vector ra, const real theta,
                                                               tensor &RM) const {
    if (ra.x() != 0.0 && ra.y() == 0.0 && ra.z() == 0.0) {
        RM = tensor(1.0, 0.0, 0.0, 0.0, cos(theta), -sin(theta), 0.0, sin(theta), cos(theta));
    } else if (ra.x() == 0.0 && ra.y() != 0.0 && ra.z() == 0.0) {
        RM = tensor(cos(theta), 0.0, sin(theta), 0.0, 1.0, 0.0, -sin(theta), 0.0, cos(theta));
    } else if (ra.x() == 0.0 && ra.y() == 0.0 && ra.z() != 0.0) {
        RM = tensor(cos(theta), -sin(theta), 0.0, sin(theta), cos(theta), 0.0, 0.0, 0.0, 1.0);
    } else if (ra.x() == 0.0 && ra.y() == 0.0 && ra.z() == 0.0) {
        LFatal("The rotation axis must be set");
    } else {
        real cosa = ra.z() / sqrt(sqr(ra.y()) + sqr(ra.z()));
        real sina = ra.y() / sqrt(sqr(ra.y()) + sqr(ra.z()));

        real cosb = sqrt(sqr(ra.y()) + sqr(ra.z())) / sqrt(sqr(ra.x()) + sqr(ra.y()) + sqr(ra.z()));
        real sinb = ra.x() / sqrt(sqr(ra.x()) + sqr(ra.y()) + sqr(ra.z()));

        tensor Rx(1.0, 0.0, 0.0, 0.0, cosa, -sina, 0.0, sina, cosa);

        tensor Rxma(1.0, 0.0, 0.0, 0.0, cosa, sina, 0.0, -sina, cosa);

        tensor Ry(cosb, 0.0, sinb, 0.0, 1.0, 0.0, -sinb, 0.0, cosb);

        tensor Rymb(cosb, 0.0, -sinb, 0.0, 1.0, 0.0, sinb, 0.0, cosb);

        tensor Rz(cos(theta), -sin(theta), 0.0, sin(theta), cos(theta), 0.0, 0.0, 0.0, 1.0);

        RM = Rx * Ry * Rz * Rymb * Rxma;
    }
}

void OpenHurricane::decomposingMeshByMetis::faceToProcessor(integerArrayArray &faceOfProc_) {
    if (HurMPI::master()) {
        /*!\brief Reorder faceZone(interiorzone,cutzone,bndzone,interfacezone).*/
        integer czSize = 0;
        for (integer zoneI = 0; zoneI < originMesh_.faceZones().size(); ++zoneI) {
            if (originMesh_.faceZones()[zoneI].isInterior() &&
                (!originMesh_.faceZones()[zoneI].isInterface())) {
                czSize++;
            }
        }
        faceZones_.resize(czSize + originMesh_.faceZones().size());

        integer count = Zero;

        for (integer zoneI = 0; zoneI < originMesh_.faceZones().size(); zoneI++) {
            if (originMesh_.faceZones()[zoneI].isInterior() &&
                (!originMesh_.faceZones()[zoneI].isInterface())) {
                faceZones_[count++] = originMesh_.faceZones()[zoneI];
            }
        }
        for (integer zoneI = 0; zoneI < originMesh_.faceZones().size(); zoneI++) {
            if (originMesh_.faceZones()[zoneI].isBnd()) {
                faceZones_[count++] = originMesh_.faceZones()[zoneI];
            }
        }
        for (integer zoneI = 0; zoneI < originMesh_.faceZones().size(); zoneI++) {
            if (originMesh_.faceZones()[zoneI].isInterior() &&
                (!originMesh_.faceZones()[zoneI].isInterface())) {
                faceZones_[count] = originMesh_.faceZones()[zoneI];
                std::string name = "cut_" + originMesh_.faceZones()[zoneI].name();
                faceZones_[count].resetName(name);
                faceZones_[count].setIndex(-originMesh_.faceZones()[zoneI].index());
                faceZones_[count].setBcType(1);
                count++;
            }
        }
        for (integer zoneI = 0; zoneI < originMesh_.faceZones().size(); zoneI++) {
            if (originMesh_.faceZones()[zoneI].isInterface()) {
                faceZones_[count++] = originMesh_.faceZones()[zoneI];
            }
        }
    }
    bcastZone(faceZones_);

    //integerArrayArray faceOfProc_(nparts_);
    integerArrayArrayArray zoneFL(faceZones_.size());
    integerArrayArray cutColor(nparts_);
    for (integer ip = 0; ip < nparts_; ip++) {
        cutColor[ip].resize(nparts_);
    }
    for (integer i = 0; i < faceZones_.size(); i++) {
        zoneFL[i].resize(nparts_);
        for (integer ip = 0; ip < nparts_; ip++) {
            zoneFL[i][ip].resize(2);
        }
    }

    faceAlloc(faceOfProc_, cutColor);

    integerArray size(nparts_, OpenHurricane::Zero);

    faceScatter(faceOfProc_, cutColor, zoneFL, size);

    checkPartition(zoneFL);

    faceShare(faceOfProc_, size, originMesh_.faces().size());
}

void OpenHurricane::decomposingMeshByMetis::perFaceToProcessor(integerArray &globalColor,
                                                               integerArrayArray &faceOfProc) {
    integerArrayArray pzColor(nparts_);
    for (integer ip = 0; ip < nparts_; ip++) {
        pzColor[ip].resize(nparts_);
    }

    perZoneAlloc(pzColor, globalColor);
    perZoneScatter(pzColor, faceOfProc, globalColor);
}

void OpenHurricane::decomposingMeshByMetis::perZoneAlloc(integerArrayArray &pzColor,
                                                         integerArray &globalColor) {
    integer ppSize = 0;
    if (HurMPI::master()) {
        ppSize = originMesh_.periodicPairZones().size();
    }
    HurMPI::bcast(&ppSize, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

    periodicPairSize_ = ppSize;

    for (integer zoneI = 0; zoneI < ppSize; ++zoneI) {
        integerArrayArray perSize(nparts_);
        integerArrayArray sdwSize(nparts_);
        for (integer ip = 0; ip < nparts_; ip++) {
            perSize[ip].resize(nparts_);
            sdwSize[ip].resize(nparts_);
            perSize[ip] = 0;
            sdwSize[ip] = 0;
        }
        integer perType = -1;
        if (HurMPI::master()) {
            perType = (short)originMesh_.periodicPairZones()[zoneI].type();
            perFaceDistribute(originMesh_.periodicPairZones()[zoneI], globalColor, perSize, sdwSize,
                              pzColor, true);
        }
        HurMPI::bcast(&perType, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        integerArrayArray szColor(nparts_);
        for (integer ip = 0; ip < nparts_; ip++) {
            szColor[ip].resize(nparts_, -1);
            pzColor[ip] = -1;
        }
        integer pz = 0;
        if (HurMPI::master()) {
            for (integer ip1 = 0; ip1 < nparts_; ip1++) {
                for (integer ip2 = 0; ip2 < nparts_; ip2++) {
                    if (perSize[ip1][ip2]) {
                        pzColor[ip1][ip2] = pz++;
                    }
                    if (sdwSize[ip1][ip2]) {
                        szColor[ip1][ip2] = pz++;
                    }
                }
            }
        }
        HurMPI::bcast(&pz, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        perZoneList newPerZones(pz);

        if (HurMPI::master()) {
            for (integer ip1 = 0; ip1 < nparts_; ip1++) {
                for (integer ip2 = 0; ip2 < nparts_; ip2++) {
                    if (perSize[ip1][ip2]) {
                        newPerZones[pzColor[ip1][ip2]].setType(short(perType));
                        newPerZones[pzColor[ip1][ip2]].setPeriodic();
                        newPerZones[pzColor[ip1][ip2]].setSendProc(ip1);
                        newPerZones[pzColor[ip1][ip2]].setRecvProc(ip2);
                        newPerZones[pzColor[ip1][ip2]].sor().resize(perSize[ip1][ip2]);
                        newPerZones[pzColor[ip1][ip2]].des().resize(perSize[ip1][ip2]);
                        newPerZones[pzColor[ip1][ip2]].sorFace().resize(perSize[ip1][ip2]);
                        newPerZones[pzColor[ip1][ip2]].desFace().resize(perSize[ip1][ip2]);
                    }
                    if (sdwSize[ip1][ip2]) {
                        newPerZones[szColor[ip1][ip2]].setType(short(perType));
                        newPerZones[szColor[ip1][ip2]].setShadow();
                        newPerZones[szColor[ip1][ip2]].setSendProc(ip1);
                        newPerZones[szColor[ip1][ip2]].setRecvProc(ip2);
                        newPerZones[szColor[ip1][ip2]].sor().resize(sdwSize[ip1][ip2]);
                        newPerZones[szColor[ip1][ip2]].des().resize(sdwSize[ip1][ip2]);
                        newPerZones[szColor[ip1][ip2]].sorFace().resize(sdwSize[ip1][ip2]);
                        newPerZones[szColor[ip1][ip2]].desFace().resize(sdwSize[ip1][ip2]);
                    }
                }
            }

            for (integer ip = 0; ip < nparts_; ip++) {
                perSize[ip] = 0;
                sdwSize[ip] = 0;
            }
            createPerZone(originMesh_.periodicPairZones()[zoneI], globalColor, perSize, sdwSize,
                          pzColor, szColor, newPerZones);
        }
        vector tranfV;
        tensor rotationMatrix;

        if (HurMPI::master()) {
            if (perType == periodicTypes::TRANSLATIONAL) {
                getPeriodicTransformVector(originMesh_.periodicPairZones()[zoneI], tranfV);
            } else {
                getPeriodicRotationalMatrix(originMesh_.periodicPairZones()[zoneI], rotationMatrix);
            }
        }
        if (perType == periodicTypes::TRANSLATIONAL) {
            for (int i = 0; i < vector::nElements_; ++i) {
                HurMPI::bcast(&tranfV[i], 1, feature<real>::MPIType, HurMPI::masterNo(),
                              HurMPI::getComm());
            }
        } else {
            for (int i = 0; i < tensor::nElements_; ++i) {
                HurMPI::bcast(&rotationMatrix[i], 1, feature<real>::MPIType, HurMPI::masterNo(),
                              HurMPI::getComm());
            }
        }
        for (integer i = 0; i < pz; i++) {
            HurMPI::bcast(&newPerZones[i].id(), 1, feature<integer>::MPIType, HurMPI::masterNo(),
                          HurMPI::getComm());
            HurMPI::bcast(&newPerZones[i].pairId(), 1, feature<integer>::MPIType,
                          HurMPI::masterNo(), HurMPI::getComm());
            integer isp = newPerZones[i].sendProc();
            integer irp = newPerZones[i].receivProc();
            integer isize = newPerZones[i].perSize();
            bool isShadow = newPerZones[i].isShadow();

            HurMPI::bcast(&isp, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                          HurMPI::getComm());
            HurMPI::bcast(&irp, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                          HurMPI::getComm());
            HurMPI::bcast(&isize, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                          HurMPI::getComm());
            HurMPI::bcast(&isShadow, 1, feature<bool>::MPIType, HurMPI::masterNo(),
                          HurMPI::getComm());

            if (!HurMPI::master()) {
                newPerZones[i].setType(short(perType));
                newPerZones[i].setSendProc(isp);
                newPerZones[i].setRecvProc(irp);
                newPerZones[i].sor().resize(isize);
                newPerZones[i].des().resize(isize);
                newPerZones[i].sorFace().resize(isize);
                newPerZones[i].desFace().resize(isize);
                if (isShadow) {
                    newPerZones[i].setShadow();
                } else {
                    newPerZones[i].setPeriodic();
                }
            }
            HurMPI::bcastList(newPerZones[i].des(), HurMPI::masterNo(), HurMPI::getComm());
            HurMPI::bcastList(newPerZones[i].sor(), HurMPI::masterNo(), HurMPI::getComm());
            HurMPI::bcastList(newPerZones[i].desFace(), HurMPI::masterNo(), HurMPI::getComm());
            HurMPI::bcastList(newPerZones[i].sorFace(), HurMPI::masterNo(), HurMPI::getComm());

            if (perType == periodicTypes::TRANSLATIONAL) {
                if (isShadow) {
                    newPerZones[i].setTranDistance(-tranfV);
                } else {
                    newPerZones[i].setTranDistance(tranfV);
                }
            } else {
                if (isShadow) {
                    newPerZones[i].setRotionalMatrix(inv(rotationMatrix));
                } else {
                    newPerZones[i].setRotionalMatrix(rotationMatrix);
                }
            }
        }
        perZones_.append(newPerZones);
    }
}

void OpenHurricane::decomposingMeshByMetis::perFaceDistribute(
    const periodicPair &pp, integerArray &globalColor, integerArrayArray &persize,
    integerArrayArray &sdwsize, const integerArrayArray &pzColor, const bool count) {
    integer perZ = pp.periodicZone();
    integer sdwZ = pp.shadowZone();
    integer firstIndex = pp.firstIndex();
    integer lastIndex = pp.lastIndex();

    for (integer i = firstIndex; i < lastIndex + 1; i++) {
        integer pF1 = pp(i).x();
        integer pF2 = pp(i).y();

        integer cell1 = originMesh_.faces()[pF1].leftCell();
        integer cell2 = originMesh_.faces()[pF2].leftCell();

        integer ip1 = partOfProcs_[cell1];
        integer ip2 = partOfProcs_[cell2];

        if (!count) {
            perZones_[pzColor[ip1][ip2]].sor()[persize[ip1][ip2]] = globalColor[pF1];
            perZones_[pzColor[ip1][ip2]].des()[persize[ip1][ip2]] = globalColor[pF2];
            perZones_[pzColor[ip2][ip1]].sor()[persize[ip2][ip1]] = globalColor[pF2];
            perZones_[pzColor[ip2][ip1]].des()[persize[ip2][ip1]] = globalColor[pF1];
            perZones_[pzColor[ip1][ip2]].setId(perZ);
            perZones_[pzColor[ip1][ip2]].setPairId(sdwZ);
            perZones_[pzColor[ip2][ip1]].setId(sdwZ);
            perZones_[pzColor[ip2][ip1]].setPairId(perZ);
        }
        persize[ip1][ip2]++;
        sdwsize[ip2][ip1]++;
    }
}

void OpenHurricane::decomposingMeshByMetis::createPerZone(
    const periodicPair &pp, integerArray &globalColor, integerArrayArray &persize,
    integerArrayArray &sdwsize, const integerArrayArray &pzColor, const integerArrayArray &szColor,
    perZoneList &newPerZone) {
    integer perZ = pp.periodicZone();
    integer sdwZ = pp.shadowZone();
    integer firstIndex = pp.firstIndex();
    integer lastIndex = pp.lastIndex();

    for (integer i = firstIndex; i < lastIndex + 1; i++) {
        integer pF1 = pp(i).x();
        integer pF2 = pp(i).y();

        integer cell1 = originMesh_.faces()[pF1].leftCell();
        integer cell2 = originMesh_.faces()[pF2].leftCell();

        integer ip1 = partOfProcs_[cell1];
        integer ip2 = partOfProcs_[cell2];

        if (pzColor[ip1][ip2] != -1) {
            newPerZone[pzColor[ip1][ip2]].sor()[persize[ip1][ip2]] = globalColor[pF1];
            newPerZone[pzColor[ip1][ip2]].des()[persize[ip1][ip2]] = globalColor[pF2];
            newPerZone[pzColor[ip1][ip2]].sorFace()[persize[ip1][ip2]] = globalColor[pF1];
            newPerZone[pzColor[ip1][ip2]].desFace()[persize[ip1][ip2]] = globalColor[pF2];
            newPerZone[pzColor[ip1][ip2]].setId(perZ);
            newPerZone[pzColor[ip1][ip2]].setPairId(sdwZ);
        }
        if (szColor[ip2][ip1] != -1) {
            newPerZone[szColor[ip2][ip1]].sor()[sdwsize[ip2][ip1]] = globalColor[pF2];
            newPerZone[szColor[ip2][ip1]].des()[sdwsize[ip2][ip1]] = globalColor[pF1];
            newPerZone[szColor[ip2][ip1]].sorFace()[sdwsize[ip2][ip1]] = globalColor[pF2];
            newPerZone[szColor[ip2][ip1]].desFace()[sdwsize[ip2][ip1]] = globalColor[pF1];
            newPerZone[szColor[ip2][ip1]].setId(sdwZ);
            newPerZone[szColor[ip2][ip1]].setPairId(perZ);
        }
        persize[ip1][ip2]++;
        sdwsize[ip2][ip1]++;
    }
}

void OpenHurricane::decomposingMeshByMetis::perZoneScatter(integerArrayArray &pzColor,
                                                           const integerArrayArray &faceOfProc,
                                                           integerArray &globalColor) {
    integer size = 0;
    integerArray sCount(nparts_, Zero);

    for (integer zoneI = 0; zoneI < perZones_.size(); ++zoneI) {
        size += perZones_[zoneI].perSize();
        integer sendP = perZones_[zoneI].sendProc();
        sCount[sendP] += 4 * perZones_[zoneI].perSize();
    }

    integerArray displ(nparts_);
    integerArray sBuf(4 * size);
    displ[0] = 0;
    for (integer ip = 1; ip < nparts_; ip++) {
        displ[ip] = displ[ip - 1] + sCount[ip - 1];
    }
    if (HurMPI::master()) {
        integerArray count(HurMPI::getProcSize(), Zero);
        for (integer zoneI = 0; zoneI < perZones_.size(); ++zoneI) {
            integer sendP = perZones_[zoneI].sendProc();
            integer start = displ[sendP] + 4 * count[sendP];
            for (integer i = 0; i < perZones_[zoneI].perSize(); i++) {
                integer globalFace =
                    faceOfProc[perZones_[zoneI].receivProc()][perZones_[zoneI].des()[i]];
                sBuf[start + 4 * i] = perZones_[zoneI].sor()[i];
                sBuf[start + 4 * i + 1] = perZones_[zoneI].receivProc();
                sBuf[start + 4 * i + 2] = perZones_[zoneI].des()[i];
                sBuf[start + 4 * i + 3] = color_[originMesh_.faces()[globalFace].leftCell()];
                count[sendP]++;
            }
        }
    }

    HurMPI::bcastList(sCount, HurMPI::masterNo(), HurMPI::getComm());

    integer rsize = sCount[HurMPI::getProcRank()];

    integerArray rBuf(rsize, Zero);
    HurMPI::scatterv(sBuf.data(), sCount.data(), displ.data(), feature<integer>::MPIType,
                     rBuf.data(), sCount[HurMPI::getProcRank()], feature<integer>::MPIType,
                     HurMPI::masterNo(), HurMPI::getComm());

    for (integer i = 0; i < sCount[HurMPI::getProcRank()] / 4; i++) {
        integerArray pair(3);
        pair[0] = rBuf[4 * i + 1];
        pair[1] = rBuf[4 * i + 2];
        pair[2] = rBuf[4 * i + 3];
        perPairMaps_.emplace(rBuf[4 * i], pair);
    }
}

void OpenHurricane::decomposingMeshByMetis::faceAlloc(integerArrayArray &faceOfProc,
                                                      integerArrayArray &cutColor) {
    integerArray nullColor;
    integerArrayArray null;
    integerArrayArrayArray nullFL;
    integerArrayArray cutSize(nparts_);
    for (integer ip = 0; ip < nparts_; ip++) {
        cutSize[ip].resize(nparts_);
        cutSize[ip] = 0;
    }
    integerArray size(nparts_, OpenHurricane::Zero);
    integer cz = 0;

    if (HurMPI::master()) {
        for (integer i = 0; i < faceZones_.size(); i++) {
            if (!faceZones_[i].isCutFace()) {
                faceDistribute(nullColor, null, cutColor, cutSize, nullFL, size, true, true, i);
            }
        }
        for (integer ip1 = 0; ip1 < nparts_; ip1++) {
            faceOfProc[ip1].resize(size[ip1]);
            for (integer ip2 = 0; ip2 < nparts_; ip2++) {
                if (cutSize[ip1][ip2] != 0) {
                    cutColor[ip1][ip2] = cz++;
                }
            }
        }
    }
    HurMPI::bcast(&cz, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
    cutZones_.resize(cz);
    integer i = 0;
    integerArray info(3 * cz);
    if (HurMPI::master()) {
        for (integer ip1 = 0; ip1 < nparts_; ip1++) {
            for (integer ip2 = 0; ip2 < nparts_; ip2++) {
                if (cutSize[ip1][ip2]) {
                    info[i++] = ip1;
                    info[i++] = ip2;
                    info[i++] = cutSize[ip1][ip2];
                }
            }
        }
    }
    HurMPI::bcastList(info, HurMPI::masterNo(), HurMPI::getComm());
    for (integer i = 0; i < cz; i++) {
        cutZones_[i].setSendProc(info[3 * i]);
        cutZones_[i].setRecvProc(info[3 * i + 1]);
        if (HurMPI::master()) {
            cutZones_[i].sor().resize(info[3 * i + 2]);
            cutZones_[i].des().resize(info[3 * i + 2]);
        }
    }

    integer localSize = size[HurMPI::masterNo()];
    HurMPI::barrierDefaultComm();

    HurMPI::scatter(size.data(), 1, feature<integer>::MPIType, &localSize, 1,
                    feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

    faces_.resize(localSize);
}

void OpenHurricane::decomposingMeshByMetis::faceDistribute(
    integerArray &globalColor, integerArrayArray &faceOfProc, integerArrayArray &ctColor,
    integerArrayArray &cutSize, integerArrayArrayArray &zoneFL, integerArray &size,
    const bool count, const bool cut, const integer id) {
    integer fZIndex = -1;
    for (integer zoneI = 0; zoneI < originMesh_.faceZones().size(); ++zoneI) {
        if (originMesh_.faceZones()[zoneI].index() == faceZones_[id].index()) {
            fZIndex = zoneI;
            break;
        }
    }

    const faceZone &fZ = originMesh_.faceZones()[fZIndex];
    integer firstIndex = fZ.firstIndex();
    integer lastIndex = fZ.lastIndex();

    if (fZ.isInterior() && (!fZ.isInterface())) {
        integer cutId = -1;
        if (!count) {
            if (cut) {
                for (integer zoneI = 0; zoneI < faceZones_.size(); ++zoneI) {
                    if (faceZones_[zoneI].index() == -fZ.index()) {
                        cutId = zoneI;
                        break;
                    }
                }
                for (integer ip = 0; ip < nparts_; ip++) {
                    zoneFL[cutId][ip] = size[ip];
                }
            } else {
                for (integer ip = 0; ip < nparts_; ip++) {
                    zoneFL[id][ip] = size[ip];
                }
            }
        }

        for (integer i = firstIndex; i < lastIndex + 1; i++) {
            const auto &leftCell = originMesh_.faces()[i].leftCell();
            const auto &rightCell = originMesh_.faces()[i].rightCell();
            integer ip1 = partOfProcs_[leftCell];
            integer ip2 = partOfProcs_[rightCell];

            if (ip1 != ip2) {
                if (!count) {
                    if (cut) {
                        cutZones_[ctColor[ip1][ip2]].sor()[cutSize[ip1][ip2]] = size[ip1];
                        cutZones_[ctColor[ip1][ip2]].des()[cutSize[ip1][ip2]] = size[ip2];
                        cutZones_[ctColor[ip2][ip1]].sor()[cutSize[ip2][ip1]] = size[ip2];
                        cutZones_[ctColor[ip2][ip1]].des()[cutSize[ip2][ip1]] = size[ip1];
                        faceOfProc[ip1][size[ip1]] = i;
                        faceOfProc[ip2][size[ip2]] = i;
                        zoneFL[cutId][ip1][1]++;
                        zoneFL[cutId][ip2][1]++;
                    } else {
                        break;
                    }
                }
                size[ip1]++;
                size[ip2]++;
                cutSize[ip1][ip2]++;
                cutSize[ip2][ip1]++;
            }
        }
        for (integer i = firstIndex; i < lastIndex + 1; i++) {
            const auto &leftCell = originMesh_.faces()[i].leftCell();
            const auto &rightCell = originMesh_.faces()[i].rightCell();
            integer ip1 = partOfProcs_[leftCell];
            integer ip2 = partOfProcs_[rightCell];

            if (ip1 == ip2) {
                if (!count) {
                    if (!cut) {
                        faceOfProc[ip1][size[ip1]] = i;
                        zoneFL[id][ip1][1]++;
                    } else {
                        break;
                    }
                }
                size[ip1]++;
            }
        }
    } else {
        if (!count) {
            for (integer ip = 0; ip < nparts_; ip++) {
                zoneFL[id][ip] = size[ip];
            }
        }
        for (integer i = firstIndex; i < lastIndex + 1; i++) {
            integer leftCell = originMesh_.faces()[i].leftCell();

            if (leftCell < 0) {
                errorAbortStr(("Cannot find the left cell index of boundary "
                               "face in face zone: " +
                               fZ.name() + ". The left cell index is: " + toString(leftCell)));
            }

            integer ip1 = partOfProcs_[leftCell];
            if (!count) {
                globalColor[i] = size[ip1];
                faceOfProc[ip1][size[ip1]] = i;
                zoneFL[id][ip1][1]++;
            }
            size[ip1]++;
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::faceScatter(integerArrayArray &faceOfProc,
                                                        integerArrayArray &ctColor,
                                                        integerArrayArrayArray &zoneFL,
                                                        integerArray &size) {
    integerArrayArray cutSize(nparts_);
    integerArray globalColor(originMesh_.faces().size(), -1);
    for (integer ip = 0; ip < nparts_; ip++) {
        cutSize[ip].resize(nparts_);
        cutSize[ip] = 0;
    }
    if (HurMPI::master()) {
        for (integer i = 0; i < faceZones_.size(); i++) {
            if (faceZones_[i].isInterior() && (!faceZones_[i].isInterface())) {
                faceDistribute(globalColor, faceOfProc, ctColor, cutSize, zoneFL, size, false,
                               false, i);
            }
        }
        for (integer i = 0; i < faceZones_.size(); i++) {
            if (faceZones_[i].isInterior() && (!faceZones_[i].isInterface())) {
                faceDistribute(globalColor, faceOfProc, ctColor, cutSize, zoneFL, size, false, true,
                               i);
            }
        }
        for (integer i = 0; i < faceZones_.size(); i++) {
            if (!(faceZones_[i].isInterior() || faceZones_[i].isCutFace())) {
                faceDistribute(globalColor, faceOfProc, ctColor, cutSize, zoneFL, size, false,
                               false, i);
            }
        }
        for (integer i = 0; i < faceZones_.size(); i++) {
            if (faceZones_[i].isInterface()) {
                faceDistribute(globalColor, faceOfProc, ctColor, cutSize, zoneFL, size, false,
                               false, i);
            }
        }
        createFaceShareZone(globalColor);
    }
    for (integer zoneI = 0; zoneI < faceZones_.size(); ++zoneI) {
        if (faceZones_[zoneI].isPeriodic() || faceZones_[zoneI].isPeriodicShadow()) {
            perFaceToProcessor(globalColor, faceOfProc);
            break;
        }
    }

    integerArray info(nparts_ + 1);
    integerArray zFL(2 * nparts_);
    integerArray fl(2);
    integerArray flCount(nparts_, 2);
    integerArray disp(nparts_);
    for (integer ip = 0; ip < nparts_; ip++) {
        disp[ip] = 2 * ip;
    }

    integerArray sBufl;
    List<short> sBufs;
    integerArray displ;
    integerArray disps;
    integerArray sCountl(nparts_);
    integerArray sCounts(nparts_);
    for (integer i = 0; i < faceZones_.size(); i++) {
        if (HurMPI::master()) {
            faceMember(i, faceOfProc, zoneFL, sBufl, sBufs, displ, disps, info);
            for (integer ip = 0; ip < nparts_; ip++) {
                zFL(2 * ip) = zoneFL[i][ip][0];
                zFL(2 * ip + 1) = zoneFL[i][ip][1];
            }
        }
        HurMPI::bcastList(info, HurMPI::masterNo(), HurMPI::getComm());
        integerArray rBufl(3 * info[nparts_]);
        List<short> rBufs(info[nparts_]);
        for (integer ip = 0; ip < nparts_; ip++) {
            sCountl[ip] = 3 * info[ip];
            sCounts[ip] = info[ip];
        }

        HurMPI::scatterv(zFL.data(), flCount.data(), disp.data(), feature<integer>::MPIType,
                         fl.data(), flCount[HurMPI::getProcRank()], feature<integer>::MPIType,
                         HurMPI::masterNo(), HurMPI::getComm());
        HurMPI::scatterv(sBufl.data(), sCountl.data(), displ.data(), feature<integer>::MPIType,
                         rBufl.data(), sCountl[HurMPI::getProcRank()], feature<integer>::MPIType,
                         HurMPI::masterNo(), HurMPI::getComm());
        HurMPI::scatterv(sBufs.data(), sCounts.data(), disps.data(), feature<short>::MPIType,
                         rBufs.data(), sCounts[HurMPI::getProcRank()], feature<short>::MPIType,
                         HurMPI::masterNo(), HurMPI::getComm());
        HurMPI::barrierDefaultComm();
        faceReconstr(fl, rBufl, rBufs, i);
    }
}

void OpenHurricane::decomposingMeshByMetis::faceMember(const integer zoneI,
                                                       const integerArrayArray &faceOfProc,
                                                       integerArrayArrayArray &zoneFL,
                                                       integerArray &sBufl, List<short> &sBufs,
                                                       integerArray &displ, integerArray &disps,
                                                       integerArray &info) {
    integer total = 0;
    for (integer ip = 0; ip < nparts_; ip++) {
        integer size = zoneFL[zoneI][ip][1] - zoneFL[zoneI][ip][0];
        info[ip] = size;
        total += size;
    }
    sBufl.resize(3 * total);
    sBufs.resize(total);

    displ.resize(nparts_);
    disps.resize(nparts_);
    displ[0] = 0;
    disps[0] = 0;
    for (integer ip = 1; ip < nparts_; ip++) {
        displ[ip] = 3 * info[ip - 1] + displ[ip - 1];
        disps[ip] = info[ip - 1] + disps[ip - 1];
    }
    integer recvSize = *std::max_element(info.data(), info.data() + nparts_);
    info[nparts_] = recvSize;

    if (faceZones_[zoneI].isCutFace()) {
        for (integer ip = 0; ip < nparts_; ip++) {
            integer startl = displ[ip];
            integer starts = disps[ip];
            for (integer i = zoneFL[zoneI][ip][0]; i < zoneFL[zoneI][ip][1]; i++) {
                integer j = i - zoneFL[zoneI][ip][0];
                integer id = faceOfProc[ip][i];
                integer ip1 = partOfProcs_[originMesh_.faces()[id].leftCell()];
                if (ip1 == ip) {
                    sBufl[3 * j + startl] = color_[originMesh_.faces()[id].leftCell()];
                } else {
                    sBufl[3 * j + startl] = color_[originMesh_.faces()[id].rightCell()];
                }
                sBufl[3 * j + 1 + startl] = -1;
                sBufl[3 * j + 2 + startl] = originMesh_.faces()[id].size();
                sBufs[j + starts] = 1;
            }
        }
    } else {
        for (integer ip = 0; ip < nparts_; ip++) {
            integer startl = displ[ip];
            integer starts = disps[ip];
            for (integer i = zoneFL[zoneI][ip][0]; i < zoneFL[zoneI][ip][1]; i++) {
                integer j = i - zoneFL[zoneI][ip][0];
                integer id = faceOfProc[ip][i];
                sBufl[3 * j + startl] = color_[originMesh_.faces()[id].leftCell()];
                if (originMesh_.faces()[id].rightCell() == -1) {
                    sBufl[3 * j + 1 + startl] = -1;
                } else {
                    sBufl[3 * j + 1 + startl] = color_[originMesh_.faces()[id].rightCell()];
                }
                sBufl[3 * j + 2 + startl] = originMesh_.faces()[id].size();
                sBufs[j + starts] = originMesh_.faces()[id].type();
            }
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::faceReconstr(const integerArray &fl,
                                                         const integerArray &rBufl,
                                                         const List<short> &rBufs,
                                                         const integer zoneI) {
    integer ip = HurMPI::getProcRank();
    faceZones_[zoneI].setFirstIndex(fl[0]);
    faceZones_[zoneI].setLastIndex(fl[1] - 1);
    integer j = 0;
    for (integer i = fl[0]; i < fl[1]; i++) {
        faces_[i].leftCell() = rBufl[3 * j];
        faces_[i].rightCell() = rBufl[3 * j + 1];
        faces_[i].resize(rBufl[3 * j + 2]);
        faces_[i].setBCType(rBufs[j]);
        j++;
    }
}

void OpenHurricane::decomposingMeshByMetis::checkPartition(const integerArrayArrayArray &zoneFL) {
    if (HurMPI::master()) {
        integerArray zoneFace(faceZones_.size(), OpenHurricane::Zero);
        integer total = 0;
        for (integer i = 0; i < faceZones_.size(); i++) {
            for (integer ip = 0; ip < nparts_; ip++) {
                zoneFace[i] += zoneFL[i][ip][1] - zoneFL[i][ip][0];
            }
            total += zoneFace[i];
        }

        integer cutCount = 0;
        for (integer zoneI = 0; zoneI < originMesh_.faceZones().size(); ++zoneI) {
            if (originMesh_.faceZones()[zoneI].isInterior() &&
                (!originMesh_.faceZones()[zoneI].isInterface())) {
                for (integer i = 0; i < faceZones_.size(); i++) {
                    if (faceZones_[i].index() == -originMesh_.faceZones()[zoneI].index()) {
                        cutCount += zoneFace[i];
                    }
                }
            } else {
                for (integer i = 0; i < faceZones_.size(); i++) {
                    if (faceZones_[i].index() == originMesh_.faceZones()[zoneI].index()) {
                        if (originMesh_.faceZones()[zoneI].size() != zoneFace[i]) {
                            std::string errMsg = std::to_string(faceZones_[i].bcType());
                            errMsg += ": error";
                            errorAbortStr(errMsg);
                        }
                    }
                }
            }
        }
        if (total - cutCount / 2 != originMesh_.faces().size()) {
            std::cout << total << " " << cutCount << " " << originMesh_.faces().size() << std::endl;
            std::string errMsg = ": interior error";
            errorAbortStr(errMsg);
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::isUnStructedMesh() {
    for (integer zoneI = 0; zoneI < originMesh_.cellZones().size(); ++zoneI) {
        if (originMesh_.cellZones()[zoneI].shapeType() != cellShapeType::shapeTypes::hexahedral) {
            unStructed_ = true;
        } else {
            hasStructed_ = true;
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::PointToProcessor(const integerArrayArray &faceOfProc) {
    if (HurMPI::master()) {
        pointZones_.resize(originMesh_.pointZones().size());
        pointZones_ = originMesh_.pointZones();
    }
    bcastZone(pointZones_);

    integerArray size(nparts_, OpenHurricane::Zero);

    integerArrayArray color(nparts_);
    integer total = 0;
    if (HurMPI::master()) {
        for (integer ip = 0; ip < nparts_; ip++) {
            color[ip].resize(originMesh_.points().size());
            color[ip] = -1;
        }
        total = originMesh_.faces().size();
    }

    pointScatter(faceOfProc, color, size);

    pointShare(faceOfProc, size, color, total);
}

void OpenHurricane::decomposingMeshByMetis::PointToProcessorFix(
    const integerArrayArray &faceOfProc) {
    if (HurMPI::master()) {
        pointZones_.resize(originMesh_.pointZones().size());
        pointZones_ = originMesh_.pointZones();
    }
    bcastZone(pointZones_);

    integerArray size(nparts_, OpenHurricane::Zero);
    List<List<Vector2D<integer>>> mapOP;
    integer nFace = 0;
    if (HurMPI::master()) {
        // Set size of mapOP
        mapOP.resize(originMesh_.points().size());
        nFace = originMesh_.faces().size();
    }

    pointScatterFix(faceOfProc, mapOP, size);

    pointShareFix(faceOfProc, size, mapOP, nFace);
}

void OpenHurricane::decomposingMeshByMetis::pointScatter(const integerArrayArray &faceOfProc,
                                                         integerArrayArray &color,
                                                         integerArray &size) {
    for (integer ip = 0; ip < nparts_; ip++) {
        integer nodes = 0;
        realList localNode;
        integerArray faceNode;
        integerArray pzSize(pointZones_.size(), Zero);
        if (HurMPI::master()) {
            for (integer faceI = 0; faceI < faceOfProc[ip].size(); faceI++) {
                nodes += originMesh_.faces()[faceOfProc[ip][faceI]].size();
            }
            for (integer zoneI = 0; zoneI < originMesh_.pointZones().size(); zoneI++) {
                integer firstIndex = originMesh_.pointZones()[zoneI].firstIndex();
                integer lastIndex = originMesh_.pointZones()[zoneI].lastIndex();
                for (integer faceI = 0; faceI < faceOfProc[ip].size(); faceI++) {
                    for (integer i = 0; i < originMesh_.faces()[faceOfProc[ip][faceI]].size();
                         i++) {
                        integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                        if (nodeI >= firstIndex && nodeI <= lastIndex) {
                            if (color[ip][nodeI] == -1) {
                                color[ip][nodeI] = size[ip]++;
                                pzSize[zoneI]++;
                            }
                        }
                    }
                }
            }

            if (ip == HurMPI::masterNo()) {
                nodes_.resize(size[ip]);
                for (integer i = 0; i < originMesh_.points().size(); ++i) {
                    integer id = color[ip][i];
                    if (id != -1) {
                        nodes_[id] = originMesh_.points()[i];
                    }
                }

                for (integer faceI = 0; faceI < faceOfProc[ip].size(); ++faceI) {
                    if (faces_[faceI].isCutFace()) {
                        integer ip1 =
                            partOfProcs_[originMesh_.faces()[faceOfProc[ip][faceI]].leftCell()];
                        if (ip1 != ip) {
                            for (integer i = 0;
                                 i < originMesh_.faces()[faceOfProc[ip][faceI]].size(); ++i) {
                                integer n = originMesh_.faces()[faceOfProc[ip][faceI]].size();
                                integer nodeI =
                                    originMesh_.faces()[faceOfProc[ip][faceI]][n - 1 - i];
                                faces_[faceI][i] = color[ip][nodeI];
                            }
                        } else {
                            for (integer i = 0;
                                 i < originMesh_.faces()[faceOfProc[ip][faceI]].size(); ++i) {
                                integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                                faces_[faceI][i] = color[ip][nodeI];
                            }
                        }
                    } else {
                        for (integer i = 0; i < originMesh_.faces()[faceOfProc[ip][faceI]].size();
                             ++i) {
                            integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                            faces_[faceI][i] = color[ip][nodeI];
                        }
                    }
                }
                pointZones_[0].setFirstIndex(0);
                pointZones_[0].setLastIndex(pzSize[0] - 1);
                for (integer zoneI = 1; zoneI < pointZones_.size(); zoneI++) {
                    integer firstIndex = pointZones_[zoneI - 1].lastIndex() + 1;
                    pointZones_[zoneI].setFirstIndex(firstIndex);
                    pointZones_[zoneI].setLastIndex(firstIndex + pzSize[zoneI] - 1);
                }
            } else {
                localNode.resize(3 * size[ip]);
                faceNode.resize(nodes);
                nodes = 0;

                for (integer i = 0; i < originMesh_.points().size(); ++i) {
                    integer id = color[ip][i];
                    if (id != -1) {
                        localNode[3 * id] = originMesh_.points()[i].x();
                        localNode[3 * id + 1] = originMesh_.points()[i].y();
                        localNode[3 * id + 2] = originMesh_.points()[i].z();
                    }
                }
                for (integer faceI = 0; faceI < faceOfProc[ip].size(); ++faceI) {
                    integer ip1 =
                        partOfProcs_[originMesh_.faces()[faceOfProc[ip][faceI]].leftCell()];
                    if (originMesh_.faces()[faceOfProc[ip][faceI]].rightCell() != -1 &&
                        ip1 !=
                            partOfProcs_[originMesh_.faces()[faceOfProc[ip][faceI]].rightCell()]) {
                        if (ip1 != ip) {
                            for (integer i = 0;
                                 i < originMesh_.faces()[faceOfProc[ip][faceI]].size(); ++i) {
                                integer n = originMesh_.faces()[faceOfProc[ip][faceI]].size();
                                integer nodeI =
                                    originMesh_.faces()[faceOfProc[ip][faceI]][n - 1 - i];
                                faceNode[nodes++] = color[ip][nodeI];
                            }
                        } else {
                            for (integer i = 0;
                                 i < originMesh_.faces()[faceOfProc[ip][faceI]].size(); ++i) {
                                integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                                faceNode[nodes++] = color[ip][nodeI];
                            }
                        }
                    } else {
                        for (integer i = 0; i < originMesh_.faces()[faceOfProc[ip][faceI]].size();
                             ++i) {
                            integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                            faceNode[nodes++] = color[ip][nodeI];
                        }
                    }
                }
                HurMPI::send(&size[ip], 1, feature<integer>::MPIType, ip, 0, HurMPI::getComm());
                HurMPI::send(&nodes, 1, feature<integer>::MPIType, ip, 1, HurMPI::getComm());
                HurMPI::sendList(localNode, ip, 2, HurMPI::getComm());
                HurMPI::sendList(faceNode, ip, 3, HurMPI::getComm());
                HurMPI::sendList(pzSize, ip, 4, HurMPI::getComm());
            }
        }
        if (ip != HurMPI::masterNo()) {
            if (HurMPI::isThisProc(ip)) {
                integer nodes;
                integer sizeN;
                HurMPI::Status status;
                HurMPI::recv(&sizeN, 1, feature<integer>::MPIType, HurMPI::masterNo(), 0,
                             HurMPI::getComm(), &status);
                HurMPI::recv(&nodes, 1, feature<integer>::MPIType, HurMPI::masterNo(), 1,
                             HurMPI::getComm(), &status);
                integerArray fN(nodes);
                realList lN(3 * sizeN);

                HurMPI::recvList(lN, HurMPI::masterNo(), 2, HurMPI::getComm(), &status);
                HurMPI::recvList(fN, HurMPI::masterNo(), 3, HurMPI::getComm(), &status);
                HurMPI::recvList(pzSize, HurMPI::masterNo(), 4, HurMPI::getComm(), &status);

                integer count = 0;
                for (integer fI = 0; fI < faces_.size(); ++fI) {
                    for (integer i = 0; i < faces_[fI].size(); ++i) {
                        faces_[fI][i] = fN[count++];
                    }
                }
                nodes_.resize(lN.size() / 3);
                for (integer i = 0; i < nodes_.size(); ++i) {
                    nodes_[i].x() = lN(3 * i);
                    nodes_[i].y() = lN(3 * i + 1);
                    nodes_[i].z() = lN(3 * i + 2);
                }

                pointZones_[0].setFirstIndex(0);
                pointZones_[0].setLastIndex(pzSize[0] - 1);
                for (integer zoneI = 1; zoneI < pointZones_.size(); zoneI++) {
                    integer firstIndex = pointZones_[zoneI - 1].lastIndex() + 1;
                    pointZones_[zoneI].setFirstIndex(firstIndex);
                    pointZones_[zoneI].setLastIndex(firstIndex + pzSize[zoneI] - 1);
                }
            }
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::pointScatterFix(const integerArrayArray &faceOfProc,
                                                            List<List<Vector2D<integer>>> &mapOP,
                                                            integerArray &sizeLP) {
    integerArray color;
    for (integer ip = 0; ip < nparts_; ip++) {
        integer nodes = 0;
        realList localNode;
        integerArray faceNode;
        integerArray pzSize(pointZones_.size(), Zero);
        if (HurMPI::master()) {
            if (color.size() == 0) {
                color.resize(originMesh_.points().size());
            }
            color = -1;

            for (integer faceI = 0; faceI < faceOfProc[ip].size(); faceI++) {
                nodes += originMesh_.faces()[faceOfProc[ip][faceI]].size();
            }

            for (integer zoneI = 0; zoneI < originMesh_.pointZones().size(); zoneI++) {
                integer firstIndex = originMesh_.pointZones()[zoneI].firstIndex();
                integer lastIndex = originMesh_.pointZones()[zoneI].lastIndex();
                for (integer faceI = 0; faceI < faceOfProc[ip].size(); faceI++) {
                    for (integer i = 0; i < originMesh_.faces()[faceOfProc[ip][faceI]].size();
                         i++) {
                        integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                        if (nodeI >= firstIndex && nodeI <= lastIndex) {
                            if (color[nodeI] == -1) {
                                mapOP[nodeI].append(integerVector2D(ip, sizeLP[ip]));
                                color[nodeI] = sizeLP[ip]++;
                                pzSize[zoneI]++;
                            }
                        }
                    }
                }
            }

            if (ip == HurMPI::masterNo()) {
                nodes_.resize(sizeLP[ip]);
                for (integer i = 0; i < originMesh_.points().size(); ++i) {
                    integer id = color[i];
                    if (id != -1) {
                        nodes_[id] = originMesh_.points()[i];
                    }
                }

                for (integer faceI = 0; faceI < faceOfProc[ip].size(); ++faceI) {
                    if (faces_[faceI].isCutFace()) {
                        integer ip1 =
                            partOfProcs_[originMesh_.faces()[faceOfProc[ip][faceI]].leftCell()];
                        if (ip1 != ip) {
                            for (integer i = 0;
                                 i < originMesh_.faces()[faceOfProc[ip][faceI]].size(); ++i) {
                                integer n = originMesh_.faces()[faceOfProc[ip][faceI]].size();
                                integer nodeI =
                                    originMesh_.faces()[faceOfProc[ip][faceI]][n - 1 - i];
                                faces_[faceI][i] = color[nodeI];
                            }
                        } else {
                            for (integer i = 0;
                                 i < originMesh_.faces()[faceOfProc[ip][faceI]].size(); ++i) {
                                integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                                faces_[faceI][i] = color[nodeI];
                            }
                        }
                    } else {
                        for (integer i = 0; i < originMesh_.faces()[faceOfProc[ip][faceI]].size();
                             ++i) {
                            integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                            faces_[faceI][i] = color[nodeI];
                        }
                    }
                }
                pointZones_[0].setFirstIndex(0);
                pointZones_[0].setLastIndex(pzSize[0] - 1);
                for (integer zoneI = 1; zoneI < pointZones_.size(); zoneI++) {
                    integer firstIndex = pointZones_[zoneI - 1].lastIndex() + 1;
                    pointZones_[zoneI].setFirstIndex(firstIndex);
                    pointZones_[zoneI].setLastIndex(firstIndex + pzSize[zoneI] - 1);
                }
            } else {
                localNode.resize(3 * sizeLP[ip]);
                faceNode.resize(nodes);
                nodes = 0;

                for (integer i = 0; i < originMesh_.points().size(); ++i) {
                    integer id = color[i];
                    if (id != -1) {
                        localNode[3 * id] = originMesh_.points()[i].x();
                        localNode[3 * id + 1] = originMesh_.points()[i].y();
                        localNode[3 * id + 2] = originMesh_.points()[i].z();
                    }
                }
                for (integer faceI = 0; faceI < faceOfProc[ip].size(); ++faceI) {
                    //if (faces_[faceI].isCutFace())//mark error
                    integer ip1 =
                        partOfProcs_[originMesh_.faces()[faceOfProc[ip][faceI]].leftCell()];
                    if (originMesh_.faces()[faceOfProc[ip][faceI]].rightCell() != -1 &&
                        ip1 !=
                            partOfProcs_[originMesh_.faces()[faceOfProc[ip][faceI]].rightCell()]) {
                        if (ip1 != ip) {
                            for (integer i = 0;
                                 i < originMesh_.faces()[faceOfProc[ip][faceI]].size(); ++i) {
                                integer n = originMesh_.faces()[faceOfProc[ip][faceI]].size();
                                integer nodeI =
                                    originMesh_.faces()[faceOfProc[ip][faceI]][n - 1 - i];
                                faceNode[nodes++] = color[nodeI];
                            }
                        } else {
                            for (integer i = 0;
                                 i < originMesh_.faces()[faceOfProc[ip][faceI]].size(); ++i) {
                                integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                                faceNode[nodes++] = color[nodeI];
                            }
                        }
                    } else {
                        for (integer i = 0; i < originMesh_.faces()[faceOfProc[ip][faceI]].size();
                             ++i) {
                            integer nodeI = originMesh_.faces()[faceOfProc[ip][faceI]][i];
                            faceNode[nodes++] = color[nodeI];
                        }
                    }
                }
                HurMPI::send(&sizeLP[ip], 1, feature<integer>::MPIType, ip, 0, HurMPI::getComm());
                HurMPI::send(&nodes, 1, feature<integer>::MPIType, ip, 1, HurMPI::getComm());
                HurMPI::sendList(localNode, ip, 2, HurMPI::getComm());
                HurMPI::sendList(faceNode, ip, 3, HurMPI::getComm());
                HurMPI::sendList(pzSize, ip, 4, HurMPI::getComm());
            }
        }

        if (ip != HurMPI::masterNo()) {
            if (HurMPI::isThisProc(ip)) {
                integer nodes;
                integer sizeN;
                HurMPI::Status status;
                HurMPI::recv(&sizeN, 1, feature<integer>::MPIType, HurMPI::masterNo(), 0,
                             HurMPI::getComm(), &status);
                HurMPI::recv(&nodes, 1, feature<integer>::MPIType, HurMPI::masterNo(), 1,
                             HurMPI::getComm(), &status);
                integerArray fN(nodes);
                realList lN(3 * sizeN);

                HurMPI::recvList(lN, HurMPI::masterNo(), 2, HurMPI::getComm(), &status);
                HurMPI::recvList(fN, HurMPI::masterNo(), 3, HurMPI::getComm(), &status);
                HurMPI::recvList(pzSize, HurMPI::masterNo(), 4, HurMPI::getComm(), &status);

                integer count = 0;
                for (integer fI = 0; fI < faces_.size(); ++fI) {
                    for (integer i = 0; i < faces_[fI].size(); ++i) {
                        faces_[fI][i] = fN[count++];
                    }
                }
                nodes_.resize(lN.size() / 3);
                for (integer i = 0; i < nodes_.size(); ++i) {
                    nodes_[i].x() = lN(3 * i);
                    nodes_[i].y() = lN(3 * i + 1);
                    nodes_[i].z() = lN(3 * i + 2);
                }

                pointZones_[0].setFirstIndex(0);
                pointZones_[0].setLastIndex(pzSize[0] - 1);
                for (integer zoneI = 1; zoneI < pointZones_.size(); zoneI++) {
                    integer firstIndex = pointZones_[zoneI - 1].lastIndex() + 1;
                    pointZones_[zoneI].setFirstIndex(firstIndex);
                    pointZones_[zoneI].setLastIndex(firstIndex + pzSize[zoneI] - 1);
                }
            }
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::pointShareFix(
    const integerArrayArray &faceOfProc, integerArray &sizeLP,
    const List<List<Vector2D<integer>>> &mapOP, const integer nFace) {
    integerArrayArray cZone(nparts_);
    integerArrayArray globalColor(nparts_, OpenHurricane::Zero);
    shareFix(faceShareZones_, faceOfProc, sizeLP, mapOP, nFace, globalColor, cZone, true);

    List<List<std::map<integer, integer>>> mapColor(nparts_);
    integerArray mapSize(nparts_, OpenHurricane::Zero);
    integerArray sBufId;
    integerArray sBufMap;
    integerArray dispI(nparts_);
    integerArray dispM(nparts_);
    integerArray sCountI(nparts_, OpenHurricane::Zero);
    integerArray sCountM(nparts_, OpenHurricane::Zero);
    integerArray info(2 * nparts_ + 2);
    if (HurMPI::master()) {
        for (integer ip = 0; ip < nparts_; ip++) {
            mapColor[ip].resize(sizeLP[ip]);
        }
        pointMapFix(mapColor, mapSize, faceOfProc, mapOP);

        for (integer ip = 0; ip < nparts_; ip++) {
            for (integer i = 0; i < sizeLP[ip]; i++) {
                if (mapColor[ip][i].size() != 0) {
                    mapSize[ip]++;
                }
            }
        }
        dispI[0] = 0;
        sCountI[0] = 2 * mapSize[0];
        integer sum = sCountI[0];
        for (integer ip = 1; ip < nparts_; ip++) {
            dispI[ip] = dispI[ip - 1] + sCountI[ip - 1];
            sCountI[ip] = 2 * mapSize[ip];
            sum += sCountI[ip];
        }
        sBufId.resize(sum);
        sum = 0;
        integer count = 0;
        for (integer ip = 0; ip < nparts_; ip++) {
            for (integer i = 0; i < sizeLP[ip]; i++) {
                if (mapColor[ip][i].size() != 0) {
                    sCountM[ip] += integer(2 * mapColor[ip][i].size());
                    sBufId[count++] = integer(mapColor[ip][i].size());
                    sBufId[count++] = i;
                }
            }
        }
        dispM[0] = 0;
        for (integer ip = 1; ip < nparts_; ip++) {
            dispM[ip] = dispM[ip - 1] + sCountM[ip - 1];
        }
        for (integer ip = 0; ip < nparts_; ip++) //mark
        {
            sum += sCountM[ip];
        }
        sBufMap.resize(sum);
        count = 0;
        std::map<integer, integer>::const_iterator iter;
        for (integer ip = 0; ip < nparts_; ip++) {
            for (integer i = 0; i < sizeLP[ip]; i++) {
                if (mapColor[ip][i].size() != 0) {
                    for (iter = mapColor[ip][i].begin(); iter != mapColor[ip][i].end(); iter++) {
                        sBufMap[count++] = iter->first;
                        sBufMap[count++] = iter->second;
                    }
                }
            }
        } //mark

        for (integer ip = 0; ip < nparts_; ip++) {
            info[ip] = sCountI[ip];
            info[nparts_ + ip] = sCountM[ip];
        }
        info[2 * nparts_] = *std::max_element(sCountI.data(), sCountI.data() + nparts_);
        info[2 * nparts_ + 1] = *std::max_element(sCountM.data(), sCountM.data() + nparts_);
    }
    HurMPI::bcastList(info, HurMPI::masterNo(), HurMPI::getComm());
    HurMPI::bcastList(sCountI, HurMPI::masterNo(), HurMPI::getComm());
    HurMPI::bcastList(sCountM, HurMPI::masterNo(), HurMPI::getComm());
    integerArray rBufId(info[2 * nparts_]);
    integerArray rBufMap(info[2 * nparts_ + 1]);
    for (integer ip = 0; ip < nparts_; ip++) {
        sCountI[ip] = info[ip];
        sCountM[ip] = info[nparts_ + ip];
    }
    HurMPI::scatterv(sBufId.data(), sCountI.data(), dispI.data(), feature<integer>::MPIType,
                     rBufId.data(), sCountI[HurMPI::getProcRank()], feature<integer>::MPIType,
                     HurMPI::masterNo(), HurMPI::getComm());
    HurMPI::scatterv(sBufMap.data(), sCountM.data(), dispM.data(), feature<integer>::MPIType,
                     rBufMap.data(), sCountM[HurMPI::getProcRank()], feature<integer>::MPIType,
                     HurMPI::masterNo(), HurMPI::getComm());

    integer ip = HurMPI::getProcRank();
    integer count = 0;
    for (integer i = 0; i < sCountI[ip] / 2; i++) {
        integerVectorList pairList(rBufId[2 * i]);
        for (integer p = 0; p < rBufId[2 * i]; p++) {
            integer ip = rBufMap[count++];
            integer id = rBufMap[count++];
            integerVector pair(ip, id, -1);
            pairList[p] = pair;
        }
        pointPairMaps_.emplace(rBufId[2 * i + 1], pairList);
    }
}

void OpenHurricane::decomposingMeshByMetis::share(
    cutZoneList &cutZones_, const integerArrayArray &faceOfProc, const integerArray &size,
    const integerArrayArray &color, const integer total, integerArrayArray &globalColor,
    integerArrayArray &cZone, const bool pShare) {
    for (integer ip = 0; ip < nparts_; ip++) {
        integer n = 0;

        for (integer i = 0; i < cutZones_.size(); ++i) {
            if (cutZones_[i].sendProc() == ip) {
                n++;
            }
        }
        cZone[ip].resize(n);
        n = 0;
        for (integer i = 0; i < cutZones_.size(); ++i) {
            if (cutZones_[i].sendProc() == ip) {
                cZone[ip][n++] = i;
            }
        }
    }
    if (HurMPI::master()) {
        integerArray mark(total, 0);
        integerArray count(nparts_, OpenHurricane::Zero);
        for (integer ip = 0; ip < nparts_; ip++) {
            globalColor[ip].resize(size[ip]);
            globalColor[ip] = -3;
        }
        for (integer ip = 0; ip < nparts_; ip++) {

            for (integer i = 0; i < cZone[ip].size(); ++i) {
                integer index = cZone[ip][i];
                integer recvProc = cutZones_[index].receivProc();

                if (pShare) {
                    if (ip < recvProc) {
                        for (integer fI = 0; fI < cutZones_[index].sor().size(); ++fI) {
                            integer localFace = cutZones_[index].sor()[fI];
                            integer globalFace = faceOfProc[ip][localFace];

                            for (integer nI = 0; nI < originMesh_.faces()[globalFace].size();
                                 ++nI) {
                                integer globalNode = originMesh_.faces()[globalFace][nI];
                                //integer localNode = color[ip][globalNode];
                                integer localNode = color[recvProc][globalNode];

                                //mark[globalNode] = -1;
                                globalColor[recvProc][localNode] = -2;
                            }
                        }
                    }
                } else {
                    if (ip > recvProc) {
                        for (integer fI = 0; fI < cutZones_[index].sor().size(); ++fI) {
                            integer localFace = cutZones_[index].sor()[fI];
                            integer globalFace = faceOfProc[ip][localFace];

                            globalColor[ip][localFace] = -1;
                        }
                    }
                }
            }

            if (pShare) {
                for (integer i = 0; i < size[ip]; i++) {
                    //if (globalColor[ip][i] != -1)
                    if (globalColor[ip][i] != -2) {
                        globalColor[ip][i] = count[ip]++;
                    }
                }
            } else {
                for (integer i = 0; i < size[ip]; i++) {
                    //if (globalColor[ip][i] != -1)
                    if (globalColor[ip][i] != -1) {
                        globalColor[ip][i] = count[ip]++;
                    }
                }
            }
        }
        for (integer ip = 1; ip < nparts_; ip++) {
            if (pShare) {
                pointOffSet_[ip] = pointOffSet_[ip - 1] + count[ip - 1];
            } else {
                faceOffSet_[ip] = faceOffSet_[ip - 1] + count[ip - 1];
            }
        }
    }
    if (pShare) {
        HurMPI::bcastList(pointOffSet_, HurMPI::masterNo(), HurMPI::getComm());
    } else {
        HurMPI::bcastList(faceOffSet_, HurMPI::masterNo(), HurMPI::getComm());
    }
}

void OpenHurricane::decomposingMeshByMetis::shareFix(
    cutZoneList &cutZones, const integerArrayArray &faceOfProc, const integerArray &sizeLP,
    const List<List<Vector2D<integer>>> &mapOP, const integer nFace, integerArrayArray &globalColor,
    integerArrayArray &cZone, const bool pShare) {
    for (integer ip = 0; ip < nparts_; ip++) {
        integer n = 0;

        for (integer i = 0; i < cutZones.size(); ++i) {
            if (cutZones[i].sendProc() == ip) {
                n++;
            }
        }
        cZone[ip].resize(n);
        n = 0;
        for (integer i = 0; i < cutZones.size(); ++i) {
            if (cutZones[i].sendProc() == ip) {
                cZone[ip][n++] = i;
            }
        }
    }
    if (HurMPI::master()) {
        integerArray mark(nFace, 0);
        integerArray count(nparts_, OpenHurricane::Zero);
        for (integer ip = 0; ip < nparts_; ip++) {
            globalColor[ip].resize(sizeLP[ip]);
            globalColor[ip] = -3;
        }
        for (integer ip = 0; ip < nparts_; ip++) {

            for (integer i = 0; i < cZone[ip].size(); ++i) {
                integer index = cZone[ip][i];
                integer recvProc = cutZones[index].receivProc();

                if (pShare) {
                    if (ip < recvProc) {
                        for (integer fI = 0; fI < cutZones[index].sor().size(); ++fI) {
                            integer localFace = cutZones[index].sor()[fI];
                            integer globalFace = faceOfProc[ip][localFace];

                            for (integer nI = 0; nI < originMesh_.faces()[globalFace].size();
                                 ++nI) {
                                integer globalNode = originMesh_.faces()[globalFace][nI];

                                integer localNode = -1;
                                for (integer kk = 0; kk < mapOP[globalNode].size(); ++kk) {
                                    if (mapOP[globalNode][kk].x() == recvProc) {
                                        localNode = mapOP[globalNode][kk].y();
                                        break;
                                    }
                                }
                                if (localNode == -1) {
                                    errorAbortStr(("Cannot find the local index of point for ip: " +
                                                   toString(recvProc)));
                                }
                                //mark[globalNode] = -1;
                                globalColor[recvProc][localNode] = -2;
                            }
                        }
                    }
                } else {
                    if (ip > recvProc) {
                        for (integer fI = 0; fI < cutZones[index].sor().size(); ++fI) {
                            integer localFace = cutZones[index].sor()[fI];
                            integer globalFace = faceOfProc[ip][localFace];

                            globalColor[ip][localFace] = -1;
                        }
                    }
                }
            }

            if (pShare) {
                for (integer i = 0; i < sizeLP[ip]; i++) {
                    //if (globalColor[ip][i] != -1)
                    if (globalColor[ip][i] != -2) {
                        globalColor[ip][i] = count[ip]++;
                    }
                }
            } else {
                for (integer i = 0; i < sizeLP[ip]; i++) {
                    //if (globalColor[ip][i] != -1)
                    if (globalColor[ip][i] != -1) {
                        globalColor[ip][i] = count[ip]++;
                    }
                }
            }
        }
        for (integer ip = 1; ip < nparts_; ip++) {
            if (pShare) {
                pointOffSet_[ip] = pointOffSet_[ip - 1] + count[ip - 1];
            } else {
                faceOffSet_[ip] = faceOffSet_[ip - 1] + count[ip - 1];
            }
        }
    }
    if (pShare) {
        HurMPI::bcastList(pointOffSet_, HurMPI::masterNo(), HurMPI::getComm());
    } else {
        HurMPI::bcastList(faceOffSet_, HurMPI::masterNo(), HurMPI::getComm());
    }
}

void OpenHurricane::decomposingMeshByMetis::pointShare(const integerArrayArray &faceOfProc,
                                                       integerArray &size,
                                                       const integerArrayArray &color,
                                                       const integer total) {
    integerArrayArray cZone(nparts_);
    integerArrayArray globalColor(nparts_, OpenHurricane::Zero);
    share(faceShareZones_, faceOfProc, size, color, total, globalColor, cZone, true);

    List<List<std::map<integer, integer>>> mapColor(nparts_);
    integerArray mapSize(nparts_, OpenHurricane::Zero);
    //integerArray paraSize(nparts_, OpenHurricane::Zero);
    integerArray sBufId;
    integerArray sBufMap;
    integerArray dispI(nparts_);
    integerArray dispM(nparts_);
    integerArray sCountI(nparts_, OpenHurricane::Zero);
    integerArray sCountM(nparts_, OpenHurricane::Zero);
    integerArray info(2 * nparts_ + 2);
    if (HurMPI::master()) {
        for (integer ip = 0; ip < nparts_; ip++) {
            mapColor[ip].resize(size[ip]);
        }
        pointMap(mapColor, mapSize, faceOfProc, color);

        for (integer ip = 0; ip < nparts_; ip++) {
            for (integer i = 0; i < size[ip]; i++) {
                if (mapColor[ip][i].size() != 0) {
                    //paraSize[ip]+= mapColor[ip][i].size();
                    mapSize[ip]++;
                }
            }
        }
        dispI[0] = 0;
        sCountI[0] = 2 * mapSize[0];
        integer sum = sCountI[0];
        for (integer ip = 1; ip < nparts_; ip++) {
            dispI[ip] = dispI[ip - 1] + sCountI[ip - 1];
            sCountI[ip] = 2 * mapSize[ip];
            sum += sCountI[ip];
        }
        sBufId.resize(sum);
        sum = 0;
        integer count = 0;
        for (integer ip = 0; ip < nparts_; ip++) {
            for (integer i = 0; i < size[ip]; i++) {
                if (mapColor[ip][i].size() != 0) {
                    sCountM[ip] += integer(2 * mapColor[ip][i].size());
                    sBufId[count++] = integer(mapColor[ip][i].size());
                    sBufId[count++] = i;
                }
            }
        }
        dispM[0] = 0;
        for (integer ip = 1; ip < nparts_; ip++) {
            dispM[ip] = dispM[ip - 1] + sCountM[ip - 1];
        }
        for (integer ip = 0; ip < nparts_; ip++) //mark
        {
            sum += sCountM[ip];
        }
        sBufMap.resize(sum);
        count = 0;
        std::map<integer, integer>::const_iterator iter;
        for (integer ip = 0; ip < nparts_; ip++) {
            for (integer i = 0; i < size[ip]; i++) {
                if (mapColor[ip][i].size() != 0) {
                    for (iter = mapColor[ip][i].begin(); iter != mapColor[ip][i].end(); iter++) {
                        sBufMap[count++] = iter->first;
                        //sBufMap[count++] = globalColor[iter->first][iter->second];
                        sBufMap[count++] = iter->second;
                    }
                }
            }
        } //mark

        for (integer ip = 0; ip < nparts_; ip++) {
            info[ip] = sCountI[ip];
            info[nparts_ + ip] = sCountM[ip];
        }
        info[2 * nparts_] = *std::max_element(sCountI.data(), sCountI.data() + nparts_);
        info[2 * nparts_ + 1] = *std::max_element(sCountM.data(), sCountM.data() + nparts_);
    }
    HurMPI::bcastList(info, HurMPI::masterNo(), HurMPI::getComm());
    HurMPI::bcastList(sCountI, HurMPI::masterNo(), HurMPI::getComm());
    HurMPI::bcastList(sCountM, HurMPI::masterNo(), HurMPI::getComm());
    integerArray rBufId(info[2 * nparts_]);
    integerArray rBufMap(info[2 * nparts_ + 1]);
    for (integer ip = 0; ip < nparts_; ip++) {
        sCountI[ip] = info[ip];
        sCountM[ip] = info[nparts_ + ip];
    }
    HurMPI::scatterv(sBufId.data(), sCountI.data(), dispI.data(), feature<integer>::MPIType,
                     rBufId.data(), sCountI[HurMPI::getProcRank()], feature<integer>::MPIType,
                     HurMPI::masterNo(), HurMPI::getComm());
    HurMPI::scatterv(sBufMap.data(), sCountM.data(), dispM.data(), feature<integer>::MPIType,
                     rBufMap.data(), sCountM[HurMPI::getProcRank()], feature<integer>::MPIType,
                     HurMPI::masterNo(), HurMPI::getComm());

    integer ip = HurMPI::getProcRank();
    integer count = 0;
    for (integer i = 0; i < sCountI[ip] / 2; i++) {
        integerVectorList pairList(rBufId[2 * i]);
        for (integer p = 0; p < rBufId[2 * i]; p++) {
            integer ip = rBufMap[count++];
            integer id = rBufMap[count++];
            integerVector pair(ip, id, -1);
            pairList[p] = pair;
        }
        pointPairMaps_.emplace(rBufId[2 * i + 1], pairList);
    }
}

void OpenHurricane::decomposingMeshByMetis::pointMap(
    List<List<std::map<integer, integer>>> &mapColor, integerArray &mapSize,
    const integerArrayArray &faceOfProc, const integerArrayArray &color) {
    List<std::map<integer, integer>> globalMap(originMesh_.points().size());
    for (integer zoneI = 0; zoneI < faceShareZones_.size(); ++zoneI) {
        integer sendP = faceShareZones_[zoneI].sendProc();
        for (integer fI = 0; fI < faceShareZones_[zoneI].sor().size(); ++fI) {

            integer localFace1 = faceShareZones_[zoneI].sor()[fI];
            //integer localFace2 = faceShareZones_[zoneI].des()[fI];
            integer globalFace1 = faceOfProc[sendP][localFace1];
            //integer globalFace2 = faceOfProc[recvP][localFace2];

            for (integer nI = 0; nI < originMesh_.faces()[globalFace1].size(); ++nI) {
                integer globalNode1 = originMesh_.faces()[globalFace1][nI];
                //integer globalNode2 = originMesh_.faces()[globalFace2][nI];
                integer localNode1 = color[sendP][globalNode1];
                //integer localNode2 = color[recvP][globalNode2];
                std::map<integer, integer>::const_iterator iter =
                    globalMap[globalNode1].find(sendP);
                if (iter == globalMap[globalNode1].end()) {
                    globalMap[globalNode1][sendP] = localNode1;
                }
            }
        }
    }

    for (integer globalNode = 0; globalNode < globalMap.size(); globalNode++) {
        for (std::map<integer, integer>::const_iterator iter1 = globalMap[globalNode].begin();
             iter1 != globalMap[globalNode].end(); iter1++) {
            integer ip1 = iter1->first;
            integer localNode1 = iter1->second;

            for (std::map<integer, integer>::const_iterator iter2 = globalMap[globalNode].begin();
                 iter2 != globalMap[globalNode].end(); iter2++) {
                integer ip2 = iter2->first;

                if (ip1 != ip2) {
                    integer localNode2 = iter2->second;
                    mapColor[ip1][localNode1][ip2] = localNode2;
                }
            }
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::pointMapFix(
    List<List<std::map<integer, integer>>> &mapColor, integerArray &mapSize,
    const integerArrayArray &faceOfProc, const List<List<Vector2D<integer>>> &mapOP) {
    List<std::map<integer, integer>> globalMap(originMesh_.points().size());
    for (integer zoneI = 0; zoneI < faceShareZones_.size(); ++zoneI) {
        integer sendP = faceShareZones_[zoneI].sendProc();
        for (integer fI = 0; fI < faceShareZones_[zoneI].sor().size(); ++fI) {
            integer localFace1 = faceShareZones_[zoneI].sor()[fI];
            integer globalFace1 = faceOfProc[sendP][localFace1];

            for (integer nI = 0; nI < originMesh_.faces()[globalFace1].size(); ++nI) {
                integer globalNode1 = originMesh_.faces()[globalFace1][nI];

                integer localNode1 = -1;
                for (integer kk = 0; kk < mapOP[globalNode1].size(); ++kk) {
                    if (mapOP[globalNode1][kk].x() == sendP) {
                        localNode1 = mapOP[globalNode1][kk].y();
                        break;
                    }
                }
                if (localNode1 == -1) {
                    errorAbortStr(
                        ("Cannot find the local index of point for ip: " + toString(sendP)));
                }

                std::map<integer, integer>::const_iterator iter =
                    globalMap[globalNode1].find(sendP);
                if (iter == globalMap[globalNode1].end()) {
                    globalMap[globalNode1][sendP] = localNode1;
                }
            }
        }
    }

    for (integer globalNode = 0; globalNode < globalMap.size(); globalNode++) {
        for (std::map<integer, integer>::const_iterator iter1 = globalMap[globalNode].begin();
             iter1 != globalMap[globalNode].end(); iter1++) {
            integer ip1 = iter1->first;
            integer localNode1 = iter1->second;

            for (std::map<integer, integer>::const_iterator iter2 = globalMap[globalNode].begin();
                 iter2 != globalMap[globalNode].end(); iter2++) {
                integer ip2 = iter2->first;

                if (ip1 != ip2) {
                    integer localNode2 = iter2->second;
                    mapColor[ip1][localNode1][ip2] = localNode2;
                }
            }
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::faceShare(const integerArrayArray &faceOfProc,
                                                      const integerArray &size,
                                                      const integer total) {
    integerArrayArray cZone(nparts_);
    integerArrayArray globalColor(nparts_);
    integerArrayArray nullColor;
    share(cutZones_, faceOfProc, size, nullColor, total, globalColor, cZone, false);
    if (HurMPI::master()) {
        //List<std::map<integer, integerVector>> map(nparts_);
        List<std::map<integer, integerArray>> map(nparts_);

        for (integer zoneI = 0; zoneI < cutZones_.size(); ++zoneI) {
            integer sendP = cutZones_[zoneI].sendProc();

            if (HurMPI::masterNo() == sendP) {
                faceMap(zoneI, globalColor, faceOfProc, facePairMaps_);
            } else {
                faceMap(zoneI, globalColor, faceOfProc, map[sendP]);
            }
        }

        for (integer ip = 0; ip < nparts_; ip++) {
            if (ip != HurMPI::masterNo()) {
                integer n = 0;
                integerArray buf(5 * static_cast<integer>(map[ip].size()));
                std::map<integer, integerArray>::const_iterator iter;
                for (iter = map[ip].begin(); iter != map[ip].end(); iter++) {
                    buf[n++] = iter->first;
                    buf[n++] = (iter->second)[0];
                    buf[n++] = (iter->second)[1];
                    buf[n++] = (iter->second)[2];
                    buf[n++] = (iter->second)[3];
                }
                HurMPI::send(&n, 1, feature<integer>::MPIType, ip, 1, HurMPI::getComm());
                HurMPI::sendList(buf, ip, 2, HurMPI::getComm());
            }
        }
    } else {
        integer bufSize = -1;
        HurMPI::Status status;
        HurMPI::recv(&bufSize, 1, feature<integer>::MPIType, HurMPI::masterNo(), 1,
                     HurMPI::getComm(), &status);
        integerArray buf(bufSize);
        HurMPI::recvList(buf, HurMPI::masterNo(), 2, HurMPI::getComm(), &status);

        integer n = 0;
        for (integer i = 0; i < bufSize / 5; i++) {
            //integerVector pair;
            integer id = buf[n++];
            integerArray pair(4);
            pair[0] = buf[n++];
            pair[1] = buf[n++];
            pair[2] = buf[n++];
            pair[3] = buf[n++];
            facePairMaps_[id] = pair;
        }
    }
}

void OpenHurricane::decomposingMeshByMetis::faceMap(
    const integer zoneI, const integerArrayArray &globalColor, const integerArrayArray &faceOfProc,
    std::map<integer, OpenHurricane::integerArray> &facePairMaps) {
    integer sendP = cutZones_[zoneI].sendProc();
    integer recvP = cutZones_[zoneI].receivProc();

    for (integer i = 0; i < cutZones_[zoneI].sor().size(); i++) {
        integer index1 = cutZones_[zoneI].sor()[i];
        integer index2 = cutZones_[zoneI].des()[i];
        integer gC = originMesh_.faces()[faceOfProc[sendP][index1]].leftCell();
        //integerVector pair(recvP, globalColor[recvP][index2], -1);
        integerArray pair(4);
        pair[0] = recvP;
        pair[1] = globalColor[recvP][index2];
        pair[2] = index2;
        if (sendP == partOfProcs_[gC]) {
            //pair.z() = color_[originMesh_.faces()[faceOfProc[sendP][index1]].rightCell()];
            pair[3] = color_[originMesh_.faces()[faceOfProc[sendP][index1]].rightCell()];
        } else {
            //pair.z() = color_[gC];
            pair[3] = color_[gC];
        }
        facePairMaps[index1] = pair;
    }
}

bool OpenHurricane::decomposingMeshByMetis::isTwoLayer() {
    for (integer zoneI = 0; zoneI < faceZones_.size(); ++zoneI) {
        if ((!faceZones_[zoneI].isInterior()) || faceZones_[zoneI].isInterface()) {
            integer firstIndex = faceZones_[zoneI].firstIndex();
            integer lastIndex = faceZones_[zoneI].lastIndex();
            for (integer i = firstIndex; i < lastIndex + 1; i++) {
                integer celli = faces_[i].leftCell();
                integer oppsite = cells_[celli].oppositeFace(i);
                if (faces_[oppsite].leftCell() == -1 || faces_[oppsite].rightCell() == -1) {
                    twoGhostLayer_ = false;
                    return false;
                }
            }
        }
    }
    twoGhostLayer_ = true;
    return true;
}

OpenHurricane::integer OpenHurricane::decomposingMeshByMetis::nInteriorFaces() const {
    integer nf = 0;
    for (integer faceZI = 0; faceZI < faceZones_.size(); ++faceZI) {
        if (faceZones_[faceZI].isInterior()) {
            nf += faceZones_[faceZI].size();
        }
    }
    return nf;
}

void OpenHurricane::decomposingMeshByMetis::createFaceShareZone(const integerArray &globalColor) {
    faceShareZones_.resize(cutZones_.size());
    faceShareZones_ = cutZones_;
    if (HurMPI::parRun() && originMesh_.interiorWallFaceMap().size() != 0) {
        integerArrayArray size(nparts_);
        for (integer ip = 0; ip < nparts_; ip++) {
            size[ip].resize(nparts_);
            size[ip] = 0;
        }

        const std::map<integer, integer> &iwMap = originMesh_.interiorWallFaceMap();
        for (std::map<integer, integer>::const_iterator iter = iwMap.begin(); iter != iwMap.end();
             iter++) {
            integer cl = originMesh_.faces()[iter->first].leftCell();
            integer cr = originMesh_.faces()[iter->second].leftCell();

            integer ip1 = partOfProcs_[cl];
            integer ip2 = partOfProcs_[cr];

            if (ip1 != ip2) {
                size[ip1][ip2]++;
                size[ip2][ip1]++;
            }
        }

        integerArrayArrayArray sor(nparts_);
        integerArrayArrayArray des(nparts_);

        for (integer ip1 = 0; ip1 < nparts_; ip1++) {
            sor[ip1].resize(nparts_);
            des[ip1].resize(nparts_);
            for (integer ip2 = 0; ip2 < nparts_; ip2++) {
                if (size[ip1][ip2] != 0) {
                    sor[ip1][ip2].resize(size[ip1][ip2]);
                    des[ip1][ip2].resize(size[ip1][ip2]);
                }
            }
            size[ip1] = 0;
        }

        for (std::map<integer, integer>::const_iterator iter = iwMap.begin(); iter != iwMap.end();
             iter++) {
            integer cl = originMesh_.faces()[iter->first].leftCell();
            integer cr = originMesh_.faces()[iter->second].leftCell();

            integer ip1 = partOfProcs_[cl];
            integer ip2 = partOfProcs_[cr];

            if (ip1 != ip2) {
                sor[ip1][ip2][size[ip1][ip2]] = globalColor[iter->first];
                des[ip1][ip2][size[ip1][ip2]] = globalColor[iter->second];
                sor[ip2][ip1][size[ip2][ip1]] = globalColor[iter->second];
                des[ip2][ip1][size[ip2][ip1]] = globalColor[iter->first];
                size[ip1][ip2]++;
                size[ip2][ip1]++;
            }
        }

        for (integer ip = 0; ip < nparts_; ip++) {
            size[ip] = -1;
        }

        for (integer i = 0; i < faceShareZones_.size(); i++) {
            size[faceShareZones_[i].sendProc()][faceShareZones_[i].receivProc()] = i;
        }

        integer count = 0;
        for (integer ip1 = 0; ip1 < nparts_; ip1++) {
            for (integer ip2 = 0; ip2 < nparts_; ip2++) {
                integer czi = size[ip1][ip2];
                if (czi != -1) {
                    if (sor[ip1][ip2].size() != 0) {
                        faceShareZones_[czi].sor().append(sor[ip1][ip2]);
                        faceShareZones_[czi].des().append(des[ip1][ip2]);
                    }
                } else {
                    if (sor[ip1][ip2].size() != 0) {
                        count++;
                    }
                }
            }
        }

        if (count != 0) {
            cutZoneList appendCZ(count);
            count = 0;

            for (integer ip1 = 0; ip1 < nparts_; ip1++) {
                for (integer ip2 = 0; ip2 < nparts_; ip2++) {
                    integer czi = size[ip1][ip2];
                    if (czi == -1) {
                        if (sor[ip1][ip2].size() != 0) {
                            appendCZ[count].setSendProc(ip1);
                            appendCZ[count].setRecvProc(ip2);
                            appendCZ[count].sor().append(sor[ip1][ip2]);
                            appendCZ[count++].des().append(des[ip1][ip2]);
                        }
                    }
                }
            }
            faceShareZones_.append(appendCZ);
        }
    }
}