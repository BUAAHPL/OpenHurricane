/*!
 * \file geometryMesh.cpp
 * \brief Main subroutines for geometry mesh.
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

#include "geometryMesh.hpp"
#include "FluentMeshRead.hpp"
#include "cellArrays.hpp"
#include "decomposingByMetis.hpp"
#include "globalMesh.hpp"
#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "iteration.hpp"
#include "kNNGrid.hpp"
#include "object.hpp"
#include "originMeshRead.hpp"
#include "relayScatter.hpp"
#include "reordering.hpp"
#include "tecplotFormat.hpp"

OpenHurricane::integer OpenHurricane::geometryMesh::getPeriodicPairSize(const perZoneList &pZL,
                                                                        integerArray &pi) const {
    for (integer i = 0; i < pZL.size(); ++i) {
        if (!pZL[i].isShadow()) {
            if (pi.size() == 0) {
                pi.append(pZL[i].id());
            } else {
                bool found = false;
                for (integer j = 0; j < pi.size(); ++j) {
                    if (pi[j] == pZL[i].id()) {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    pi.append(pZL[i].id());
                }
            }
        }
    }
    return pi.size();
}

void OpenHurricane::geometryMesh::getPeriodicPairList(const perZoneList &pZL,
                                                      periodicPairList &ppl) const {
    integerArray pi;
    getPeriodicPairSize(pZL, pi);

    ppl.resize(periodicPairSize_);

    integerArray nSizeL(HurMPI::getProcSize(), Zero);
    integerArray displs;
    displs.resize(HurMPI::getProcSize(), Zero);
    nSizeL[HurMPI::getProcRank()] = pi.size();
    HurMPI::allGatherList(nSizeL, HurMPI::getComm());
    for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
        displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
    }

    integer ssi = 0;
    for (integer i = 0; i < nSizeL.size(); ++i) {
        ssi += nSizeL[i];
    }
    integerArray allPi(ssi);

    HurMPI::Request request;
    HurMPI::iallGatherv(pi.data(), pi.size(), feature<integer>::MPIType, allPi.data(),
                        nSizeL.data(), displs.data(), feature<integer>::MPIType, HurMPI::getComm(),
                        &request);
    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

    integerArray ppi;
    for (integer i = 0; i < allPi.size(); ++i) {
        if (ppi.size() == 0) {
            ppi.append(allPi[i]);
        } else {
            bool found = false;
            for (integer j = 0; j < ppi.size(); ++j) {
                if (ppi[j] == allPi[i]) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                ppi.append(allPi[i]);
            }
        }
    }

    periodicPairList temPPL(ppi.size());
    for (integer i = 0; i < ppi.size(); ++i) {
        for (integer j = 0; j < pZL.size(); ++j) {
            if (pZL[j].id() == ppi[i] && pZL[j].isSendFromThisProc()) {
                getPeriodicPair(pZL[j], temPPL[i]);
            }
        }
    }

    HurMPI::barrier(HurMPI::getComm());
    for (integer i = 0; i < ppi.size(); ++i) {
        periodicPair::gatherPeriodicPair(temPPL[i], ppl[i], HurMPI::masterNo());
    }
}

void OpenHurricane::geometryMesh::getPeriodicPair(const perZone &pZ, periodicPair &ppl) const {
    if (pZ.isShadow()) {
        return;
    }
    ppl.resize(pZ.sorFace().size());
    const auto pid = pZ.id();
    const auto sid = pZ.pairId();
    const auto sendP = pZ.sendProc();
    const auto recvP = pZ.receivProc();
    ppl.setPeriodicZone(pid);
    ppl.setShadowZone(sid);
    ppl.setType(pZ.type());
    integer whichSZi = -1;
    integer whichPZi = -1;

    for (integer i = 0; i < faceZones_.size(); ++i) {
        if (pid == faceZones_[i].index()) {
            if (whichPZi >= 0) {
                errorAbortStr(("The index of periodic zone: " + toString(pid) + " is duplicated"));
            }
            whichPZi = i;
        }
        if (sid == faceZones_[i].index()) {
            if (whichSZi >= 0) {
                errorAbortStr(("The index of shadow zone: " + toString(sid) + " is duplicated"));
            }
            whichSZi = i;
        }

        if (whichPZi >= 0 && whichSZi >= 0) {
            break;
        }
    }

    if (HurMPI::parRun()) {
        for (integer i = 0; i < pZ.sorFace().size(); ++i) {
            ppl[i].x() = globalMeshInfoPtr_->globalFaceIndeces().toGlobalIndex(
                whichPZi, pZ.sorFace()[i], sendP);
            ppl[i].y() = globalMeshInfoPtr_->globalFaceIndeces().toGlobalIndex(
                whichSZi, pZ.desFace()[i], recvP);
        }
    } else {
        for (integer i = 0; i < pZ.sorFace().size(); ++i) {
            ppl[i].x() = pZ.sorFace()[i];
            ppl[i].y() = pZ.desFace()[i];
        }
    }
}

void OpenHurricane::geometryMesh::getCellOfProc() {
    if (HurMPI::master()) {
        cellOfProc_.resize(HurMPI::getProcSize());
        integerArray nsize(HurMPI::getProcSize(), Zero);
        for (integer i = 0; i < partOfProcs_.size(); ++i) {
            nsize[partOfProcs_[i]]++;
        }

        for (integer ip = 0; ip < cellOfProc_.size(); ++ip) {
            cellOfProc_[ip].resize(nsize[ip]);
        }

        nsize = Zero;
        for (integer i = 0; i < partOfProcs_.size(); ++i) {
            cellOfProc_[partOfProcs_[i]][nsize[partOfProcs_[i]]] = i;
            nsize[partOfProcs_[i]]++;
        }
    }
}

void OpenHurricane::geometryMesh::transferFromDCM(const controller &cont,
                                                  decomposingMeshByMetis &dcM) {
    setBaseMesh(dcM.points().size(), dcM.faces().size(), dcM.nInteriorFaces(), dcM.cells().size(),
                dcM.twoGhostLayer());
    partOfProcs_.transfer(dcM.partOfProcs());
    getCellOfProc();
    periodicPairSize_ = dcM.periodicPairSize();
    points_.transfer(dcM.points());
    faces_.transfer(dcM.faces());
    cells_.transfer(dcM.cells());
    pointZones_.transfer(dcM.pointZones());
    faceZones_.transfer(dcM.faceZones());
    cutZones_.transfer(dcM.cutZones());
    perZones_.transfer(dcM.perZones());
    cellZones_.transfer(dcM.cellZones());
    pairProcessSize_ = cutZones_.size();
    pairPeriodicProcessSize_ = perZones_.size();
    perFaceMaps_ = dcM.perPairMaps();

    point vo = dcM.origin();
    point ax = dcM.axis();

    setOriginAndAxis(vo, ax);

    if (HurMPI::parRun()) {
        globalMeshInfoPtr_.reset(new globalMesh(*this, dcM));
    }
}

hur_nodiscard OpenHurricane::real
OpenHurricane::geometryMesh::meshScaleFactor(const string &meshUnit) const {
    if (meshUnit == "mm") {
        return real(0.001);
    } else if (meshUnit == "cm") {
        return real(0.01);
    } else if (meshUnit == "dm") {
        return real(0.1);
    } else if (meshUnit == "km") {
        return real(1000.0);
    } else if (meshUnit == "um") {
        return real(1.0e-6);
    }

    return real(1.0);
}

OpenHurricane::geometryMesh::geometryMesh(const object &ob)
    : baseMesh(), registerTable(ob), geometryMeshCore(), globalMeshInfoPtr_(nullptr),
      globalFaceZoneL_(), perFaceMaps_(), pairProcessSize_(0), pairPeriodicProcessSize_(0),
      periodicPairSize_(0), cellOrthoPtr_(nullptr), faceCellIntersectionPtr_(nullptr),
      fIntSectFCentrePtr_(nullptr), faceSkewnessWeightPtr_(nullptr), skewFacePtr_(nullptr),
      isSkewFacePtr_(nullptr), LDUFaceMapPtr_(nullptr), processCutPtr_(nullptr),
      processPerPtr_(nullptr), sorCellCentrePtr_(nullptr), sorKNN_(nullptr), tarOfProc_(nullptr),
      cellLoadWeightsPtr_(nullptr) {}

OpenHurricane::geometryMesh::geometryMesh(object &&ob)
    : baseMesh(), registerTable(std::move(ob)), geometryMeshCore(), globalMeshInfoPtr_(nullptr),
      globalFaceZoneL_(), perFaceMaps_(), pairProcessSize_(0), pairPeriodicProcessSize_(0),
      periodicPairSize_(0), cellOrthoPtr_(nullptr), faceCellIntersectionPtr_(nullptr),
      fIntSectFCentrePtr_(nullptr), faceSkewnessWeightPtr_(nullptr), skewFacePtr_(nullptr),
      isSkewFacePtr_(nullptr), LDUFaceMapPtr_(nullptr), processCutPtr_(nullptr),
      processPerPtr_(nullptr), sorCellCentrePtr_(nullptr), sorKNN_(nullptr), tarOfProc_(nullptr),
      cellLoadWeightsPtr_(nullptr) {}

OpenHurricane::geometryMesh::geometryMesh(const object &ob, const controller &cont,
                                          pointField &meshPoints, faceList &meshFaces,
                                          const integer nInternFaces, cellList &meshCells,
                                          pointZoneList &pointZ, faceZoneList &faceZ,
                                          cutZoneList &cutZ, perZoneList &perZ, cellZoneList &cellZ,
                                          const bool twoLayer)
    : baseMesh(cont, meshPoints.size(), meshFaces.size(), nInternFaces, meshCells.size(), twoLayer),
      registerTable(ob), geometryMeshCore(), globalMeshInfoPtr_(nullptr), globalFaceZoneL_(),
      perFaceMaps_(), pairProcessSize_(0), pairPeriodicProcessSize_(0), periodicPairSize_(0),
      cellOrthoPtr_(nullptr), faceCellIntersectionPtr_(nullptr), fIntSectFCentrePtr_(nullptr),
      faceSkewnessWeightPtr_(nullptr), skewFacePtr_(nullptr), isSkewFacePtr_(nullptr),
      LDUFaceMapPtr_(nullptr), processCutPtr_(nullptr), processPerPtr_(nullptr),
      sorCellCentrePtr_(nullptr), sorKNN_(nullptr), tarOfProc_(nullptr),
      cellLoadWeightsPtr_(nullptr) {
    points_.transfer(meshPoints);
    faces_.transfer(meshFaces);
    cells_.transfer(meshCells);
    pointZones_.transfer(pointZ);
    faceZones_.transfer(faceZ);
    cutZones_.transfer(cutZ);
    cellZones_.transfer(cellZ);
    perZones_.transfer(perZ);
    pairProcessSize_ = cutZones_.size();
    pairPeriodicProcessSize_ = perZones_.size();

    setareAllHexOrQuadri();
    calcBashMesh();

    if (HurMPI::parRun()) {
        allMeshCellNumber_ = nCells();
        HurMPI::allReduce(allMeshCellNumber_, MPI_SUM);
    } else {
        allMeshCellNumber_ = nCells();
    }
}

OpenHurricane::geometryMesh::geometryMesh(object &&ob, const controller &cont,
                                          pointField &meshPoints, faceList &meshFaces,
                                          const integer nInternFaces, cellList &meshCells,
                                          pointZoneList &pointZ, faceZoneList &faceZ,
                                          cutZoneList &cutZ, perZoneList &perZ, cellZoneList &cellZ,
                                          const bool twoLayer)
    : baseMesh(cont, meshPoints.size(), meshFaces.size(), nInternFaces, meshCells.size(), twoLayer),
      registerTable(std::move(ob)), geometryMeshCore(), globalMeshInfoPtr_(nullptr),
      globalFaceZoneL_(), perFaceMaps_(), pairProcessSize_(0), pairPeriodicProcessSize_(0),
      periodicPairSize_(0), cellOrthoPtr_(nullptr), faceCellIntersectionPtr_(nullptr),
      fIntSectFCentrePtr_(nullptr), faceSkewnessWeightPtr_(nullptr), skewFacePtr_(nullptr),
      isSkewFacePtr_(nullptr), LDUFaceMapPtr_(nullptr), processCutPtr_(nullptr),
      processPerPtr_(nullptr), sorCellCentrePtr_(nullptr), sorKNN_(nullptr), tarOfProc_(nullptr),
      cellLoadWeightsPtr_(nullptr) {
    points_.transfer(meshPoints);
    faces_.transfer(meshFaces);
    cells_.transfer(meshCells);
    pointZones_.transfer(pointZ);
    faceZones_.transfer(faceZ);
    cutZones_.transfer(cutZ);
    cellZones_.transfer(cellZ);
    perZones_.transfer(perZ);
    pairProcessSize_ = cutZones_.size();
    pairPeriodicProcessSize_ = perZones_.size();

    setareAllHexOrQuadri();
    calcBashMesh();
    if (HurMPI::parRun()) {
        allMeshCellNumber_ = nCells();
        HurMPI::allReduce(allMeshCellNumber_, MPI_SUM);
    } else {
        allMeshCellNumber_ = nCells();
    }
    //setWriteControl();
}

OpenHurricane::geometryMesh::geometryMesh(const object &ob, const controller &cont)
    : baseMesh(cont), registerTable(ob), geometryMeshCore(), globalMeshInfoPtr_(nullptr),
      globalFaceZoneL_(), perFaceMaps_(), pairProcessSize_(0), pairPeriodicProcessSize_(0),
      periodicPairSize_(0), cellOrthoPtr_(nullptr), faceCellIntersectionPtr_(nullptr),
      fIntSectFCentrePtr_(nullptr), faceSkewnessWeightPtr_(nullptr), skewFacePtr_(nullptr),
      isSkewFacePtr_(nullptr), LDUFaceMapPtr_(nullptr), processCutPtr_(nullptr),
      processPerPtr_(nullptr), sorCellCentrePtr_(nullptr), sorKNN_(nullptr), tarOfProc_(nullptr),
      cellLoadWeightsPtr_(nullptr) {
    readMesh(cont);
}

OpenHurricane::geometryMesh::geometryMesh(object &&ob, const controller &cont)
    : baseMesh(cont), registerTable(std::move(ob)), geometryMeshCore(), globalMeshInfoPtr_(nullptr),
      globalFaceZoneL_(), perFaceMaps_(), pairProcessSize_(0), pairPeriodicProcessSize_(0),
      periodicPairSize_(0), cellOrthoPtr_(nullptr), faceCellIntersectionPtr_(nullptr),
      fIntSectFCentrePtr_(nullptr), faceSkewnessWeightPtr_(nullptr), skewFacePtr_(nullptr),
      isSkewFacePtr_(nullptr), LDUFaceMapPtr_(nullptr), processCutPtr_(nullptr),
      processPerPtr_(nullptr), sorCellCentrePtr_(nullptr), sorKNN_(nullptr), tarOfProc_(nullptr),
      cellLoadWeightsPtr_(nullptr) {
    readMesh(cont);
}

OpenHurricane::geometryMesh::geometryMesh(object &&ob, const controller &cont,
                                          const std::string &meshStr)
    : baseMesh(cont), registerTable(std::move(ob)), geometryMeshCore(), globalMeshInfoPtr_(nullptr),
      globalFaceZoneL_(), perFaceMaps_(), pairProcessSize_(0), pairPeriodicProcessSize_(0),
      periodicPairSize_(0), cellOrthoPtr_(nullptr), faceCellIntersectionPtr_(nullptr),
      fIntSectFCentrePtr_(nullptr), faceSkewnessWeightPtr_(nullptr), skewFacePtr_(nullptr),
      isSkewFacePtr_(nullptr), LDUFaceMapPtr_(nullptr), processCutPtr_(nullptr),
      processPerPtr_(nullptr), sorCellCentrePtr_(nullptr), sorKNN_(nullptr), tarOfProc_(nullptr),
      cellLoadWeightsPtr_(nullptr) {
    real scaleFactor = 1.0;
    const auto &meshSetCont = cont.subController("meshSet");

    // Get mesh unit from controller file.
    string unitw = meshSetCont.findWordOrDefault("unit", "m");

    originMeshRead *gmRead = new FluentMeshRead(meshStr, HurMPI::getProcSize());
    (*gmRead).read(unitw);

    if (gmRead->checkCoupledWall()) {
        setCoupledWall(cont, *gmRead);
    }

    // Get scale factor for mesh.
    if (HurMPI::master()) {
        originCellIndex_.transfer((*gmRead).originCellIndex());
        scaleFactor = meshScaleFactor(unitw);
    }
    HurMPI::bcast(&scaleFactor, 1, feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

    if (HurMPI::master()) {
        allMeshNodeNumber_ = gmRead->globalNodes();
        allMeshFaceNumber_ = gmRead->faces().size();
    }
    HurMPI::bcast(&allMeshNodeNumber_, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                  HurMPI::getComm());
    HurMPI::bcast(&allMeshFaceNumber_, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                  HurMPI::getComm());

    decomposingMeshByMetis dcM(*gmRead, HurMPI::getProcSize(), meshSetCont);
    dcM.decomposing();
    dcM.meshToProcessor();

    auto reord = reordering::creator(dcM, meshSetCont);
    HurMPI::barrier();
    reord->reorder();
    perm_.transfer(reord->perm());
    iperm_.transfer(reord->iperm());
    gmRead->clear();
    transferFromDCM(cont, dcM);
    setareAllHexOrQuadri();
    HurDelete(gmRead);

    // Scale mesh with a scalar factor
    scaleMesh(scaleFactor);

    if (meshSetCont.found("scaleMesh")) {
        // Scale mesh with a vector factor.
        auto vScale = meshSetCont.findType<vector>("scaleMesh", vector(0.0));
        scaleMesh(vScale);
    }

    calcBashMesh();

    if (HurMPI::parRun()) {
        allMeshCellNumber_ = nCells();
        HurMPI::allReduce(allMeshCellNumber_, MPI_SUM);
    } else {
        allMeshCellNumber_ = nCells();
    }
}

hur_nodiscard OpenHurricane::fileName
OpenHurricane::geometryMesh::getMeshFileName(const controller &cont) const {
    fileName meshFile;
    if (cont.found("meshFile")) {
        meshFile = cont.findWord("meshFile");
    } else {
        meshFile = cont.subController("meshSet").findWord("meshFile", true);
    }
    if (cont.found("meshStr")) {
        const_cast<controller &>(cont).remove("meshStr");
    }
    if (!meshFile.isAbsolute()) {
        meshFile = Iteration().inputPath() / meshFile;
    }
    return meshFile;
}

void OpenHurricane::geometryMesh::readMesh(const controller &cont) {
    real scaleFactor = 1.0;
    if (!cont.found("meshSet")) {
        LFatal("Cannot find meshSet section in: %s", cont.name().c_str());
    }
    const auto &meshSetCont = cont.subController("meshSet");

    if ((cont.found("meshSet") && meshSetCont.found("meshFile", true)) || cont.found("meshFile")) {
        // Get mesh unit from controller file.
        string unitw = meshSetCont.findWordOrDefault("unit", "m");

        auto gmRead =
            originMeshRead::creator(getMeshFileName(cont), HurMPI::getProcSize(), meshSetCont);

        gmRead->read(unitw);

        if (gmRead->checkCoupledWall()) {
            setCoupledWall(cont, *gmRead);
        }

        // Get scale factor for mesh.
        if (HurMPI::master()) {
            originCellIndex_.transfer(gmRead->originCellIndex());
            scaleFactor = meshScaleFactor(unitw);
        }
        HurMPI::bcast(&scaleFactor, 1, feature<real>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());

        if (HurMPI::master()) {
            allMeshNodeNumber_ = gmRead->globalNodes();
            allMeshFaceNumber_ = gmRead->faces().size();
        }
        HurMPI::bcast(&allMeshNodeNumber_, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        HurMPI::bcast(&allMeshFaceNumber_, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());

        decomposingMeshByMetis dcM(*gmRead, HurMPI::getProcSize(), meshSetCont);
        dcM.decomposing();
        dcM.meshToProcessor();

        auto reord = reordering::creator(dcM, meshSetCont);
        HurMPI::barrier();
        reord->reorder();
        perm_.transfer(reord->perm());
        iperm_.transfer(reord->iperm());
        gmRead->clear();
        transferFromDCM(cont, dcM);
        setareAllHexOrQuadri();

        // Scale mesh with a scalar factor
        scaleMesh(scaleFactor);

        if (meshSetCont.found("scaleMesh")) {
            // Scale mesh with a vector factor.
            auto vScale = meshSetCont.findType<vector>("scaleMesh", vector(0.0));
            scaleMesh(vScale);
        }

        calcBashMesh();

        if (HurMPI::parRun()) {
            allMeshCellNumber_ = nCells();
            HurMPI::allReduce(allMeshCellNumber_, MPI_SUM);
        } else {
            allMeshCellNumber_ = nCells();
        }

        if (!gmRead->isHurricaneMesh()) {
            if (!Iteration().isGUIMode()) {
                const auto convertMeshName = Iteration().meshFileName();
                Iteration().writeMesh(convertMeshName);
            }
        }
    } else {
        LFatal("No mesh information");
    }
}

void OpenHurricane::geometryMesh::setCoupledWall(const controller &cont, originMeshRead &gmRead) {
    if (!cont.found("boundaryCondition")) {
        LFatal("Cannot find boundary condition in %s", cont.name().c_str());
    }

    controller &subCont = const_cast<controller &>(cont).subController("boundaryCondition");

    for (integer i = 0; i < gmRead.faceZones().size(); i++) {

        if (subCont.found(gmRead.faceZones()[i].name() + "_SHADOW")) {
            controller &wallCont = subCont.subController(gmRead.faceZones()[i].name());
            if (wallCont.findWord("bcType") == "wall") {
                controller &sdwCont =
                    subCont.subController(gmRead.faceZones()[i].name() + "_SHADOW");
                if (!wallCont.found("adjacentCellZone") || !sdwCont.found("adjacentCellZone")) {
                    errorAbortStr(("Cannot find specification of adjacent cellzone in " +
                                   gmRead.faceZones()[i].name() + " and its shadow zone."));
                }

                faceList &faces = const_cast<faceList &>(gmRead.faces());

                bool flag = true;
                for (integer fi = gmRead.faceZones()[i].firstIndex();
                     fi <= gmRead.faceZones()[i].lastIndex(); fi++) {
                    integer cl = faces[fi].leftCell();
                    integer cr = faces[fi].rightCell();
                    if (cr != 0) {
                        for (integer czi = 0; czi < gmRead.cellZones().size(); czi++) {
                            const cellZone &CZ = gmRead.cellZones()[czi];
                            if (cl >= CZ.firstIndex() && cl <= CZ.lastIndex()) {
                                if (CZ.name() != wallCont.findWord("adjacentcellzone") &&
                                    CZ.name() != sdwCont.findWord("adjacentcellzone")) {
                                    flag = false;
                                }
                            } else if (cr >= CZ.firstIndex() && cl <= CZ.lastIndex()) {
                                if (CZ.name() != wallCont.findWord("adjacentcellzone") &&
                                    CZ.name() != sdwCont.findWord("adjacentcellzone")) {
                                    flag = false;
                                }
                            }
                        }
                    }
                }
                if (!flag) {
                    std::string errMsg("Error in adjacent cellzone specification of ");
                    errMsg += gmRead.faceZones()[i].name() + " and its shadow.";
                    errorAbortStr(errMsg);
                }

                for (integer czi = 0; czi < gmRead.cellZones().size(); czi++) {
                    const cellZone &CZ = gmRead.cellZones()[czi];
                    if (CZ.name() == wallCont.findWord("adjacentcellzone")) {
                        //flag = true;
                        cellList &cells = const_cast<cellList &>(gmRead.cells());
                        faceList sdwfaces(gmRead.faceZones()[i].size());
                        integer count = 0;
                        for (integer fi = gmRead.faceZones()[i].firstIndex();
                             fi <= gmRead.faceZones()[i].lastIndex(); fi++) {
                            integer cl = faces[fi].leftCell();
                            integer cr = faces[fi].rightCell();
                            if (cr != 0) {
                                if (cl >= CZ.firstIndex() && cl <= CZ.lastIndex()) {
                                    sdwfaces[count] = faces[fi];
                                    faces[fi].rightCell() = -1;
                                    for (integer pi = 0; pi < faces[fi].size(); pi++) {
                                        sdwfaces[count][pi] = faces[fi][faces[fi].size() - pi - 1];
                                    }
                                    sdwfaces[count].leftCell() = sdwfaces[count].rightCell();
                                    sdwfaces[count].rightCell() = -1;
                                    integer local = cells[sdwfaces[count].leftCell()].findFace(fi);
                                    cells[sdwfaces[count].leftCell()].facesList()[local] =
                                        faces.size() + count;
                                } else {
                                    sdwfaces[count] = faces[fi];
                                    sdwfaces[count].rightCell() = -1;
                                    for (integer pi = 0; pi < faces[fi].size(); pi++) {
                                        faces[fi][pi] = sdwfaces[count][faces[fi].size() - pi - 1];
                                    }
                                    faces[fi].leftCell() = faces[fi].rightCell();
                                    faces[fi].rightCell() = -1;
                                    integer local = cells[faces[fi].leftCell()].findFace(fi);
                                    cells[faces[fi].leftCell()].facesList()[local] =
                                        faces.size() + count;
                                }
                                gmRead.interiorWallFaceMap()[fi] = faces.size() + count;
                                count++;
                            }
                        }
                        faceZone sdw = gmRead.faceZones()[i];
                        sdw.resetName(gmRead.faceZones()[i].name() + "_SHADOW");
                        sdw.setIndex(-gmRead.faceZones()[i].index());
                        sdw.setFirstIndex(gmRead.faces().size());
                        sdw.setLastIndex(gmRead.faces().size() + count - 1);
                        const_cast<faceZoneList &>(gmRead.faceZones()).append(sdw);
                        sdwfaces.resize(count);
                        faces.append(sdwfaces);
                    }
                }
            }
        }
    }
}

OpenHurricane::geometryMesh::~geometryMesh() noexcept {
    clearGeoMesh();
}

void OpenHurricane::geometryMesh::scaleMesh(const real scaleFactor) {
    Pout << "    Info: scale the mesh with the factor: " << scaleFactor << std::endl;
    points_ *= max(tiny, scaleFactor);
    if (perZones_.size() != 0) {
        for (integer i = 0; i < perZones_.size(); ++i) {
            if (perZones_[i].isTranslational()) {
                perZones_[i].setTranDistance(perZones_[i].tranDistance() * max(tiny, scaleFactor));
            }
        }
    }
    origin_ *= max(tiny, scaleFactor);
    axis_ *= max(tiny, scaleFactor);
}

void OpenHurricane::geometryMesh::scaleMesh(const vector &scaleFactor) {
    Pout << "    Info: scale the mesh with the factor: " << toString(scaleFactor) << std::endl;
    for (integer i = 0; i < points_.size(); ++i) {
        points_[i] = componentMultiply(points_[i], scaleFactor);
    }
    if (perZones_.size() != 0) {
        for (integer i = 0; i < perZones_.size(); ++i) {
            if (perZones_[i].isTranslational()) {
                perZones_[i].setTranDistance(
                    componentMultiply(perZones_[i].tranDistance(), scaleFactor));
            }
        }
    }
    origin_ = componentMultiply(origin_, scaleFactor);
    axis_ = componentMultiply(axis_, scaleFactor);
}

void OpenHurricane::geometryMesh::calcNumberOfCutFaces() const {
    if (!totalCutFacesPtr_) {
        integer n = 0;
        for (integer faceI = 0; faceI < faceZones_.size(); ++faceI) {
            if (faceZones_[faceI].isCutFace()) {
                n += faceZones_[faceI].size();
            }
        }
        totalCutFacesPtr_.reset(new integer(n));
    }
}

/*!\brief Ghost cell index. Note ghost cell ID for the per and cut cell have been asigned
               in the pretreat code following the same rule, i.e.
               GhostID=ID+totalInternalCellNumber.*/
void OpenHurricane::geometryMesh::createGhostCell() {
    if (!isGhostCellCreated()) {
#ifdef HUR_DEBUG
        Pout << "geometryMesh::createGhostCell() : creating ghost cells" << std::endl;
#endif // HUR_DEBUG

        cellList &cs = cells_;
        faceList &fs = faces_;
        const faceZoneList &fz = faceZones_;

        // Creating first ghost cell layer.
        integer n = nCells();
        //cellList ghostCell(nBoundFaces());
        cellList ghostCell(nBoundFaces() * integer(ghostCellsType()));
        integer ghostCellI = 0;
        for (integer fzI = 0; fzI < fz.size(); fzI++) {
            if (!fz[fzI].isInterior()) {
                for (integer i = 0; i < fz[fzI].size(); i++) {
                    integer fi = fz[fzI].firstIndex() + i;

                    fs[fi].rightCell() = n++;
                    //cs.resize(n + 1);

                    ghostCell[ghostCellI].facesList().resize(1, fi);
                    //cs[n].facesList().resize(1, fi);

                    integer nnodes = fs[fi].size();

                    ghostCell[ghostCellI].nodesList().resize(nnodes);
                    //cs[n].nodesList().resize(nnodes);

                    ghostCell[ghostCellI].nodesList() = fs[fi];
                    //cs[n].nodesList() = fs[fi];
                    ghostCellI++;
                }
            }
        }

        cs.append(ghostCell);

        if (cs.size() != nTotalCells()) {
            std::string errMsg = " The number of ghost cells isn't equal to "
                                 "total boundary faces";
            errMsg += "\n";
            errMsg += " Ghost cell size: " + std::to_string(cs.size() - nCells());
            errMsg += "\n";
            errMsg += " Total boundary face size: " + std::to_string(nBoundFaces());
            errorAbortStr(errMsg);
        }
        setIsGhostCellCreated(true);

        if (HurMPI::parRun()) {
            const std::map<integer, integerArray> &facePairMap =
                globalMeshInfoPtr_->globalFaceIndeces().facePairMap();
            createFirstLayerCutZone(facePairMap);
            createFirstLayerPerZone(perFaceMaps_);
            for (integer layer = 2; layer <= ghostCellsType(); layer++) {
                createSecondGhostCellLayer(layer, facePairMap, perFaceMaps_, cutZones_, perZones_);
            }
        } else {
            std::map<integer, integerArray> null;
            createFirstLayerPerZone(perFaceMaps_);
            for (integer layer = 2; layer <= ghostCellsType(); layer++) {
                createSecondGhostCellLayer(layer, null, perFaceMaps_, cutZones_, perZones_);
            }
        }
    }
}

hur_nodiscard OpenHurricane::integer OpenHurricane::geometryMesh::nFaceZoneExceptCutZone() const {
    integer nfz = 0;

    for (integer i = 0; i < faceZones_.size(); ++i) {
        if (!faceZones_[i].isCutFace()) {
            nfz++;
        }
    }
    return nfz;
}

void OpenHurricane::geometryMesh::createFirstLayerCutZone(
    const std::map<integer, integerArray> &_faceMap) {
#ifdef HUR_DEBUG
    Pout << "geometryMesh::createFirstLayerCutZone() : Creating ghost cells "
            "for decomposing face zones..."
         << std::endl;
#endif // HUR_DEBUG
    // It is no need for serial running to run this subroutine.
    if (!HurMPI::parRun()) {
        return;
    }

    // Face list
    faceList &fs = faces_;

    integerArray cutIndex;

    for (integer czi = 0; czi < cutZones_.size(); ++czi) {
        if (cutZones_[czi].isThisProcReceiv()) {
            cutIndex.append(czi);
        }
        cutZones_[czi].clearList();
    }

    integerArrayArray sors(cutIndex.size());
    integerArrayArray dess(cutIndex.size());
    integerArray sorsSize(cutIndex.size(), Zero);

    for (integer i = 0; i < cutIndex.size(); ++i) {
        sors[i].resize(integer(_faceMap.size()));
        dess[i].resize(integer(_faceMap.size()));
    }

    // The information storaged in _faceMap is as follows
    // iter->first: the local index of the cut face of this process.
    // iter->second: A list se[4]:
    //               se[0]: The remote process id.
    //               se[1]: The output index of the cut face.
    //               se[2]: The local index of the cut face in the remote process.
    //               se[3]: The left cell index of the cut face in the remote process.
    for (std::map<integer, integerArray>::const_iterator iter = _faceMap.begin();
         iter != _faceMap.end(); ++iter) {
        for (integer czi = 0; czi < cutIndex.size(); ++czi) {
            if (iter->second[0] == cutZones_[cutIndex[czi]].sendProc()) {
                dess[czi][sorsSize[czi]] = fs[iter->first].rightCell();
                sors[czi][sorsSize[czi]] = iter->second[3];
                sorsSize[czi]++;
                break;
            }
        }
    }

    // Setting cutZoneList
    for (integer czi = 0; czi < cutIndex.size(); ++czi) {
        cutZones_[cutIndex[czi]].sor().resize(sorsSize[czi]);
        cutZones_[cutIndex[czi]].des().resize(sorsSize[czi]);
        for (integer i = 0; i < sorsSize[czi]; ++i) {
            cutZones_[cutIndex[czi]].sor()[i] = sors[czi][i];
            cutZones_[cutIndex[czi]].des()[i] = dess[czi][i];
        }
        sors[czi].clear();
        dess[czi].clear();
    }

    // Broadcast cutzone information.
    for (integer czi = 0; czi < cutZones_.size(); ++czi) {
        // Firstly, broadcast size of sorList (or desList, they are the same as each other).
        integer desSize = cutZones_[czi].sor().size();
        HurMPI::bcast(&desSize, 1, feature<integer>::MPIType, cutZones_[czi].receivProc(),
                      HurMPI::getComm());

        // Secondly, allocate storage for sorList and desList if the process id is not the same with cutZones_[czi].receivProc()
        if (!cutZones_[czi].isThisProcReceiv()) {
            cutZones_[czi].sor().resize(desSize);
            cutZones_[czi].des().resize(desSize);
        }
        // Barrier for all processes
        HurMPI::barrier(HurMPI::getComm());

        // Broadcasting
        HurMPI::bcastList(cutZones_[czi].sor(), cutZones_[czi].receivProc(), HurMPI::getComm());
        HurMPI::bcastList(cutZones_[czi].des(), cutZones_[czi].receivProc(), HurMPI::getComm());
    }
#ifdef HUR_DEBUG
    Pout << "geometryMesh::createFirstLayerCutZone() : Finish creating ghost "
            "cells for decomposing face zones..."
         << std::endl;
#endif // HUR_DEBUG
}

void OpenHurricane::geometryMesh::createFirstLayerPerZone(
    const std::map<integer, integerArray> &_faceMap) {
    if (perZones_.size() == 0) {
        return;
    }
#ifdef HUR_DEBUG
    Pout << "geometryMesh::createFirstLayerPerZone() : Creating ghost cells "
            "for periodic face zones..."
         << std::endl;
#endif // HUR_DEBUG
    // Face list
    faceList &fs = faces_;

    for (integer pzi = 0; pzi < perZones_.size(); ++pzi) {
        perZones_[pzi].desFace() = perZones_[pzi].des();
        if (perZones_[pzi].isThisProcReceiv()) {
            for (integer i = 0; i < perZones_[pzi].des().size(); ++i) {
                const integer fi = perZones_[pzi].des()[i];
                perZones_[pzi].des()[i] = fs[fi].rightCell();
            }
        }
        if (perZones_[pzi].isSendFromThisProc()) {
            perZones_[pzi].sorFace() = perZones_[pzi].sor();
            for (integer i = 0; i < perZones_[pzi].des().size(); ++i) {
                const integer fi = perZones_[pzi].sor()[i];
                perZones_[pzi].sor()[i] = fs[fi].leftCell();
            }
        }
    }
#ifdef HUR_DEBUG
    Pout << "geometryMesh::createFirstLayerPerZone() : Finish creating ghost "
            "cells for periodic face zones..."
         << std::endl;
#endif // HUR_DEBUG
}

/* Create second or more layers ghost cells.*/
void OpenHurricane::geometryMesh::createSecondGhostCellLayer(
    const short nLayer, const std::map<integer, integerArray> &_faceMap,
    const std::map<integer, integerArray> &_facePeriMap, cutZoneList &cZ, perZoneList &pZ) {
    if (!SNCFPtr_) {
        LFatal("Second neighbour cell of face pointer is null!");
    }

    if (HurMPI::parRun()) {
        cutZoneList newCZ(pairProcessSize_);

        for (integer czi = 0; czi < pairProcessSize_; ++czi) {
            newCZ[czi].setSendProc(cZ[czi].sendProc());
            newCZ[czi].setRecvProc(cZ[czi].receivProc());
            newCZ[czi].sor().resize(cZ[czi].cutSize());
            newCZ[czi].des().resize(cZ[czi].cutSize());
        }

        integerArray cutIndex;
        integerArray cutIndexRecv;

        for (integer czi = 0; czi < pairProcessSize_; ++czi) {
            if (cutZones_[czi].isSendFromThisProc()) {
                cutIndex.append(czi);
            } else if (cutZones_[czi].isThisProcReceiv()) {
                cutIndexRecv.append(czi);
            }
        }

        integerArrayArray sencondSorList(cutIndex.size());
        integerArrayArray sencondSorFaceList(cutIndex.size());
        integerArray sorsSize(cutIndex.size(), Zero);

        for (integer ctii = 0; ctii < cutIndex.size(); ctii++) {
            sencondSorList[ctii].resize(cutZones_[cutIndex[ctii]].sor().size());
            sencondSorFaceList[ctii].resize(cutZones_[cutIndex[ctii]].sor().size());
        }

        integerArrayArray sencondDesList(cutIndexRecv.size());
        integerArrayArray sencondDesFaceList(cutIndexRecv.size());
        for (integer ctii = 0; ctii < cutIndexRecv.size(); ctii++) {
            sencondDesList[ctii].resize(cutZones_[cutIndexRecv[ctii]].des().size());
            sencondDesFaceList[ctii].resize(cutZones_[cutIndexRecv[ctii]].des().size());
        }

        // The information storaged in _faceMap is as follows
        // iter->first: the local index of the cut face of this process.
        // iter->second: A list se[4]:
        //               se[0]: The remote process id.
        //               se[1]: The output index of the cut face.
        //               se[2]: The local index of the cut face in the remote process.
        //               se[3]: The left cell index of the cut face in the remote process.
        for (std::map<integer, integerArray>::const_iterator iter = _faceMap.begin();
             iter != _faceMap.end(); ++iter) {
            for (integer ctii = 0; ctii < cutIndex.size(); ++ctii) {
                if (iter->second[0] == cutZones_[cutIndex[ctii]].receivProc()) {
                    sencondSorList[ctii][sorsSize[ctii]] =
                        (*SNCFPtr_)[iter->first][2 * (nLayer - 1) - 2];
                    sencondSorFaceList[ctii][sorsSize[ctii]] = iter->second[2];
                    sorsSize[ctii]++;
                    break;
                }
            }
        }

        HurMPI::Status status;
        for (integer czi = 0; czi < pairProcessSize_; ++czi) {
            if (cutZones_[czi].isSendFromThisProc()) {
                for (integer ctii = 0; ctii < cutIndex.size(); ++ctii) {
                    if (cutIndex[ctii] == czi) {
                        HurMPI::sendList(sencondSorList[ctii], cutZones_[czi].receivProc(),
                                         cutZones_[czi].sendProc() + 1, HurMPI::getComm());
                        HurMPI::sendList(sencondSorFaceList[ctii], cutZones_[czi].receivProc(),
                                         cutZones_[czi].sendProc() + 2, HurMPI::getComm());
                    }
                }
            } else if (cutZones_[czi].isThisProcReceiv()) {
                for (integer ctii = 0; ctii < cutIndexRecv.size(); ++ctii) {
                    if (cutIndexRecv[ctii] == czi) {
                        HurMPI::recvList(sencondDesList[ctii], cutZones_[czi].sendProc(),
                                         cutZones_[czi].sendProc() + 1, HurMPI::getComm(), &status);
                        HurMPI::recvList(sencondDesFaceList[ctii], cutZones_[czi].sendProc(),
                                         cutZones_[czi].sendProc() + 2, HurMPI::getComm(), &status);
                    }
                }
            }
            HurMPI::barrier(HurMPI::getComm());
        }

        for (integer czi = 0; czi < pairProcessSize_; ++czi) {
            for (integer ctii = 0; ctii < cutIndexRecv.size(); ++ctii) {
                if (cutIndexRecv[ctii] == czi) {
                    for (integer fi = 0; fi < sencondDesList[ctii].size(); ++fi) {
                        newCZ[czi].sor()[fi] = sencondDesList[ctii][fi];
                        integer cutFaceI = sencondDesFaceList[ctii][fi];
                        newCZ[czi].des()[fi] = (*SNCFPtr_)[cutFaceI][2 * (nLayer - 1) - 1];
                    }
                }
            }
        }

        // Broadcast cutzone information.
        for (integer czi = 0; czi < newCZ.size(); ++czi) {
            // Broadcasting
            HurMPI::bcastList(newCZ[czi].sor(), newCZ[czi].receivProc(), HurMPI::getComm());
            HurMPI::bcastList(newCZ[czi].des(), newCZ[czi].receivProc(), HurMPI::getComm());
            HurMPI::barrier(HurMPI::getComm());
        }
        cZ.append(newCZ);
    }
    if (pairPeriodicProcessSize_ != 0) {
        perZoneList newPZ(pairPeriodicProcessSize_);

        for (integer pzi = 0; pzi < pairPeriodicProcessSize_; ++pzi) {
            newPZ[pzi] = perZones_[pzi];
        }
        for (integer pzi = 0; pzi < pairPeriodicProcessSize_; ++pzi) {
            if (newPZ[pzi].isThisProcReceiv()) {
                for (integer i = 0; i < perZones_[pzi].des().size(); ++i) {
                    const integer fi = perZones_[pzi].desFace()[i];
                    newPZ[pzi].des()[i] = (*SNCFPtr_)[fi][2 * (nLayer - 1) - 1];
                }
            }
            if (newPZ[pzi].isSendFromThisProc()) {
                for (integer i = 0; i < perZones_[pzi].sor().size(); ++i) {
                    const integer fi = perZones_[pzi].sorFace()[i];
                    newPZ[pzi].sor()[i] = (*SNCFPtr_)[fi][2 * (nLayer - 1) - 2];
                }
            }
        }
        pZ.append(newPZ);
    }
}

void OpenHurricane::geometryMesh::clearGeoMesh() noexcept {
    totalCutFacesPtr_.clear();
    totalPerFacesPtr_.clear();
    globalMeshInfoPtr_.clear();
    cellOrthoPtr_.clear();
    faceSkewnessWeightPtr_.clear();
    faceCellIntersectionPtr_.clear();
    fIntSectFCentrePtr_.clear();
    skewFacePtr_.clear();
    isSkewFacePtr_.clear();
    LDUFaceMapPtr_.clear();
    indexMap_.clear();

    processCutPtr_.clear();
    processPerPtr_.clear();
    sorCellCentrePtr_.clear();
    sorKNN_.clear();
    tarOfProc_.clear();

    cellLoadWeightsPtr_.clear();
    CRSAdrrPtr_.clear();
}

void OpenHurricane::geometryMesh::checkMeshFaceSkewness(real &maxSkew) const {
    realArray skews =
        meshCheckers::faceSkewness(*this, points(), faceCentre(), faceArea(), cellCentre());

    maxSkew = max(skews);

    HurMPI::allReduce(maxSkew, MPI_MAX);

    integer nSkewsLarge = 0;

    if (report) {
        for (integer i = 0; i < skews.size(); ++i) {
            if (skews[i] > baseMeshWarningThreshold::skewenessThreshold) {
                if (report) {
                    LWarning("Severe face skewness: %e in face: %d face centre: %s", skews[i], i,
                             toString(faceCentre()[i]).c_str());
                }
                nSkewsLarge++;
            }
        }

        HurMPI::allReduce(nSkewsLarge, MPI_SUM);
        if (nSkewsLarge > 0) {
            Pout << "     Warning: " << nSkewsLarge << " faces exceed high face skewness limit: "
                 << baseMeshWarningThreshold::skewenessThreshold << std::endl;
        }
    }
    Pout << "      maximum face skewness: " << maxSkew << std::endl;
}

void OpenHurricane::geometryMesh::checkMeshOrthogonality(real &minOrtho) const {
    realArray ortho =
        meshCheckers::faceOrthogonality(*this, faceArea(), faceCentre(), cellCentre());

    makeCellOrtho(ortho);
    minOrtho = min(ortho);

    HurMPI::allReduce(minOrtho, MPI_MIN);

    integer nOrthoSmall = 0;

    if (report) {
        for (integer i = 0; i < ortho.size(); ++i) {
            if (ortho[i] < baseMeshWarningThreshold::OrthogonalityThreshold) {
                if (report) {
                    LWarning("Severe face non-orthogonality: %e in face: %d face centre: %s",
                             ortho[i], i, toString(faceCentre()[i]).c_str());
                }
                nOrthoSmall++;
            }
        }

        HurMPI::allReduce(nOrthoSmall, MPI_SUM);
        if (nOrthoSmall > 0) {
            Pout << "     Warning: " << nOrthoSmall
                 << " faces exceed low face orthogonality limit: "
                 << baseMeshWarningThreshold::OrthogonalityThreshold << std::endl;
        }
    }
    Pout << "      minimum face orthogonality: " << minOrtho << std::endl;
}

void OpenHurricane::geometryMesh::checkMeshVolRatio(real &maxRatio) const {
    realArray ratio = meshCheckers::volRatio(*this, cellVolume());

    maxRatio = max(ratio);

    HurMPI::allReduce(maxRatio, MPI_MAX);

    integer nVolRatioLarge = 0;

    if (report) {
        for (integer i = 0; i < ratio.size(); ++i) {
            if (ratio[i] > baseMeshWarningThreshold::volumeRatioThreshold) {
                if (report) {
                    LWarning("Highly cell volume ratio: %e along face: %d face centre: %s",
                             ratio[i], i, toString(faceCentre()[i]).c_str());
                }
                nVolRatioLarge++;
            }
        }

        HurMPI::allReduce(nVolRatioLarge, MPI_SUM);
        if (nVolRatioLarge > 0) {
            Pout << "     Warning: " << nVolRatioLarge << " cells exceed high volume ratio limit: "
                 << baseMeshWarningThreshold::volumeRatioThreshold << std::endl;
        }
    }
    Pout << "      maximum volume ratio: " << maxRatio << std::endl;
}

void OpenHurricane::geometryMesh::checkMeshAspectRatio(real &maxRatio) const {
    const realArray &ratio = aspectRatio();

    maxRatio = max(ratio);

    HurMPI::allReduce(maxRatio, MPI_MAX);

    integer nAspectRatioLarge = 0;

    if (report) {
        for (integer i = 0; i < ratio.size(); ++i) {
            if (ratio[i] > baseMeshWarningThreshold::aspectRatioThreshold) {
                if (report) {
                    LWarning("Highly cell aspect ratio: %e along face: %d face centre: %s",
                             ratio[i], i, toString(faceCentre()[i]).c_str());
                }
                nAspectRatioLarge++;
            }
        }

        HurMPI::allReduce(nAspectRatioLarge, MPI_SUM);
        if (nAspectRatioLarge > 0) {
            Pout << "     Warning: " << nAspectRatioLarge
                 << " cells exceed high aspect ratio limit: "
                 << baseMeshWarningThreshold::aspectRatioThreshold << std::endl;
        }
    }
    Pout << "      maximum cell aspect ratio: " << maxRatio << std::endl;
}

void OpenHurricane::geometryMesh::checkAndRportMeshQuality() const {
    Pout << std::endl << "    Info: Statistics of mesh..." << std::endl;

    real minArea = large;
    real maxArea = -large;
    meshCheckers::minMaxFaceArea(faceArea(), faceZones(), minArea, maxArea);

    real minVol = large;
    real maxVol = -large;
    real totalVol;
    meshCheckers::minMaxCellVol(cellVolume(), cellZones(), minVol, maxVol, totalVol);

    HurMPI::reduce(minArea, MPI_MIN);
    HurMPI::reduce(maxArea, MPI_MAX);
    HurMPI::reduce(minVol, MPI_MIN);
    HurMPI::reduce(maxVol, MPI_MAX);
    HurMPI::reduce(totalVol, MPI_SUM);

    Pout.setScientic() << "      Face area magnitude: " << std::endl
                       << "        minimum face area: " << minArea << " m^2" << std::endl
                       << "        maximum face area: " << maxArea << " m^2" << std::endl
                       << "      Cell volume: " << std::endl
                       << "        minimum cell volume: " << minVol << " m^3" << std::endl
                       << "        maximum cell volume: " << maxVol << " m^3" << std::endl
                       << "        total cell volume: " << totalVol << " m^3" << std::endl;

    Pout.unsetScientic();

    vector minP;
    vector maxP;
    meshCheckers::minMaxPointDomain(points_, nPoints(), minP, maxP);

    HurMPI::reduceVS(minP, MPI_MIN);
    HurMPI::reduceVS(maxP, MPI_MAX);

    Pout.setScientic() << std::endl
                       << "      Mesh domain: " << std::endl
                       << "        minimum x: " << minP[0] << "m   maximum x: " << maxP[0] << "m"
                       << std::endl
                       << "        minimum y: " << minP[1] << "m   maximum y: " << maxP[1] << "m"
                       << std::endl
                       << "        minimum z: " << minP[2] << "m   maximum z: " << maxP[2] << "m"
                       << std::endl
                       << std::endl;

    vector minCC;
    vector maxCC;
    meshCheckers::minMaxCellCentre(cellCentre(), cellZones(), minCC, maxCC);

    HurMPI::reduceVS(minCC, MPI_MIN);
    HurMPI::reduceVS(maxCC, MPI_MAX);
    Pout.setScientic() << std::endl
                       << "      Mesh cell center domain: " << std::endl
                       << "        minimum x: " << minCC[0] << "m   maximum x: " << maxCC[0] << "m"
                       << std::endl
                       << "        minimum y: " << minCC[1] << "m   maximum y: " << maxCC[1] << "m"
                       << std::endl
                       << "        minimum z: " << minCC[2] << "m   maximum z: " << maxCC[2] << "m"
                       << std::endl
                       << std::endl;
    Pout.unsetScientic();

    real maxSkew;
    checkMeshFaceSkewness(maxSkew);

    real minOrtho;
    checkMeshOrthogonality(minOrtho);

    real maxVolRatio;
    checkMeshVolRatio(maxVolRatio);

    real maxAspectRatio;
    checkMeshAspectRatio(maxAspectRatio);
}

void OpenHurricane::geometryMesh::calcCellOrtho(const realArray &faceOrtho,
                                                realArray &cellOrtho) const {
    for (integer i = 0; i < nCells(); ++i) {
        cellOrtho[i] = real(1);
        for (integer j = 0; j < cells()[i].faceSize(); ++j) {
            const auto fi = cells()[i].facei(j);
            cellOrtho[i] = min(cellOrtho[i], faceOrtho[fi]);
        }
    }
}

void OpenHurricane::geometryMesh::makeCellOrtho(const realArray &faceOrtho) const {
    if (!cellOrthoPtr_) {
        cellOrthoPtr_.reset(new realArray(nCells(), real(1)));
        calcCellOrtho(faceOrtho, *cellOrthoPtr_);
    }
}

void OpenHurricane::geometryMesh::makeCellOrtho() const {
    if (!cellOrthoPtr_) {
        realArray ortho =
            meshCheckers::faceOrthogonality(*this, faceArea(), faceCentre(), cellCentre());
        makeCellOrtho(ortho);
    }
}

hur_nodiscard const OpenHurricane::realArray &OpenHurricane::geometryMesh::cellOrtho() const {
    if (!cellOrthoPtr_) {
        makeCellOrtho();
    }
    return *cellOrthoPtr_;
}

void OpenHurricane::geometryMesh::makeFaceCellIntersectionPoint() const {
    faceCellIntersectionPtr_.clear();
    faceCellIntersectionPtr_.reset(new vectorArray(nFaces()));
    for (integer iface = 0; iface < nFaces(); ++iface) {
        const auto cl = faces()[iface].leftCell();
        const auto cr = faces()[iface].rightCell();
        const auto ee = cellCentre()[cr] - cellCentre()[cl];
        const auto n = faceArea()[iface].normalized();
        const auto e = ee.normalized();

        (*faceCellIntersectionPtr_)[iface] =
            ((faceCentre()[iface] - cellCentre()[cl]) * n) / (e * n) * e + cellCentre()[cl];
    }
}

hur_nodiscard const OpenHurricane::vectorArray &
OpenHurricane::geometryMesh::faceCellIntersectionPoint() const {
    if (!faceCellIntersectionPtr_) {
        makeFaceCellIntersectionPoint();
    }
    return *faceCellIntersectionPtr_;
}

void OpenHurricane::geometryMesh::makeFIntSectFCentre() const {
    fIntSectFCentrePtr_.clear();
    fIntSectFCentrePtr_.reset(new vectorArray(nFaces()));
    const auto &fcis = faceCellIntersectionPoint();
    for (integer iface = 0; iface < nFaces(); ++iface) {
        (*fIntSectFCentrePtr_)[iface] = faceCentre()[iface] - fcis[iface];
    }
}

hur_nodiscard const OpenHurricane::vectorArray &
OpenHurricane::geometryMesh::fIntSectFCentre() const {
    if (!fIntSectFCentrePtr_) {
        makeFIntSectFCentre();
    }
    return *fIntSectFCentrePtr_;
}

void OpenHurricane::geometryMesh::makeFaceSkewnessWeight(const real minSkewness) const {
    faceSkewnessWeightPtr_.clear();

    faceSkewnessWeightPtr_.reset(new realArray(nFaces()));

    const auto &rfd = faceCellIntersectionPoint();
    const auto &skfa = skewFace();

    for (integer iface = 0; iface < nFaces(); ++iface) {
        (*faceSkewnessWeightPtr_)[iface] = faceWeight()[iface];
    }

    for (integer i = 0; i < skfa.size(); ++i) {
        const integer iface = skfa[i];
        const auto cl = faces()[iface].leftCell();
        const auto cr = faces()[iface].rightCell();
        const auto drl = dist(cellCentre()[cr], cellCentre()[cl]);
        const auto dffd = dist(faceCentre()[iface], rfd[iface]);
        if (dffd / drl > max(minSkewness, tiny)) {
            const auto drfd = dist(cellCentre()[cr], rfd[iface]);
            (*faceSkewnessWeightPtr_)[iface] = drfd / max(drl, veryTiny);
            if ((*faceSkewnessWeightPtr_)[iface] >= 1.0) {
                errorAbortStr(("face weight lager than 1, which is " +
                               toString((*faceSkewnessWeightPtr_)[iface])));
            }
        }
    }
}

hur_nodiscard const OpenHurricane::realArray &
OpenHurricane::geometryMesh::faceSkewnessWeight(const real minSkewness) const {
    if (!faceSkewnessWeightPtr_) {
        makeFaceSkewnessWeight(minSkewness);
    }
    return *faceSkewnessWeightPtr_;
}

void OpenHurricane::geometryMesh::makeSkewFace(const real minSkewness) const {
    skewFacePtr_.clear();

    const auto &rfd = faceCellIntersectionPoint();
    integer count = Zero;
    for (integer fzi = 0; fzi < faceZones().size(); ++fzi) {
        if (faceZones()[fzi].isInterior() || faceZones()[fzi].isCutFace() ||
            faceZones()[fzi].isPeriodic() || faceZones()[fzi].isPeriodicShadow()) {
            for (integer iface = faceZones()[fzi].firstIndex();
                 iface <= faceZones()[fzi].lastIndex(); ++iface) {
                const auto cl = faces()[iface].leftCell();
                const auto cr = faces()[iface].rightCell();
                const auto drl = dist(cellCentre()[cr], cellCentre()[cl]);
                const auto dffd = dist(faceCentre()[iface], rfd[iface]);
                if (dffd / drl > max(minSkewness, tiny)) {
                    count++;
                }
            }
        }
    }

    skewFacePtr_.reset(new integerArray(count));
    count = Zero;
    for (integer fzi = 0; fzi < faceZones().size(); ++fzi) {
        if (faceZones()[fzi].isInterior() || faceZones()[fzi].isCutFace() ||
            faceZones()[fzi].isPeriodic() || faceZones()[fzi].isPeriodicShadow()) {
            for (integer iface = faceZones()[fzi].firstIndex();
                 iface <= faceZones()[fzi].lastIndex(); ++iface) {
                const auto cl = faces()[iface].leftCell();
                const auto cr = faces()[iface].rightCell();
                const auto drl = dist(cellCentre()[cr], cellCentre()[cl]);
                const auto dffd = dist(faceCentre()[iface], rfd[iface]);
                if (dffd / drl > max(minSkewness, tiny)) {
                    (*skewFacePtr_)[count] = iface;
                    count++;
                }
            }
        }
    }
}

hur_nodiscard const OpenHurricane::integerArray &
OpenHurricane::geometryMesh::skewFace(const real minSkewness) const {
    if (!skewFacePtr_) {
        makeSkewFace(minSkewness);
    }
    return *skewFacePtr_;
}

void OpenHurricane::geometryMesh::makeIsSkewFace(const real minSkewness) const {
    isSkewFacePtr_.clear();
    isSkewFacePtr_.reset(new boolList(nFaces(), false));

    const auto &skf = skewFace(minSkewness);

    for (integer iskf = 0; iskf < skf.size(); ++iskf) {
        (*isSkewFacePtr_)[skf[iskf]] = true;
    }
}

hur_nodiscard const OpenHurricane::boolList &
OpenHurricane::geometryMesh::isSkewFace(const real minSkewness) const {
    if (!isSkewFacePtr_) {
        makeIsSkewFace(minSkewness);
    }
    return *isSkewFacePtr_;
}

hur_nodiscard OpenHurricane::integer
OpenHurricane::geometryMesh::findPointCell(const point &p) const {
    for (integer celli = 0; celli < nCells(); ++celli) {
        if (pointInCell(p, celli)) {
            return celli;
        }
    }
    return integer(-1);
}

hur_nodiscard bool OpenHurricane::geometryMesh::pointInCell(const point &p,
                                                            const integer celli) const {
#ifdef HUR_DEBUG
    if (celli < 0 || celli > nCells()) {
        errorAbortStr(("Invalid cell index: " + toString(celli)));
    }
#endif // HUR_DEBUG

    const auto &cel = cells()[celli];
    const auto &f = faces();
    const auto &fa = faceArea();
    const auto &fc = faceCentre();
    for (integer i = 0; i < cel.faceSize(); ++i) {
        const auto fi = cel.facei(i);
        const auto cl = f[fi].leftCell();
        const auto cr = f[fi].rightCell();
        const auto pf = p - fc[fi];
        auto normal = fa[fi];
        if (cl == celli) {
            normal = -normal;
        }
        if ((normal * pf) > 0) {
            return false;
        }
    }
    return true;
}

void OpenHurricane::geometryMesh::calcLDUFaceMap() const {
    if (!LDUFaceMapPtr_) {
        LDUFaceMapPtr_.reset(new integerArray(nInteriorFaces(), -1));
        integer nCountFace = 0;
        for (integer celli = 0; celli < nCells(); ++celli) {
            integerVector2DList cellFace(cells()[celli].faceSize());
            for (integer fi = 0; fi < cells()[celli].faceSize(); ++fi) {
                const auto facei = cells()[celli].facei(fi);
                const auto cl = faces()[facei].leftCell();
                const auto cr = faces()[facei].rightCell();
                cellFace[fi][0] = celli == cl ? cr : cl;
                cellFace[fi][1] = facei;
            }
            std::sort(
                cellFace.begin(), cellFace.end(),
                [](const integerVector2D &v1, const integerVector2D &v2) { return v1[0] < v2[0]; });

            for (integer fi = 0; fi < cells()[celli].faceSize(); ++fi) {
                if (celli < cellFace[fi][0]) {
                    if (cellFace[fi][1] < nInteriorFaces()) {
                        (*LDUFaceMapPtr_)[cellFace[fi][1]] = nCountFace++;
                    }
                }
            }
        }
        if (min(*LDUFaceMapPtr_) < 0) {
            LFatal("Not all interior faces have been set in LDUFaceMap");
        }
    }
}

hur_nodiscard const OpenHurricane::integerArray &OpenHurricane::geometryMesh::LDUFaceMap() const {
    if (!LDUFaceMapPtr_) {
        calcLDUFaceMap();
    }
    return *LDUFaceMapPtr_;
}

void OpenHurricane::geometryMesh::calcCRSAdrr() const {
    if (!CRSAdrrPtr_) {
        List<integer> leftCells(nInteriorFaces(), -1);
        List<integer> rightCells(nInteriorFaces(), -1);
        for (integer facei = 0; facei < nInteriorFaces(); ++facei) {
            leftCells[facei] = faces()[facei].leftCell();
            rightCells[facei] = faces()[facei].rightCell();
        }

        integer nIntf = 0;

        for (integer fzi = 0; fzi < faceZones_.size(); ++fzi) {
            if (faceZones_[fzi].isCutFace() || faceZones_[fzi].isPeriodic() ||
                faceZones_[fzi].isPeriodicShadow()) {
                nIntf += faceZones_[fzi].size();
            }
        }

        List<integer> intfCell(nIntf);
        nIntf = 0;
        List<integer> intfDes(nTotalCells(), -1);

        for (integer fzi = 0; fzi < faceZones_.size(); ++fzi) {
            if (faceZones_[fzi].isCutFace() || faceZones_[fzi].isPeriodic() ||
                faceZones_[fzi].isPeriodicShadow()) {
                for (integer fi = faceZones_[fzi].firstIndex(); fi <= faceZones_[fzi].lastIndex();
                     ++fi) {
                    intfCell[nIntf] = faces_[fi].leftCell();
                    intfDes[faces_[fi].leftCell()] = nIntf;
                    intfDes[faces_[fi].rightCell()] = nIntf++;
                }
            }
        }

        integer cntPerFaceZ = 0;
        for (integer fzi = 0; fzi < faceZones_.size(); ++fzi) {
            if (faceZones_[fzi].isPeriodic() || faceZones_[fzi].isPeriodicShadow()) {
                cntPerFaceZ++;
            }
        }

        integer cntCut = 0;
        for (auto &e : cutZones_) {
            if (e.isSendFromThisProc()) {
                cntCut++;
            } else if (e.isThisProcReceiv()) {
                cntCut++;
            }
        }
        List<cutZone> intfCutZone(cntCut);
        cntCut = 0;
        for (auto &e : cutZones_) {
            if (e.isSendFromThisProc()) {
                auto &icutZ = intfCutZone[cntCut++];
                icutZ = e;
            } else if (e.isThisProcReceiv()) {
                const auto &des = e.des();
                auto &icutZ = intfCutZone[cntCut++];
                icutZ = e;
                for (integer nn = 0; nn < des.size(); ++nn) {
                    icutZ.des()[nn] = intfDes[des[nn]];
                }
            }
        }

        integer cntPer = 0;
        if (cntPerFaceZ != 0) {
            for (auto &e : perZones_) {
                if (e.isSameProc()) {
                    if (e.isSendFromThisProc()) {
                        cntPer++;
                    }
                } else if (e.isSendFromThisProc()) {
                    cntPer++;
                } else if (e.isThisProcReceiv()) {
                    cntPer++;
                }
            }
        }
        List<perZone> intfPerZone(cntPer);
        cntPer = 0;
        if (cntPerFaceZ != 0) {
            for (auto &e : perZones_) {
                if (e.isSameProc()) {
                    if (e.isSendFromThisProc()) {
                        auto &iperZ = intfPerZone[cntPer++];
                        iperZ = e;
                        const auto &des = e.des();
                        for (integer nn = 0; nn < des.size(); ++nn) {
                            iperZ.des()[nn] = intfDes[des[nn]];
                        }
                    }
                } else if (e.isSendFromThisProc()) {
                    auto &iperZ = intfPerZone[cntPer++];
                    iperZ = e;
                } else if (e.isThisProcReceiv()) {
                    const auto &des = e.des();
                    auto &iperZ = intfPerZone[cntPer++];
                    iperZ = e;
                    for (integer nn = 0; nn < des.size(); ++nn) {
                        iperZ.des()[nn] = intfDes[des[nn]];
                    }
                }
            }
        }
        CRSMatrixAddrCutInterface CRSIntf(intfCutZone, intfPerZone, intfCell, true);
        CRSAdrrPtr_.reset(new CRSMatrixAddressing(nCells(), leftCells, rightCells, CRSIntf, true));
    }
}

void OpenHurricane::geometryMesh::makeProcessCutTopology() const {
    if (!processCutPtr_) {
        processCutPtr_.reset(new processCutTopology(this->cutZones_, this->ghostCellsType()));
    }
}

void OpenHurricane::geometryMesh::makeProcessPerTopology() const {
    if (!processPerPtr_) {
        processPerPtr_.reset(new processPerTopology(this->perZones_, this->ghostCellsType()));
    }
}

void OpenHurricane::geometryMesh::writeASCIIToBinary(const char *str, fileOsstream &fos) {
    int value = 0;
    while ((*str) != '\0') {
        value = (int)*str;
        fos.write(reinterpret_cast<const char *>(&value), sizeof(int));
        str++;
    }

    char null_char[] = "";
    value = (int)*null_char;
    fos.write(reinterpret_cast<const char *>(&value), sizeof(int));
}

void OpenHurricane::geometryMesh::writeMeshToTecplot(fileOsstream &fos) const {
    if (HurMPI::master()) {
        // Tecplot file header
        writeTecplotHeader(fos, tecplotFormat::GRID);
    }
    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (HurMPI::isThisProc(ip)) {
            // Zone header.
            writeTecplotZoneHeader(fos, tecplotFormat::GRID, HurMPI::parRun());

            // Points
            writePointToTecplot(fos);

            // Cells
            writeCellToTecplot(fos);

            // Face connectivity if necessarily
            writeFaceConnectTecplot(fos);
        }

        // Wait
        HurMPI::barrier(HurMPI::getComm());
    }
}

void OpenHurricane::geometryMesh::writeMeshToTecplot(fileOsstream &fos, const integer fzid) const {
    const auto &name = globalFaceZoneInfo(fzid).fZ().name();

    if (HurMPI::master()) {
        //Pout << "write face zone: " << name.c_str() << std::endl;
        // Tecplot file header
        writeTecplotHeader(fos, fzid, tecplotFormat::GRID);

        //Pout << "write zone header" << std::endl;
        // Zone header.
        writeTecplotZoneHeaderMaster(fos, fzid, tecplotFormat::GRID);
    }
    // Wait
    HurMPI::barrier(HurMPI::getComm());
    // Points
    writePointToTecplot(fos, fzid);

    // Faces
    writeFaceToTecplot(fos, fzid);
}

void OpenHurricane::geometryMesh::writeTecplotHeader(fileOsstream &fos, const short fType) const {
    // Note: It is case-insensitive.
    /* For example:
    ------------------------------------------------------------------
                    Title="OpenHurricane"
                    variables="x","y","z","rho","u","v","w","p","t","qs","Yp"
                    FileType=FULL
    ------------------------------------------------------------------*/

    const iteration &iter = this->Iteration();
    string caseName = iter.configName().name(true);
    std::string titleStr;
    titleStr = "TITLE=\"";
    titleStr += caseName;

    std::string vStr;
    vStr = "variables=\"x\",\"y\",\"z\"";

    std::string fileType;
    fileType = "FILETYPE = ";
    if (fType == tecplotFormat::GRID) {
        titleStr += "-Grid\"";
        fileType += "GRID";
    } else if (fType == tecplotFormat::SOLUTION) {
        titleStr += "-SOLUTION\"";
        fileType += "SOLUTION";
        vStr += ",";
        vStr += outputTitleNameDoc();
    } else {
        titleStr += "\"";
        fileType += "FULL";
        vStr += outputTitleNameDoc();
    }
    fos.os() << titleStr.c_str() << std::endl;
    fos.os() << vStr.c_str() << std::endl;
    fos.os() << fileType.c_str() << std::endl;
}

void OpenHurricane::geometryMesh::writeTecplotHeader(fileOsstream &fos,
                                                     const stringList &outVarName,
                                                     const short fType) const {
    // Note: It is case-insensitive.
    /* For example:
    ------------------------------------------------------------------
                    Title="OpenHurricane"
                    variables="x","y","z","rho","u","v","w","p","t","qs","Yp"
                    FileType=FULL
    ------------------------------------------------------------------*/

    const iteration &iter = this->Iteration();
    string caseName = iter.configName().name(true);
    std::string titleStr;
    titleStr = "TITLE=\"";
    titleStr += caseName;

    std::string vStr;
    vStr = "variables=\"x\",\"y\",\"z\"";

    std::string fileType;
    fileType = "FILETYPE = ";
    if (fType == tecplotFormat::GRID) {
        titleStr += "-Grid\"";
        fileType += "GRID";
    } else if (fType == tecplotFormat::SOLUTION) {
        titleStr += "-SOLUTION\"";
        fileType += "SOLUTION";
        vStr += ",";
        vStr += outputTitleNameDoc(outVarName);
    } else {
        titleStr += "\"";
        fileType += "FULL";
        vStr += outputTitleNameDoc(outVarName);
    }
    fos.os() << titleStr.c_str() << std::endl;
    fos.os() << vStr.c_str() << std::endl;
    fos.os() << fileType.c_str() << std::endl;
}

void OpenHurricane::geometryMesh::writeTecplotHeader(fileOsstream &fos, const integer fzid,
                                                     const short fType) const {
    // Note: It is case-insensitive.
    /* For example:
    ------------------------------------------------------------------
                    Title="OpenHurricane"
                    variables="x","y","z","rho","u","v","w","p","t","qs","Yp"
                    FileType=FULL
    ------------------------------------------------------------------*/

    const iteration &iter = this->Iteration();
    string caseName = iter.configName().name(true);
    std::string titleStr;
    titleStr = "TITLE=\"";
    titleStr += caseName;
    titleStr += "-";
    titleStr += globalFaceZoneInfo(fzid).fZ().name();

    std::string vStr;
    vStr = "variables=\"x\",\"y\",\"z\"";

    std::string fileType;
    fileType = "FILETYPE = ";
    if (fType == tecplotFormat::GRID) {
        titleStr += "-Grid\"";
        fileType += "GRID";
    } else if (fType == tecplotFormat::SOLUTION) {
        titleStr += "-SOLUTION\"";
        fileType += "SOLUTION";
        vStr += ",";
        vStr += outputTitleNameDoc();
    } else {
        titleStr += "\"";
        fileType += "FULL";
        vStr += outputTitleNameDoc();
    }
    fos.os() << titleStr.c_str() << std::endl;
    fos.os() << vStr.c_str() << std::endl;
    fos.os() << fileType.c_str() << std::endl;
}

void OpenHurricane::geometryMesh::writeTecplotHeader(fileOsstream &fos, const integer fzid,
                                                     const stringList &outVarName,
                                                     const short fType) const {
    // Note: It is case-insensitive.
    /* For example:
    ------------------------------------------------------------------
                    Title="OpenHurricane"
                    variables="x","y","z","rho","u","v","w","p","t","qs","Yp"
                    FileType=FULL
    ------------------------------------------------------------------*/

    const iteration &iter = this->Iteration();
    string caseName = iter.configName().name(true);
    std::string titleStr;
    titleStr = "TITLE=\"";
    titleStr += caseName;
    titleStr += "-";
    titleStr += globalFaceZoneInfo(fzid).fZ().name();

    std::string vStr;
    vStr = "variables=\"x\",\"y\",\"z\"";

    std::string fileType;
    fileType = "FILETYPE = ";
    if (fType == tecplotFormat::GRID) {
        titleStr += "-Grid\"";
        fileType += "GRID";
    } else if (fType == tecplotFormat::SOLUTION) {
        titleStr += "-SOLUTION\"";
        fileType += "SOLUTION";
        vStr += ",";
        vStr += outputTitleNameDoc(outVarName);
    } else {
        titleStr += "\"";
        fileType += "FULL";
        vStr += outputTitleNameDoc(outVarName);
    }
    fos.os() << titleStr.c_str() << std::endl;
    fos.os() << vStr.c_str() << std::endl;
    fos.os() << fileType.c_str() << std::endl;
}

void OpenHurricane::geometryMesh::writeTecplotZoneHeader(fileOsstream &fos, const short fileType,
                                                         bool faceConnect) const {
    // Note: It is case-insensitive.
    /* For example:
    ------------------------------------------------------------------
                    Title="OpenHurricane"
                    variables="x","y","z","rho","u","v","w","p","t","qs","Yp"
                    FileType=FULL
                    Zone T="0",nodes=1208,elements=601,datapacking=block,zonetype=fequadrilateral,varlocation=([4-11]=cellcentered),FACENEIGHBORMODE=GLOBALONETOONE,FACENEIGHBORCONNECTIONS=4
    ------------------------------------------------------------------*/

    std::string zoneStr;
    zoneStr = "ZONE T=\"";
    zoneStr += std::to_string(HurMPI::getProcRank());
    zoneStr += "\",NODES=";
    zoneStr += toString(nPoints());
    zoneStr += ",elements=";
    zoneStr += toString(nCells());
    zoneStr += ",datapacking=block,zonetype=febrick";

    if (fileType == tecplotFormat::SOLUTION) {
        zoneStr += ",varlocation = ([";
        if (outputTitleSize() == 1) {
            zoneStr += "1";
        } else {
            zoneStr += "1-";
            zoneStr += std::to_string(outputTitleSize());
        }
        zoneStr += "]=cellcentered)";
    } else if (fileType == tecplotFormat::FULL) {
        zoneStr += ",varlocation = ([";
        if (outputTitleSize() == 1) {
            zoneStr += "4";
        } else {
            zoneStr += "4-";
            zoneStr += std::to_string(outputTitleSize() + 3);
        }
        zoneStr += "]=cellcentered)";
    }

    if (faceConnect) {
        zoneStr += ",FACENEIGHBORMODE=GLOBALONETOONE,";
        zoneStr += "FACENEIGHBORCONNECTIONS=";
        integer connectSize = globalMeshInfoPtr_->globalFaceIndeces().sharedFaceSize();
        zoneStr += toString(connectSize);
    }

    fos.os() << zoneStr.c_str() << std::endl;
}

void OpenHurricane::geometryMesh::writeTecplotZoneHeader(fileOsstream &fos,
                                                         const stringList &outVarName,
                                                         const short fileType,
                                                         bool faceConnect) const {
    // Note: It is case-insensitive.
    /* For example:
    ------------------------------------------------------------------
                    Title="OpenHurricane"
                    variables="x","y","z","rho","u","v","w","p","t","qs","Yp"
                    FileType=FULL
                    Zone T="0",nodes=1208,elements=601,datapacking=block,zonetype=fequadrilateral,varlocation=([4-11]=cellcentered),FACENEIGHBORMODE=GLOBALONETOONE,FACENEIGHBORCONNECTIONS=4
    ------------------------------------------------------------------*/

    std::string zoneStr;
    zoneStr = "ZONE T=\"";
    zoneStr += std::to_string(HurMPI::getProcRank());
    zoneStr += "\",NODES=";
    zoneStr += toString(nPoints());
    zoneStr += ",elements=";
    zoneStr += toString(nCells());
    zoneStr += ",datapacking=block,zonetype=febrick";

    const integer outSize = outputTitleSize(outVarName);
    if (fileType == tecplotFormat::SOLUTION) {
        zoneStr += ",varlocation = ([";
        if (outSize == 1) {
            zoneStr += "1";
        } else {
            zoneStr += "1-";
            zoneStr += std::to_string(outSize);
        }
        zoneStr += "]=cellcentered)";
    } else if (fileType == tecplotFormat::FULL) {
        zoneStr += ",varlocation = ([";
        if (outSize == 1) {
            zoneStr += "4";
        } else {
            zoneStr += "4-";
            zoneStr += std::to_string(outSize + 3);
        }
        zoneStr += "]=cellcentered)";
    }

    if (faceConnect) {
        zoneStr += ",FACENEIGHBORMODE=GLOBALONETOONE,";
        zoneStr += "FACENEIGHBORCONNECTIONS=";
        integer connectSize = globalMeshInfoPtr_->globalFaceIndeces().sharedFaceSize();
        zoneStr += toString(connectSize);
    }

    fos.os() << zoneStr.c_str() << std::endl;
}

void OpenHurricane::geometryMesh::writeTecplotZoneHeaderMaster(fileOsstream &fos,
                                                               const integer fzid,
                                                               const short fileType) const {
    // Note: It is case-insensitive.
    /* For example:
    ------------------------------------------------------------------
                    Title="OpenHurricane"
                    variables="x","y","z","rho","u","v","w","p","t","qs","Yp"
                    FileType=FULL
                    Zone T="0",nodes=1208,elements=601,datapacking=block,zonetype=fequadrilateral,varlocation=([4-11]=cellcentered),FACENEIGHBORMODE=GLOBALONETOONE,FACENEIGHBORCONNECTIONS=4
    ------------------------------------------------------------------*/

    std::string zoneStr;
    zoneStr = "ZONE T=\"";
    zoneStr += std::to_string(fzid);
    zoneStr += "\",NODES=";
    zoneStr += toString(globalFaceZoneInfo(fzid).totalNodes());
    zoneStr += ",elements=";
    zoneStr += toString(globalFaceZoneInfo(fzid).totalFaces());
    zoneStr += ",datapacking=block,zonetype=fequadrilateral";

    if (fileType == tecplotFormat::SOLUTION) {
        zoneStr += ",varlocation = ([";
        if (outputTitleSize() == 1) {
            zoneStr += "1";
        } else {
            zoneStr += "1-";
            zoneStr += std::to_string(outputTitleSize());
        }
        zoneStr += "]=cellcentered)";
    } else if (fileType == tecplotFormat::FULL) {
        zoneStr += ",varlocation = ([";
        if (outputTitleSize() == 1) {
            zoneStr += "4";
        } else {
            zoneStr += "4-";
            zoneStr += std::to_string(outputTitleSize() + 3);
        }
        zoneStr += "]=cellcentered)";
    }

    fos.os() << zoneStr.c_str() << std::endl;
}

void OpenHurricane::geometryMesh::writeTecplotZoneHeaderMaster(fileOsstream &fos,
                                                               const integer fzid,
                                                               const stringList &outVarName,
                                                               const short fileType) const {
    // Note: It is case-insensitive.
    /* For example:
    ------------------------------------------------------------------
                    Title="OpenHurricane"
                    variables="x","y","z","rho","u","v","w","p","t","qs","Yp"
                    FileType=FULL
                    Zone T="0",nodes=1208,elements=601,datapacking=block,zonetype=fequadrilateral,varlocation=([4-11]=cellcentered),FACENEIGHBORMODE=GLOBALONETOONE,FACENEIGHBORCONNECTIONS=4
    ------------------------------------------------------------------*/

    std::string zoneStr;
    zoneStr = "ZONE T=\"";
    zoneStr += std::to_string(fzid);
    zoneStr += "\",NODES=";
    zoneStr += toString(globalFaceZoneInfo(fzid).totalNodes());
    zoneStr += ",elements=";
    zoneStr += toString(globalFaceZoneInfo(fzid).totalFaces());
    zoneStr += ",datapacking=block,zonetype=fequadrilateral";

    const integer outSize = outputTitleSize(outVarName);
    if (fileType == tecplotFormat::SOLUTION) {
        zoneStr += ",varlocation = ([";
        if (outSize == 1) {
            zoneStr += "1";
        } else {
            zoneStr += "1-";
            zoneStr += std::to_string(outSize);
        }
        zoneStr += "]=cellcentered)";
    } else if (fileType == tecplotFormat::FULL) {
        zoneStr += ",varlocation = ([";
        if (outSize == 1) {
            zoneStr += "4";
        } else {
            zoneStr += "4-";
            zoneStr += std::to_string(outSize + 3);
        }
        zoneStr += "]=cellcentered)";
    }

    fos.os() << zoneStr.c_str() << std::endl;
}

void OpenHurricane::geometryMesh::writePointToTecplot(fileOsstream &fos) const {
    points_.writeToStreamWithFactor(fos);
}

void OpenHurricane::geometryMesh::writePointToTecplot(fileOsstream &fos, const integer fzid) const {
    if (!HurMPI::master()) {
        return;
    }

    const auto &p = points_;
    std::stringstream sstr;
    const auto &nodel = globalFaceZoneInfo(fzid).facePoints();
    //const auto& nodeMapl = globalFaceZoneInfo(fzid).writeNodeList();

    //std::cout << "69 = " << nodel[68][0] << " 71 = " << nodel[70][0] << std::endl;

    for (int j = 0; j < point::nElements_; ++j) {
        sstr.clear();
        sstr.str("");
        integer i = 0;
        while (i < nodel.size()) {

            //integer pi = nodeMapl[i++];
            sstr << std::setprecision(feature<real>::precision) << nodel[i++][j] << " ";
            //sstr << std::setprecision(feature<real>::precision) << p[pi][j] * f << " ";
            if (i < nodel.size()) {
                sstr << std::setprecision(feature<real>::precision) << nodel[i++][j] << " ";
                //pi = nodeMapl[i++];
                //sstr << std::setprecision(feature<real>::precision) << p[pi][j] * f << " ";
            }
            if (i < nodel.size()) {
                sstr << std::setprecision(feature<real>::precision) << nodel[i++][j] << " ";
                //pi = nodeMapl[i++];
                //sstr << std::setprecision(feature<real>::precision) << p[pi][j] * f << " ";
            }
            sstr << "\n";
        }
        fos.os() << sstr.str().c_str();
    }
}

void OpenHurricane::geometryMesh::writeCellToTecplot(fileOsstream &fos) const {
    for (integer n = 0; n < nCells(); ++n) {
        const integerList &nl = cells_[n].nodesList();
        if (cells_[n].isTetrahedral()) {
            // For tetrahedral cell the cell connectivity is: N1,N2,N3,N3,N4,N4,N4,N4
            fos.os() << nl[0] + 1 << " " << nl[1] + 1 << " " << nl[2] + 1 << " " << nl[2] + 1 << " "
                     << nl[3] + 1 << " " << nl[3] + 1 << " " << nl[3] + 1 << " " << nl[3] + 1
                     << std::endl;
        } else if (cells_[n].isPyramid()) {
            // For pyramid cell: N1,N2,N3,N4,N5,N5,N5,N5
            fos.os() << nl[0] + 1 << " " << nl[1] + 1 << " " << nl[2] + 1 << " " << nl[3] + 1 << " "
                     << nl[4] + 1 << " " << nl[4] + 1 << " " << nl[4] + 1 << " " << nl[4] + 1
                     << std::endl;
        } else if (cells_[n].isWedge()) {
            // For wedge cell: N1,N2,N3,N3,N4,N5,N6,N6
            fos.os() << nl[0] + 1 << " " << nl[1] + 1 << " " << nl[2] + 1 << " " << nl[2] + 1 << " "
                     << nl[3] + 1 << " " << nl[4] + 1 << " " << nl[5] + 1 << " " << nl[5] + 1
                     << std::endl;
        } else if (cells_[n].isHexahedral()) {
            // For hexahedral cell: N1,N2,N3,N4,N5,N6,N7,N8
            fos.os() << nl[0] + 1 << " " << nl[1] + 1 << " " << nl[2] + 1 << " " << nl[3] + 1 << " "
                     << nl[4] + 1 << " " << nl[5] + 1 << " " << nl[6] + 1 << " " << nl[7] + 1
                     << std::endl;
        }
    }
}

void OpenHurricane::geometryMesh::writeFaceToTecplot(fileOsstream &fos, const integer fzid) const {
    const auto &gfz = globalFaceZoneInfo(fzid);
    const auto &pm = gfz.antiOrderMap();
    const auto &fz = gfz.fZ();
    const auto &f = faces();
    std::string str;
    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (HurMPI::isThisProc(ip)) {
            std::stringstream sstr;
            for (integer fi = fz.firstIndex(); fi <= fz.lastIndex(); ++fi) {
                if (f[fi].isTriangular()) {
                    /*fos.os() << pm.at(f[fi][0]) + 1 << " " << pm.at(f[fi][1]) + 1 << " "
                            << pm.at(f[fi][2]) + 1 << " " << pm.at(f[fi][2]) + 1 << std::endl;*/
                    sstr << pm.at(f[fi][0]) + 1 << " " << pm.at(f[fi][1]) + 1 << " "
                         << pm.at(f[fi][2]) + 1 << " " << pm.at(f[fi][2]) + 1 << std::endl;
                } else if (f[fi].isQuadrilateral()) {
                    /*fos.os() << pm.at(f[fi][0]) + 1 << " " << pm.at(f[fi][1]) + 1 << " "
                            << pm.at(f[fi][2]) + 1 << " " << pm.at(f[fi][3]) + 1 << std::endl;*/
                    sstr << pm.at(f[fi][0]) + 1 << " " << pm.at(f[fi][1]) + 1 << " "
                         << pm.at(f[fi][2]) + 1 << " " << pm.at(f[fi][3]) + 1 << std::endl;
                    /*sstr << pm.at(f[fi][0]) + 1 <<" "<< f[fi][0] << " " << pm.at(f[fi][1]) + 1 << " "<< f[fi][1]<<" "
                            << pm.at(f[fi][2]) + 1 <<" "<< f[fi][2] << " " << pm.at(f[fi][3]) + 1 <<" "<< f[fi][3] << std::endl;*/
                } else {
                    LFatal("Other type of face not supported yet");
                }
            }
            str = sstr.str();
        }
        //HurMPI::barrier(HurMPI::getComm());
    }
    HurMPI::gatherString(str, HurMPI::masterNo(), HurMPI::getComm());
    if (HurMPI::master()) {
        fos.os() << str.c_str();
    }
    HurMPI::barrier(HurMPI::getComm());
}

void OpenHurricane::geometryMesh::writeFaceConnectTecplot(fileOsstream &fos) const {
    if (HurMPI::parRun()) {
        // Face connectivity information(for global-one-to-one):
        //       cz, fz, zr, cr
        //       cz - the cell number in the current zone
        //       fz - the number of the cell face in the current zone
        //       zr - the remote zone number
        //       cr - the cell number of the neighboring cell in the remote zone

        const globalFaceIndex &gfi = globalMeshInfoPtr_->globalFaceIndeces();
        const typename globalFaceIndex::facePairMapType &gfiMap = gfi.facePairMap();
        for (typename globalFaceIndex::facePairMapType::const_iterator iter = gfiMap.begin();
             iter != gfiMap.end(); ++iter) {
            integer fi = iter->first;
            integer cz = faces_[fi].leftCell();
            integer fz = -1;
            for (integer i = 0; i < cells_[cz].faceSize(); i++) {
                if (fi == cells_[cz].facei(i)) {
                    fz = i;
                    break;
                }
            }

            if (fz == -1) {
                LFatal("Shared face map is not right!");
            }

            fz++; // fz+1

            if (cells_[cz].isHexahedral()) {
                // nothing to be done
            } else if (cells_[cz].isTetrahedral()) {
                if (fz == 1) // f1 -> f5
                {
                    fz = 5;
                } else if (fz == 2) // f2 -> f3
                {
                    fz = 3;
                } else if (fz == 3) // f3 -> f1
                {
                    fz = 1;
                } else if (fz == 4) // f4 -> f2
                {
                    fz = 2;
                }
            } else if (cells_[cz].isPyramid()) {
                // nothing to be done
                // Because the face number order of cell is the same as that of hexahedral cell.
            } else if (cells_[cz].isWedge()) {
                if (fz > 3) {
                    fz++; // f4 -> f5; f5 -> f6
                }
            }

            // Because the start of array in C/C++ is zero.
            integer zr = iter->second[0] + 1;
            integer cr = iter->second[3] + 1;

            // Writting cz, fz, zr, cr.
            fos.os() << cz + 1 << " " << fz << " " << zr << " " << cr << std::endl;
        }
    }
}

void OpenHurricane::geometryMesh::writeOutputToTecplot(fileOsstream &fos,
                                                       const short fileType) const {
    if (HurMPI::master()) {
        // Header for tecplot file. Only once.
        writeTecplotHeader(fos, fileType);
    }
    HurMPI::barrier(HurMPI::getComm());
    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (HurMPI::isThisProc(ip)) {
            // Zone header. Only once for every process.
            writeTecplotZoneHeader(fos, fileType, HurMPI::parRun());

            // If the option is a FULL or GRID , write out the points of the mesh.
            if (fileType != tecplotFormat::SOLUTION) {
                writePointToTecplot(fos);
            }

            // Write out the solution.
            if (fileType != tecplotFormat::GRID) {
                writeOutput(fos);
            }

            // Write out cell connectivity and face connectivity if necessarily.
            if (fileType != tecplotFormat::SOLUTION) {
                writeCellToTecplot(fos);
                writeFaceConnectTecplot(fos);
            }
            fos.os().flush();
        }

        // Wait for the previous process finishing the write operation.
        HurMPI::barrier(HurMPI::getComm());
    }
}

void OpenHurricane::geometryMesh::writeOutputToTecplot(fileOsstream &fos,
                                                       const stringList &outVarName,
                                                       const short fileType) const {
    if (HurMPI::master()) {
        // Header for tecplot file. Only once.
        writeTecplotHeader(fos, outVarName, fileType);
    }
    HurMPI::barrier(HurMPI::getComm());
    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (HurMPI::isThisProc(ip)) {
            // Zone header. Only once for every process.
            writeTecplotZoneHeader(fos, outVarName, fileType, HurMPI::parRun());

            // If the option is a FULL or GRID , write out the points of the mesh.
            if (fileType != tecplotFormat::SOLUTION) {
                writePointToTecplot(fos);
            }

            // Write out the solution.
            if (fileType != tecplotFormat::GRID) {
                writeOutput(fos, outVarName);
            }

            // Write out cell connectivity and face connectivity if necessarily.
            if (fileType != tecplotFormat::SOLUTION) {
                writeCellToTecplot(fos);
                writeFaceConnectTecplot(fos);
            }
            fos.os().flush();
        }

        // Wait for the previous process finishing the write operation.
        HurMPI::barrier(HurMPI::getComm());
    }
}

void OpenHurricane::geometryMesh::writeOutputToTecplot(fileOsstream &fos, const integer fzid,
                                                       const short fileType) const {
    const auto &name = globalFaceZoneInfo(fzid).fZ().name();
    if (HurMPI::master()) {
        // Header for tecplot file. Only once.
        writeTecplotHeader(fos, fzid, fileType);

        // Zone header.
        writeTecplotZoneHeaderMaster(fos, fzid, fileType);
    }
    HurMPI::barrier(HurMPI::getComm());

    // If the option is a FULL or GRID , write out the points of the mesh.
    if (fileType != tecplotFormat::SOLUTION) {
        writePointToTecplot(fos, fzid);
    }

    // Write out the solution.
    if (fileType != tecplotFormat::GRID) {
        writeOutput(fos, fzid);
    }

    // Write out cell connectivity and face connectivity if necessarily.
    if (fileType != tecplotFormat::SOLUTION) {

        // Faces
        writeFaceToTecplot(fos, fzid);
    }
}

void OpenHurricane::geometryMesh::writeOutputToTecplot(fileOsstream &fos, const integer fzid,
                                                       const stringList &outVarName,
                                                       const short fileType) const {
    const auto &name = globalFaceZoneInfo(fzid).fZ().name();
    if (HurMPI::master()) {
        // Header for tecplot file. Only once.
        writeTecplotHeader(fos, fzid, outVarName, fileType);

        // Zone header.
        writeTecplotZoneHeaderMaster(fos, fzid, outVarName, fileType);
    }
    HurMPI::barrier(HurMPI::getComm());

    // If the option is a FULL or GRID , write out the points of the mesh.
    if (fileType != tecplotFormat::SOLUTION) {
        writePointToTecplot(fos, fzid);
    }

    // Write out the solution.
    if (fileType != tecplotFormat::GRID) {
        writeOutput(fos, fzid, outVarName);
    }

    // Write out cell connectivity and face connectivity if necessarily.
    if (fileType != tecplotFormat::SOLUTION) {
        // Faces
        writeFaceToTecplot(fos, fzid);
    }
}

void OpenHurricane::geometryMesh::writeMesh(hdf5O &fos, const string &gridGroupName) const {
    if (HurMPI::master()) {
        fos.createGroup(gridGroupName);
        // The dimension of the computation program.
        fos.writeIntegerAttributeToGroup(DIMENSIONSET, string("nDimension"), gridGroupName);

        //	The number of nodes (points) in the complete mesh.
        fos.writeIntegerAttributeToGroup(allMeshNodeNumber_, string("nPoints"), gridGroupName);
        // The number of the zones of nodes (points).
        fos.writeIntegerAttributeToGroup(pointZones().size(), string("nPointZones"), gridGroupName);

        // The number of faces in the complete mesh.
        fos.writeIntegerAttributeToGroup(allMeshFaceNumber_, string("nFaces"), gridGroupName);
        // The number of zones (except cut zones) of faces in the complete mesh.
        fos.writeIntegerAttributeToGroup(nFaceZoneExceptCutZone(), string("nFaceZones"),
                                         gridGroupName);

        // The number of cells in the complete mesh.
        fos.writeIntegerAttributeToGroup(allMeshCellNumber_, string("nCells"), gridGroupName);
        // The number of zones of cells in the complete mesh.
        fos.writeIntegerAttributeToGroup(cellZones().size(), string("nCellZones"), gridGroupName);

        fos.writeStringAttributeToGroup("m", "unit", gridGroupName);

        fos.writeIntegerAttributeToGroup(periodicPairSize_, string("periodicPairSize"),
                                         gridGroupName);

        fos.writeVectorAttributeToGroup(origin_, string("origin"), gridGroupName);
        fos.writeVectorAttributeToGroup(axis_, string("axis"), gridGroupName);
    }
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLDebug("    Calling writePoint(fos, gridGroupName).");
    }
#endif // HUR_FULL_LOGGER
    writePoint(fos, gridGroupName);
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLDebug("    Calling writePoint(fos, gridGroupName). Done.");
    }
#endif // HUR_FULL_LOGGER

    writeCell(fos, gridGroupName);
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLDebug("    Calling writeCell(fos, gridGroupName). Done");
    }
#endif // HUR_FULL_LOGGER

    writeFace(fos, gridGroupName);
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLDebug("    Calling writeFace(fos, gridGroupName). Done");
    }
#endif // HUR_FULL_LOGGER

    writePeriodicPairList(fos, gridGroupName);
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLDebug("    Calling writePeriodicPairList(fos, gridGroupName). Done");
    }
#endif // HUR_FULL_LOGGER

    if (cellLoadWeightsPtr_) {
        writeCellLoadWeights(fos, gridGroupName);
    }
}

void OpenHurricane::geometryMesh::writePoint(hdf5O &fos, const string &gridGroupName) const {
    if (!HurMPI::parRun()) {
        integerList pointZId(pointZones_.size());
        if (pointZones_.size() == 1) {
            string dataName = "pointZone";
            dataName += toString(pointZones_[0].index());
            writePointZoneAttribute(fos, pointZones_[0], points_, gridGroupName, dataName);
            pointZId[0] = pointZones_[0].index();
        } else {
            for (integer pzid = 0; pzid < pointZones_.size(); ++pzid) {
                pointField tpp(pointZones_[pzid].size());
                integer pi = 0;
                for (integer i = pointZones_[pzid].firstIndex(); i <= pointZones_[pzid].lastIndex();
                     ++i) {
                    tpp[pi++] = points_[i];
                }
                string dataName = "pointZone";
                dataName += toString(pointZones_[pzid].index());
                writePointZoneAttribute(fos, pointZones_[pzid], tpp, gridGroupName, dataName);
                pointZId[pzid] = pointZones_[pzid].index();
            }
        }
        fos.write(pointZId, gridGroupName, "pointZoneId");
    } else {
        integerList pointZId(pointZones_.size());
        integer offset = 0;
        for (integer pzid = 0; pzid < pointZones_.size(); ++pzid) {
            pointZone tpz;
            tpz = pointZones_[pzid];
            pointField tpp;

            globalMeshInfoPtr_->globalPointIndeces().getNodeList(pzid, tpp);

            if (HurMPI::master()) {
                tpz.setFirstIndex(offset);
                tpz.setLastIndex(offset + tpp.size() - 1);
                offset += tpp.size();
                string dataName = "pointZone";
                dataName += toString(pointZones_[pzid].index());
                writePointZoneAttribute(fos, tpz, tpp, gridGroupName, dataName);
            }
            pointZId[pzid] = pointZones_[pzid].index();
        }
        if (HurMPI::master()) {
            fos.write(pointZId, gridGroupName, "pointZoneId");
        }
    }
}

void OpenHurricane::geometryMesh::writeCell(hdf5O &fos, const string &gridGroupName) const {
    if (!HurMPI::parRun()) {
        integerList cellZId(cellZones_.size());
        for (integer czid = 0; czid < cellZones_.size(); ++czid) {
            string dataName = "cellZone";
            dataName += toString(cellZones_[czid].index());
            integerList cellTypes;

            string numName = "cellOriginIndex";
            numName += toString(cellZones_[czid].index());
            integerList cellOriNo;
            cellOriNo.resize(cellZones_[czid].size());
            integer ci = 0;
            for (integer i = cellZones_[czid].firstIndex(); i <= cellZones_[czid].lastIndex();
                 ++i) {
                cellOriNo[ci++] = cells_[i].orignIndex();
            }

            if (cellZones_[czid].isMixed()) {
                cellTypes.resize(cellZones_[czid].size());
                ci = 0;
                for (integer i = cellZones_[czid].firstIndex(); i <= cellZones_[czid].lastIndex();
                     ++i) {
                    cellTypes[ci++] = cells_[i].shapeType();
                }
            } else {
                cellTypes.resize(1, cellZones_[czid].shapeType());
            }
            fos.write(cellTypes, gridGroupName, dataName);
            writeCellZoneAttribute(fos, cellZones_[czid], gridGroupName, dataName);
            fos.write(cellOriNo, gridGroupName, numName);
            cellZId[czid] = cellZones_[czid].index();
        }
        fos.write(cellZId, gridGroupName, "cellZoneId");
    } else {
        integer offset = 0;
        integerList cellZId(cellZones_.size());

        // Create the decomposing cell group of the grid.
        string decomposingGroup = gridGroupName + "Decompose";
        if (HurMPI::master()) {
            fos.createGroup(decomposingGroup);
        }

        for (integer czid = 0; czid < cellZones_.size(); ++czid) {
            cellZone tcz;
            tcz = cellZones_[czid];

            // Get the total number of cells in this cell zone from all processes.
            integer csize = cellZones_[czid].size();
            HurMPI::reduce(csize, MPI_SUM);

            string dataName = "cellZone";
            if (HurMPI::master()) {
                tcz.setFirstIndex(offset);
                tcz.setLastIndex(offset + csize - 1);
                offset += csize;
                dataName += toString(cellZones_[czid].index());
            }

            string numName = "cellOriginIndex";
            if (HurMPI::master()) {
                numName += toString(cellZones_[czid].index());
            }
            integerList cellTypes;

            integerList cellOriNo;
            integerList cellOriginNum(cellZones_[czid].size());

            integer ci = 0;
            for (integer i = cellZones_[czid].firstIndex(); i <= cellZones_[czid].lastIndex();
                 ++i) {
                cellOriginNum[ci++] = cells_[i].orignIndex();
            }
            integerList nSizeL(HurMPI::getProcSize(), Zero);
            integerList displs;
            if (HurMPI::master()) {
                displs.resize(HurMPI::getProcSize(), Zero);
            }
            nSizeL[HurMPI::getProcRank()] = cellOriginNum.size();
            HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());
            if (HurMPI::master()) {
                cellOriNo.resize(csize);
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
                }
            }
            /*HurMPI::gatherv
            (
                    cellOriginNum.data(),
                    cellOriginNum.size(),
                    feature<integer>::MPIType,
                    cellOriNo.data(),
                    nSizeL.data(),
                    displs.data(),
                    feature<integer>::MPIType,
                    HurMPI::masterNo(),
                    HurMPI::getComm()
            );*/

            HurMPI::Request request;
            HurMPI::igatherv(cellOriginNum.data(), cellOriginNum.size(), feature<integer>::MPIType,
                             cellOriNo.data(), nSizeL.data(), displs.data(),
                             feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm(),
                             &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (cellZones_[czid].isMixed()) {
                integerList tmpCellTypes(cellZones_[czid].size());

                ci = 0;
                for (integer i = cellZones_[czid].firstIndex(); i <= cellZones_[czid].lastIndex();
                     ++i) {
                    tmpCellTypes[ci++] = cells_[i].shapeType();
                }
                if (HurMPI::master()) {
                    cellTypes.resize(csize);
                }
                /*HurMPI::gatherv
                (
                        tmpCellTypes.data(),
                        tmpCellTypes.size(),
                        feature<integer>::MPIType,
                        cellTypes.data(),
                        nSizeL.data(),
                        displs.data(),
                        feature<integer>::MPIType,
                        HurMPI::masterNo(),
                        HurMPI::getComm()
                );*/

                HurMPI::Request request;
                HurMPI::igatherv(tmpCellTypes.data(), tmpCellTypes.size(),
                                 feature<integer>::MPIType, cellTypes.data(), nSizeL.data(),
                                 displs.data(), feature<integer>::MPIType, HurMPI::masterNo(),
                                 HurMPI::getComm(), &request);
                HurMPI::wait(&request, MPI_STATUSES_IGNORE);
            } else {
                cellTypes.resize(1, cellZones_[czid].shapeType());
            }

            if (HurMPI::master()) {
                fos.write(cellTypes, gridGroupName, dataName);
                writeCellZoneAttribute(fos, tcz, gridGroupName, dataName);
                fos.write(cellOriNo, gridGroupName, numName);

                fos.write(nSizeL, decomposingGroup, dataName);
            }
            cellZId[czid] = cellZones_[czid].index();
        }
        if (HurMPI::master()) {
            fos.write(cellZId, gridGroupName, "cellZoneId");
        }
    }
}

void OpenHurricane::geometryMesh::writeFace(hdf5O &fos, const string &gridGroupName) const {
    if (!HurMPI::parRun()) {
        integerList faceZId;
        for (integer fzid = 0; fzid < faceZones_.size(); ++fzid) {
            if (faceZones_[fzid].isCutFace()) {
                continue;
            }
            faceZId.append(faceZones_[fzid].index());
            string dataName = "faceZone";
            dataName += toString(faceZones_[fzid].index());
            integerArrayArray faceConnect(faceZones_[fzid].size());
            bool isMix = false;
            if (faceZones_[fzid].isMixed()) {
                isMix = true;
            }
            integer fi = 0;
            for (integer i = faceZones_[fzid].firstIndex(); i <= faceZones_[fzid].lastIndex();
                 ++i) {
                integer offset = 0;
                if (isMix) {
                    faceConnect[fi].resize(faces_[i].size() + 1 + 2);
                    faceConnect[fi][0] = faces_[i].size();
                    offset = 1;
                } else {
                    faceConnect[fi].resize(faces_[i].size() + 2);
                }
                for (integer j = 0; j < faces_[i].size(); ++j) {
                    faceConnect[fi][j + offset] = faces_[i][j];
                }
                faceConnect[fi][faces_[i].size() + offset] = faces_[i].leftCell();
                if (faces_[i].rightCell() >= nCells()) {
                    faceConnect[fi][faces_[i].size() + offset + 1] = -1;
                } else {
                    faceConnect[fi][faces_[i].size() + offset + 1] = faces_[i].rightCell();
                }
                fi++;
            }

            fos.writeArrayArray(faceConnect, gridGroupName, dataName);
            writeFaceZoneAttribute(fos, faceZones_[fzid], gridGroupName, dataName);
        }
        fos.write(faceZId, gridGroupName, "faceZoneId");
    } else {
        //auto cutZL = getCutZonesIdArrays();
        integer _offset = 0;
        integerList faceZId;
        for (integer fzid = 0; fzid < faceZones_.size(); ++fzid) {
            if (faceZones_[fzid].isCutFace()) {
                continue;
            }
            faceZId.append(faceZones_[fzid].index());
            string dataName;
            if (HurMPI::master()) {
                dataName = "faceZone";
                dataName += toString(faceZones_[fzid].index());
            }
            integerArrayArray faceConnect;
            faceZone tfz;
            globalMeshInfoPtr_->gatherFaceZone(fzid, tfz, faceConnect, _offset);

            if (HurMPI::master()) {
                fos.writeArrayArray(faceConnect, gridGroupName, dataName);
                writeFaceZoneAttribute(fos, tfz, gridGroupName, dataName);
            }
        }
        if (HurMPI::master()) {
            fos.write(faceZId, gridGroupName, "faceZoneId");
        }
    }
}

void OpenHurricane::geometryMesh::writePeriodicPairList(hdf5O &fos,
                                                        const string &gridGroupName) const {
    if (periodicPairSize_ == 0) {
        return;
    }
#ifdef HUR_DEBUG
    Pout << "geometryMesh::writePeriodicPairList : Writting periodic pair list" << std::endl;
#endif // HUR_DEBUG

    periodicPairList ppl;
    getPeriodicPairList(perZones_, ppl);

    if (HurMPI::master()) {
        integerList pplpZ(ppl.size());
        for (integer i = 0; i < ppl.size(); ++i) {
            string dataName;
            dataName = "periodicPair";
            dataName += toString(ppl[i].periodicZone());
            fos.write(ppl[i], gridGroupName, dataName);
            writePeriodicPairAttribute(fos, ppl[i], gridGroupName, dataName);
            pplpZ[i] = ppl[i].periodicZone();
        }
        fos.write(pplpZ, gridGroupName, "periodicPairList");
    }
#ifdef HUR_DEBUG
    Pout << "geometryMesh::writePeriodicPairList : Finish writting periodic "
            "pair list"
         << std::endl;
#endif // HUR_DEBUG
}

void OpenHurricane::geometryMesh::writePointZoneAttribute(hdf5O &fos, const pointZone &pz,
                                                          const pointField &p,
                                                          const string &gridGroupName,
                                                          const string &dataName) const {
    fos.write(p, p.size(), gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(pz.index(), "index", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(pz.firstIndex(), "firstIndex", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(pz.lastIndex(), "lastIndex", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(pz.type(), "type", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(pz.ND(), "nDimension", gridGroupName, dataName);

    fos.writeStringAttributeToDataset(pz.name(), "name", gridGroupName, dataName);
}

void OpenHurricane::geometryMesh::writeCellZoneAttribute(hdf5O &fos, const cellZone &cz,
                                                         const string &gridGroupName,
                                                         const string &dataName) const {
    fos.writeIntegerAttributeToDataset(cz.index(), "index", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(cz.firstIndex(), "firstIndex", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(cz.lastIndex(), "lastIndex", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(cz.type(), "type", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(cz.shapeType(), "cellType", gridGroupName, dataName);

    fos.writeStringAttributeToDataset(cz.name(), "name", gridGroupName, dataName);
}

void OpenHurricane::geometryMesh::writeFaceZoneAttribute(hdf5O &fos, const faceZone &fz,
                                                         const string &gridGroupName,
                                                         const string &dataName) const {
    fos.writeIntegerAttributeToDataset(fz.index(), "index", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(fz.firstIndex(), "firstIndex", gridGroupName, dataName);
    fos.writeIntegerAttributeToDataset(fz.lastIndex(), "lastIndex", gridGroupName, dataName);

    fos.writeIntegerAttributeToDataset(fz.bcType(), "bcType", gridGroupName, dataName);

    fos.writeIntegerAttributeToDataset(fz.faceType(), "faceType", gridGroupName, dataName);

    fos.writeStringAttributeToDataset(fz.name(), "name", gridGroupName, dataName);
}

void OpenHurricane::geometryMesh::writePeriodicPairAttribute(hdf5O &fos, const periodicPair &pp,
                                                             const string &gridGroupName,
                                                             const string &dataName) const {
    fos.writeIntegerAttributeToDataset(pp.periodicZone(), "periodicZone", gridGroupName, dataName);

    fos.writeIntegerAttributeToDataset(pp.shadowZone(), "shadowZone", gridGroupName, dataName);

    fos.writeIntegerAttributeToDataset((integer)pp.type(), "type", gridGroupName, dataName);
}

hur_nodiscard OpenHurricane::integerArrayArray
OpenHurricane::geometryMesh::getCutZonesIdArrays() const {
    integerArrayArray cutId(faceZones_.size());
    for (integer fzid = 0; fzid < faceZones_.size(); ++fzid) {
        if (faceZones_[fzid].isCutFace()) {
            for (integer i = 0; i < faceZones_.size(); ++i) {
                if (faceZones_[i].index() == -faceZones_[fzid].index()) {
                    cutId[i].append(fzid);
                }
            }
        }
    }
    return cutId;
}

hur_nodiscard OpenHurricane::pointField OpenHurricane::geometryMesh::gatherAllPoints() const {
    if (HurMPI::parRun()) {
        pointField pp;
        for (integer pzid = 0; pzid < pointZones_.size(); ++pzid) {
            pointField tpp;
            globalMeshInfo().globalPointIndeces().getNodeList(pzid, tpp);
            if (HurMPI::master()) {
                pp.append(tpp);
            }
        }
        return pp;
    } else {
        if (pointZones_.size() == 1) {
            return points_;
        } else {
            pointField pp;
            for (integer pzid = 0; pzid < pointZones_.size(); ++pzid) {
                pointField tpp;
                globalMeshInfo().globalPointIndeces().getNodeList(pzid, tpp);
                pp.append(tpp);
            }
            return pp;
        }
    }
    return pointField();
}

void OpenHurricane::geometryMesh::writeCellLoadWeights(hdf5O &fos,
                                                       const string &gridGroupName) const {
    if (!cellLoadWeightsPtr_) {
        return;
    }

    if (!HurMPI::parRun()) {
        for (integer czid = 0; czid < cellZones_.size(); ++czid) {
            string numName = "cellLoadWeights";
            numName += toString(cellZones_[czid].index());
            integerList cellLoadWgt;
            cellLoadWgt.resize(cellZones_[czid].size());
            integer ci = 0;
            for (integer i = cellZones_[czid].firstIndex(); i <= cellZones_[czid].lastIndex();
                 ++i) {
                cellLoadWgt[ci++] = cellLoadWeights()[i];
            }
            fos.write(cellLoadWgt, gridGroupName, numName);
        }

    } else {
        for (integer czid = 0; czid < cellZones_.size(); ++czid) {
            // Get the total number of cells in this cell zone from all processes.
            integer csize = cellZones_[czid].size();
            HurMPI::reduce(csize, MPI_SUM);

            string numName = "cellLoadWeights";
            if (HurMPI::master()) {
                numName += toString(cellZones_[czid].index());
            }

            integerList cellLoadWgt;
            integerList cellLoadWgtNum(cellZones_[czid].size());

            integer ci = 0;
            for (integer i = cellZones_[czid].firstIndex(); i <= cellZones_[czid].lastIndex();
                 ++i) {
                cellLoadWgtNum[ci++] = cellLoadWeights()[i];
            }
            integerList nSizeL(HurMPI::getProcSize(), Zero);
            integerList displs;
            if (HurMPI::master()) {
                displs.resize(HurMPI::getProcSize(), Zero);
            }
            nSizeL[HurMPI::getProcRank()] = cellLoadWgtNum.size();
            HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());
            if (HurMPI::master()) {
                cellLoadWgt.resize(csize);
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
                }
            }

            HurMPI::Request request;
            HurMPI::igatherv(cellLoadWgtNum.data(), cellLoadWgtNum.size(),
                             feature<integer>::MPIType, cellLoadWgt.data(), nSizeL.data(),
                             displs.data(), feature<integer>::MPIType, HurMPI::masterNo(),
                             HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                fos.write(cellLoadWgt, gridGroupName, numName);
            }
        }
    }
}

void OpenHurricane::geometryMesh::calcNumberOfPerFaces() const {
    if (!totalPerFacesPtr_) {
        integer n = 0;
        for (integer faceI = 0; faceI < faceZones_.size(); ++faceI) {
            if (faceZones_[faceI].isPeriodic() || faceZones_[faceI].isPeriodicShadow()) {
                n += faceZones_[faceI].size();
            }
        }
        totalPerFacesPtr_.reset(new integer(n));
    }
}

hur_nodiscard const OpenHurricane::globalMesh &OpenHurricane::geometryMesh::globalMeshInfo() const {
    if (!globalMeshInfoPtr_) {
        LFatal("Attempt to access a null global mesh pointer.");
    }
    return *globalMeshInfoPtr_;
}

hur_nodiscard const OpenHurricane::globalFaceZoneIndex &
OpenHurricane::geometryMesh::globalFaceZoneInfo(const integer fzid) const {
    if (globalFaceZoneL_.size() == 0) {
        globalFaceZoneL_.append(new globalFaceZoneIndex(*this, fzid));
        return globalFaceZoneL_[0];
    } else {
        for (integer i = 0; i < globalFaceZoneL_.size(); ++i) {
            if (fzid == globalFaceZoneL_[i].id()) {
                return globalFaceZoneL_[i];
            }
        }
        globalFaceZoneL_.append(new globalFaceZoneIndex(*this, fzid));
        return globalFaceZoneL_[globalFaceZoneL_.size() - 1];
    }
}

hur_nodiscard OpenHurricane::vectorArray OpenHurricane::geometryMesh::allCellCentre() const {
    integerArray nSizeL;
    integerArray displs;
    integer allSize;
    vectorArray rootV;
    realArray rootF;
    realArray cmpF;
    if (HurMPI::parRun()) {
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
        if (HurMPI::master()) {
            displs.resize(HurMPI::getProcSize(), Zero);
        }
        HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());
        allSize = 0;
        if (HurMPI::master()) {
            for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
            }
            for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                allSize += nSizeL[ip];
            }
        }
        cmpF.resize(this->internalArraySize());
        if (HurMPI::master()) {
            rootF.resize(allSize);
            rootV.resize(allSize);
        }
        HurMPI::barrier(HurMPI::getComm());

        for (integer i = 0; i < feature<vector>::nElements_; ++i) {
            for (integer j = 0; j < this->internalArraySize(); ++j) {
                cmpF[j] = this->cellCentre()[j][i];
            }

            HurMPI::Request request;
            HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                             rootF.data(), nSizeL.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            for (integer j = 0; j < rootF.size(); ++j) {
                rootV[j][i] = rootF[j];
            }
        }
    } else {
        rootV.resize(this->internalArraySize());
        for (integer i = 0; i < this->internalArraySize(); i++) {
            rootV[i] = this->cellCentre()[i];
        }
    }
    return rootV;
}

bool OpenHurricane::geometryMesh::getInterpolationSource(const hdf5I &fos) {
    bool isExist = true;
    if (HurMPI::master()) {
        if (!fos.exist("cellCentre")) {
            isExist = false;
        }
    }
    HurMPI::bcast(&isExist, 1, feature<bool>::MPIType);
    if (!isExist) {
        return false;
    }

    vectorArray sor;
    fos.read(sor, "cellCentre");

    const controller &subCont = cont().subController("iteration");

    integer k_ = -1;
    if (subCont.found("nearNodes")) {
        k_ = subCont.findType<integer>("nearNodes", k_);
    } else {
        k_ = 1;
    }

    boundBox ipBB(this->cellCentre());
    vector span = ipBB.span();
    realArray minMax(2 * feature<vector>::nElements_);
    for (integer i = 0; i < feature<vector>::nElements_; i++) {
        minMax[2 * i] = ipBB.min()[i] - 0.1 * span[i];
        minMax[2 * i + 1] = ipBB.max()[i] + 0.1 * span[i];
    }

    realArray bbGather(HurMPI::getProcSize() * 2 * feature<vector>::nElements_);
    HurMPI::gather(minMax.data(), minMax.size(), feature<real>::MPIType, bbGather.data(),
                   minMax.size(), feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

    List<boundBox> allBB(HurMPI::getProcSize());
    integerArrayArray contains(HurMPI::getProcSize());

    if (HurMPI::master()) {
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            integer displ = ip * 2 * feature<vector>::nElements_;
            for (integer i = 0; i < feature<vector>::nElements_; i++) {
                allBB[ip].min()[i] = bbGather[2 * i + displ];
                allBB[ip].max()[i] = bbGather[2 * i + 1 + displ];
            }
        }

        integerArray count(HurMPI::getProcSize(), Zero);
        for (integer i = 0; i < sor.size(); i++) {
            for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
                if (allBB[ip].contains(sor[i])) {
                    count[ip]++;
                    //break;
                }
            }
        }
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            contains[ip].resize(count[ip]);
            count[ip] = 0;
        }
        for (integer i = 0; i < sor.size(); i++) {
            for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
                if (allBB[ip].contains(sor[i])) {
                    contains[ip][count[ip]++] = i;
                    //break;
                }
            }
        }
    }

    vectorArray sorCC;
    relayScatterFunc::relayScatter(contains, sor, sorCC);

    vectorArray tarCC(nCells());
    for (integer n = 0; n < nCells(); n++) {
        tarCC[n] = this->cellCentre()[n];
    }

    kNNGrid knn(subCont, sorCC, tarCC);

    knn.nearestNbr();

    integer use = 0;
    integerArray color(sorCC.size(), -1);
    for (integer n = 0; n < knn.nearestNbr().size(); n++) {
        for (integer k = 0; k < knn.k(); k++) {
            integer index = knn.nearestNbr()[n][k];
            if (color[index] == -1) {
                color[index] = use++;
            }
        }
    }

    sorCellCentrePtr_.reset(new vectorArray());
    tarOfProc_.reset(new integerArrayArray(HurMPI::getProcSize()));
    sorKNN_.reset(new integerArrayArray(this->nCells()));

    sorCellCentrePtr_->transfer(sorCC);
    for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
        (*tarOfProc_)[ip].transfer(contains[ip]);
    }
    for (integer i = 0; i < this->nCells(); i++) {
        (*sorKNN_)[i].transfer(knn.nearestNbr()[i]);
    }

    return true;
}

void OpenHurricane::meshCheckers::minMaxFaceArea(const vectorArray &faceArea,
                                                 const faceZoneList &fZL, real &minArea,
                                                 real &maxArea) {
#ifdef HUR_DEBUG
    Pout << "    Info: Checking face area..." << std::endl;
#endif // HUR_DEBUG

    realArray faceMagAreaF(mag(faceArea));
    minArea = large;
    maxArea = -large;

    for (integer fzI = 0; fzI < fZL.size(); ++fzI) {
        const faceZone &fzs = fZL[fzI];
        for (integer fI = fzs.firstIndex(); fI <= fzs.lastIndex(); ++fI) {
            if (faceMagAreaF[fI] < veryTiny) {
                if (report) {
                    LInfo(" Zero of negative face area detected in face zone: %s. Normalized face "
                          "area: %e",
                          fzs.name().c_str(), faceMagAreaF[fI]);
                }
            }
            minArea = min(minArea, faceMagAreaF[fI]);
            maxArea = max(maxArea, faceMagAreaF[fI]);
        }
    }

    if (minArea < veryTiny) {
        LFatal("Zero of negative face area are detected.");
    }
}

void OpenHurricane::meshCheckers::minMaxCellVol(const realArray &cellVol, const cellZoneList &cZL,
                                                real &minVol, real &maxVol, real &totalVol) {
#ifdef HUR_DEBUG
    Pout << "    Info: Checking cell volume..." << std::endl;
#endif // HUR_DEBUG

    minVol = large;
    maxVol = -large;
    totalVol = Zero;

    for (integer czI = 0; czI < cZL.size(); ++czI) {
        const cellZone &czs = cZL[czI];
        for (integer cI = czs.firstIndex(); cI <= czs.lastIndex(); ++cI) {
            if (cellVol[cI] < veryTiny) {
                if (report) {
                    LInfo("Zero of negative cell volume detected in cell zone: %s. Normalized cell "
                          "volume: %e",
                          czs.name().c_str(), cellVol[cI]);
                }
            }
            minVol = min(minVol, cellVol[cI]);
            maxVol = max(maxVol, cellVol[cI]);
            totalVol += cellVol[cI];
        }
    }

    if (minVol < veryTiny) {
        LFatal("Zero of negative cell volumes are detected.");
    }
}

void OpenHurricane::meshCheckers::minMaxCellCentre(const pointField &cellCntr,
                                                   const cellZoneList &cZL, vector &mincCComponents,
                                                   vector &maxcCComponents) {
#ifdef HUR_DEBUG
    Pout << "    Info: Checking cell centre..." << std::endl;
#endif // HUR_DEBUG

    mincCComponents = vector(large);
    maxcCComponents = vector(-large);

    for (integer cI = 0; cI < cZL.size(); ++cI) {
        {
            mincCComponents = componentMin(mincCComponents, cellCntr[cI]);
            maxcCComponents = componentMax(maxcCComponents, cellCntr[cI]);
        }
    }
}

void OpenHurricane::meshCheckers::minMaxPointDomain(const pointField &point, const integer nPoints,
                                                    vector &minPointCompo, vector &maxPointCompo) {
#ifdef HUR_DEBUG
    Pout << "    Info: Checking mesh domain..." << std::endl;
#endif // HUR_DEBUG
    minPointCompo = vector(large);
    maxPointCompo = vector(-large);
    for (integer pI = 0; pI < nPoints; ++pI) {
        minPointCompo = componentMin(minPointCompo, point[pI]);
        maxPointCompo = componentMax(maxPointCompo, point[pI]);
    }
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::meshCheckers::faceSkewness(
    const geometryMesh &mesh, const pointField &points, const vectorArray &faceCentres,
    const vectorArray &faceAreas, const vectorArray &cellCentres) {
    const faceList &faces = mesh.faces();
    realArray fSkew(mesh.nFaces());

    for (integer faceI = 0; faceI < mesh.nFaces(); ++faceI) {
        const auto &cl = faces[faceI].leftCell();
        const auto &cr = faces[faceI].rightCell();
        vector Cpf = faceCentres[faceI] - cellCentres[cl];
        vector d = cellCentres[cr] - cellCentres[cl];

        vector sv = Cpf - ((faceAreas[faceI] * Cpf) / ((faceAreas[faceI] * d) + rootVeryTiny)) * d;
        vector svHat = sv / (mag(sv) + rootVeryTiny);

        real fd = 0.2 * mag(d) + rootVeryTiny;

        const face &f = faces[faceI];
        for (integer pointI = 0; pointI < f.size(); ++pointI) {
            fd = max(fd, mag(svHat * (points[f[pointI]] - faceCentres[faceI])));
        }
        fSkew[faceI] = mag(sv) / fd;
    }
    return fSkew;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::meshCheckers::volRatio(const geometryMesh &mesh, const realArray &volume) {
    const faceList &faces = mesh.faces();

    realArray ratio(mesh.nFaces(), real(1.0));

    for (integer faceI = 0; faceI < mesh.nFaces(); ++faceI) {
        const auto &cl = faces[faceI].leftCell();
        const auto &cr = faces[faceI].rightCell();

        ratio[faceI] = OpenHurricane::max(volume[cl], volume[cr]) /
                       (OpenHurricane::min(volume[cl], volume[cr]) + veryTiny);
    }

    return ratio;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::meshCheckers::faceOrthogonality(
    const geometryMesh &mesh, const vectorArray &faceAreas, const vectorArray &faceCentres,
    const vectorArray &cellCentres) {
    const faceList &faces = mesh.faces();
    realArray ortho(mesh.nFaces(), 1.0);

    for (integer faceI = 0; faceI < mesh.nFaces(); ++faceI) {
        const auto &cl = faces[faceI].leftCell();
        const auto &cr = faces[faceI].rightCell();

        vector c = cellCentres[cl] - cellCentres[cr];

        vector f = cellCentres[cl] - faceCentres[faceI];

        real or1 = (c * faceAreas[faceI]) / (mag(c) * mag(faceAreas[faceI]) + rootVeryTiny);
        real or2 = (f * faceAreas[faceI]) / (mag(f) * mag(faceAreas[faceI]) + rootVeryTiny);

        ortho[faceI] = min(or1, or2);
    }

    return ortho;
}
