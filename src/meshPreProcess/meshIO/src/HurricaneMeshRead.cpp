/*!
 * \file HurricaneMeshRead.cpp
 * \brief Main subroutines for reading origin mesh.
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

#include "HurricaneMeshRead.hpp"
#include "hdf5I.hpp"

namespace OpenHurricane {
    createClassNameStr(HurricaneMeshRead, "OpenHurricane");
    registerObjFty(originMeshRead, HurricaneMeshRead, controller);
} // namespace OpenHurricane

void OpenHurricane::HurricaneMeshRead::readOriginAndAxis(const hdf5I &fos,
                                                         const string &gridGroupName) {
    if (!fos.exist(gridGroupName, "origin")) {
        fos.readVectorAttributeFromGroup(origin_, "origin", gridGroupName);
    }
    if (!fos.exist(gridGroupName, "axis")) {
        fos.readVectorAttributeFromGroup(axis_, "axis", gridGroupName);
    }
}

void OpenHurricane::HurricaneMeshRead::readCellZonesHurricane(const hdf5I &fos,
                                                              const string &gridGroupName,
                                                              const integer nCZ) {
    cellZones_.resize(nCZ);
    originCellIndex_.resize(nCZ);
    integerList cellZId;
    fos.read(cellZId, gridGroupName, "cellZoneId");

    for (integer i = 0; i < nCZ; ++i) {
        string dataName = "cellZone";
        dataName += toString(cellZId[i]);

        integer id;
        fos.readIntegerAttributeFromDataset(id, "index", gridGroupName, dataName);
        cellZones_[i].setIndex(id);

        integer fi;
        fos.readIntegerAttributeFromDataset(fi, "firstIndex", gridGroupName, dataName);
        cellZones_[i].setFirstIndex(fi);

        integer li;
        fos.readIntegerAttributeFromDataset(li, "lastIndex", gridGroupName, dataName);
        cellZones_[i].setLastIndex(li);

        integer t;
        fos.readIntegerAttributeFromDataset(t, "type", gridGroupName, dataName);
        cellZones_[i].setType(t);

        integer ct;
        fos.readIntegerAttributeFromDataset(ct, "cellType", gridGroupName, dataName);
        cellZones_[i].setCellType(ct);

        std::string n;
        fos.readStringAttributeFromDataset(n, "name", gridGroupName, dataName);

        originCellIndex_[i].resize(cellZones_[i].size());
        string dataNameOri = "cellOriginIndex";
        dataNameOri += toString(cellZId[i]);
        fos.read(originCellIndex_[i], gridGroupName, dataNameOri);

        cellZones_[i].resetName(n);
        if (cellZones_[i].isMixed()) {
            integerList typeList;
            fos.read(typeList, gridGroupName, dataName);
            for (integer j = cellZones_[i].firstIndex(); j <= cellZones_[i].lastIndex(); ++j) {
                integer nn = j - cellZones_[i].firstIndex();
                cells_[j].setType(typeList[j]);
                cells_[j].setOrignIndex(originCellIndex_[i][nn]);
            }
        } else {
            for (integer j = cellZones_[i].firstIndex(); j <= cellZones_[i].lastIndex(); ++j) {
                integer nn = j - cellZones_[i].firstIndex();
                cells_[j].setType(cellZones_[i].shapeType());
                cells_[j].setOrignIndex(originCellIndex_[i][nn]);
            }
        }
    }
}

void OpenHurricane::HurricaneMeshRead::setCellsHurricane(const hdf5I &fos,
                                                         const string &gridGroupName,
                                                         const integer nC) {
    cells_.resize(nC);
}

void OpenHurricane::HurricaneMeshRead::readFaceZonesHurricane(const hdf5I &fos,
                                                              const string &gridGroupName,
                                                              const integer nFZ) {
    faceZones_.resize(nFZ);
    integerList faceZId;
    fos.read(faceZId, gridGroupName, "faceZoneId");
    globalInterior_ = 0;
    for (integer i = 0; i < nFZ; ++i) {
        string dataName = "faceZone";
        dataName += toString(faceZId[i]);

        integer id;
        fos.readIntegerAttributeFromDataset(id, "index", gridGroupName, dataName);
        faceZones_[i].setIndex(id);

        integer fi;
        fos.readIntegerAttributeFromDataset(fi, "firstIndex", gridGroupName, dataName);
        faceZones_[i].setFirstIndex(fi);

        integer li;
        fos.readIntegerAttributeFromDataset(li, "lastIndex", gridGroupName, dataName);
        faceZones_[i].setLastIndex(li);

        integer bcT;
        fos.readIntegerAttributeFromDataset(bcT, "bcType", gridGroupName, dataName);
        faceZones_[i].setBcType(bcT);

        integer fT;
        fos.readIntegerAttributeFromDataset(fT, "faceType", gridGroupName, dataName);
        faceZones_[i].setFaceType(fT);

        std::string n;
        fos.readStringAttributeFromDataset(n, "name", gridGroupName, dataName);
        faceZones_[i].resetName(n);

        integerArrayArray faceConnect;
        fos.readArrayArray(faceConnect, gridGroupName, dataName);
        readFacesHurricane(fos, gridGroupName, faceZones_[i], faceConnect);
        if (faceZones_[i].isInterior()) {
            globalInterior_ += faceZones_[i].size();
        }
    }
}

void OpenHurricane::HurricaneMeshRead::setFacesHurricane(const hdf5I &fos,
                                                         const string &gridGroupName,
                                                         const integer nF) {
    faces_.resize(nF);
}

void OpenHurricane::HurricaneMeshRead::readFacesHurricane(const hdf5I &fos,
                                                          const string &gridGroupName,
                                                          const faceZone &fz,
                                                          const integerArrayArray &faceConnect) {
    if (fz.isMixed() || fz.isPolygonal()) {
        integer fi = 0;
        for (integer j = fz.firstIndex(); j <= fz.lastIndex(); ++j) {
            faces_[j].resize(faceConnect[fi][0]);
            for (integer i = 0; i < faces_[j].size(); ++i) {
                faces_[j][i] = faceConnect[fi][i + 1];
            }
            faces_[j].leftCell() = faceConnect[fi][faces_[j].size() + 1];
            faces_[j].rightCell() = faceConnect[fi][faces_[j].size() + 2];
            faces_[j].setBCType(fz.bcType());
            fi++;
        }
    } else if (fz.isLinear()) {
        integer fi = 0;
        for (integer j = fz.firstIndex(); j <= fz.lastIndex(); ++j) {
            faces_[j].resize(2);
            for (integer i = 0; i < faces_[j].size(); ++i) {
                faces_[j][i] = faceConnect[fi][i];
            }
            faces_[j].leftCell() = faceConnect[fi][faces_[j].size()];
            faces_[j].rightCell() = faceConnect[fi][faces_[j].size() + 1];
            faces_[j].setBCType(fz.bcType());
            fi++;
        }

    } else if (fz.isTriangular()) {
        integer fi = 0;
        for (integer j = fz.firstIndex(); j <= fz.lastIndex(); ++j) {
            faces_[j].resize(3);
            for (integer i = 0; i < faces_[j].size(); ++i) {
                faces_[j][i] = faceConnect[fi][i];
            }
            faces_[j].leftCell() = faceConnect[fi][faces_[j].size()];
            faces_[j].rightCell() = faceConnect[fi][faces_[j].size() + 1];
            faces_[j].setBCType(fz.bcType());
            fi++;
        }
    } else if (fz.isQuadrilateral()) {
        integer fi = 0;
        for (integer j = fz.firstIndex(); j <= fz.lastIndex(); ++j) {
            faces_[j].resize(4);
            for (integer i = 0; i < faces_[j].size(); ++i) {
                faces_[j][i] = faceConnect[fi][i];
            }
            faces_[j].leftCell() = faceConnect[fi][faces_[j].size()];
            faces_[j].rightCell() = faceConnect[fi][faces_[j].size() + 1];
            faces_[j].setBCType(fz.bcType());
            fi++;
        }
    }
}

void OpenHurricane::HurricaneMeshRead::readPointZonesHurricane(const hdf5I &fos,
                                                               const string &gridGroupName,
                                                               const integer nPZ) {
    pointZones_.resize(nPZ);
    integerList pointZId;
    fos.read(pointZId, gridGroupName, "pointZoneId");
    for (integer i = 0; i < nPZ; ++i) {
        string dataName = "pointZone";
        dataName += toString(pointZId[i]);

        integer id;
        fos.readIntegerAttributeFromDataset(id, "index", gridGroupName, dataName);
        pointZones_[i].setIndex(id);

        integer fi;
        fos.readIntegerAttributeFromDataset(fi, "firstIndex", gridGroupName, dataName);
        pointZones_[i].setFirstIndex(fi);

        integer li;
        fos.readIntegerAttributeFromDataset(li, "lastIndex", gridGroupName, dataName);
        pointZones_[i].setLastIndex(li);

        integer t;
        fos.readIntegerAttributeFromDataset(t, "type", gridGroupName, dataName);
        pointZones_[i].setPointType(t);

        integer nd;
        fos.readIntegerAttributeFromDataset(nd, "nDimension", gridGroupName, dataName);
        pointZones_[i].setND(nd);

        std::string n;
        fos.readStringAttributeFromDataset(n, "name", gridGroupName, dataName);
        pointZones_[i].resetName(n);

        vectorArray tp;
        fos.read(tp, gridGroupName, dataName);

        readPointsHurricane(fos, gridGroupName, pointZones_[i], tp);
    }
}

void OpenHurricane::HurricaneMeshRead::setPointsHurricane(const hdf5I &fos,
                                                          const string &gridGroupName,
                                                          const integer nP) {
    points_.resize(nP);
}

void OpenHurricane::HurricaneMeshRead::readPointsHurricane(const hdf5I &fos,
                                                           const string &gridGroupName,
                                                           const pointZone &pz,
                                                           const vectorArray &tp) {
    integer pi = 0;
    for (integer j = pz.firstIndex(); j <= pz.lastIndex(); ++j) {
        points_[j] = tp[pi];
        pi++;
    }
}

void OpenHurricane::HurricaneMeshRead::setPeriodicPairListSizeHurricane(const hdf5I &fos,
                                                                        const string &gridGroupName,
                                                                        const integer nPP) {
    periodicPairZone_.resize(nPP);
}

void OpenHurricane::HurricaneMeshRead::readPeriodicPairListHurricane(const hdf5I &fos,
                                                                     const string &gridGroupName) {
    integerList pplpZ;
    fos.read(pplpZ, gridGroupName, "periodicPairList");
    integer offset = 0;
    for (integer i = 0; i < pplpZ.size(); ++i) {
        string dataName;
        dataName = "periodicPair";
        dataName += toString(pplpZ[i]);

        integer pzid;
        fos.readIntegerAttributeFromDataset(pzid, "periodicZone", gridGroupName, dataName);
        periodicPairZone_[i].setPeriodicZone(pzid);

        integer szid;
        fos.readIntegerAttributeFromDataset(szid, "shadowZone", gridGroupName, dataName);
        periodicPairZone_[i].setShadowZone(szid);

        integer t;
        fos.readIntegerAttributeFromDataset(t, "type", gridGroupName, dataName);
        periodicPairZone_[i].setType(t);

        fos.read(periodicPairZone_[i], gridGroupName, dataName);
        periodicPairZone_[i].setFirstIndex(1 + offset);
        periodicPairZone_[i].setLastIndex(periodicPairZone_[i].size() + offset);
    }
}

OpenHurricane::integer OpenHurricane::HurricaneMeshRead::readGridIntAttrHurricane(
    const hdf5I &fos, const string &gridGroupName, const string &attrName) {
    integer ci;
    fos.readIntegerAttributeFromGroup(ci, attrName, gridGroupName);
    return ci;
}

void OpenHurricane::HurricaneMeshRead::readHurricane(const hdf5I &fos, const string &gridGroupName,
                                                     string &gridUnit) {
    if (HurMPI::master()) {
        if (hasBeenRead_) {
            LFatal("Attempt to read the origin mesh file again!");
        }

        integer nd;
        fos.readIntegerAttributeFromGroup(nd, "nDimension", gridGroupName);
        //string gridUnit;
        fos.readStringAttributeFromGroup(gridUnit, "unit", gridGroupName);

        const auto np = readGridIntAttrHurricane(fos, gridGroupName, "nPoints");
        globalNodes_ = np;
        const auto npz = readGridIntAttrHurricane(fos, gridGroupName, "nPointZones");

        const auto nc = readGridIntAttrHurricane(fos, gridGroupName, "nCells");
        globalCells_ = nc;
        const auto ncz = readGridIntAttrHurricane(fos, gridGroupName, "nCellZones");

        const auto nf = readGridIntAttrHurricane(fos, gridGroupName, "nFaces");
        globalFaces_ = nf;
        const auto nfz = readGridIntAttrHurricane(fos, gridGroupName, "nFaceZones");

        const auto npp = readGridIntAttrHurricane(fos, gridGroupName, "periodicPairSize");

        setCellsHurricane(fos, gridGroupName, nc);
        readCellZonesHurricane(fos, gridGroupName, ncz);
        readCellLoadWeightHurricane(fos, gridGroupName, nc, ncz);

        setFacesHurricane(fos, gridGroupName, nf);
        readFaceZonesHurricane(fos, gridGroupName, nfz);

        setPointsHurricane(fos, gridGroupName, np);
        readPointZonesHurricane(fos, gridGroupName, npz);

        readOriginAndAxis(fos, gridGroupName);

        setPeriodicPairListSizeHurricane(fos, gridGroupName, npp);
        if (npp != 0) {
            readPeriodicPairListHurricane(fos, gridGroupName);
        }

        originMeshRead::formingFaces(cells_.size(), faces_, cells_);
        originMeshRead::formingNodes(cells_.size(), faces_, cells_);
        hasBeenRead_ = true;

        printMeshInfo();

        // Read the number of the processoers used in the mesh file.
        int nProcessor;
        fos.readIntegerAttributeFromFile(nProcessor, "nProcessor");
        if (nProcessor != 1) {
            originMeshDecompSize_ = nProcessor;
            decomposeList_.resize(cellZones_.size());
            string decomposingGroup = gridGroupName + "Decompose";
            for (integer i = 0; i < cellZones_.size(); ++i) {
                string dataName = "cellZone";
                dataName += toString(cellZones_[i].index());
                fos.read(decomposeList_[i], decomposingGroup, dataName);
            }
        }
    }

    HurMPI::bcast(&origin_[0], 3, feature<real>::MPIType, HurMPI::masterNo());
    HurMPI::bcast(&axis_[0], 3, feature<real>::MPIType, HurMPI::masterNo());
}

OpenHurricane::HurricaneMeshRead::HurricaneMeshRead() : originMeshRead() {}

OpenHurricane::HurricaneMeshRead::HurricaneMeshRead(const fileName &fN, const int nP)
    : originMeshRead(fN, nP) {}

OpenHurricane::HurricaneMeshRead::HurricaneMeshRead(const std::string &str, const int nP)
    : originMeshRead(str, nP) {}

void OpenHurricane::HurricaneMeshRead::reading(string &gridUnit) {
    hdf5I myh5(fileName_);
    if (!myh5.isHDF5File()) {
        LFatal("Cannot read mesh: %s. Because it is not a HURRICANE mesh. Please check!",
               fileName_.c_str());
    }
    myh5.open();
    readHurricane(myh5, "grid", gridUnit);
    myh5.close();
}

void OpenHurricane::HurricaneMeshRead::readCellLoadWeightHurricane(const hdf5I &fos,
                                                                   const string &gridGroupName,
                                                                   const integer nC,
                                                                   const integer nCZ) {
    cellLoadWeights_.resize(nC, Zero);

    integerList cellZId;
    fos.read(cellZId, gridGroupName, "cellZoneId");
    bool found = true;
    for (integer i = 0; i < nCZ; ++i) {
        integerArray cellLoadWgt(cellZones_[i].size());
        string dataNameOri = "cellLoadWeights";
        dataNameOri += toString(cellZId[i]);
        if (!fos.exist(gridGroupName, dataNameOri)) {
            found = false;
            break;
        }
        fos.read(cellLoadWgt, gridGroupName, dataNameOri);
        for (integer j = cellZones_[i].firstIndex(); j <= cellZones_[i].lastIndex(); ++j) {
            integer nn = j - cellZones_[i].firstIndex();
            cellLoadWeights_[j] = cellLoadWgt[nn];
        }
    }

    if (!found) {
        cellLoadWeights_.clear();
    }
}
