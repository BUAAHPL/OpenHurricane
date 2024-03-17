/*!
 * \file CGNSMeshRead.cpp
 * \brief Main subroutines for reading origin mesh.
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

#include "CGNSMeshRead.hpp"

#ifdef USES_CGNS
#include "cgnsIO.hpp"
#include "formFaceFromCell.hpp"

namespace OpenHurricane {
    createClassNameStr(CGNSMeshRead, "CGNS");
    registerObjFty(originMeshRead, CGNSMeshRead, controller);
} // namespace OpenHurricane

hur_nodiscard int OpenHurricane::CGNSMeshRead::readBase(const cgnsIO &cff) const {
    return cff.readNBases();
    ;
}

void OpenHurricane::CGNSMeshRead::checkNumBase(const int nB) const {
#ifdef HUR_DEBUG
    Pout << "    Info: the number of base is " << nB << std::endl;
#endif // HUR_DEBUG
    if (nB != 1) {
        LFatal("The number of base must be 1 in %s at current program", fileName_.c_str());
    }
}

void OpenHurricane::CGNSMeshRead::checkNumZone(const int nZone) const {
#ifdef HUR_DEBUG
    Pout << "    Info: the number of base is " << nZone << std::endl;
#endif // HUR_DEBUG
    if (nZone != 1) {
        LFatal("The number of zone must be 1 in %s at current program", fileName_.c_str());
    }
}

hur_nodiscard int OpenHurricane::CGNSMeshRead::readNumZone(const cgnsIO &cff,
                                                           const string &baseName,
                                                           const int iBase) const {
#ifdef HUR_DEBUG
    const auto nz = cff.readNZones(iBase + 1);
    Pout << "    Info: the number of zones is " << nz << " in base \"" << baseName << "\""
         << std::endl;
    return nz;
#else  // HUR_DEBUG
    return cff.readNZones(iBase + 1);
#endif // HUR_DEBUG
}

void OpenHurricane::CGNSMeshRead::setZoneSize(const cgnsIO &cff, const int nBase) {
    checkNumBase(nBase);

    integer tCells = 0;
    integer tPoints = 0;
    integer tCellZone = 0;
    integer tPointZone = 0;
    integer tFaceZone = 0;
    for (int iBase = 0; iBase < nBase; ++iBase) {
        string baseName;
        integer cellDim = 0;
        integer physcDim = 0;
        cff.readBase(iBase + 1, baseName, cellDim, physcDim);
#ifdef HUR_DEBUG
        Pout << "    Info: reading base: " << iBase + 1 << " of which the name is \"" << baseName
             << "\", the cell dimension is " << cellDim << ", and the physical dimension is "
             << physcDim << std::endl;
#endif // HUR_DEBUG

        const auto nZone = readNumZone(cff, baseName, iBase);
        checkNumZone(nZone);

        for (int iZone = 0; iZone < nZone; ++iZone) {
            tPointZone++;

            integer nPoints = 0, nCells = 0, nBoundPoints = 0;
            string zoneName;
            cff.readZone(iBase + 1, iZone + 1, zoneName, nPoints, nCells, nBoundPoints);
            tCells += nCells;
            tPoints += (nPoints + nBoundPoints);
#ifdef HUR_DEBUG
            Pout << "    Info: the zone: " << iZone + 1 << " of which the name is \"" << zoneName
                 << "\", and has " << nPoints << " vertices, " << nCells << " cells and "
                 << nBoundPoints << " boundary vertices" << std::endl;
#endif // HUR_DEBUG
            const auto nSection = cff.readNSections(iBase + 1, iZone + 1);

            tFaceZone += nSection;
#ifdef HUR_DEBUG
            const auto nBoundCondition = cff.readNBndCon(iBase + 1, iZone + 1);
            Pout << "    Info: " << nSection << " sections and " << nBoundCondition
                 << " boundary conditions in " << zoneName << std::endl;
#endif // HUR_DEBUG
        }
    }

    cells_.resize(tCells);

    pointZones_.resize(tPointZone);
    //cellZones_.resize(tCellZone);
    faceZones_.resize(tFaceZone);

    tPointZone = 0;
    for (int iBase = 0; iBase < nBase; ++iBase) {
        string baseName;
        integer cellDim = 0;
        integer physcDim = 0;
        cff.readBase(iBase + 1, baseName, cellDim, physcDim);

        const auto nZone = readNumZone(cff, baseName, iBase);
        for (int iZone = 0; iZone < nZone; ++iZone) {
            integer nPoints = 0, nCells = 0, nBoundPoints = 0;
            string zoneName;
            cff.readZone(iBase + 1, iZone + 1, zoneName, nPoints, nCells, nBoundPoints);
            pointZones_[tPointZone].resetName(zoneName);
            pointZones_[tPointZone].setIndex(iZone);
        }
    }
}

void OpenHurricane::CGNSMeshRead::setPoints(vectorArray &coord) {
    points_.transfer(coord);
    globalNodes_ = points_.size();
    pointZones_[0].setFirstIndex(0);
    pointZones_[0].setLastIndex(globalNodes_ - 1);
}

void OpenHurricane::CGNSMeshRead::setCells(integerListList &cellConn, const integer start,
                                           const integer end) {
    integer count = 0;
    for (integer icell = start; icell <= end; ++icell) {
        cells_[icell].nodesList().transfer(cellConn[count++]);
    }
}

void OpenHurricane::CGNSMeshRead::setNGlobalFace() {
    globalFaces_ = faces_.size();
    integer nintf = 0;
    for (const auto &e : faceZones_) {
        if (e.isInterior()) {
            nintf += e.size();
        }
    }
    globalInterior_ = nintf;

#ifdef HUR_DEBUG
    integer nf = 0;
    for (const auto &e : faceZones_) {
        nf += e.size();
    }
    if (nf != globalFaces_) {
        LFatal(
            "The number of global face is not equal to the summation of faces of all face zones");
    }
#endif // HUR_DEBUG

    for (const auto &e : faceZones_) {
        for (integer fi = e.firstIndex(); fi <= e.lastIndex(); ++fi) {
            faces_[fi].setBCType(e.bcType());
        }
    }
}

void OpenHurricane::CGNSMeshRead::setCellOriginalIndex() {
    for (integer ci = 0; ci < cells_.size(); ++ci) {
        cells_[ci].setOrignIndex(ci);
    }

    originCellIndex_.resize(cellZones_.size());
    for (integer czi = 0; czi < cellZones_.size(); ++czi) {
        originCellIndex_[czi].resize(cellZones_[czi].size());

        integer countCi = 0;
        for (integer ci = cellZones_[czi].firstIndex(); ci <= cellZones_[czi].lastIndex(); ++ci) {
            originCellIndex_[czi][countCi] = ci;
            countCi++;
        }
    }
}

void OpenHurricane::CGNSMeshRead::formingFaces(const integerListListList &faceZoneEleConn,
                                               const integer faceTableCapacity) {
    formFaceFromCell ffc(cells_, faces_, cellZones_, faceZones_, faceZoneEleConn,
                         2 * faceTableCapacity);
    ffc();
}

void OpenHurricane::CGNSMeshRead::getFaceZoneEleConn(const cgnsIO &cff,
                                                     integerListListList &faceZoneEleConn) {
    const auto nBase = readBase(cff);
    for (int iBase = 0; iBase < nBase; ++iBase) {
        string baseName;
        integer cellDim = 0;
        integer physcDim = 0;
        cff.readBase(iBase + 1, baseName, cellDim, physcDim);
        const auto nZone = readNumZone(cff, baseName, iBase);

        for (int iZone = 0; iZone < nZone; ++iZone) {
            const auto nSection = cff.readNSections(iBase + 1, iZone + 1);
            const auto nBoundCondition = cff.readNBndCon(iBase + 1, iZone + 1);

            for (integer is = 0; is < nSection; ++is) {
                string sectn;
                ElementType_t tpptr;
                integer start, end;
                int nbdnry;
                int parentFlag;
                cff.readSection(iBase + 1, iZone + 1, is + 1, sectn, &tpptr, start, end, &nbdnry,
                                &parentFlag);

                if (start <= cells_.size() && end <= cells_.size()) {
                    // Do nothing
                } else {
                    const auto eleSize = end - start + 1;
                    const auto eleDataSize = cff.readElementDataSize(iBase + 1, iZone + 1, is + 1);
                    cgsize_t *ele = new cgsize_t[eleDataSize];
                    cgsize_t *pare = nullptr;
                    cff.readElement(iBase + 1, iZone + 1, is + 1, ele, pare);
                    faceZoneEleConn[is - cellZones_.size()] =
                        cff.parsingElementConn(ele, eleSize, eleDataSize, tpptr);
                    HurDeleteDynArray(ele);
                }
            }
        }
    }
}

void OpenHurricane::CGNSMeshRead::getGridConnectivity(const cgnsIO &cff, const integer nBase,
                                                      const integer nZone) {
    for (int iBase = 0; iBase < nBase; ++iBase) {
        for (int iZone = 0; iZone < nZone; ++iZone) {
            const auto nconns = cff.readNConns(iBase + 1, iZone + 1);
            if (nconns != 0) {
                string connName;
                for (integer iC = 0; iC < nconns; ++iC) {
                    cff.readConnInfo(iBase + 1, iZone + 1, iC + 1, connName);
                }
            }
        }
    }
}

OpenHurricane::CGNSMeshRead::CGNSMeshRead() : originMeshRead() {}

OpenHurricane::CGNSMeshRead::CGNSMeshRead(const fileName &fN, const int nP)
    : originMeshRead(fN, nP) {}

OpenHurricane::CGNSMeshRead::CGNSMeshRead(const std::string &str, const int nP)
    : originMeshRead(str, nP) {}

void OpenHurricane::CGNSMeshRead::reading(string &gridUnit) {
    if (HurMPI::master()) {
        if (hasBeenRead_) {
            LFatal("Attempt to read the origin mesh file again!");
        }
        cgnsIO cff(fileName_, cgnsIO::READ);
        cff.open();

#ifdef HUR_DEBUG
        Pout << "    Info: the CGNS file version is " << cff.cgVersion()
             << ". And the precision is " << cff.cgPrecision() << std::endl;
#endif // HUR_DEBUG

        integer faceTableCapacity = 0;

        const auto nBase = readBase(cff);
        setZoneSize(cff, nBase);
        integer nZone = 0;
        for (int iBase = 0; iBase < nBase; ++iBase) {
            string baseName;
            integer cellDim = 0;
            integer physcDim = 0;
            cff.readBase(iBase + 1, baseName, cellDim, physcDim);
            nZone = readNumZone(cff, baseName, iBase);

            for (int iZone = 0; iZone < nZone; ++iZone) {
                integer nPoints = 0, nCells = 0, nBoundPoints = 0;
                string zoneName;
                cff.readZone(iBase + 1, iZone + 1, zoneName, nPoints, nCells, nBoundPoints);
                globalCells_ = nCells;
                faceTableCapacity += nCells;
                vectorArray coord;
                cff.readCoord(iBase + 1, iZone + 1, nPoints + nBoundPoints, coord);
                setPoints(coord);
                const auto nSection = cff.readNSections(iBase + 1, iZone + 1);
                const auto nBoundCondition = cff.readNBndCon(iBase + 1, iZone + 1);
                integer cellZS = 0;
                for (integer is = 0; is < nSection; ++is) {
                    string sectn;
                    ElementType_t tpptr;
                    integer start, end;
                    int nbdnry;
                    int parentFlag;
                    cff.readSection(iBase + 1, iZone + 1, is + 1, sectn, &tpptr, start, end,
                                    &nbdnry, &parentFlag);
                    const auto eleSize = end - start + 1;
                    const auto eleDataSize = cff.readElementDataSize(iBase + 1, iZone + 1, is + 1);
#ifdef FULL_DEBUG
                    Pout << sectn << " " << ElementTypeName[tpptr] << " " << eleDataSize << " "
                         << parentFlag << std::endl;
#endif // FULL_DEBUG

                    if (start <= nCells && end <= nCells) {
                        cgsize_t *ele = new cgsize_t[eleDataSize];
                        cgsize_t *pare = nullptr;
                        cff.readElement(iBase + 1, iZone + 1, is + 1, ele, pare);
                        cellZones_.append(cellZone(sectn, is, start - 1, end - 1,
                                                   cellZoneTypes::FLUID,
                                                   cff.checkElementType(tpptr)));
                        cff.parsingCellElementConn(ele, start - 1, end - 1, eleDataSize, tpptr,
                                                   cells_);
                        HurDeleteDynArray(ele);
                        faceZones_[is].resetName("int_" + sectn);
                        faceZones_[is].setIndex(is + 1);
                    } else {
                        faceZones_[is].resetName(sectn);
                        faceZones_[is].setIndex(is + 1);
                        faceZones_[is].setFaceType(cff.checkFaceElementType(tpptr));
                        faceTableCapacity += eleSize;
                    }
                }
                for (integer ib = 0; ib < nBoundCondition; ++ib) {
                    cff.readBndConInfo(iBase + 1, iZone + 1, ib + 1, faceZones_);
                }
            }
        }

        integerListListList faceZoneEleConn(faceZones_.size() - cellZones_.size());
        getFaceZoneEleConn(cff, faceZoneEleConn);
        formingFaces(faceZoneEleConn, faceTableCapacity);

        getGridConnectivity(cff, nBase, nZone);

        cff.close();

        hasBeenRead_ = true;
        setNGlobalFace();
        setCellOriginalIndex();

        printMeshInfo();
    }
}
#endif // USES_CGNS