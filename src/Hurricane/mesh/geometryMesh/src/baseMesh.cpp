/*!
 * \file baseMesh.cpp
 * \brief The subroutines and functions of base mesh
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

#include "baseMesh.hpp"
#include "commonInclude.hpp"

OpenHurricane::real OpenHurricane::baseMeshWarningThreshold::OrthogonalityThreshold =
    std::cos(OpenHurricane::degToRad(60.0));

namespace OpenHurricane {
    const std::string baseMesh::className_ = "baseMesh";

    void OpenHurricane::baseMesh::createCellNeighbourCellTypeMap() {
        cellNeighbourCellTypeMap_.emplace("FACE", cellNeighbourCellType::faceNeighbour);
        cellNeighbourCellTypeMap_.emplace("FACENEIGHBOUR", cellNeighbourCellType::faceNeighbour);
        cellNeighbourCellTypeMap_.emplace("TWOFACE", cellNeighbourCellType::twoFaceNeighbour);
        cellNeighbourCellTypeMap_.emplace("TWOFACENEIGHBOUR",
                                          cellNeighbourCellType::twoFaceNeighbour);
        cellNeighbourCellTypeMap_.emplace("SECONDFACE", cellNeighbourCellType::twoFaceNeighbour);
        cellNeighbourCellTypeMap_.emplace("SECONDFACENEIGHBOUR",
                                          cellNeighbourCellType::twoFaceNeighbour);
        cellNeighbourCellTypeMap_.emplace("NODE", cellNeighbourCellType::nodeNeighbour);
        cellNeighbourCellTypeMap_.emplace("NODENEIGHBOUR", cellNeighbourCellType::nodeNeighbour);
    }

    baseMesh::baseMesh()
        : baseMeshCore(), cont_(controller::null),
          ghostCellType_(Options::ghostCellLayers::NO_GHOST_CELL), isGhostCellCreated_(false),
          twoLayer_(false), cellNeighbourCell_(cellNeighbourCellType::faceNeighbour),
          cellNeighbourCellTypeMap_() {}

    OpenHurricane::baseMesh::baseMesh(const controller &cont, const integer nPoint,
                                      const integer nFace, const integer nInteriorFace,
                                      const integer nCell, const bool twoLayer)
        : baseMeshCore(), cont_(cont), ghostCellType_(Options::ghostCellLayers::NO_GHOST_CELL),
          isGhostCellCreated_(false), twoLayer_(twoLayer),
          cellNeighbourCell_(cellNeighbourCellType::faceNeighbour), cellNeighbourCellTypeMap_() {

        setGhostCellsType(1);

        if (ghostCellType_ > 1) {
            setSecondNeighbourCellsSize(short(2) * (ghostCellType_ - 1));
        }

        nTotalCells_ = nCells_ + nGhostCells_;

        string cellNeighbourCellW = "face";
        if (cont_.found("spatialScheme")) {
            const auto &spsCont = cont_.subController("spatialScheme");
            if (spsCont.found("cellNeighbourCell")) {
                cellNeighbourCellW = spsCont.findWord("cellNeighbourCell");
            }

            if (spsCont.found("gradient")) {
                string grdw = spsCont.findWord("gradient");
                if (grdw == "leastSquareGrad") {
                    if (spsCont.found("weightType")) {
                        string wtw = spsCont.findWord("weightType");

                        // For WLSQG type cell-neighbour-cell must be face
                        if (wtw == "WLSQG") {
                            cellNeighbourCellW = "face";
                        }
                    }
                }
            }
        }

        createCellNeighbourCellTypeMap();

        stringToUpperCase(cellNeighbourCellW);

        auto iter = cellNeighbourCellTypeMap_.find(cellNeighbourCellW);
        if (iter != cellNeighbourCellTypeMap_.end()) {
            cellNeighbourCell_ = iter->second;
        }
        cellNeighbourCellTypeMap_.clear();
    }

    baseMesh::baseMesh(const controller &cont)
        : baseMeshCore(), cont_(cont), ghostCellType_(Options::ghostCellLayers::NO_GHOST_CELL),
          isGhostCellCreated_(false), twoLayer_(false),
          cellNeighbourCell_(cellNeighbourCellType::faceNeighbour), cellNeighbourCellTypeMap_() {

        setGhostCellsType(1);

        if (ghostCellType_ > 1) {
            setSecondNeighbourCellsSize(short(2) * (ghostCellType_ - 1));
        }

        string cellNeighbourCellW = "face";
        if (cont_.found("spatialScheme")) {
            const auto &spsCont = cont_.subController("spatialScheme");
            if (spsCont.found("cellNeighbourCell")) {
                cellNeighbourCellW = spsCont.findWord("cellNeighbourCell");
            }

            if (spsCont.found("gradient")) {
                string grdw = spsCont.findWord("gradient");
                if (grdw == "leastSquareGrad") {
                    if (spsCont.found("weightType")) {
                        string wtw = spsCont.findWord("weightType");

                        // For WLSQG type cell-neighbour-cell must be face
                        if (wtw == "WLSQG") {
                            cellNeighbourCellW = "face";
                        }
                    }
                }
            }
        }

        createCellNeighbourCellTypeMap();
        stringToUpperCase(cellNeighbourCellW);

        auto iter = cellNeighbourCellTypeMap_.find(cellNeighbourCellW);
        if (iter != cellNeighbourCellTypeMap_.end()) {
            cellNeighbourCell_ = iter->second;
        }

        cellNeighbourCellTypeMap_.clear();
    }

    baseMesh::~baseMesh() noexcept {
        clearMesh();
    }

    void OpenHurricane::baseMesh::setBaseMesh(const integer nPoint, const integer nFace,
                                              const integer nInteriorFace, const integer nCell,
                                              const bool twoLayer) {
        nPoints_ = nPoint;
        nFaces_ = nFace;
        nInteriorFaces_ = nInteriorFace;
        nCells_ = nCell;
        nBoundFaces_ = nFaces_ - nInteriorFaces_;

        nGhostCells_ = ghostCellType_ * nBoundFaces_;

        nTotalCells_ = nCells_ + nGhostCells_;
    }

    void OpenHurricane::baseMesh::clearMesh() noexcept {
        baseMeshCore::clearMesh();
    }

} // namespace OpenHurricane

void OpenHurricane::baseMesh::calcCellAspectRatio() const {
    if (!aspectRatioPtr_) {
        aspectRatioPtr_.reset(new realArray(nCells()));
        realArray &cAR = *aspectRatioPtr_;

        const cellList &cl = cells();

        const vectorArray &fCtrs = faceCentre();
        const vectorArray &fAreas = faceArea();

        const vectorArray &cCtrs = cellCentre();
        const realArray &cVols = cellVolume();

        for (integer ci = 0; ci < nCells(); ++ci) {
            real smax = -large;
            real dmax = -large;

            for (integer fLI = 0; fLI < cl[ci].facesList().size(); ++fLI) {
                integer fI = cl[ci].facei(fLI);

                smax = max(smax, mag(fAreas[fI]));
                dmax = max(dmax, mag(fCtrs[fI] - cCtrs[ci]));
            }
            cAR[ci] = 2.0 * dmax * smax / cVols[ci];
        }
    }
}

void OpenHurricane::baseMesh::calcCellDeltaMax() const {
    if (!deltaMaxPtr_) {
        deltaMaxPtr_.reset(new realArray(nCells()));
        realArray &dM = *deltaMaxPtr_;

        const cellList &cl = cells();
        const pointField &ps = points();

        for (integer i = 0; i < nCells_; ++i) {
            integer j = cl[i].nodei(0);
            real xmin = ps[j].x();
            real ymin = ps[j].y();
            real zmin = ps[j].z();
            real xmax = ps[j].x();
            real ymax = ps[j].y();
            real zmax = ps[j].z();

            for (integer k = 1; k < cl[i].nodeSize(); ++k) {
                j = cl[i].nodei(k);
                xmin = min(xmin, ps[j].x());
                ymin = min(ymin, ps[j].y());
                zmin = min(zmin, ps[j].z());
                xmax = max(xmax, ps[j].x());
                ymax = max(ymax, ps[j].y());
                zmax = max(zmax, ps[j].z());
            }

            const real dx = xmax - xmin;
            const real dy = ymax - ymin;
            const real dz = zmax - zmin;
            dM[i] = max(max(dx, dy), dz);
        }
    }
}

void OpenHurricane::baseMesh::calcCellCentreAndVol() const {
    if (!cellCentrePtr_ && !cellVolumePtr_) {
        cellCentrePtr_.reset(new vectorArray(nTotalCells_));
        vectorArray &cCtrs = *cellCentrePtr_;

        cellVolumePtr_.reset(new realArray(nTotalCells_));
        realArray &cVols = *cellVolumePtr_;

        makeCellCentreAndVol(points(), faceCentre(), faceArea(), cCtrs, cVols);

        if (!isGhostCellCreated_) {
            LFatal("Not creating ghost cells");
        }
        makeGhostCellCentreAndVol(faceCentre(), faceArea(), cCtrs, cVols);
    }
}

void OpenHurricane::baseMesh::makeCellCentreAndVol(const pointField &p, const vectorArray &fCtrs,
                                                   const vectorArray &fAreas, vectorArray &cellCtrs,
                                                   realArray &cellVols) const {
    const cellList &cs = cells();
    const faceList &fs = faces();
    real volmin = large;
    for (integer celli = 0; celli < nCells_; celli++) {
        const cell &c = cs[celli];
        integer nCPoints = c.nodeSize();

        vector cCtrs = Zero;

        cCtrs = p[c.nodei(0)];
        for (integer pi = 1; pi < nCPoints; pi++) {
            cCtrs += p[c.nodei(pi)];
        }

        cCtrs /= real(nCPoints);

        cellCtrs[celli] = cCtrs;

        //Cell volume Gauss intergration based on Eq.s(4.14)~(4.15), Eq.(5.13) in Ref.[1]:
        //[1] Blazek J., Computational fluid Dynamics - Principles and Applications, 2001
        real cVol = real(0.0);
        for (integer fi = 0; fi < c.faceSize(); fi++) {
            real cVolf = real(0.0);
            cVolf = fAreas[c.facei(fi)] * (fCtrs[c.facei(fi)] - p[c.nodei(0)]);
            integer cl = fs[c.facei(fi)].leftCell();
            if (celli == cl) {
                cVolf = -cVolf;
            }
            cVol += cVolf;
        }
        cVol *= real(1.0 / 3.0);

        cellVols[celli] = cVol;
        volmin = min(cVol, volmin);
    }
    const_cast<real &>(minSize_) = pow(volmin, real(1.0 / 3.0));
}

void OpenHurricane::baseMesh::makeGhostCellCentreAndVol(const vectorArray &faceCtrs,
                                                        const vectorArray &faceAreas,
                                                        vectorArray &cellCtrs,
                                                        realArray &cellVols) const {

    if (!isGhostCellCreated()) {
        LFatal("The ghost cells hadn't been created!");
    }
    const faceZoneList &fz = faceZones();
    const faceList &fs = faces();

    for (integer zoneI = 0; zoneI < fz.size(); zoneI++) {
        if (fz[zoneI].isInterior()) {
            continue;
        } else if (fz[zoneI].isCutFace()) {
            continue;
        } else if (fz[zoneI].isPeriodic() || fz[zoneI].isPeriodicShadow()) {
            continue;
        } else if (fz[zoneI].isSymmetric()) {
            for (integer faceI = fz[zoneI].firstIndex(); faceI <= fz[zoneI].lastIndex(); faceI++) {
                integer cl = fs[faceI].leftCell();
                integer cr = fs[faceI].rightCell();

                const vector faceN = faceAreas[faceI].normalized();
                real rr = (cellCtrs[cl] - faceCtrs[faceI]) * faceN;

                cellCtrs[cr] = cellCtrs[cl] - real(2.0) * rr * faceN;
                cellVols[cr] = cellVols[cl];
            }
        } else {
            for (integer faceI = fz[zoneI].firstIndex(); faceI <= fz[zoneI].lastIndex(); faceI++) {
                integer cl = fs[faceI].leftCell();
                integer cr = fs[faceI].rightCell();

                cellCtrs[cr] = real(2.0) * faceCtrs[faceI] - cellCtrs[cl];
                cellVols[cr] = cellVols[cl];
            }
        }
    }

    const cutZoneList &cutZ = cutZones();
    for (integer cutI = 0; cutI < cutZ.size(); cutI++) {
        cutZ[cutI].transfer(cellVols);
        cutZ[cutI].transferVS(cellCtrs);
    }

    const perZoneList &perZ = perZones();
    for (integer perI = 0; perI < perZ.size(); perI++) {
        perZ[perI].transfer(cellVols);
        perZ[perI].transferPoint(cellCtrs, origin_);
    }
}

void OpenHurricane::baseMesh::calcAdjoinnCellCentre() const {
    if (!adjoinCellCtrPtr_) {
        adjoinCellCtrPtr_.reset(new vectorArray(nFaces_));
        vectorArray &acCtrs = *adjoinCellCtrPtr_;

        const auto &cellCtrs = cellCentre();
        const auto &fs = faces();

        for (integer fi = 0; fi < fs.size(); ++fi) {
            const auto &cl = fs[fi].leftCell();
            const auto &cr = fs[fi].rightCell();

            acCtrs[fi] = cellCtrs[cr] - cellCtrs[cl];
        }
    }
}

void OpenHurricane::baseMesh::calcCellNeighbourCells() const {
#ifdef HUR_DEBUG
    Pout << "baseMesh::calcCellNeighbourCells() : calculating cell neighbour "
            "cells"
         << std::endl;
#endif // HUR_DEBUG

    if (!CNCPtr_) {
        // Create the storage
        CNCPtr_.reset(new integerListList(nCells()));

        // The total number for the neighbour cells for every cell
        integerList NCPC;

        const cellList &cs = cells();
        const faceList &fs = faces();
        for (integer cellI = 0; cellI < nCells(); cellI++) {
            integer nnb = 0;
            NCPC.clear();
            for (integer nf1 = 0; nf1 < cs[cellI].faceSize(); nf1++) {
                integer idxf = cs[cellI].facei(nf1);
                const integer &cl = fs[idxf].leftCell();
                const integer &cr = fs[idxf].rightCell();
                integer nbc = cl + cr - cellI;
                bool hasIncluded = false;

                for (integer m = 0; m < nnb; m++) {
                    if (nbc == NCPC[m]) {
                        hasIncluded = true;
                        break;
                    }
                }
                if (hasIncluded) {
                    continue;
                }

                nnb++;
                NCPC.append(nbc);
            }
            // Second layer face-neighbour cells
            if (cellNeighbourCell_ == cellNeighbourCellType::twoFaceNeighbour) {
                for (integer nf1 = 0; nf1 < cs[cellI].faceSize(); nf1++) {
                    integer idxf = cs[cellI].facei(nf1);
                    const integer &CL = fs[idxf].leftCell();
                    const integer &CR = fs[idxf].rightCell();
                    integer nbc = CL + CR - cellI;
                    // Second layer face-neighbour cells
                    for (integer nf2 = 0; nf2 < cs[nbc].faceSize(); nf2++) {
                        idxf = cs[nbc].facei(nf2);
                        const auto &cl = fs[idxf].leftCell();
                        const auto &cr = fs[idxf].rightCell();
                        integer nbc1 = cl + cr - nbc;

                        if (nbc1 == nbc || nbc1 == cellI) {
                            continue;
                        }
                        //To check if  nb1 cell is already in the neighbour list.
                        bool hasIncluded = false;
                        for (integer m = 0; m < NCPC.size(); m++) {
                            if (nbc1 == NCPC[m]) {
                                hasIncluded = true;
                                break;
                            }
                        }
                        if (hasIncluded) {
                            continue;
                        }

                        //Form the node-neighbour stencil.
                        //To check if they share a node. If NOT, exclude the cell in the neighbour list.
                        if (cellNeighbourCell_ == cellNeighbourCellType::nodeNeighbour) {
                            bool isShared = false;
                            for (integer nn1 = 0; nn1 < cs[cellI].nodeSize(); nn1++) {
                                for (integer nn2 = 0; nn1 < cs[nbc1].nodeSize(); nn2++) {
                                    if (cs[cellI].nodei(nn1) == cs[nbc1].nodei(nn2)) {
                                        isShared = true;
                                        break;
                                    }
                                }
                                if (isShared) {
                                    break;
                                }
                            }
                            if (!isShared) {
                                continue;
                            }
                        }
                        NCPC.append(nbc1);
                    }
                }
            }

            CNCPtr_->operator[](cellI).append(NCPC);
        }

    } else {
        if (report) {
            LWarning("cell neighbour cells has already been calculated");
        }
    }
}

void OpenHurricane::baseMesh::calcCellWeight(const string &wgtType) const {
#ifdef HUR_DEBUG
    Pout << "baseMesh::calcCellWeight() : calculating weights of WLSQ stencil "
            "cells"
         << std::endl;
#endif // HUR_DEBUG
    calcCellNeighbourCells();
    if (!cellWeightPtr_) {
        cellWeightPtr_.reset(new vectorArrayArray(nCells_)); //nTotalCells_ instead of nCells_
        vectorArrayArray &cWgt = *cellWeightPtr_;

        for (integer celli = 0; celli < (*CNCPtr_).size(); celli++) {
            realArrayArray r(4);
            for (integer i = 0; i < 4; i++) {
                r[i].resize(4);
                r[i] = real(0.0);
            }
            realArray weight((*CNCPtr_)[celli].size());
            for (integer i = 0; i < (*CNCPtr_)[celli].size(); i++) {
                cWgt[celli].resize((*CNCPtr_)[celli].size());
                integer nbr = (*CNCPtr_)[celli][i];
                vector disV = (*cellCentrePtr_)[nbr] - (*cellCentrePtr_)[celli];
                real dist = disV.magnitude();
                if (wgtType == "WLSQ0") {
                    weight[i] = real(1.0);
                } else if (wgtType == "WLSQ1") {
                    weight[i] = real(1.0) / dist;
                } else if (wgtType == "WLSQ2") {
                    weight[i] = real(1.0) / (dist * dist);
                } else if (wgtType == "WLSQ3") {
                    weight[i] = real(1.0) / (dist * dist * dist);
                } else if (wgtType == "WLSQG") {
                    integer facei = cells()[celli].facesList()[i];
                    vector normal = (*faceAreaPtr_)[facei].normalized();
                    vector fDistV = (*faceCenterPtr_)[facei] - (*cellCentrePtr_)[celli];

                    real s = (*faceAreaPtr_)[facei].magnitude();
                    real dist0 = mag(normal * disV);
                    real l = mag(fDistV * normal);
                    weight[i] = real(4.0) * l * l * s / (dist0 * dist0 * dist);
                } else {
                    LFatal("Unknown reordering method: %s\nValid reordering method of current "
                           "program are:\nWLSQ0\nWLSQ1\nWLSQ2\nWLSQ3\nWLSQG\n",
                           wgtType.c_str());
                }
                r[1][1] += disV.x() * disV.x() * weight[i];
                r[1][2] += disV.x() * disV.y() * weight[i];
                r[1][3] += disV.x() * disV.z() * weight[i];
                r[2][2] += disV.y() * disV.y() * weight[i];
                r[2][3] += disV.y() * disV.z() * weight[i];
                r[3][3] += disV.z() * disV.z() * weight[i];
            }
            r[1][1] = std::sqrt(r[1][1]);
            r[1][2] = r[1][2] / r[1][1];
            r[1][3] = r[1][3] / r[1][1];
            r[2][2] = std::sqrt(r[2][2] - r[1][2] * r[1][2]);
            r[2][3] = (r[2][3] - r[1][2] * r[1][3]) / r[2][2];
            r[3][3] = std::sqrt(r[3][3] - (r[1][3] * r[1][3] + r[2][3] * r[2][3]));
            real b = (r[1][2] * r[2][3] - r[1][3] * r[2][2]) / (r[1][1] * r[2][2]);
            for (integer i = 0; i < (*CNCPtr_)[celli].size(); i++) {
                integer nbr = (*CNCPtr_)[celli][i];
                vector disV = (*cellCentrePtr_)[nbr] - (*cellCentrePtr_)[celli];
                real r12or11 = r[1][2] / r[1][1];
                real r23or22 = r[2][3] / r[2][2];
                real a1 = disV.x() / (r[1][1] * r[1][1]);
                real a2 = (disV.y() - r12or11 * disV.x()) / (r[2][2] * r[2][2]);
                real a3 = (disV.z() - r23or22 * disV.y() + b * disV.x()) / (r[3][3] * r[3][3]);
                cWgt[celli][i].x() = a1 - r12or11 * a2 + b * a3;
                cWgt[celli][i].y() = a2 - r23or22 * a3;
                cWgt[celli][i].z() = a3;
                cWgt[celli][i] *= weight[i];
            }
        }
    }
}

void OpenHurricane::baseMesh::calcFaceCentreAndArea() const {
    if (!faceAreaPtr_ && !faceCenterPtr_) {
        faceCenterPtr_.reset(new vectorArray(faces().size()));
        vectorArray &fCtrs = *faceCenterPtr_;

        faceAreaPtr_.reset(new vectorArray(faces().size()));
        vectorArray &fAreas = *faceAreaPtr_;

        const auto &p = points();
        const faceList &fs = faces();

        for (integer facei = 0; facei < fs.size(); ++facei) {
            const face &f = fs[facei];
            integer nFPoints = f.size();

            if (nFPoints == 3) {
                fCtrs[facei] = (real(1.0 / 3.0)) * (p[f[0]] + p[f[1]] + p[f[2]]);
                fAreas[facei] = (real(0.5)) * ((p[f[1]] - p[f[0]]) ^ (p[f[2]] - p[f[0]]));
            } else if (nFPoints == 4) {
                fCtrs[facei] = (real(0.25)) * (p[f[0]] + p[f[1]] + p[f[2]] + p[f[3]]);
                fAreas[facei] = (real(0.5) * ((p[f[2]] - p[f[0]]) ^ (p[f[3]] - p[f[1]])));
            } else {
                vector sumN = Zero;
                real sumA = real(0.0);
                vector sumAc = Zero;

                point fCentre = p[f[0]];
                for (integer pi = 1; pi < nFPoints; pi++) {
                    fCentre += p[f[pi]];
                }
                fCentre /= real(nFPoints);

                for (integer pi = 0; pi < nFPoints; pi++) {
                    const point &nextPoint = p[f[(pi + 1) % nFPoints]];

                    vector c = p[f[pi]] + nextPoint + fCentre;
                    vector n = (nextPoint - p[f[pi]]) ^ (fCentre - p[f[pi]]);
                    real a = mag(n);

                    sumN += n;
                    sumA += a;
                    sumAc += a * c;
                }

                if (sumA < rootVeryTiny) {
                    fCtrs[facei] = fCentre;
                    fAreas[facei] = Zero;
                } else {
                    fCtrs[facei] = real(1.0 / 3.0) * sumAc / sumA;
                    fAreas[facei] = real(0.5) * sumN;
                }
            }
        }
    }
}

void OpenHurricane::baseMesh::calcFaceCtrToCellCtr() const {
    if (!faceCtrToLeftCellCtrPtr_ && !faceCtrToRightCellCtrPtr_) {
        faceCtrToLeftCellCtrPtr_.reset(new vectorArray(faces().size()));
        auto &fCtrToLCCtrs = *faceCtrToLeftCellCtrPtr_;
        faceCtrToRightCellCtrPtr_.reset(new vectorArray(faces().size()));
        auto &fCtrToRCCtrs = *faceCtrToRightCellCtrPtr_;

        const auto &fCtrs = faceCentre();
        const auto &cellCtrs = cellCentre();
        const faceList &fs = faces();

        for (integer facei = 0; facei < fs.size(); ++facei) {
            const auto &cl = fs[facei].leftCell();
            const auto &cr = fs[facei].rightCell();
            fCtrToLCCtrs[facei] = fCtrs[facei] - cellCtrs[cl];
            fCtrToRCCtrs[facei] = fCtrs[facei] - cellCtrs[cr];
        }
    }
}

void OpenHurricane::baseMesh::calcFaceWeight() const {
    if (!faceWeightPtr_) {
        faceWeightPtr_.reset(new realArray(nFaces_));
        realArray &fWgt = *faceWeightPtr_;

        const faceList &fs = faces();
        const cellList &cs = cells();

        const vectorArray &fCtrs = faceCentre();
        const vectorArray &cCtrs = cellCentre();

        for (integer faceI = 0; faceI < fs.size(); ++faceI) {
            const integer cl = fs[faceI].leftCell();
            const integer cr = fs[faceI].rightCell();
            if (cr == -1) {
                fWgt[faceI] = real(1.0);
            } else {
                vector fToCl = fCtrs[faceI] - cCtrs[cl];
                vector fToCr = fCtrs[faceI] - cCtrs[cr];
                fWgt[faceI] = mag(fToCr) / (mag(fToCl) + mag(fToCr));
            }
        }
    }
}

void OpenHurricane::baseMesh::calcNodeCellWeight() const {
#ifdef HUR_DEBUG
    Pout << "baseMesh::calcNodeCellWeight() : calculating weights of point "
            "neighbour cells"
         << std::endl;
#endif // HUR_DEBUG
    if (!nodeWeightPtr_) {

        nodeWeightPtr_.reset(new realArrayArray(nPoints_));
        realArrayArray &nCWgt = *nodeWeightPtr_;

        const pointField &p = points();
        const cellList &cs = cells();
        const vectorArray &cCtrs = cellCentre();
        const integerListList &PNC = pointNeighbourCells();

        for (integer pi = 0; pi < p.size(); ++pi) {
            const integerList &PINC = PNC[pi];
            const point &p0 = p[pi];
            nCWgt[pi].resize(PINC.size());
            tensor I = Zero;
            vector Rs = Zero;

            for (integer i = 0; i < PINC.size(); ++i) {
                vector tempV = (cCtrs[PINC[i]] - p0);
                I += inv(tempV.magSqr()) * (tempV & tempV);
                Rs += tempV;
            }
            vector lamb = Rs / I;

            for (integer i = 0; i < PINC.size(); ++i) {
                vector tempV = (cCtrs[PINC[i]] - p0);
                nCWgt[pi][i] = real(1.0) + inv(tempV.magSqr()) * lamb * tempV;
                nCWgt[pi][i] = min(real(2), max(nCWgt[pi][i], real(0)));
            }
        }
    }
}

void OpenHurricane::baseMesh::calcPointNeighbourCells() const {
#ifdef HUR_DEBUG
    Pout << "baseMesh::calcPointNeighbourCells() : calculating point neighbour "
            "cells"
         << std::endl;
#endif // HUR_DEBUG

    if (!isGhostCellCreated()) {
        LFatal("Must create dummy cell first.");
    }

    if (!PNCPtr_) {
        // Create the storage
        PNCPtr_.reset(new integerListList(nPoints()));

        // The total number for the neighbour cells for every point
        integerList NPPC(nPoints(), integer(0));

        const cellList &cs = cells();

        for (integer cellI = 0; cellI < nCells(); cellI++) {
            for (integer pi = 0; pi < cs[cellI].nodeSize(); pi++) {
                const integer nodeI = cs[cellI].nodei(pi);
                NPPC[nodeI]++;
            }
        }

        const auto &fzl = faceZones();
        const auto &fl = faces();

        for (integer i = 0; i < fzl.size(); ++i) {
            if (fzl[i].isBnd()) {
                for (integer fi = fzl[i].firstIndex(); fi <= fzl[i].lastIndex(); ++fi) {
                    for (integer ni = 0; ni < fl[fi].size(); ++ni) {
                        const integer nodeI = fl[fi][ni];
                        NPPC[nodeI]++;
                    }
                }
            }
        }

        integerListList &pointCellsAddr = *PNCPtr_;

        for (integer nodeI = 0; nodeI < nPoints(); nodeI++) {
            pointCellsAddr[nodeI].resize(NPPC[nodeI]);
        }

        NPPC = 0;

        for (integer cellI = 0; cellI < nCells(); cellI++) {
            for (integer pi = 0; pi < cs[cellI].nodeSize(); pi++) {
                const integer nodeI = cs[cellI].nodei(pi);
                pointCellsAddr[nodeI][NPPC[nodeI]++] = cellI;
            }
        }
        for (integer i = 0; i < fzl.size(); ++i) {
            if (fzl[i].isBnd()) {
                for (integer fi = fzl[i].firstIndex(); fi <= fzl[i].lastIndex(); ++fi) {
                    for (integer ni = 0; ni < fl[fi].size(); ++ni) {
                        const integer nodeI = fl[fi][ni];
                        const integer cr = fl[fi].rightCell();
                        pointCellsAddr[nodeI][NPPC[nodeI]++] = cr;
                    }
                }
            }
        }
    } else {
        if (report) {
            LWarning("point neighbour cells has already been calculated");
        }
    }
}

void OpenHurricane::baseMesh::setareAllHexOrQuadri() {
    const cellZoneList &cz = cellZones();
    for (integer i = 0; i < cz.size(); i++) {
        if (cz[i].shapeType() != cellShapeType::shapeTypes::hexahedral &&
            cz[i].shapeType() != cellShapeType::shapeTypes::quadrilateral &&
            cz[i].shapeType() != cellShapeType::shapeTypes::mixed) {
            areAllHexOrQuadri_ = false;
            return;
        }
    }
    for (integer i = 0; i < cz.size(); i++) {
        if (cz[i].shapeType() == cellShapeType::shapeTypes::mixed) {
            const cellList &cs = cells();
            for (integer cellI = cz[i].firstIndex(); cellI < cz[i].lastIndex() + 1; cellI++) {
                if (cs[cellI].shapeType() != cellShapeType::shapeTypes::hexahedral &&
                    cs[cellI].shapeType() != cellShapeType::shapeTypes::quadrilateral) {
                    areAllHexOrQuadri_ = false;
                    return;
                }
            }
        }
    }
    areAllHexOrQuadri_ = true;
    return;
}

void OpenHurricane::baseMesh::modifyFaceZone() {
    faceZoneList &fZL = const_cast<faceZoneList &>(faceZones());
    const controller &interCont = cont_.subController("boundaryCondition");
    const auto &mapEle = interCont.mapEntries();

    for (auto &e : mapEle) {
        if (interCont.isController(e.first)) {
            if (interCont.subController(e.first).findWord("bcType") == "symmetry") {
                for (integer i = 0; i < fZL.size(); i++) {
                    if (fZL[i].name() == e.first) {
                        if (fZL[i].bcType() != faceBCType::bcTypes::SYMMETRY) {
                            const_cast<faceZone &>(fZL[i]).setBcType(faceBCType::bcTypes::SYMMETRY);
                        }
                    }
                }
            }
        }
    }
}
