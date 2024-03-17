/*!
 * \file Venkatakrishnan.cpp
 * \brief Main subroutines for Venkatakrishnan's limiter for piecewise linear reconstruction.
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

#include "Venkatakrishnan.hpp"

#include "meshElements.hpp"
namespace OpenHurricane {
    createClassNameStr(Venkatakrishnan, "Venk");
    registerObjFty(limitersForLinear, Venkatakrishnan, controller);
} // namespace OpenHurricane

void OpenHurricane::Venkatakrishnan::calcLimiters(
    const geometryArray<real, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &grad,
    geometryArray<real, cellMesh> &limiters) const {
    const auto &mesh = cellQ.mesh();

    limiters.setComponent(feature<real>::elementType(1));

    const integer nCells = mesh.nCells();
    const auto &cells = mesh.cells();
    const auto &faces = mesh.faces();

    // Face-midpoint
    //const auto& fCtr = mesh.faceCentre();

    // Cell-centroid
    const auto &cCtr = mesh.cellCentre();

    const auto &cV = mesh.cellVolume();
    const auto &fW = mesh.faceWgt();
    //const auto& cnc = mesh.cellNeighbourCells();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    for (integer cellI = 0; cellI < nCells; ++cellI) {
        // Umax and Umin stand for the minimum/maximum values of all surrounding
        // cells j and including the cell i itself.

        real Umax = cellQ[cellI];
        real Umin = cellQ[cellI];

        for (integer i = 0; i < cells[cellI].faceSize(); ++i) {
            const integer fI = cells[cellI].facei(i);
            const auto &cl = faces[fI].leftCell();
            const auto &cr = faces[fI].rightCell();
            const integer m = (cellI == cl) ? cr : cl;
            //const auto& m = cnc[cellI][i];
            //const integer m = faces[fI].leftCell() + faces[fI].rightCell() - cellI;
            //const integer m = faces[fI].oppositeCell(cellI);

            Umax = max(Umax, cellQ[m]);
            Umin = min(Umin, cellQ[m]);
        }

        for (integer i = 0; i < cells[cellI].faceSize(); ++i) {
            const integer fI = cells[cellI].facei(i);

            // Vector pointing from cell centroid to face midpoint.
            /*const vector dFC = fCtr[fI] - cCtr[cellI];

            real delta2 = grad[cellI] * dFC;*/
            //real delta2 = grad[cellI] * (fCtr[fI] - cCtr[cellI]);
            real delta2;
            if (faces[fI].leftCell() == cellI) {
                delta2 = grad[cellI] * faceLeftCellCtr[fI];
                //const auto& dFC = faceLeftCellCtr[fI];
                //delta2 = grad[cellI].x() * dFC.x() + grad[cellI].y() * dFC.y() + grad[cellI].z() * dFC.z();
            } else {
                delta2 = grad[cellI] * faceRightCellCtr[fI];
                //const auto& dFC = faceLeftCellCtr[fI];
                //delta2 = grad[cellI].x() * dFC.x() + grad[cellI].y() * dFC.y() + grad[cellI].z() * dFC.z();
            }

            const real delta = 1.0 / fW[fI];
            real delta2s =
                componentMultiply(componentSign(delta2), componentAdd(componentMag(delta2), tiny));

            real limiterj;
            //const real eps2 = pow3(K_) * cV[cellI];
            const real eps2 = K3_ * cV[cellI];

            if (delta2 > tiny) {
                real delta1Max = Umax - cellQ[cellI];

                limiterj = 1.0 / delta2s *
                           (((sqr(delta1Max) + eps2) * delta2s + delta * sqr(delta2s) * delta1Max) /
                            (sqr(delta1Max) + delta * sqr(delta2s) + delta1Max * delta2s + eps2));

            } else if (delta2 < -tiny) {
                real delta1Min = Umin - cellQ[cellI];

                limiterj = 1.0 / delta2s *
                           (((sqr(delta1Min) + eps2) * delta2s + delta * sqr(delta2s) * delta1Min) /
                            (sqr(delta1Min) + delta * sqr(delta2s) + delta1Min * delta2s + eps2));
            } else {
                limiterj = 1.0;
            }

            limiters[cellI] = min(limiters[cellI], limiterj);
        }
    }
}

void OpenHurricane::Venkatakrishnan::calcLimiters(
    const geometryArray<vector, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &grad,
    geometryArray<vector, cellMesh> &limiters) const {
    const auto &mesh = cellQ.mesh();

    limiters.setComponent(feature<vector>::elementType(1));

    const integer nCells = mesh.nCells();
    const auto &cells = mesh.cells();
    const auto &faces = mesh.faces();

    // Face-midpoint
    //const auto& fCtr = mesh.faceCentre();

    // Cell-centroid
    //const auto& cCtr = mesh.cellCentre();

    const auto &cV = mesh.cellVolume();
    const auto &fW = mesh.faceWgt();

    //const auto& cnc = mesh.cellNeighbourCells();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    for (integer cellI = 0; cellI < nCells; ++cellI) {
        vector Umax = cellQ[cellI];
        vector Umin = cellQ[cellI];

        for (integer i = 0; i < cells[cellI].faceSize(); ++i) {
            const integer fI = cells[cellI].facei(i);
            const integer m = faces[fI].leftCell() + faces[fI].rightCell() - cellI;
            //const integer m = (cellI == faces[fI].leftCell()) ? faces[fI].rightCell() : faces[fI].leftCell();
            /*const auto& cl = faces[fI].leftCell();
            const auto& cr = faces[fI].rightCell();
            const integer m = (cellI == cl) ? cr : cl;*/
            //const integer& m = cnc[cellI][i];
            Umax = componentMax(Umax, cellQ[m]);
            Umin = componentMin(Umin, cellQ[m]);
        }

        for (integer i = 0; i < cells[cellI].faceSize(); ++i) {
            const integer fI = cells[cellI].facei(i);

            // Vector pointing from cell centroid to face midpoint.
            /*vector dFC = fCtr[fI] - cCtr[cellI];

            vector delta2 = grad[cellI] * dFC;*/

            //vector delta2 = grad[cellI] * (fCtr[fI] - cCtr[cellI]);
            vector delta2;
            if (faces[fI].leftCell() == cellI) {
                delta2 = grad[cellI] * faceLeftCellCtr[fI];
            } else {
                delta2 = grad[cellI] * faceRightCellCtr[fI];
            }
            vector delta2s =
                componentMultiply(componentSign(delta2), componentAdd(componentMag(delta2), tiny));

            const real delta = 1.0 / fW[fI];
            vector limiterj;

            //real eps2 = pow3(K_) * cV[cellI];
            const real eps2 = K3_ * cV[cellI];

            for (int d = 0; d < feature<vector>::nElements_; ++d) {
                if (delta2[d] > tiny) {
                    real delta1Max = Umax[d] - cellQ[cellI][d];
                    /*limiterj[d] = 1.0 / delta2[d] * (
                            ((sqr(delta1Max) + eps2) * delta2[d] + 2.0 * sqr(delta2[d]) * delta1Max)
                            / (sqr(delta1Max) + 2.0 * sqr(delta2[d]) + delta1Max * delta2[d] + eps2)
                            );*/
                    limiterj[d] = 1.0 / delta2s[d] *
                                  (((sqr(delta1Max) + eps2) * delta2s[d] +
                                    delta * sqr(delta2s[d]) * delta1Max) /
                                   (sqr(delta1Max) + delta * sqr(delta2s[d]) +
                                    delta1Max * delta2s[d] + eps2));
                } else if (delta2[d] < -tiny) {
                    real delta1Min = Umin[d] - cellQ[cellI][d];

                    limiterj[d] = 1.0 / delta2s[d] *
                                  (((sqr(delta1Min) + eps2) * delta2s[d] +
                                    delta * sqr(delta2s[d]) * delta1Min) /
                                   (sqr(delta1Min) + delta * sqr(delta2s[d]) +
                                    delta1Min * delta2s[d] + eps2));
                } else {
                    limiterj[d] = 1.0;
                }

                limiters[cellI][d] = min(limiters[cellI][d], limiterj[d]);
            }
        }
    }
}