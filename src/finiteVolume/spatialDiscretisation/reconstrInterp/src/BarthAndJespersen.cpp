/*!
 * \file BarthAndJespersen.cpp
 * \brief Main subroutines for Barth and Jespersen limiter for piecewise linear reconstruction.
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

#include "BarthAndJespersen.hpp"

namespace OpenHurricane {
    createClassNameStr(BarthAndJespersen, "Barth");
    registerObjFty(limitersForLinear, BarthAndJespersen, controller);
} // namespace OpenHurricane

void OpenHurricane::BarthAndJespersen::calcLimiters(
    const geometryArray<real, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &grad,
    geometryArray<real, cellMesh> &limiters) const {
    // Mesh
    const auto &mesh = cellQ.mesh();

    // Initialize to 1.0
    limiters.setComponent(feature<real>::elementType(1));

    const integer nCells = mesh.nCells();
    const auto &cells = mesh.cells();
    const auto &faces = mesh.faces();

    // Face-midpoint
    //const auto& fCtr = mesh.faceCentre();

    // Cell-centroid
    //const auto& cCtr = mesh.cellCentre();
    //const auto& cnc = mesh.cellNeighbourCells();

    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    for (integer cellI = 0; cellI < nCells; ++cellI) {
        real Umax = cellQ[cellI];
        real Umin = cellQ[cellI];

        for (integer i = 0; i < cells[cellI].faceSize(); ++i) {
            const integer fI = cells[cellI].facei(i);
            const integer m = faces[fI].leftCell() + faces[fI].rightCell() - cellI;
            //const auto& m = cnc[cellI][i];
            Umax = componentMax(Umax, cellQ[m]);
            Umin = componentMin(Umin, cellQ[m]);
        }

        for (integer i = 0; i < cells[cellI].faceSize(); ++i) {
            const integer fI = cells[cellI].facei(i);

            // Vector pointing from cell centroid to face midpoint.
            /*vector dFC = fCtr[fI] - cCtr[cellI];

            real delta2 = grad[cellI] * dFC;*/

            //real delta2 = grad[cellI] * (fCtr[fI] - cCtr[cellI]);
            real delta2;
            if (faces[fI].leftCell() == cellI) {
                delta2 = grad[cellI] * faceLeftCellCtr[fI];
            } else {
                delta2 = grad[cellI] * faceRightCellCtr[fI];
            }

            real delta2s =
                componentMultiply(componentSign(delta2), componentAdd(componentMag(delta2), tiny));

            real limiterj;

            if (delta2 > tiny) {
                limiterj = min(real(1.0), (Umax - cellQ[cellI]) / delta2s);
            } else if (delta2 < -tiny) {
                limiterj = min(real(1.0), (Umin - cellQ[cellI]) / delta2s);
            } else {
                limiterj = 1.0;
            }

            limiters[cellI] = min(limiters[cellI], limiterj);
        }
    }
}

void OpenHurricane::BarthAndJespersen::calcLimiters(
    const geometryArray<vector, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &grad,
    geometryArray<vector, cellMesh> &limiters) const {
    // Mesh
    const auto &mesh = cellQ.mesh();

    // Initialize all components to 1.0
    limiters.setComponent(feature<vector>::elementType(1));

    const integer nCells = mesh.nCells();
    const auto &cells = mesh.cells();
    const auto &faces = mesh.faces();

    // Face-midpoint
    //const auto& fCtr = mesh.faceCentre();

    // Cell-centroid
    //const auto& cCtr = mesh.cellCentre();
    //const auto& cnc = mesh.cellNeighbourCells();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    for (integer cellI = 0; cellI < nCells; ++cellI) {
        vector Umax = cellQ[cellI];
        vector Umin = cellQ[cellI];

        for (integer i = 0; i < cells[cellI].faceSize(); ++i) {
            const integer fI = cells[cellI].facei(i);
            const integer m = faces[fI].leftCell() + faces[fI].rightCell() - cellI;
            //const auto& m = cnc[cellI][i];

            Umax = componentMax(Umax, cellQ[m]);
            Umin = componentMin(Umin, cellQ[m]);
        }

        for (integer i = 0; i < cells[cellI].faceSize(); ++i) {
            const integer fI = cells[cellI].facei(i);

            // Vector pointing from cell centroid to face midpoint.
            /*const vector dFC = fCtr[fI] - cCtr[cellI];

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

            vector limiterj;

            for (int d = 0; d < feature<vector>::nElements_; ++d) {

                if (delta2[d] > tiny) {
                    limiterj[d] = min(real(1.0), (Umax[d] - cellQ[cellI][d]) / delta2s[d]);
                } else if (delta2[d] < -tiny) {
                    limiterj[d] = min(real(1.0), (Umin[d] - cellQ[cellI][d]) / delta2s[d]);
                } else {
                    limiterj[d] = 1.0;
                }

                limiters[cellI][d] = min(limiters[cellI][d], limiterj[d]);
            }
        }
    }
}