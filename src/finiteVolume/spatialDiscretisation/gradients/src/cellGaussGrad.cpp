/*!
 * \file cellBasedGreenGauss.cpp
 * \brief Main subroutines for cell-based Green-Gauss gradient approach.
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

#include "cellGaussGrad.hpp"

namespace OpenHurricane {
    createClassNameStr(cellGaussGrad, "cellGaussGrad");
    registerObjFty(gradient, cellGaussGrad, controller);
} // namespace OpenHurricane

void OpenHurricane::cellGaussGrad::gaussGrad(
    const geometryArray<real, cellMesh> &cellQ,
    geometryArray<typename outerProduct<vector, real>::type, cellMesh> &grad) const {

    const auto &mesh = cellQ.mesh();
    const integer nCells = mesh.nCells();
    for (integer i = 0; i < nCells; ++i) {
        grad[i] = Zero;
    }

    const integer nFaces = mesh.nFaces();

    const auto &fl = mesh.faces();
    const auto &fw = mesh.faceWeight();
    const auto &fA = mesh.faceArea();
    for (integer i = 0; i < nFaces; ++i) {
        const auto &cl = fl[i].leftCell();
        const auto &cr = fl[i].rightCell();

        const real wl = fw[i];
        const real wr = real(1.0) - wl;
        const real qf = wl * cellQ[cl] + wr * cellQ[cr];

        grad[cl] -= qf * fA[i];

        if (cr < nCells) {
            grad[cr] += qf * fA[i];
        }
    }

    const auto &vol = mesh.cellVolume();
    for (integer i = 0; i < nCells; ++i) {
        grad[i] /= vol[i];
    }
}

void OpenHurricane::cellGaussGrad::gaussGrad(
    const geometryArray<vector, cellMesh> &cellQ,
    geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &grad) const {
    const auto &mesh = cellQ.mesh();
    const integer nCells = mesh.nCells();
    for (integer i = 0; i < nCells; ++i) {
        grad[i] = Zero;
    }
    const integer nFaces = mesh.nFaces();

    const auto &fl = mesh.faces();
    const auto &fw = mesh.faceWeight();
    const auto &fA = mesh.faceArea();
    for (integer i = 0; i < nFaces; ++i) {
        const auto &cl = fl[i].leftCell();
        const auto &cr = fl[i].rightCell();

        const real wl = fw[i];
        const real wr = real(1.0) - wl;
        const vector qf = wl * cellQ[cl] + wr * cellQ[cr];

        grad[cl] -= qf & fA[i];

        if (cr < nCells) {
            grad[cr] += qf & fA[i];
        }
    }

    const auto &vol = mesh.cellVolume();
    for (integer i = 0; i < nCells; ++i) {
        grad[i] /= vol[i];
    }
}