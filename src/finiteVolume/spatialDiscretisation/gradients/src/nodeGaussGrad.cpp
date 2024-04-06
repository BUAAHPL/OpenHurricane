/*!
 * \file nodeGaussGrad.cpp
 * \brief Main subroutines for node-based Green-Gauss gradient approach.
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

#include "nodeGaussGrad.hpp"

namespace OpenHurricane {
    createClassNameStr(nodeGaussGrad, "nodeGaussGrad");
    registerObjFty(gradient, nodeGaussGrad, controller);
} // namespace OpenHurricane

template <class Type>
void OpenHurricane::nodeGaussGrad::gaussGrad(
    const geometryArray<Type, cellMesh> &cellQ,
    geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &grad) const {
    const runtimeMesh &mesh = cellQ.mesh();
    for (integer i = 0; i < mesh.nCells(); ++i) {
        grad[i] = Zero;
    }
    const realArrayArray &nodeWeight = mesh.NCWgt();
    const integerListList &PNC = mesh.pointNeighbourCells();

    const integer nFaces = mesh.nFaces();
    const integer nPoints = mesh.nPoints();

    const faceList &fl = mesh.faces();

    // Value in the nodes.
    Array<Type> pV(nPoints, Zero);
    for (integer pi = 0; pi < nPoints; ++pi) {
        real weightTotal = real(0.0);
        for (integer i = 0; i < PNC[pi].size(); ++i) {
            pV[pi] += (nodeWeight[pi][i] * cellQ[PNC[pi][i]]);
            weightTotal += nodeWeight[pi][i];
        }
        pV[pi] /= weightTotal;       
    }

    // Value in the centre of the faces.
    Array<Type> fV(nFaces, Zero);
    for (integer fi = 0; fi < nFaces; ++fi) {
        for (integer i = 0; i < fl[fi].size(); ++i) {
            fV[fi] += pV[fl[fi][i]];
        }
        fV[fi] /= real(fl[fi].size());
    }
    pV.clear();

    const auto &fA = mesh.faceArea();
    for (integer i = 0; i < nFaces; ++i) {
        const auto &cl = fl[i].leftCell();
        const auto &cr = fl[i].rightCell();
        grad[cl] -= fV[i] & fA[i];
      
        if (cr < mesh.nCells()) {
            grad[cr] += fV[i] & fA[i];
        }
    }

    const auto &vol = mesh.cellVolume();
    for (integer i = 0; i < mesh.nCells(); ++i) {
        grad[i] /= vol[i];
    }
}
void OpenHurricane::nodeGaussGrad::calcGrad(const geometryArray<real, cellMesh> &cellQ,
                                            geometryArray<vector, cellMesh> &grad) const {
    gaussGrad(cellQ, grad);
}

void OpenHurricane::nodeGaussGrad::calcGrad(const geometryArray<vector, cellMesh> &cellQ,
                                            geometryArray<tensor, cellMesh> &grad) const {
    gaussGrad(cellQ, grad);
}