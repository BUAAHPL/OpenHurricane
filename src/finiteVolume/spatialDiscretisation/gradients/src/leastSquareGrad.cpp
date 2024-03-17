/*!
 * \file leastSquareGrad.cpp
 * \brief Main subroutines for least square gradient approach.
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

#include "leastSquareGrad.hpp"

namespace OpenHurricane {
    createClassNameStr(leastSquareGrad, "leastSquareGrad");
    registerObjFty(gradient, leastSquareGrad, controller);
} // namespace OpenHurricane

template <class Type>
void OpenHurricane::leastSquareGrad::WLSGrad(
    const geometryArray<Type, cellMesh> &cellQ,
    geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &grad) const {
    const runtimeMesh &mesh = cellQ.mesh();
    const integer nCells = mesh.nCells();
    for (integer celli = 0; celli < nCells; celli++) {
        grad[celli] = Zero;
    }
    const vectorArrayArray &cellWeight = mesh.cWgt(key_);
    const integerListList &CNC = mesh.cellNeighbourCells();

    for (integer celli = 0; celli < nCells; celli++) {
        for (integer i = 0; i < CNC[celli].size(); i++) {
            const Type delta = cellQ[CNC[celli][i]] - cellQ[celli];
            grad[celli] += delta * cellWeight[celli][i];
        }
    }
}

void OpenHurricane::leastSquareGrad::calcGrad(const geometryArray<real, cellMesh> &cellQ,
                                              geometryArray<vector, cellMesh> &grad) const {
    WLSGrad(cellQ, grad);
}

void OpenHurricane::leastSquareGrad::calcGrad(const geometryArray<vector, cellMesh> &cellQ,
                                              geometryArray<tensor, cellMesh> &grad) const {    
    const runtimeMesh &mesh = cellQ.mesh();
    const integer nCells = mesh.nCells();
    for (integer celli = 0; celli < nCells; celli++) {
        grad[celli] = Zero;
    }
    const vectorArrayArray &cellWeight = mesh.cWgt(key_);
    const integerListList &CNC = mesh.cellNeighbourCells();
    for (integer celli = 0; celli < nCells; celli++) {
        for (integer i = 0; i < CNC[celli].size(); i++) {
            const vector delta = cellQ[CNC[celli][i]] - cellQ[celli];
            grad[celli] += delta & cellWeight[celli][i];
        }
    }
}