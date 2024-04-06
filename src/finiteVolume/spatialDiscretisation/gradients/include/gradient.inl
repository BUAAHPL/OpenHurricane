#include "gradient.hpp"
/*!
 * \file gradient.inl
 * \brief In-Line subroutines of the <i>gradient.hpp</i> file.
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
#pragma once

template <class Type>
inline void OpenHurricane::gradient::grad(
    const geometryArray<Type, cellMesh> &cellQ,
    geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &grad) const {
    calcGrad(cellQ, grad);
    updateBoundaryField(grad);
}

template <class Type>
inline OpenHurricane::geometryArray<
    typename OpenHurricane::outerProduct<OpenHurricane::vector, Type>::type,
    OpenHurricane::cellMesh>
OpenHurricane::gradient::grad(const geometryArray<Type, cellMesh> &cellQ,
                              const string &name) const {
    const runtimeMesh &rMesh = cellQ.mesh();
    geometryArray<typename outerProduct<vector, Type>::type, cellMesh> gradients(
        object(name, rMesh, object::NOT_WRITE, object::TEMPORARY), rMesh);
    grad(cellQ, gradients);
    return gradients;
}

template <class Type>
inline void OpenHurricane::gradient::grad(geometryArray<Type, cellMesh> &cellQ) const {
    calcGrad(cellQ, cellQ.grad());
    updateBoundaryField(cellQ.grad());
}
