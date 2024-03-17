#include "faceGrad.hpp"
/*!
 * \file faceGrad.inl
 * \brief In-Line subroutines of the <i>faceGrad.hpp</i> file.
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
#pragma once

inline OpenHurricane::geometryArray<
    typename OpenHurricane::outerProduct<OpenHurricane::vector, OpenHurricane::real>::type,
    OpenHurricane::faceMesh>
OpenHurricane::faceGrad::grad(
    const geometryArray<real, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &cellGrad) {
    const runtimeMesh &mesh = cellQ.mesh();
    const integer nFaces = mesh.nFaces();

    geometryArray<typename outerProduct<vector, real>::type, faceMesh> fG(
        object(cellQ.name() + "_faceGrad", mesh, object::NOT_WRITE, object::TEMPORARY), mesh);

    grad(cellQ, cellGrad, fG);

    return fG;
}

inline OpenHurricane::geometryArray<
    typename OpenHurricane::outerProduct<OpenHurricane::vector, OpenHurricane::vector>::type,
    OpenHurricane::faceMesh>
OpenHurricane::faceGrad::grad(
    const geometryArray<vector, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &cellGrad) {
    const runtimeMesh &mesh = cellQ.mesh();
    const integer nFaces = mesh.nFaces();

    geometryArray<typename outerProduct<vector, vector>::type, faceMesh> fG(
        object(cellQ.name() + "_faceGrad", mesh, object::NOT_WRITE, object::TEMPORARY), mesh);

    grad(cellQ, cellGrad, fG);
    return fG;
}