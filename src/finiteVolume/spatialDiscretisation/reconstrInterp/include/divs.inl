#include "divs.hpp"
/*!
 * \file divs.inl
 * \brief In-Line subroutines of the <i>divs.hpp</i> file.
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

inline void OpenHurricane::fv::div(const geometryArray<vector, cellMesh> &cellQ,
                                   geometryArray<real, cellMesh> &divQ) {
    GaussDiv::calcDiv(cellQ, divQ);
}

inline void OpenHurricane::fv::div(const geometryArray<tensor, cellMesh> &cellQ,
                                   geometryArray<vector, cellMesh> &divQ) {
    GaussDiv::calcDiv(cellQ, divQ);
}

inline void OpenHurricane::fv::div(const geometryArray<vector, cellMesh> &cellQ,
                                   const geometryArray<vector, faceMesh> &faceQ,
                                   geometryArray<real, cellMesh> &divQ) {
    GaussDiv::calcDiv(cellQ, faceQ, divQ);
}

inline void OpenHurricane::fv::div(const geometryArray<tensor, cellMesh> &cellQ,
                                   const geometryArray<tensor, faceMesh> &faceQ,
                                   geometryArray<vector, cellMesh> &divQ) {
    GaussDiv::calcDiv(cellQ, faceQ, divQ);
}

inline OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>
OpenHurricane::fv::div(const geometryArray<vector, cellMesh> &cellQ) {
    return GaussDiv::div(cellQ);
}

inline OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::cellMesh>
OpenHurricane::fv::div(const geometryArray<tensor, cellMesh> &cellQ) {
    return GaussDiv::div(cellQ);
}

inline OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>
OpenHurricane::fv::div(const geometryArray<vector, cellMesh> &cellQ,
                       const geometryArray<vector, faceMesh> &faceQ) {
    return GaussDiv::div(cellQ, faceQ);
}

inline OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::cellMesh>
OpenHurricane::fv::div(const geometryArray<tensor, cellMesh> &cellQ,
                       const geometryArray<tensor, faceMesh> &faceQ) {
    return GaussDiv::div(cellQ, faceQ);
}