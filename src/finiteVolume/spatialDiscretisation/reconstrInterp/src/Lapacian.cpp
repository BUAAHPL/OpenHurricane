/*!
 * \file Lapacian.cpp
 * \brief Main subroutines for Lapacian.
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

#include "Lapacian.hpp"
#include "cellMesh.hpp"
#include "faceMesh.hpp"
#include "gradients.hpp"

namespace OpenHurricane {
    createClassNameStr(Lapacian, "Lapacian");
}

void OpenHurricane::Lapacian::calcLapacian(const geometryArray<real, cellMesh> &cellQ,
                                           geometryArray<real, cellMesh> &lap) {
    const auto &fgQ = fv::gradf(cellQ);

    GaussDiv::calcDiv(cellQ.grad(), fgQ, lap);
}

void OpenHurricane::Lapacian::calcLapacian(const geometryArray<vector, cellMesh> &cellQ,
                                           geometryArray<vector, cellMesh> &lap) {
    const auto &fgQ = fv::gradf(cellQ);

    GaussDiv::calcDiv(cellQ.grad(), fgQ, lap);
}

OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>
OpenHurricane::Lapacian::lapacian(const geometryArray<real, cellMesh> &cellQ) {
    geometryArray<real, cellMesh> ff(object("lapacian(" + cellQ.name() + ")", cellQ.mesh(),
                                            object::NOT_WRITE, object::TEMPORARY),
                                     cellQ.mesh());
    calcLapacian(cellQ, ff);
    return ff;
}

OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::cellMesh>
OpenHurricane::Lapacian::lapacian(const geometryArray<vector, cellMesh> &cellQ) {
    geometryArray<vector, cellMesh> ff(object("lapacian(" + cellQ.name() + ")", cellQ.mesh(),
                                              object::NOT_WRITE, object::TEMPORARY),
                                       cellQ.mesh());
    calcLapacian(cellQ, ff);
    return ff;
}
