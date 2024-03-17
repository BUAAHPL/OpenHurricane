/*!
 * \file viscousFlux.cpp
 * \brief Main subroutines of viscous flux.
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

#include "viscousFlux.hpp"

void OpenHurricane::visFluxPerZone(const geometryArray<OpenHurricane::symmTensor, faceMesh> &gf,
                                   const integer faceZoneId, Array<OpenHurricane::vector> &rhs) {
    const auto &mesh = gf.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fA = mesh.faceArea();

    for (integer fi = fzl[faceZoneId].firstIndex(); fi <= fzl[faceZoneId].lastIndex(); ++fi) {
        const auto &cl = fl[fi].leftCell();
        const auto &cr = fl[fi].rightCell();
        const vector flux = gf[fi] * fA[fi];

        rhs[cl] -= flux;
        if (cr < mesh.nCells()) {
            rhs[cr] += flux;
        }
    }
}

OpenHurricane::Array<OpenHurricane::vector> OpenHurricane::visFlux(
    const OpenHurricane::geometryArray<OpenHurricane::symmTensor, OpenHurricane::faceMesh> &gf) {
    const auto &mesh = gf.mesh();
    Array<vector> rhs(mesh.nTotalCells(), Zero);

    const auto &fzl = mesh.faceZones();

    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        visFluxPerZone(gf, fzi, rhs);
    }
    return rhs;
}

void OpenHurricane::visFlux(
    const OpenHurricane::geometryArray<OpenHurricane::symmTensor, OpenHurricane::faceMesh> &gf,
    geometryArray<vector, cellMesh> &cellQ) {
    const auto &mesh = gf.mesh();
    auto &rhs = cellQ.rhs();

    const auto &fzl = mesh.faceZones();

    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        visFluxPerZone(gf, fzi, rhs);
    }
}