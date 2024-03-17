/*!
 * \file flowRate.cpp
 * \brief Main subroutines for  computing flow rate of a quantity through a surface.
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

#include "flowRate.hpp"
#include "fVArraysInclude.hpp"

namespace OpenHurricane {
    createClassName(flowRate);
    registerObjFty(surfaceIntegrals, flowRate, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real
OpenHurricane::flowRate::getFlowRate(const cellRealArray &phi) const {
    const auto &faces = mesh().faces();
    const auto &fW = mesh().faceWeight();
    const auto &fA = mesh().faceArea();
    const auto &fzl = mesh().faceZones();
    const auto &rho = mesh().findObject<cellRealArray>("rho");
    const auto &v = mesh().findObject<cellVectorArray>("v");
    real fR = Zero;
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const auto fzi = zoneIdList_[i];

        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            const auto wl = fW[fi];

            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const auto rhol = rho[cl] * wl + rho[cr] * (1.0 - wl);
            const auto vl = v[cl] * wl + v[cr] * (1.0 - wl);
            const auto phif = phi[cl] * wl + phi[cr] * (1.0 - wl);
            fR += phif * rhol * (vl * fA[fi]);
        }
    }
    HurMPI::allReduce(fR, MPI_SUM);
    return fR;
}

hur_nodiscard OpenHurricane::vector
OpenHurricane::flowRate::getFlowRate(const cellVectorArray &phi) const {
    const auto &faces = mesh().faces();
    const auto &fW = mesh().faceWeight();
    const auto &fA = mesh().faceArea();
    const auto &fzl = mesh().faceZones();
    const auto &rho = mesh().findObject<cellRealArray>("rho");
    const auto &v = mesh().findObject<cellVectorArray>("v");
    vector fR = Zero;
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const auto fzi = zoneIdList_[i];

        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            const auto wl = fW[fi];

            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const auto rhol = rho[cl] * wl + rho[cr] * (1.0 - wl);
            const auto vl = v[cl] * wl + v[cr] * (1.0 - wl);
            const auto phif = phi[cl] * wl + phi[cr] * (1.0 - wl);
            fR += phif * (rhol * (vl * fA[fi]));
        }
    }
    HurMPI::allReduceVS(fR, MPI_SUM);

    return fR;
}

OpenHurricane::flowRate::flowRate(const iteration &iter, const runtimeMesh &mesh,
                                  const controller &cont, const integerList &zd)
    : surfaceIntegrals(iter, mesh, cont, zd) {}
