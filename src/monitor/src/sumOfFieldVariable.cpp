/*!
 * \file sumOfFieldVariable.cpp
 * \brief Main subroutines for computing the sum of a specified field variable on a surface.
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

#include "sumOfFieldVariable.hpp"
#include "fVArraysInclude.hpp"
namespace OpenHurricane {
    createClassName(sumOfFieldVariable);
    registerObjFty(surfaceIntegrals, sumOfFieldVariable, controller);
} // namespace OpenHurricane

OpenHurricane::sumOfFieldVariable::sumOfFieldVariable(const iteration &iter,
                                                      const runtimeMesh &mesh,
                                                      const controller &cont, const integerList &zd)
    : surfaceIntegrals(iter, mesh, cont, zd) {}

hur_nodiscard OpenHurricane::real
OpenHurricane::sumOfFieldVariable::surIntegral(const cellRealArray &phi) const {
    const auto &faces = mesh().faces();
    const auto &fW = mesh().faceWeight();
    const auto &fA = mesh().faceArea();
    const auto &fzl = mesh().faceZones();
    real fR = Zero;
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const auto fzi = zoneIdList_[i];

        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            const auto wl = fW[fi];

            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const auto phif = phi[cl] * wl + phi[cr] * (1.0 - wl);
            fR += phif;
        }
    }
    HurMPI::allReduce(fR, MPI_SUM);

    return fR;
}

hur_nodiscard OpenHurricane::vector
OpenHurricane::sumOfFieldVariable::surIntegral(const cellVectorArray &phi) const {
    const auto &faces = mesh().faces();
    const auto &fW = mesh().faceWeight();
    const auto &fA = mesh().faceArea();
    const auto &fzl = mesh().faceZones();
    vector fR = Zero;
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const auto fzi = zoneIdList_[i];

        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            const auto wl = fW[fi];

            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const auto phif = phi[cl] * wl + phi[cr] * (1.0 - wl);
            fR += phif;
        }
    }
    HurMPI::allReduceVS(fR, MPI_SUM);

    return fR;
}
