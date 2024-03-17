/*!
 * \file areaWeightedUniformityIndex.cpp
 * \brief Main subroutines for computing the area-weighted uniformity index on a surface.
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

#include "areaWeightedUniformityIndex.hpp"
#include "fVArraysInclude.hpp"

namespace OpenHurricane {
    createClassName(areaWeightedUniformityIndex);
    registerObjFty(surfaceIntegrals, areaWeightedUniformityIndex, controller);
} // namespace OpenHurricane

OpenHurricane::areaWeightedUniformityIndex::areaWeightedUniformityIndex(const iteration &iter,
                                                                        const runtimeMesh &mesh,
                                                                        const controller &cont,
                                                                        const integerList &zd)
    : areaWeightedAverage(iter, mesh, cont, zd) {}

hur_nodiscard OpenHurricane::real
OpenHurricane::areaWeightedUniformityIndex::surIntegral(const cellRealArray &phi) const {
    const auto &faces = mesh().faces();
    const auto &fW = mesh().faceWeight();
    const auto &fA = mesh().faceArea();
    const auto &fzl = mesh().faceZones();
    const real phia = areaWeightedAverage::surIntegral(phi);
    real fR = Zero;
    real A = Zero;
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const auto fzi = zoneIdList_[i];

        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            const auto wl = fW[fi];

            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();
            const auto phif = phi[cl] * wl + phi[cr] * (1.0 - wl);
            fR += mag(phif - phia) * fA[fi].magnitude();
            A += fA[fi].magnitude();
        }
    }
    HurMPI::allReduce(fR, MPI_SUM);
    HurMPI::allReduce(A, MPI_SUM);

    return 1 - fR / max(tiny, real(2) * mag(phia) * A);
}

hur_nodiscard OpenHurricane::vector
OpenHurricane::areaWeightedUniformityIndex::surIntegral(const cellVectorArray &phi) const {
    const auto &faces = mesh().faces();
    const auto &fW = mesh().faceWeight();
    const auto &fA = mesh().faceArea();
    const auto &fzl = mesh().faceZones();
    const vector phia = areaWeightedAverage::surIntegral(phi);
    real A = Zero;
    vector fR = Zero;
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const auto fzi = zoneIdList_[i];

        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            const auto wl = fW[fi];

            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();
            const auto phif = phi[cl] * wl + phi[cr] * (1.0 - wl);
            fR += componentMag(phif - phia) * fA[fi].magnitude();
            A += fA[fi].magnitude();
        }
    }
    HurMPI::allReduceVS(fR, MPI_SUM);
    HurMPI::allReduce(A, MPI_SUM);

    return vector(1) -
           componentDivide(fR, componentMax(vector(tiny), real(2) * componentMag(phia) * A));
}
