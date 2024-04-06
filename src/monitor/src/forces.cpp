/*!
 * \file forces.cpp
 * \brief Main subroutines for base class of computing forces.
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

#include "forces.hpp"
#include "controllerSwitch.hpp"
#include "gradients.hpp"

namespace OpenHurricane {
    createClassName(forces);
    createObjFty(forces, controller);
} // namespace OpenHurricane

void OpenHurricane::forces::getForces(vector &Fp, vector &Fv) const {
    Fp = Zero;
    Fv = Zero;

    const auto &faces = mesh().faces();
    const auto &fW = mesh().faceWeight();
    const auto &fA = mesh().faceArea();
    const auto &fzl = mesh().faceZones();
    const auto &v = mesh().findObject<cellVectorArray>("v");
    const auto &mu = mesh().findObject<cellRealArray>("mu");
    const auto &p = mesh().findObject<cellRealArray>("p");
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const auto fzi = zoneIdList_[i];

        const auto gradVf = fv::gradf(v, fzi);

        integer counti = Zero;
        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            const auto wl = fW[fi];

            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const auto vl = v[cl] * wl + v[cr] * (1.0 - wl);
            const auto muf = mu[cl] * wl + mu[cr] * (1.0 - wl);
            const auto pf = p[cl] * wl + p[cr] * (1.0 - wl);
            const auto tauf = muf * twoSymm(gradVf[counti]) -
                              real(2.0 / 3.0) * (div(diagToVector(gradVf[counti])) * I);

            Fv += (tauf * fA[fi]);
            Fp -= (pf * fA[fi]);
            counti++;
        }
    }
    HurMPI::allReduceVS(Fp, MPI_SUM);
    HurMPI::allReduceVS(Fv, MPI_SUM);
}

void OpenHurricane::forces::setArraySize() {
    integer count = Zero;
    count++;
    if (reportCoordinateComponent_) {
        count++;
        count++;
        count++;
    }
    if (reportPressureComponent_) {
        count++;
        if (reportCoordinateComponent_) {
            count++;
            count++;
            count++;
        }
    }
    if (reportViscousComponent_) {
        count++;
        if (reportCoordinateComponent_) {
            count++;
            count++;
            count++;
        }
    }
    reportArray_.resize(count);
    reportName_.resize(count);
}

OpenHurricane::real OpenHurricane::forces::pDynFree() const {
    real pdyn =
        0.5 * mesh().Iteration().refValues().rho() * sqr(mesh().Iteration().refValues().vMag());

    return pdyn;
}

OpenHurricane::forces::forces(const iteration &iter, const runtimeMesh &mesh,
                              const controller &cont, const integerList &zoneIdList)
    : iter_(iter), mesh_(mesh), zoneIdList_(zoneIdList), repTp_(reportType::force),
      reportPressureComponent_(false), reportViscousComponent_(false),
      reportCoordinateComponent_(false), reportArray_(), reportName_() {
    if (cont.found("reportType")) {
        string rt = cont.findWord("reportType");
        trim(rt);
        if (rt == "force") {
            repTp_ = reportType::force;
        } else if (rt == "coefficient") {
            repTp_ = reportType::coefficient;
        } else {
            errorAbortStr(("Unknown report type: " + rt + " in " + cont.name()));
        }
    }
    if (cont.found("options")) {
        const auto &optCont = cont.subController("options");
        reportPressureComponent_ =
            controllerSwitch(optCont)("reportPressureComponent", reportPressureComponent_);

        reportViscousComponent_ =
            controllerSwitch(optCont)("reportViscousComponent", reportViscousComponent_);

        reportCoordinateComponent_ =
            controllerSwitch(optCont)("reportCoordinateComponent", reportCoordinateComponent_);
    }
    setArraySize();
}

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::forces>
OpenHurricane::forces::creator(const iteration &iter, const runtimeMesh &mesh,
                               const controller &cont, const integerList &zd) {
    string forceType = cont.findWord("forceType");
    defineInObjCreator(forces, forceType, controller, (iter, mesh, cont, zd));
}