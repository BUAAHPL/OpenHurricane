/*!
 * \file moment.cpp
 * \brief Main subroutines for moment.
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

#include "moment.hpp"
#include "fVArraysInclude.hpp"
#include "gradients.hpp"

namespace OpenHurricane {
    createClassName(moment);
    registerObjFty(forces, moment, controller);
} // namespace OpenHurricane

void OpenHurricane::moment::setReportArray() const {
    vector MpV, MvV;
    getMoment(MpV, MpV);
    const real Mp = momentAxis_ * MpV;
    const real Mv = momentAxis_ * MvV;
    const real Mtotal = Mp + Mv;

    integer count = Zero;
    reportArray_[count++] = Mtotal;
    if (reportCoordinateComponent_) {
        reportArray_[count++] = MpV[0] + MvV[0];
        reportArray_[count++] = MpV[1] + MvV[1];
        reportArray_[count++] = MpV[2] + MvV[2];
    }
    if (reportPressureComponent_) {
        reportArray_[count++] = Mp;
        if (reportCoordinateComponent_) {
            reportArray_[count++] = MpV[0];
            reportArray_[count++] = MpV[1];
            reportArray_[count++] = MpV[2];
        }
    }
    if (reportViscousComponent_) {
        reportArray_[count++] = Mv;
        if (reportCoordinateComponent_) {
            reportArray_[count++] = MvV[0];
            reportArray_[count++] = MvV[1];
            reportArray_[count++] = MvV[2];
        }
    }
    if (repTp_ == reportType::coefficient) {
        reportArray_ /= (pDynFree() * mesh().Iteration().refValues().area() *
                         mesh().Iteration().refValues().length());
    }
}

void OpenHurricane::moment::setReportName() {
    integer count = Zero;
    reportName_[count++] = "Moment";
    if (reportCoordinateComponent_) {
        reportName_[count++] = "Mx";
        reportName_[count++] = "My";
        reportName_[count++] = "Mz";
    }
    if (reportPressureComponent_) {
        reportName_[count++] = "Mp";
        if (reportCoordinateComponent_) {
            reportName_[count++] = "Mpx";
            reportName_[count++] = "Mpy";
            reportName_[count++] = "Mpz";
        }
    }
    if (reportViscousComponent_) {
        reportName_[count++] = "Mv";
        if (reportCoordinateComponent_) {
            reportName_[count++] = "Mvx";
            reportName_[count++] = "Mvy";
            reportName_[count++] = "Mvz";
        }
    }
}

void OpenHurricane::moment::getMoment(vector &Mp, vector &Mv) const {
    Mp = Zero;
    Mv = Zero;

    const auto &faces = mesh().faces();
    const auto &fW = mesh().faceWeight();
    const auto &fA = mesh().faceArea();
    const auto &fC = mesh().faceCentre();
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

            const vector Fv = (tauf * fA[fi]);
            const vector Fp = -(pf * fA[fi]);

            Mp += ((fC[fi] - momentCenter_) ^ Fp);
            Mv += ((fC[fi] - momentCenter_) ^ Fv);
            counti++;
        }
    }
    HurMPI::allReduceVS(Mp, MPI_SUM);
    HurMPI::allReduceVS(Mv, MPI_SUM);
}

OpenHurricane::moment::moment(const iteration &iter, const runtimeMesh &mesh,
                              const controller &cont, const integerList &zd)
    : forces(iter, mesh, cont, zd), momentCenter_(0, 0, 0), momentAxis_(0, 0, 1) {
    if (cont.found("momentCenter")) {
        momentCenter_ = cont.findType<vector>("momentCenter", momentCenter_);
        momentCenter_ = momentCenter_.normalized();
    }
    if (cont.found("momentAxis")) {
        momentAxis_ = cont.findType<vector>("momentAxis", momentAxis_);
        momentAxis_ = momentAxis_.normalized();
    }

    setReportName();
}
