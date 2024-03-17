/*!
 * \file lift.cpp
 * \brief Main subroutines for lift force.
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

#include "lift.hpp"
#include "controllerSwitch.hpp"
#include "fVArraysInclude.hpp"
#include "gradients.hpp"

namespace OpenHurricane {
    createClassName(lift);
    registerObjFty(forces, lift, controller);
} // namespace OpenHurricane

void OpenHurricane::lift::setReportArray() const {
    vector FpV, FvV;
    real MLx;
    getForces(FpV, FvV, MLx);
    const real Fp = forceVector_ * FpV;
    const real Fv = forceVector_ * FvV;
    const real Ftotal = Fp + Fv;

    integer count = Zero;
    reportArray_[count++] = Ftotal;
    if (reportCoordinateComponent_) {
        reportArray_[count++] = FpV[0] + FvV[0];
        reportArray_[count++] = FpV[1] + FvV[1];
        reportArray_[count++] = FpV[2] + FvV[2];
    }
    if (reportPressureComponent_) {
        reportArray_[count++] = Fp;
        if (reportCoordinateComponent_) {
            reportArray_[count++] = FpV[0];
            reportArray_[count++] = FpV[1];
            reportArray_[count++] = FpV[2];
        }
    }
    if (reportViscousComponent_) {
        reportArray_[count++] = Fv;
        if (reportCoordinateComponent_) {
            reportArray_[count++] = FvV[0];
            reportArray_[count++] = FvV[1];
            reportArray_[count++] = FvV[2];
        }
    }
    if (repTp_ == reportType::coefficient) {
        reportArray_ /= (pDynFree() * mesh().Iteration().refValues().area());
    }

    if (reportPressureCenter_) {
        reportArray_[count++] = MLx / Ftotal;
    }
}

void OpenHurricane::lift::setReportName() {
    integer count = Zero;
    reportName_[count++] = "Lift";
    if (reportCoordinateComponent_) {
        reportName_[count++] = "Lx";
        reportName_[count++] = "Ly";
        reportName_[count++] = "Lz";
    }
    if (reportPressureComponent_) {
        reportName_[count++] = "Lp";
        if (reportCoordinateComponent_) {
            reportName_[count++] = "Lpx";
            reportName_[count++] = "Lpy";
            reportName_[count++] = "Lpz";
        }
    }
    if (reportViscousComponent_) {
        reportName_[count++] = "Lv";
        if (reportCoordinateComponent_) {
            reportName_[count++] = "Lvx";
            reportName_[count++] = "Lvy";
            reportName_[count++] = "Lvz";
        }
    }
    if (reportPressureCenter_) {
        reportName_[count++] = "PreCntr";
    }
}

void OpenHurricane::lift::getForces(vector &Fp, vector &Fv, real &MLx) const {
    Fp = Zero;
    Fv = Zero;
    MLx = Zero;

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

            const auto vl = v[cl] * wl + v[cr] * (real(1.0) - wl);
            const auto muf = mu[cl] * wl + mu[cr] * (real(1.0) - wl);
            const auto pf = p[cl] * wl + p[cr] * (1.0 - wl);
            const auto tauf = muf * twoSymm(gradVf[counti]) -
                              real(2.0 / 3.0) * (div(diagToVector(gradVf[counti])) * I);

            Fv += (tauf * fA[fi]);
            Fp -= (pf * fA[fi]);

            MLx += forceVector_ * (Fp + Fv) * fC[fi].x();
            counti++;
        }
    }
    HurMPI::allReduceVS(Fp, MPI_SUM);
    HurMPI::allReduceVS(Fv, MPI_SUM);
    HurMPI::allReduce(MLx, MPI_SUM);
}

OpenHurricane::lift::lift(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                          const integerList &zd)
    : forces(iter, mesh, cont, zd), forceVector_(0, 1, 0), reportPressureCenter_(false) {
    if (cont.found("forceVector")) {
        forceVector_ = cont.findType<vector>("forceVector", forceVector_);
        forceVector_ = forceVector_.normalized();
    }
    if (cont.found("options")) {
        const auto &optCont = cont.subController("options");
        reportPressureCenter_ =
            controllerSwitch(optCont)("reportPressureCenter", reportPressureCenter_);
        /*if (optCont.found("reportPressureCenter"))
        {
                string pw = optCont.findWord("reportPressureCenter");
                trim(pw);
                stringToUpperCase(pw);
                if (pw == "ON")
                {
                        reportPressureCenter_ = true;
                }
                else if (pw == "OFF")
                {
                        reportPressureCenter_ = false;
                }
                else
                {
                        errorAbortStr("Unknown type: " + pw + " in " + optCont.name(), HUR_FUNCTION);
                }
        }*/
    }
    if (reportPressureCenter_) {
        reportArray_.resize(reportName_.size() + 1);
        reportName_.resize(reportName_.size() + 1);
    }
    setReportName();
}
