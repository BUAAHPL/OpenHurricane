/*!
 * \file force.cpp
 * \brief Main subroutines for force.
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

#include "force.hpp"
#include "fVArraysInclude.hpp"

namespace OpenHurricane {
    createClassName(force);
    registerObjFty(forces, force, controller);
} // namespace OpenHurricane

OpenHurricane::real OpenHurricane::force::getForce() const {
    vector Fp, Fv;
    getForces(Fp, Fv);
    return forceVector_ * (Fp + Fv);
}

OpenHurricane::real OpenHurricane::force::getForce(real &Fp, real &Fv) const {
    vector FpV, FvV;
    getForces(FpV, FvV);
    Fp = forceVector_ * FpV;
    Fv = forceVector_ * FvV;
    return Fp + Fv;
}

void OpenHurricane::force::setReportArray() const {
    vector FpV, FvV;
    getForces(FpV, FvV);
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
}

void OpenHurricane::force::setReportName() {
    integer count = Zero;
    reportName_[count++] = "Force";
    if (reportCoordinateComponent_) {
        reportName_[count++] = "Fx";
        reportName_[count++] = "Fy";
        reportName_[count++] = "Fz";
    }
    if (reportPressureComponent_) {
        reportName_[count++] = "Fp";
        if (reportCoordinateComponent_) {
            reportName_[count++] = "Fpx";
            reportName_[count++] = "Fpy";
            reportName_[count++] = "Fpz";
        }
    }
    if (reportViscousComponent_) {
        reportName_[count++] = "Fv";
        if (reportCoordinateComponent_) {
            reportName_[count++] = "Fvx";
            reportName_[count++] = "Fvy";
            reportName_[count++] = "Fvz";
        }
    }
}

OpenHurricane::force::force(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                            const integerList &zd)
    : forces(iter, mesh, cont, zd), forceVector_(1, 0, 0) {
    if (cont.found("forceVector")) {
        forceVector_ = cont.findType<vector>("forceVector", forceVector_);
        forceVector_ = forceVector_.normalized();
    }
    setReportName();
}
