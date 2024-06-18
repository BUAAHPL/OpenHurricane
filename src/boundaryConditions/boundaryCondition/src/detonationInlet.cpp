/*!
 * \file detonationInlet.cpp
 * \brief Main subroutines for detonationInlet.
 * \author Chen Zhenyi
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
#include "detonationInlet.hpp"
#include "parsingDirection.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::detonationInlet::className_ = "detonationInlet";

namespace OpenHurricane {
    registerObjFty(vectorBoundary, detonationInlet, controller);
} // namespace OpenHurricane

void OpenHurricane::detonationInlet::getValue(real &value, const std::string &name,
                                              const controller &cont) const {
    if (cont.found(name)) {
        value = cont.findType<real>(name, value);
    } else {
        LFatal("Boundary %s does not specify %s value.\n Please check.",
               boundaryZone_.name().c_str(), name.c_str());
    }
}

OpenHurricane::detonationInlet::detonationInlet(const faceZone &fZ,
                                                geometryArray<vector, cellMesh> &gf,
                                                const controller &cont)
    : Base(fZ, gf, cont), directionType_(normalToBoundary), inletType_() {
    Base::setSpecified();
    getValue(Pt_, std::string("totalPressure"), cont);

    getValue(Tt_, std::string("totalTemperature"), cont);

    const auto dtw = parsingDirection::getDirection(const_cast<controller &>(cont),
                                                    varArray_.mesh(), fZ, direct_);
    if (dtw != "normalToBoundary") {
        directionType_ = givenDirect;
    }
    const auto &tb = this->varField().tb();
    const auto &thTable = tb.findObject<thermoList>("thermoList");
    yiBnd_.resize(thTable.species().size(), real(1));
    for (integer i = 0; i < thTable.species().size(); ++i) {
        auto &spn = thTable.species()[i].name();
        yiBnd_[i] = cont.findOrDefault<real>(spn, 0.0);
    }
}

void OpenHurricane::detonationInlet::updateBoundary() {
    const auto &tb = this->varField().tb();
    auto &v = this->varField();
    auto &rho = tb.findObjectRef<cellRealArray>("rho");
    auto &p = tb.findObjectRef<cellRealArray>("p");
    auto &T = tb.findObjectRef<cellRealArray>("T");
    const auto &g = tb.findObject<cellRealArray>("gamma");

    const auto &thTable = tb.findObject<thermoList>("thermoList");
    const auto &fL = this->mesh().faces();
    const auto &fA = this->mesh().faceArea();
    const auto cstep = p.mesh().Iteration().cStep();
    integer subCstep = -1;
    if (p.mesh().Iteration().hasSubIteration()) {
        subCstep = p.mesh().Iteration().subIter().cSubStep();
    }
    if (subCstep > 1) {
        return;
    }

    real Ma = 1.0;
    real g0 = 1.4;
    real Tin = Tt_ / (1 + 0.5 * (g0 - 1) * sqr(Ma));
    real pin = Pt_ / pow(1 + 0.5 * (g0 - 1) * sqr(Ma), g0 / (g0 - 1));
    real gin = thTable.gamma(pin, Tin, yiBnd_);
    integer count = 0;

    vector vbnd = Zero;
    real pbnd = Zero;
    real rhobnd = Zero;
    real Tbnd = Zero;

    while (fabs(gin - g0) / g0 > 1e-4) {
        g0 = gin;
        Tin = Tt_ / (1 + 0.5 * (g0 - 1) * sqr(Ma));
        pin = Pt_ / pow(1 + 0.5 * (g0 - 1) * sqr(Ma), g0 / (g0 - 1));
        gin = thTable.gamma(pin, Tin, yiBnd_);

        if (count++ > 1000) {
            PLWarning(
                "Maximum number of iterations exceeded: 1000 in detonationInlet undateBoundary");
            break;
        }
    }
    const real Pcr = Pt_ / pow(1 + 0.5 * (gin - 1) * sqr(Ma), gin / (gin - 1));
    const real Rg = thTable.species().Rm(yiBnd_);

    if (inletType_.size() != boundaryZone_.size()) {
        inletType_.resize(boundaryZone_.size());
    }

    for (integer fi = boundaryZone_.firstIndex(); fi < boundaryZone_.lastIndex() + 1; fi++) {
        const integer celli = fL[fi].leftCell();
        const integer ghostId = fL[fi].rightCell();
        vector normal = fA[fi].normalized();
        bool negativeGhostDensity = false;
        real pi0 = p[celli];
        if (cstep == 0 && subCstep > 0) {
            pi0 = p.lastArray()[celli];
        }

        auto id = fi - boundaryZone_.firstIndex();
        if (subCstep <= 1) {
            if (pi0 >= Pt_) {
                inletType_[id] = inletStateType::blocking;
            } else {
                inletType_[id] = inletStateType::supersonic;
            }
        }

        if (inletType_[id] != inletStateType::blocking) {
            if (Pcr < pi0 && pi0 < Pt_) {
                inletType_[id] = inletStateType::subsonic;
            } else {
                inletType_[id] = inletStateType::supersonic;
            }
        }

        if (inletType_[id] == inletStateType::blocking) {
            vbnd = Zero;

            pbnd = p[celli];
            Tbnd = T[celli];
            rhobnd = rho[celli];

            p[ghostId] = p[celli];
            T[ghostId] = T[celli];
            rho[ghostId] = rho[celli];
        }

        else {
            real vv = 0;
            if (inletType_[id] == inletStateType::subsonic) {
                pbnd = pi0;
                vv =
                    sqrt(2.0 * gin / (gin - 1) * Rg * Tt_ * (1 - pow(pbnd / Pt_, (gin - 1) / gin)));
            } else {
                pbnd = Pcr;
                vv = sqrt(2.0 * gin / (gin + 1) * Rg * Tt_);
            }
            Tbnd = Tt_ * pow(pbnd / Pt_, (gin - 1) / gin);
           
            if (directionType_ == normalToBoundary) {
                vbnd = vv * normal;
            } else {
                vbnd = vv * direct_;
            }

            rhobnd = pbnd / Rg / Tbnd;

            if (2.0 * rhobnd - rho[celli] < 0.0) {
                negativeGhostDensity = true;
                rho[ghostId] = rhobnd;
                p[ghostId] = pbnd;
                T[ghostId] = Tbnd;
            } else {
                p[ghostId] = real(2.0) * pbnd - p[celli];
                T[ghostId] = real(2.0) * Tbnd - T[celli];
                rho[ghostId] = real(2.0) * rhobnd - rho[celli];
            }
        }
        v[ghostId] = real(2.0) * vbnd - v[celli];
    }
}
