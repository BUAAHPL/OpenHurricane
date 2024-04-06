/*!
 * \file SSTBoundarySet.cpp
 * \brief Bounadry conditions setting subroutines for the SST turbulence model.
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

#include "SST.hpp"

void OpenHurricane::SST::bndValueSetting(const controller &cont) {
    const controller &interCont = cont.subController("boundaryCondition");
    const faceZoneList &fZL = mesh().faceZones();

    for (integer i = 0; i < fZL.size(); i++) {
        real k = 0.0;
        real w = 0.0;
        bool steamb = false;
        if (interCont.found(fZL[i].name())) {
            const auto bcKey = interCont.subController(fZL[i].name()).findWord("bcType");
            if (bcKey == "pressureFarField" || bcKey == "pressure-far-field" ||
                bcKey == "pressureInlet" || bcKey == "pressure-inlet" || bcKey == "velocityInlet" ||
                bcKey == "velocity-inlet" || bcKey == "massFlowInlet" ||
                bcKey == "mass-flow-inlet" || bcKey == "subsonicInlet" ||
                bcKey == "supersonicInlet" || bcKey == "syntheticTurbInlet") {
                controller &bcCont =
                    const_cast<controller &>(interCont.subController(fZL[i].name()));

                real rho = bcCont.findType<real>("rho", rho);
                real p = bcCont.findType<real>("p", p);
                real T = bcCont.findType<real>("T", T);
                real v = bcCont.findType<real>("v", v);
                auto yi =
                    getBoundariesFromController::getSpeciesMassFractions(bcCont, flowM_.mixtures());
                real mu = flowM_.mixtures().transTable().mu(p, T, yi);

                inletBndSetting(bcCont);

                const real cmu = 0.09;
                switch (bndSpMethod_) {
                case OpenHurricane::RANSModel::viscosityRatio:
                    k = real(1.5) * sqr(v * intensity());
                    w = rho * k / (mu * getViscosityRatio());
                    break;
                case OpenHurricane::RANSModel::intensityAndLength:
                    k = real(1.5) * sqr(v * intensity());
                    w = sqrt(k) / (cmu * getLengthScale());
                    break;
                case OpenHurricane::RANSModel::intensityAndHYdraulicD:
                    k = real(1.5) * sqr(v * intensity());
                    w = sqrt(k) / (pow(cmu, real(0.25)) * real(0.07) * getHydraulicD());
                    break;
                case OpenHurricane::RANSModel::origianalTurbEquation:
                    if (bcCont.found("wFactor")) {
                        real wFactor = 0.0;
                        real comDomainLength = 0.0;
                        real kFactor = 0.0;
                        wFactor = bcCont.findType<real>("wFactor", wFactor);
                        comDomainLength = bcCont.findType<real>("comDomainLength", comDomainLength);

                        if (bcCont.found("kFactor")) {
                            kFactor = bcCont.findType<real>("kFactor", kFactor);
                            w = wFactor * v / comDomainLength;
                            k = kFactor * mu / rho * v / comDomainLength;
                        } else {
                            LFatal(
                                "Turbulent parameters specification found in boundary type: %s are "
                                "not enough.",
                                bcCont.findWord("bcType").c_str());
                        }
                    }
                    break;
                case OpenHurricane::RANSModel::givenDirectly:
                    if (bcCont.found("w")) {
                        w = bcCont.findType<real>("w", w);
                    } else {
                        LFatal("Turbulent parameters specification found in boundary type: %s are "
                               "not enough.",
                               bcCont.findWord("bcType").c_str());
                    }
                    if (bcCont.found("k")) {
                        k = bcCont.findType<real>("k", k);
                    } else {
                        LFatal("Turbulent parameters specification found in boundary type: %s are "
                               "not enough.",
                               bcCont.findWord("bcType").c_str());
                    }
                    break;
                default:
                    break;
                }

                if (!steamb) {
                    kamb_ = k;
                    wamb_ = w;
                    steamb = true;
                }
                bcCont.add<real>(std::string("kt"), k);
                bcCont.add<real>(std::string("wt"), w);
                bcCont.add(std::string("mutbcType"), string("mutKOmegaInlet"));
            } else if (bcKey == "wall") {
                controller &bcCont =
                    const_cast<controller &>(interCont.subController(fZL[i].name()));
                omegaWallFunction_ = controllerSwitch(bcCont)("wallFunction", omegaWallFunction_);

                bcCont.add(std::string("mutbcType"), string("mutLowReWallTreatment"));

                if (omegaWallFunction_) {
                    bcCont.add(std::string("ktbcType"), string("kWallTreatment"));
                    bcCont.add(std::string("wtbcType"), string("omegaNearWallTreatment"));
                    omegaWallFunctionFaceZoneList_.push_back(i);
                } else {
                    bcCont.add(std::string("ktbcType"), string("kWallTreatment"));
                    bcCont.add(std::string("wtbcType"), string("omegaNearWallTreatment"));
                }
            }
        }
    }
}
