/*!
 * \file SpalartAllmarasBoundarySet.cpp
 * \brief Bounadry conditions setting subroutines for the SpalartAllmaras turbulence model.
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

#include "SpalartAllmaras.hpp"

void OpenHurricane::SpalartAllmaras::bndValueSetting(const controller &cont) {
    const controller &interCont = cont.subController("boundaryCondition");
    const faceZoneList &fZL = mesh().faceZones();

    for (integer i = 0; i < fZL.size(); i++) {
        real nut = 0.0;
        const auto bcKey = interCont.subController(fZL[i].name()).findWord("bcType");
        if (interCont.found(fZL[i].name())) {
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
                //realArray yi = flowM_.thermo().getYi(bcCont, flowM_.mixtures(), false);
                auto yi =
                    getBoundariesFromController::getSpeciesMassFractions(bcCont, flowM_.mixtures());
                real mu = flowM_.mixtures().transTable().mu(p, T, yi);

                inletBndSetting(bcCont);
                const real cmu = 0.09;
                switch (bndSpMethod_) {
                case OpenHurricane::RANSModel::viscosityRatio:
                    nut = getViscosityRatio() * mu / rho;
                    break;
                case OpenHurricane::RANSModel::intensityAndLength:
                    nut = cmu * sqrt(real(1.5)) * v * intensity() * getLengthScale();
                    break;
                case OpenHurricane::RANSModel::intensityAndHYdraulicD:
                    nut = cmu * sqrt(real(1.5)) * v * intensity() * getHydraulicD();
                    break;
                case OpenHurricane::RANSModel::origianalTurbEquation:
                    if (bcCont.found("nutFactor")) {
                        const real nutFactor = bcCont.findType<real>("nutFactor", nutFactor);
                        nut = nutFactor * mu / rho;
                        //bcCont.add<real>(std::string("nut"), nut);
                    } else {
                        LFatal("Turbulent parameters specification found in boundary type: %s are "
                               "not enough.",
                               bcCont.findWord("bcType").c_str());
                    }
                    break;
                case OpenHurricane::RANSModel::givenDirectly:
                    if (bcCont.found("nut")) {
                        nut = bcCont.findType<real>("nut", nut);
                    } else {
                        LFatal("Turbulent parameters specification found in boundary type: %s are "
                               "not enough.",
                               bcCont.findWord("bcType").c_str());
                    }
                    break;
                default:
                    break;
                }
                bcCont.add<real>(std::string("nut"), nut);
                bcCont.add(std::string("mutbcType"), string("mutSpalartAllmarasInlet"));
            }

            else if (interCont.subController(fZL[i].name()).findWord("bcType") == "wall") {
                controller &bcCont =
                    const_cast<controller &>(interCont.subController(fZL[i].name()));
                bcCont.add(std::string("mutbcType"), string("mutLowReWallTreatment"));
                string wallFunction = bcCont.findWord("wallFunction");
                if (wallFunction == string("on")) {
                    LFatal("Wall boundary condition of faceZone named %s does not include wall "
                           "function in SA model",
                           fZL[i].name().c_str());
                } else if (wallFunction == string("off")) {
                    bcCont.add(std::string("nutbcType"), string("nutWallTreatment"));
                }
            }
        }
    }
}