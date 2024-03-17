/*!
 * \file getSyntheticTurbInlet.cpp
 * \brief Main subroutines for parsing synthetic turbulence generator inlet boundary condition from controller.
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
#include "getBoundariesFromController.hpp"
#include "rhoThermo.hpp"

void OpenHurricane::getBoundariesFromController::getSyntheticTurbInlet(const mixture &mixtures,
                                                                   const controller &cont,
                                                                   controller &bcCont,
                                                                   const faceZone &fz) {
    bcCont.add(std::string("defultType"), string("fixedValue"));
    bcCont.add(std::string("vbcType"), string("syntheticTurbInlet"));
    bcCont.add(std::string("rhobcType"), string("interior"));
    bcCont.add(std::string("pbcType"), string("interior"));
    bcCont.add(std::string("TbcType"), string("interior"));
    auto yi = getSpeciesMassFractions(bcCont, mixtures, true);
    vector direct;
    const auto directTypeW =
        parsingDirection::getDirection(bcCont, mixtures.mesh(), fz, direct, true);
    if (!bcCont.found("T")) {
        LFatal("The temperature magnitude must specified in face zone: %s", fz.name().c_str());
    }
    real T = bcCont.findType<real>("T", T);

    if (!bcCont.found("momentumGivenBy")) {
        LFatal("The momentum specification must specified in face zone: %s", fz.name().c_str());
    }
    const auto gbw = bcCont.findWord("momentumGivenBy");
    real p = constant::physicalConstant::Patm;

    real g = mixtures.thermalTable().gamma(p, T, yi);
    real mufree = mixtures.transTable().mu(p, T, yi);

    real a = sqrt(g * mixtures.species().Rm(yi) * T);

    real ma;
    real vMag = 0;
    real Re;
    real rho = Zero;
    if (gbw == "pressureAndMach") {
        if (!bcCont.found("p")) {
            LFatal("The pressure must specified in face zone: %s", fz.name().c_str());
        }
        if (!bcCont.found("ma")) {
            LFatal("The Mach number  must specified in face zone: %s", fz.name().c_str());
        }
        p = bcCont.findType<real>("p", p);
        ma = bcCont.findType<real>("ma", ma);
        g = mixtures.thermalTable().gamma(p, T, yi);
        mufree = mixtures.transTable().mu(p, T, yi);
        a = sqrt(g * mixtures.species().Rm(yi) * T);
        vMag = ma * a;
        rho = mixtures.thermalTable().eos().rhom(p, T, yi);
        Re = rho * vMag / mufree;
        bcCont.add(std::string("v"), vMag);
        bcCont.add(std::string("rho"), rho);
        bcCont.add(std::string("Re"), Re);
    } else if (gbw == "pressureAndVMag") {
        if (!bcCont.found("p")) {
            LFatal("The pressure must specified in face zone: %s", fz.name().c_str());
        }
        if (!bcCont.found("v")) {
            LFatal("The velocity magnitude must specified in face zone: %s", fz.name().c_str());
        }
        p = bcCont.findType<real>("p", p);
        vMag = bcCont.findType<real>("v", vMag);
        g = mixtures.thermalTable().gamma(p, T, yi);
        mufree = mixtures.transTable().mu(p, T, yi);
        a = sqrt(g * mixtures.species().Rm(yi) * T);
        ma = vMag / a;
        rho = mixtures.thermalTable().eos().rhom(p, T, yi);
        Re = rho * vMag / mufree;
        bcCont.add(std::string("ma"), ma);
        bcCont.add(std::string("rho"), rho);
        bcCont.add(std::string("Re"), Re);
    } else if (gbw == "ReynoldAndMach") {
        if (!bcCont.found("Re")) {
            LFatal("The Reynold number must specified in face zone: %s", fz.name().c_str());
        }
        if (!bcCont.found("ma")) {
            LFatal("The Mach number must specified in face zone: %s", fz.name().c_str());
        }
        Re = bcCont.findType<real>("Re", Re);
        ma = bcCont.findType<real>("ma", ma);
        rho = Re * mufree / (ma * a);
        p = mixtures.thermalTable().eos().pm(rho, T, yi);

        real g0 = mixtures.thermalTable().gamma(p, T, yi);
        real mufree0 = mixtures.transTable().mu(p, T, yi);
        real a0 = sqrt(g0 * mixtures.species().Rm(yi) * T);
        integer count = 0;
        while ((fabs(g - g0) / g > 1e-4) || (fabs(mufree - mufree0) / mufree > 1e-4)) {
            g = g0;
            mufree = mufree0;
            a = a0;

            rho = Re * mufree / (ma * a);
            p = mixtures.thermalTable().eos().pm(rho, T, yi);
            g0 = mixtures.thermalTable().gamma(p, T, yi);
            mufree0 = mixtures.transTable().mu(p, T, yi);
            a0 = sqrt(g0 * mixtures.species().Rm(yi) * T);
            if (count++ > 1000) {
                break;
            }
        }
        g = g0;
        mufree = mufree0;
        a = a0;
        rho = Re * mufree / (ma * a);
        p = mixtures.thermalTable().eos().pm(rho, T, yi);
        vMag = ma * a;
        bcCont.add(std::string("p"), p);
        bcCont.add(std::string("v"), vMag);
        bcCont.add(std::string("rho"), rho);
    } else if (gbw == "ReynoldAndVMag") {
        if (!bcCont.found("Re")) {
            LFatal("The Reynold number must specified in face zone: %s", fz.name().c_str());
        }
        if (!bcCont.found("v")) {
            LFatal("The velocity magnitude must specified in face zone: %s", fz.name().c_str());
        }
        Re = bcCont.findType<real>("Re", Re);
        vMag = bcCont.findType<real>("v", vMag);
        rho = Re * mufree / vMag;
        ma = vMag / a;
        p = mixtures.thermalTable().eos().pm(rho, T, yi);
        real g0 = mixtures.thermalTable().gamma(p, T, yi);
        real mufree0 = mixtures.transTable().mu(p, T, yi);
        real a0 = sqrt(g0 * mixtures.species().Rm(yi) * T);
        integer count = 0;
        while ((fabs(g - g0) / g > 1e-4) || (fabs(mufree - mufree0) / mufree > 1e-4)) {
            g = g0;
            mufree = mufree0;
            a = a0;

            rho = Re * mufree / vMag;
            p = mixtures.thermalTable().eos().pm(rho, T, yi);
            g0 = mixtures.thermalTable().gamma(p, T, yi);
            mufree0 = mixtures.transTable().mu(p, T, yi);
            a0 = sqrt(g0 * mixtures.species().Rm(yi) * T);
            if (count++ > 1000) {
                break;
            }
        }
        g = g0;
        mufree = mufree0;
        a = a0;
        rho = Re * mufree / vMag;
        p = mixtures.thermalTable().eos().pm(rho, T, yi);
        ma = vMag / a;
        bcCont.add(std::string("p"), p);
        bcCont.add(std::string("ma"), ma);
        bcCont.add(std::string("rho"), rho);
    } else {
        LFatal("Unknown momentum specification method: %s in face zone: %s", gbw.c_str(),
               fz.name().c_str());
    }

    bcCont.add(std::string("gamma"), g);
    bcCont.add(std::string("mu"), mufree);
    real E = mixtures.thermalTable().ha_p(p, T, yi) + real(0.5) * vMag * vMag - p / rho;
    bcCont.add(std::string("E"), E);

    real Rgas = mixtures.species().Rm(yi);
    bcCont.add(std::string("Rgas"), Rgas);

    if (fz.bcType() != faceBCType::bcTypes::PRESSUREFARFIELD) {
        const_cast<faceZone &>(fz).setBcType(faceBCType::bcTypes::PRESSUREFARFIELD);
    }
    real k = 0;
    real epsilon = 0;
    if (bcCont.findWord("specificationMethod") == "intensityAndLength") {
        real intensity_ = 0;
        real lengthScale_ = 0;
        if (bcCont.found("intensity")) {
            real intensity_ = 0;
            intensity_ = bcCont.findType<real>("intensity", intensity_) * real(0.01);
            if (bcCont.found("length-scale")) {
                lengthScale_ = bcCont.findType<real>("length-scale", lengthScale_);
            } else {
                std::string errMsg;
                errMsg = "Turbulent parameters specification found in boundary "
                         "type: ";
                errMsg += bcCont.findWord("bcType");
                errMsg += " are not enough.";
                errMsg += "\nPlease check.";
                errorAbortStr(errMsg);
            }
        } else {
            std::string errMsg;
            errMsg = "No turbulent parameter specification found in boundary type: ";
            errMsg += bcCont.findWord("bcType");
            errMsg += "\nPlease check.";
            errorAbortStr(errMsg);
        }
        k = 1.5 * sqr(vMag * intensity_);
        epsilon = k * sqrt(k) / max(lengthScale_, veryTiny);
    } else if (bcCont.findWord("specificationMethod") == "intensityAndHydraulicD") {
        real intensity_ = 0;
        real hydraulicD_ = 0;
        if (bcCont.found("intensity")) {
            intensity_ = bcCont.findType<real>("intensity", intensity_) * real(0.01);
            if (bcCont.found("Hydraulic-Diameter")) {
                hydraulicD_ = bcCont.findType<real>("Hydraulic-Diameter", hydraulicD_);
            } else {
                std::string errMsg;
                errMsg = "Turbulent parameters specification found in boundary "
                         "type: ";
                errMsg += bcCont.findWord("bcType");
                errMsg += " are not enough.";
                errMsg += "\nPlease check.";
                errorAbortStr(errMsg);
            }
        } else {
            std::string errMsg;
            errMsg = "No turbulent parameter specification found in boundary type: ";
            errMsg += bcCont.findWord("bcType");
            errMsg += "\nPlease check.";
            errorAbortStr(errMsg);
        }
        const real cmu = 0.09;
        k = 1.5 * sqr(vMag * intensity_);
        real lengthScale_ = 0.07 * hydraulicD_ / pow(cmu, real(0.75));
        epsilon = k * sqrt(k) / max(lengthScale_, veryTiny);
    } else if (bcCont.findWord("specificationMethod") == "viscosityRatio") {
        real intensity_ = 0;
        real viscosityRatio_ = 0;
        if (bcCont.found("viscosity-ratio")) {
            viscosityRatio_ = bcCont.findType<real>("viscosity-ratio", viscosityRatio_);
            if (bcCont.found("intensity")) {
                intensity_ = bcCont.findType<real>("intensity", intensity_) * real(0.01);
            } else {
                std::string errMsg;
                errMsg = "Turbulent parameters specification found in boundary "
                         "type: ";
                errMsg += bcCont.findWord("bcType");
                errMsg += " are not enough.";
                errMsg += "\nPlease check.";
                errorAbortStr(errMsg);
            }
        } else {
            std::string errMsg;
            errMsg = "No turbulent parameter specification found in boundary type: ";
            errMsg += bcCont.findWord("bcType");
            errMsg += "\nPlease check.";
            errorAbortStr(errMsg);
        }

        const real cmu = 0.09;
        k = 1.5 * sqr(vMag * intensity_);
        epsilon = rho * cmu * sqr(k) / (mufree * viscosityRatio_);
    } else if (bcCont.findWord("specificationMethod") == "kAndEpsilon") {
        if (bcCont.found("k")) {
            k = bcCont.findType<real>("k", k);
            if (bcCont.found("epsilon")) {
                epsilon = bcCont.findType<real>("epsilon", epsilon);
            } else {
                std::string errMsg;
                errMsg = "Turbulent parameters specification found in boundary "
                         "type: ";
                errMsg += bcCont.findWord("bcType");
                errMsg += " are not enough.";
                errMsg += "\nPlease check.";
                errorAbortStr(errMsg);
            }
        } else {
            std::string errMsg;
            errMsg = "No turbulent parameter specification found in boundary type: ";
            errMsg += bcCont.findWord("bcType");
            errMsg += "\nPlease check.";
            errorAbortStr(errMsg);
        }
    } else if (bcCont.findWord("specificationMethod") == "kAndOmega") {
        if (bcCont.found("k")) {
            k = bcCont.findType<real>("k", k);
            if (bcCont.found("w")) {
                epsilon = bcCont.findType<real>("w", epsilon);
                epsilon = 0.09 * epsilon * k;
            } else {
                std::string errMsg;
                errMsg = "Turbulent parameters specification found in boundary "
                         "type: ";
                errMsg += bcCont.findWord("bcType");
                errMsg += " are not enough.";
                errMsg += "\nPlease check.";
                errorAbortStr(errMsg);
            }
        } else {
            std::string errMsg;
            errMsg = "No turbulent parameter specification found in boundary type: ";
            errMsg += bcCont.findWord("bcType");
            errMsg += "\nPlease check.";
            errorAbortStr(errMsg);
        }
    } else {
        k = 0;
        epsilon = 0;
        std::string errMsg;
        errMsg = "Turbulent parameters specification found in boundary type: ";
        errMsg += bcCont.findWord("bcType");
        errMsg += " are not enough.";
        errMsg += "\nPlease check.";
        errorAbortStr(errMsg);
    }

    bcCont.add(std::string("k"), k);
    bcCont.add(std::string("epsilon"), epsilon);
}
