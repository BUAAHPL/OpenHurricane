/*!
 * \file getSupersonicInlet.cpp
 * \brief Main subroutines for parsing supersonic inlet boundary condition from controller.
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
#include "getBoundariesFromController.hpp"
#include "rhoThermo.hpp"

void OpenHurricane::getBoundariesFromController::getSupersonicInlet(const mixture &mixtures,
                                                                    const controller &cont,
                                                                    controller &bcCont,
                                                                    const faceZone &fz) {
    bcCont.add(std::string("defultType"), string("fixedValue"));
    bcCont.add(std::string("vbcType"), string("fixedValue"));
    bcCont.add(std::string("rhobcType"), string("fixedValue"));
    bcCont.add(std::string("pbcType"), string("fixedValue"));
    bcCont.add(std::string("TbcType"), string("fixedValue"));
    auto yi = getSpeciesMassFractions(bcCont, mixtures, true);
    vector direct;
    const auto directTypeW =
        parsingDirection::getDirection(bcCont, mixtures.mesh(), fz, direct, true);
    if (!bcCont.found("T")) {
        LFatal("The temperature must specified in face zone: %s", fz.name().c_str());
    }
    real T = bcCont.findType<real>("T", T);

    if (!bcCont.found("momentumGivenBy")) {
        LFatal("The momentum must specified in face zone: %s", fz.name().c_str());
    }
    const auto gbw = bcCont.findWord("momentumGivenBy");
    real p = constant::physicalConstant::Patm;

    real g = mixtures.thermalTable().gamma(p, T, yi);
    real mufree = 0;
    if (!mixtures.inviscous()) {
        mufree = mixtures.transTable().mu(p, T, yi);
    }

    real a = sqrt(g * mixtures.species().Rm(yi) * T);

    real ma;
    real vMag = 0;
    real Re = 0;
    real rho = Zero;
    if (gbw == "pressureAndMach") {
        if (!bcCont.found("p")) {
            LFatal("The pressure must specified in face zone: %s", fz.name().c_str());
        }
        if (!bcCont.found("ma")) {
            LFatal("The Mach number must specified in face zone: %s", fz.name().c_str());
        }
        p = bcCont.findType<real>("p", p);
        ma = bcCont.findType<real>("ma", ma);
        g = mixtures.thermalTable().gamma(p, T, yi);

        a = sqrt(g * mixtures.species().Rm(yi) * T);
        vMag = ma * a;
        rho = mixtures.thermalTable().eos().rhom(p, T, yi);

        if (!mixtures.inviscous()) {
            mufree = mixtures.transTable().mu(p, T, yi);
            Re = rho * vMag / mufree;
        }

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

        a = sqrt(g * mixtures.species().Rm(yi) * T);
        ma = vMag / a;
        rho = mixtures.thermalTable().eos().rhom(p, T, yi);

        if (!mixtures.inviscous()) {
            mufree = mixtures.transTable().mu(p, T, yi);
            Re = rho * vMag / mufree;
        }

        bcCont.add(std::string("ma"), ma);
        bcCont.add(std::string("rho"), rho);
        bcCont.add(std::string("Re"), Re);
    } else if (gbw == "ReynoldAndMach") {
        if (mixtures.inviscous()) {
            LFatal("The option: \"ReynoldAndMach\" cannot be used in inviscous mixtures in face "
                   "zone: %s",
                   fz.name().c_str());
        }
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
        if (mixtures.inviscous()) {
            LFatal("The option: \"ReynoldAndVMag\" cannot be used in inviscous mixtures in face "
                   "zone: %s",
                   fz.name().c_str());
        }
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
}
