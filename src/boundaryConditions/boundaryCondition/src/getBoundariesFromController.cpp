/*!
 * \file getBoundariesFromController.cpp
 * \brief Main subroutines for parsing boundary condition from controller.
 * \author Rao Sihang
 * \version V2.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of OpenHurricane project (Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation)
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

void OpenHurricane::getBoundariesFromController::setBoundariesController(controller &cont,
                                                                     const mixture &mixtures) {
    const auto &fZL = mixtures.mesh().faceZones();

    if (!cont.found("boundaryCondition")) {
        LFatal("Cannot find boundary condition in %s", cont.name().c_str());
    }
    controller &interCont = cont.subController("boundaryCondition");

    for (integer i = 0; i < fZL.size(); i++) {
        if (interCont.found(fZL[i].name())) {
            controller &bcCont = interCont.subController(fZL[i].name());
            const auto bcKey = bcCont.findWord("bcType");

            if (bcKey == "pressureFarField" || bcKey == "pressure-far-field") {
                getPressureFarField(mixtures, cont, bcCont,
                                    fZL[i]); // [1] Pressure far field
            } else if (bcKey == "supersonicInlet") {
                getSupersonicInlet(mixtures, cont, bcCont,
                                   fZL[i]); // [2] Supersonic inlet
            } else if (bcKey == "pressureInlet" || bcKey == "pressure-inlet") {
                getPressureInlet(mixtures, cont, bcCont,
                                 fZL[i]); // [3] Pressure inlet
            } else if (bcKey == "subsonicInlet") {
                getSubsonicInlet(mixtures, cont, bcCont,
                                 fZL[i]); // [4] subsonic inlet
            } else if (bcKey == "detonationInlet" || bcKey == "detonation-inlet") {
                getDetonationInlet(mixtures, cont, bcCont,
                                   fZL[i]); //   detonation inlet
            }
            //=====================================================
            //   Velocity inlet
            //=====================================================
            else if (bcKey == "velocityInlet" || bcKey == "velocity-inlet") {
                getVelocityInlet(mixtures, cont, bcCont, fZL[i]);
            } else if (bcKey == "massFlowInlet" || bcKey == "mass-flow-inlet") {
                getMassFlowInlet(mixtures, cont, bcCont,
                                 fZL[i]); //   Mass flow inlet
            }
            //=====================================================
            //   Mass flow outlet
            //=====================================================
            else if (bcKey == "massFlowOutlet" || bcKey == "mass-flow-outlet") {
                LFatal("The mass flow outlet is not yet supported in current program used in face "
                       "zone: %s",
                       fZL[i].name().c_str());
            } else if (bcKey == "pressureOutlet" || bcKey == "pressure-outlet") {
                getPressureOutlet(cont, bcCont, fZL[i]); //   Pressure outlet
            } else if (bcKey == "outflow" || bcKey == "out-flow" || bcKey == "subsonicOutlet" ||
                       bcKey == "supersonicOutlet") {
                getOutflow(cont, bcCont,
                           fZL[i]); //   Outflow, subsonic outlet or supersonic outlet.
            } else if (bcKey == "wall") {
                getWallCondition(mixtures, cont, bcCont, fZL[i]); //   Wall
            } else if (bcKey == "symmetry") {                     //   Symmetry
                if (fZL[i].bcType() != faceBCType::bcTypes::SYMMETRY) {
                    const_cast<faceZone &>(fZL[i]).setBcType(faceBCType::bcTypes::SYMMETRY);
                }
            } else if (bcKey == "periodic") { //   periodic
                if (!fZL[i].isPeriodic() && !fZL[i].isPeriodicShadow()) {
                    LFatal("The periodic codition must be specific in the mesh");
                }
                bcCont.add(std::string("bcType"), string("interior"), true);
            } else if (bcKey == "syntheticTurbInlet") {
                getSyntheticTurbInlet(mixtures, cont, bcCont, fZL[i]);
            } else {
                LFatal("Unknown boundary type: %s for face zone: %s", bcKey.c_str(),
                       fZL[i].name().c_str());
            }
        } else if (fZL[i].isInterior()) {
            controller addCont(fZL[i].name(), interCont);
            addCont.add(std::string("bcType"), string("interior"));
            interCont.add(fZL[i].name(), addCont);
        } else if (fZL[i].isCutFace()) {
            controller addCont(fZL[i].name(), interCont);
            addCont.add(std::string("bcType"), string("interior"));
            interCont.add(fZL[i].name(), addCont);
        } else if (fZL[i].isPeriodic() || fZL[i].isPeriodicShadow()) {
            controller addCont(fZL[i].name(), interCont);
            addCont.add(std::string("bcType"), string("interior"));
            interCont.add(fZL[i].name(), addCont);
        } else {
            LFatal("Boundary condition of faceZone[%d] named %s has not been specified.",i,
                   fZL[i].name().c_str());
        }
    }

    checkInitializationSetting(cont, mixtures.mesh());
}

OpenHurricane::realArray OpenHurricane::getBoundariesFromController::getSpeciesMassFractions(
    controller &cont, const mixture &mixtures, const bool addToCont) {
    const auto &species = mixtures.species();
    realArray yi(species.size(), Zero);
    if (species.size() == 1) {
        yi = 1.0;
    } else {
        if (cont.found("species")) {
            const controller &specCont = cont.subController("species");
            integer givenFlag = 0; // 0 = mass fraction; 1 = mole fraction.
            if (specCont.found("givenBy")) {
                const auto &flagW = specCont.findWord("givenBy");
                if (flagW == "massFraction") {
                    givenFlag = 0;
                } else if (flagW == "moleFraction") {
                    givenFlag = 1;
                } else {
                    LFatal("Unknown type: %s for given species in: %s", flagW.c_str(),
                           specCont.name().c_str());
                }
            }

            real yisum = Zero;
            for (integer i = 0; i < species.size() - 1; ++i) {
                std::string spn = species[i].name();
                yi[i] = specCont.findOrDefault<real>(spn, 0.0);
                if (yi[i] < 0) {
                    LFatal("Value must not lower than zero for species: %s in %s", spn.c_str(),
                           specCont.name().c_str());
                }
                yisum += yi[i];
            }
            if (yisum > 1.0) {
                LFatal("The summation of all species must not greater "
                       "than one in %s",
                       specCont.name().c_str());
            }

            yi[species.size() - 1] = 1.0 - yisum;

            if (givenFlag == 1) {
                species.Xi2Yi(yi, yi);
            }
        } else {
            yi = real(1) / real(species.size());
        }
    }

    if (addToCont) {
        for (integer i = 0; i < species.size(); ++i) {
            std::string spn = species[i].name();
            if (cont.found(spn)) {
                cont.remove(spn);
            }
            cont.add<real>(spn, yi[i]);
        }
    }
    return yi;
}

void OpenHurricane::getBoundariesFromController::getPressureFarField(const mixture &mixtures,
                                                                 const controller &cont,
                                                                 controller &bcCont,
                                                                 const faceZone &fz) {
    if (bcCont.found("bcType")) {
        bcCont.remove("bcType");
    }
    bcCont.add("bcType", string("pressureFarField"));

    if (!bcCont.found("defultType")) {
        bcCont.add(std::string("rhobcType"), string("Riemann"));
        bcCont.add(std::string("vbcType"), string("interior"));
        bcCont.add(std::string("pbcType"), string("interior"));
        bcCont.add(std::string("TbcType"), string("interior"));
        bcCont.add(std::string("defultType"), string("fixedValue"));
    }
    auto yi = getSpeciesMassFractions(bcCont, mixtures, true);

    vector direct;
    const auto directTypeW =
        parsingDirection::getDirection(bcCont, mixtures.mesh(), fz, direct, true);
    if (!bcCont.found("T")) {
        LFatal("The temperature must specified in face zone: %s", fz.name().c_str());
    }
    real T = bcCont.findType<real>("T", T);

    if (!bcCont.found("momentumGivenBy")) {
        LFatal("The momentum specification method must specified in face zone: %s",
               fz.name().c_str());
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
            LFatal("The Mach number must specified in face zone: %s", fz.name().c_str());
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
}

void OpenHurricane::getBoundariesFromController::getVelocityInlet(const mixture &mixtures,
                                                              const controller &cont,
                                                              controller &bcCont,
                                                              const faceZone &fz) {
    getSupersonicInlet(mixtures, cont, bcCont, fz);
    if (fz.bcType() != faceBCType::bcTypes::VELOCITYINLET) {
        const_cast<faceZone &>(fz).setBcType(faceBCType::bcTypes::VELOCITYINLET);
    }
}

void OpenHurricane::getBoundariesFromController::getSubsonicInlet(const mixture &mixtures,
                                                              const controller &cont,
                                                              controller &bcCont,
                                                              const faceZone &fz) {
    bcCont.add(std::string("rhobcType"), string("subsonicInlet"));
    bcCont.add(std::string("vbcType"), string("interior"));
    bcCont.add(std::string("pbcType"), string("interior"));
    bcCont.add(std::string("TbcType"), string("interior"));
    bcCont.add(std::string("defultType"), string("fixedValue"));

    realArray yi = getSpeciesMassFractions(bcCont, mixtures, true);

    if (!bcCont.found("totalPressure")) {
        LFatal("The total pressure must specified in face zone: %s", fz.name().c_str());
    }
    if (!bcCont.found("totalTemperature")) {
        LFatal("The total temperature must specified in face zone: %s", fz.name().c_str());
    }
    if (!bcCont.found("staticPressure")) {
        LFatal("The static pressure must specified in face zone: %s", fz.name().c_str());
    }

    vector direct;
    const auto directTypeW =
        parsingDirection::getDirection(bcCont, mixtures.mesh(), fz, direct, true);
    /*if (cont.subController("initialization").findWord("initFromBoundary") == fz.name())
    {*/
    real p0 = bcCont.findType<real>(string("totalPressure"), p0);
    real T0 = bcCont.findType<real>(string("totalTemperature"), T0);
    real p = bcCont.findType<real>(string("staticPressure"), p);

    real T = T0;
    real Ti = T0;
    for (integer i = 0; i < 20; ++i) {
        real g = mixtures.thermalTable().gamma(p, Ti, yi);
        T = T0 * pow(p0 / p, (real(1.0) - g) / g);
        if (mag(Ti - T) / max(T0, tiny) < 1e-4) {
            break;
        }
        Ti = T;
    }
    real cp = mixtures.thermalTable().cp0(T, yi);
    real v = sqrt(2 * cp * (T0 - T));
    real rho = mixtures.thermalTable().rho(p, T, yi);
    bcCont.add(std::string("rho"), rho);
    bcCont.add(std::string("v"), v);
    bcCont.add(std::string("p"), p);
    bcCont.add(std::string("T"), T);
    //real E = mixtures.thermalTable().ea_p(p0, T0, yi);
    real E = mixtures.thermalTable().ea_p(p, T, yi) + real(0.5) * v * v;
    bcCont.add(std::string("E"), E);

    real g = mixtures.thermalTable().gamma(p, T, yi);
    bcCont.add(std::string("gamma"), g);
    real mufree = mixtures.transTable().mu(p, T, yi);
    bcCont.add(std::string("mu"), mufree);

    real Rgas = mixtures.species().Rm(yi);
    bcCont.add(std::string("Rgas"), Rgas);
    //}
    if (fz.bcType() != faceBCType::bcTypes::PRESSUREINLET) {
        const_cast<faceZone &>(fz).setBcType(faceBCType::bcTypes::PRESSUREINLET);
    }
}

void OpenHurricane::getBoundariesFromController::getPressureInlet(const mixture &mixtures,
                                                              const controller &cont,
                                                              controller &bcCont,
                                                              const faceZone &fz) {
    bcCont.add(std::string("rhobcType"), string("pressureInlet"));
    bcCont.add(std::string("vbcType"), string("interior"));
    bcCont.add(std::string("pbcType"), string("interior"));
    bcCont.add(std::string("TbcType"), string("interior"));
    bcCont.add(std::string("defultType"), string("fixedValue"));

    realArray yi = getSpeciesMassFractions(bcCont, mixtures, true);

    if (!bcCont.found("totalPressure")) {
        LFatal("The total pressure must specified in face zone: %s", fz.name().c_str());
    }
    if (!bcCont.found("totalTemperature")) {
        LFatal("The total temperature must specified in face zone: %s", fz.name().c_str());
    }
    if (!bcCont.found("staticPressure")) {
        LFatal("The static pressure must specified in face zone: %s", fz.name().c_str());
    }

    vector direct;
    const auto directTypeW =
        parsingDirection::getDirection(bcCont, mixtures.mesh(), fz, direct, true);
    /*if (cont.subController("initialization").findWord("initFromBoundary") == fz.name())
    {*/
    real p0 = bcCont.findType<real>(string("totalPressure"), p0);
    real T0 = bcCont.findType<real>(string("totalTemperature"), T0);
    real p = bcCont.findType<real>(string("staticPressure"), p);
    real T = T0;
    real Ti = T0;
    for (integer i = 0; i < 20; ++i) {
        real g = mixtures.thermalTable().gamma(p, Ti, yi);
        T = T0 * pow(p0 / p, (real(1.0) - g) / g);
        if (mag(Ti - T) / max(T0, tiny) < 1e-4) {
            break;
        }
        Ti = T;
    }
    real cp = mixtures.thermalTable().cp0(T, yi);
    real v = sqrt(2 * cp * (T0 - T));
    real rho = mixtures.thermalTable().rho(p, T, yi);
    bcCont.add(std::string("rho"), rho);
    bcCont.add(std::string("v"), v);
    bcCont.add(std::string("p"), p);
    bcCont.add(std::string("T"), T);
    //real E = mixtures.thermalTable().ea_p(p0, T0, yi);
    real E = mixtures.thermalTable().ea_p(p, T, yi) + real(0.5) * v * v;
    bcCont.add(std::string("E"), E);

    real g = mixtures.thermalTable().gamma(p, T, yi);
    real mufree = mixtures.transTable().mu(p, T, yi);
    bcCont.add(std::string("gamma"), g);
    bcCont.add(std::string("mu"), mufree);

    real Rgas = mixtures.species().Rm(yi);
    bcCont.add(std::string("Rgas"), Rgas);

    /*}*/
    if (fz.bcType() != faceBCType::bcTypes::PRESSUREINLET) {
        const_cast<faceZone &>(fz).setBcType(faceBCType::bcTypes::PRESSUREINLET);
    }
}

void OpenHurricane::getBoundariesFromController::getMassFlowInlet(const mixture &mixtures,
                                                              const controller &cont,
                                                              controller &bcCont,
                                                              const faceZone &fz) {
    bcCont.add(std::string("rhobcType"), string("massFlowInlet"));
    bcCont.add(std::string("vbcType"), string("interior"));
    bcCont.add(std::string("pbcType"), string("interior"));
    bcCont.add(std::string("TbcType"), string("interior"));
    bcCont.add(std::string("defultType"), string("fixedValue"));
    const runtimeMesh &mesh = mixtures.mesh();
    if (!bcCont.found("givenBy")) {
        LFatal("The mass flow specification method must specified in face zone: %s",
               fz.name().c_str());
    }
    const auto massFlowTypeW = bcCont.findWord("givenBy");
    vector direct;
    const auto directTypeW = parsingDirection::getDirection(bcCont, mesh, fz, direct, true);
    vector normal;
    parsingDirection::getFaceNormal(mesh, fz, normal);
    if (massFlowTypeW == "massFlowRate") {
        if (bcCont.found("massFlowRate")) {
            real flux = bcCont.findType<real>("massFlowRate", flux);
            const vectorArray &fA = mesh.faceArea();
            real area = 0;
            realList recv(HurMPI::getProcSize(), Zero);
            for (integer fi = fz.firstIndex(); fi < fz.lastIndex() + 1; fi++) {
                area += fA[fi].magnitude();
            }
            HurMPI::allGather(&area, 1, feature<real>::MPIType, recv.data(), 1,
                              feature<real>::MPIType, HurMPI::getComm());
            HurMPI::barrier();
            area = 0;
            for (integer ip = 0; ip < recv.size(); ip++) {
                area += recv[ip];
            }
            if (directTypeW != "normalToBoundary") {
                area *= cos(direct, normal);
            }
            if (area == real(0)) {
                LFatal("The mass flux in face zone: %s is inf, please check!\n The direction is "
                       "set to be %s",
                       fz.name().c_str(), directTypeW.c_str());
            }
            flux /= area;
            if (bcCont.found("massFlux")) {
                bcCont.remove(std::string("massFlux"));
            }
            bcCont.add("massFlux", flux);
        } else {
            LFatal("Unknown mass flow rate unset in face zone: %s", fz.name().c_str());
        }
    } else if (massFlowTypeW == "massFlux") {
        if (!bcCont.found("massFlux")) {
            LFatal("Unknown mass flux unset in face zone: %s", fz.name().c_str());
        }
    } else {
        LFatal("Unknown mass flow specification method: %s in face zone: %s", massFlowTypeW.c_str(),
               fz.name().c_str());
    }

    if (fz.bcType() != faceBCType::bcTypes::MASSFLOWINLET) {
        const_cast<faceZone &>(fz).setBcType(faceBCType::bcTypes::MASSFLOWINLET);
    }
    realArray yi = getSpeciesMassFractions(bcCont, mixtures, true);
    if (!bcCont.found("totalTemperature")) {
        LFatal("The total temperature must specified in face zone: %s", fz.name().c_str());
    }
    real T0 = bcCont.findType<real>("totalTemperature", T0);
    if (!bcCont.found("p")) {
        LFatal("The pressure must specified in face zone: %s", fz.name().c_str());
    }
    real p = bcCont.findType<real>("p", p);

    real flux = bcCont.findType<real>("massFlux", flux);
    real Test = T0;
    /*real Tnew = 300;
    real Ttol = 300 * real(1.0e-4);*/
    real Tnew = 0.9 * T0;
    real Ttol = T0 * real(1.0e-4);
    real gnew;
    int iter = 0;
    do {
        real gold = mixtures.thermalTable().gamma(p, Test, yi);
        real rhoold = mixtures.thermalTable().rho(p, Test, yi);
        real value =
            T0 - Test * (real(1.0) + real(0.5) * (gold - 1) * flux * flux / (gold * p * rhoold));
        gnew = mixtures.thermalTable().gamma(p, Tnew, yi);
        real rhonew = mixtures.thermalTable().rho(p, Tnew, yi);
        real valuenew =
            T0 - Tnew * (real(1.0) + real(0.5) * (gnew - 1) * flux * flux / (gnew * p * rhonew));
        real ki = (Tnew - Test) / (valuenew - value);
        Test = Tnew;
        integer flag = 1;
        Tnew = mixtures.thermalTable().limit(Tnew - valuenew * ki);
        if (iter++ > 100) {
#ifdef HUR_DEBUG
            PLWarning("Maximum number of rhoThermos exceeded: %d", 100);
#endif // HUR_DEBUG
            break;
        }
    } while (mag(Tnew - Test) > Ttol);

    real T = Tnew;
    real p0 = pow(T0 / T, gnew / (gnew - 1)) * p;

    real h0 = mixtures.thermalTable().ha_p(p0, T0, yi);
    bcCont.add(std::string("totalEnthalpy"), h0);
    //if (cont.subController("initialization").findWord("initFromBoundary") == fz.name())
    //{
    //real flux = bcCont.findType<real>("massFlux", flux);
    real rho = mixtures.thermalTable().rho(p, T, yi);
    real v = flux / rho;
    real E = mixtures.thermalTable().ea_p(p, T, yi) + real(0.5) * v * v;
    bcCont.add(std::string("E"), E);
    bcCont.add(std::string("T"), T);
    bcCont.add(std::string("rho"), rho);
    bcCont.add(std::string("v"), v);
    //Pout << " flux = " << flux << " v = " << v << " rho= " << rho << " T = " << T << std::endl;
    //}

    real g = mixtures.thermalTable().gamma(p, T, yi);
    real mufree = mixtures.transTable().mu(p, T, yi);
    bcCont.add(std::string("gamma"), g);
    bcCont.add(std::string("mu"), mufree);

    real Rgas = mixtures.species().Rm(yi);
    bcCont.add(std::string("Rgas"), Rgas);
    if (mixtures.species().size() > 1) {
        for (integer isp = 0; isp < mixtures.species().size(); ++isp) {
            bcCont.add(mixtures.species()[isp].name() + "bcType",
                       string("fixedValueExtrapolate"));
        }
    }
}

void OpenHurricane::getBoundariesFromController::getDetonationInlet(const mixture &mixtures,
                                                                const controller &cont,
                                                                controller &bcCont,
                                                                const faceZone &fz) {
    bcCont.add(std::string("vbcType"), string("detonationInlet"));
    bcCont.add(std::string("rhobcType"), string("interior"));
    bcCont.add(std::string("pbcType"), string("interior"));
    bcCont.add(std::string("TbcType"), string("interior"));
    bcCont.add(std::string("defultType"), string("fixedValue"));

    const runtimeMesh &mesh = mixtures.mesh();

    vector direct;
    const auto directTypeW = parsingDirection::getDirection(bcCont, mesh, fz, direct, true);
    vector normal;
    parsingDirection::getFaceNormal(mesh, fz, normal);

    if (fz.bcType() != faceBCType::bcTypes::DETONATIONINLET) {
        const_cast<faceZone &>(fz).setBcType(faceBCType::bcTypes::DETONATIONINLET);
    }
    realArray yi = getSpeciesMassFractions(bcCont, mixtures, true);
    if (!bcCont.found("totalTemperature")) {
        LFatal("The total temperature must specified in face zone: %s", fz.name().c_str());
    }
    real T0 = bcCont.findType<real>("totalTemperature", T0);
    if (!bcCont.found("totalPressure")) {
        LFatal("The total pressure must specified in face zone: %s", fz.name().c_str());
    }
    real p0 = bcCont.findType<real>("totalPressure", p0);

    integer count = 0;
    real g0 = 1.4;
    real Ma = 1.0;
    real T = T0 / (1 + real(0.5) * (g0 - 1) * sqr(Ma));
    real p = p0 / pow(1 + real(0.5) * (g0 - 1) * sqr(Ma), g0 / (g0 - 1));
    real g = mixtures.thermalTable().gamma(p, T, yi);

    while (fabs(g - g0) / g0 > 1e-4) {
        g0 = g;
        T = T0 / (1 + real(0.5) * (g0 - 1) * sqr(Ma));
        p = p0 / pow(1 + real(0.5) * (g0 - 1) * sqr(Ma), g0 / (g0 - 1));
        g = mixtures.thermalTable().gamma(p, T, yi);

        if (count++ > 1000) {
            break;
        }
    }

    real rho = mixtures.thermalTable().rho(p, T, yi);
    real v = Ma * sqrt(g * p / rho);
    real E = mixtures.thermalTable().ea_p(p, T, yi) + real(0.5) * v * v;
    bcCont.add(std::string("E"), E);
    bcCont.add(std::string("T"), T);
    bcCont.add(std::string("rho"), rho);
    bcCont.add(std::string("v"), v);

    real mufree = mixtures.transTable().mu(p, T, yi);
    bcCont.add(std::string("gamma"), g);
    bcCont.add(std::string("mu"), mufree);

    real Rgas = mixtures.species().Rm(yi);
    bcCont.add(std::string("Rgas"), Rgas);
    if (mixtures.species().size() > 1) {
        for (integer isp = 0; isp < mixtures.species().size(); ++isp) {
            bcCont.add(mixtures.species()[isp].name() + "bcType",
                       string("fixedValue"));
        }
    }
}

void OpenHurricane::getBoundariesFromController::getOutflow(const controller &cont, controller &bcCont,
                                                        const faceZone &fz) {
    bcCont.remove("bcType");
    bcCont.add("bcType", string("outflow"));
    //controller& bcCont = interCont.subController(fZL[i].name());
    const auto bcKey = bcCont.findWord("bcType");
    if (bcCont.found("p")) {
        bcCont.add(std::string("defultType"), string("outflow"));
        bcCont.add(std::string("pbcType"), string("outflowPressure"));
    } else {
        if (bcKey == "subsonicOutlet") {
            LFatal("Pressure must be specified in the subsonic outlet.");
        }
    }
    if (fz.bcType() != faceBCType::bcTypes::OUTFLOW) {
        const_cast<faceZone &>(fz).setBcType(faceBCType::bcTypes::OUTFLOW);
    }
}

void OpenHurricane::getBoundariesFromController::getPressureOutlet(const controller &cont,
                                                               controller &bcCont,
                                                               const faceZone &fz) {
    if (!bcCont.found("p")) {
        LFatal("The pressure must specified in face zone: %s", fz.name().c_str());
    }
    //controller& bcCont = interCont.subController(fZL[i].name());
    bcCont.add(std::string("rhobcType"), string("pressureOutlet"));
    //bcCont.add(std::string("vbcType"), string("interior"));
    bcCont.add(std::string("pbcType"), string("interior"));
    //bcCont.add(std::string("TbcType"), string("interior"));
    bcCont.add(std::string("defultType"), string("zeroGradient"));
    /*bcCont.add(std::string("pbcType"), string("pressureOutlet"));
    bcCont.add(std::string("EbcType"), string("totalEnergyExtrapolate"));
    bcCont.add(std::string("defultType"), string("zeroGradient"));*/
    if (fz.bcType() != faceBCType::bcTypes::PRESSUREOUTLET) {
        const_cast<faceZone &>(fz).setBcType(faceBCType::bcTypes::PRESSUREOUTLET);
    }
}

void OpenHurricane::getBoundariesFromController::getWallCondition(const mixture &mixtures,
                                                              const controller &cont,
                                                              controller &bcCont,
                                                              const faceZone &fz) {
    if (bcCont.found("thermal")) {
        if (!bcCont.subController("thermal").found("thermalCondition")) {
            LFatal("The thermal condition must specified in face zone: %s", fz.name().c_str());
        }
        string thermoCondition = bcCont.subController("thermal").findWord("thermalCondition");
        if (thermoCondition == string("isothermal")) {
            if (!bcCont.subController("thermal").found("T")) {
                LFatal("The temperature must specified in face zone: %s", fz.name().c_str());
            }
            bcCont.add(std::string("T"),
                       bcCont.subController("thermal").findControlElePtr("T")->parameterContEle());
            bcCont.add(std::string("rhobcType"), string("interior"));
        } else if (thermoCondition == string("adiabatic")) {
            bcCont.add(std::string("rhobcType"), string("interior"));
        } else {
            LFatal("Ubknown thermal condition type: %s in face zone: %s", thermoCondition.c_str(),
                   fz.name().c_str());
        }
        bcCont.add(std::string("TbcType"), string(thermoCondition + "Wall"));
    } else {
        LFatal("The thermal condition must specified in face zone: %s", fz.name().c_str());
    }

    if (bcCont.found("momentum")) {
        if (cont.subController("flow").findWord("flowModel") == "EulerFlow") {
            bcCont.add(std::string("vbcType"), string("invSlipWall"));
        } else {
            if (!bcCont.subController("momentum").found("shearCondition")) {
                LFatal("The shear condition must specified in face zone: %s", fz.name().c_str());
            }
            string momentumCondition = bcCont.subController("momentum").findWord("shearCondition");
            if (momentumCondition == "noSlip") {
                // Doing nothing
            } else if (momentumCondition == "invSlip") {
                // Doing nothing
            } else {
                LFatal("Unknown shear condition type: %s in face zone: %s",
                       momentumCondition.c_str(), fz.name().c_str());
            }
            bcCont.add("vbcType", string(momentumCondition + "Wall"));
        }
    } else {
        LFatal("The momentum condition must specified in face zone: %s", fz.name().c_str());
    }

    if (bcCont.found("species")) {
        string speciesCondition = bcCont.subController("species").findWord("type");
        if (speciesCondition == "zeroGradient") {
            // Doing nothing
        } else if (speciesCondition == "specifiedSpecies") {
            const auto &species = mixtures.species();
            if (species.size() == 1) {
                LFatal("The specifiedSpecies in %s can only be used in multi-species systems",
                       fz.name().c_str());
            }
            getSpeciesMassFractions(bcCont, mixtures, true);
            for (integer isp = 0; isp < species.size(); ++isp) {
                bcCont.add(species[isp].name() + "bcType",
                           string("fixedValueExtrapolate"));
            }
        } else {
            LFatal("Unknown species condition type: %s in face zone: %s", speciesCondition.c_str(),
                   fz.name().c_str());
        }
    }

    bcCont.add("defultType", string("zeroGradient"));
    if (fz.bcType() != faceBCType::bcTypes::WALL) {
        const_cast<faceZone &>(fz).setBcType(faceBCType::bcTypes::WALL);
    }
}

void OpenHurricane::getBoundariesFromController::addBcTypeToController(const string &fzName,
                                                                   controller &interCont,
                                                                   std::string &bcType) {
    controller addCont(fzName, interCont);
    addCont.add(std::string("bcType"), string(bcType));
    interCont.add(fzName, addCont);
}

void OpenHurricane::getBoundariesFromController::checkInitializationSetting(const controller &cont,
                                                                        const runtimeMesh &mesh) {
    if (!cont.found("initialization")) {
        LFatal("Cannot find initialization controller in: %s", cont.name().c_str());
    }
    const auto &initCont = cont.subController("initialization");
    if (initCont.found("initFromBoundary")) {
        bool finded = false;
        const auto initW = initCont.findWord("initFromBoundary");
        for (integer i = 0; i < mesh.faceZones().size(); i++) {
            if (initW == mesh.faceZones()[i].name()) {
                if (mesh.faceZones()[i].isInterior()) {
                    LFatal("Cannot set initializing from %s because it is interior face. Please "
                           "check in %s",
                           initW.c_str(), initCont.name().c_str());
                }
                finded = true;
                break;
            }
        }
        if (!finded) {
            LFatal("Cannot find boundary \"%s\" for initialization in %s", initW.c_str(),
                   initCont.name().c_str());
        }
    }
}
