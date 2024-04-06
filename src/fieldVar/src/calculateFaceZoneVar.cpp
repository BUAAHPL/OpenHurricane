/*!
 * \file calculateFaceZoneVar.cpp
 * \brief Main subroutines for calculating Face-Zone-based Variables.
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

#include "calculateFaceZoneVar.hpp"

void OpenHurricane::calculateFaceZoneVar::addToFuncMap(
    faceZoneFieldParameterFuncMap &funcMap) const {
    funcMap.addFunc(setKeyName("Ma"), calcFaceZoneMachNumber);
    funcMap.addFunc(setKeyName("pt"), calcFaceZoneTotalPressure);
    funcMap.addFunc(setKeyName("pt_constCp"), calcFaceZoneTotalPressureCCP);
    funcMap.addFunc(setKeyName("Tt"), calcFaceZoneTotalTemperature);
    funcMap.addFunc(setKeyName("Tt_constCp"), calcFaceZoneTotalTemperatureCCp);

    funcMap.addFunc(setKeyName("Qw"), wallFaceZoneHeatFlux);
    funcMap.addFunc(setKeyName("cf"), wallFaceZoneFrictionCoefficient);
    funcMap.addFunc(setKeyName("cpAbs"), wallFaceZoneAbsPressCoefficient);
    funcMap.addFunc(setKeyName("cpRel"), wallFaceZoneRelPressCoefficient);
    funcMap.addFunc(setKeyName("uPlus"), wallFaceZoneUPlus);
    funcMap.addFunc(setKeyName("YPlus"), wallFaceZoneYPlus);
    funcMap.addFunc(setKeyName("ch"), wallFaceZoneHeatCoefficient);
    funcMap.addFunc(setKeyName("vorticity"), calcFaceZoneVorticity);
    funcMap.addFunc(setKeyName("viscousRatio"), calcViscousRatio);
    funcMap.addFunc(setKeyName("u0"), calcU0);
    funcMap.addFunc(setKeyName("v0"), calcV0);
    funcMap.addFunc(setKeyName("w0"), calcW0);
}

OpenHurricane::string OpenHurricane::calculateFaceZoneVar::setKeyName(const string &key) {
    return key + basicName();
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcFaceZoneMachNumber(const flowModel &flow,
                                                            const integer fzi) {
    const auto vf = fv::interpolate(flow.v(), fzi);
    const auto gamaf = fv::interpolate(flow.gama(), fzi);
    const auto pf = fv::interpolate(flow.p(), fzi);
    const auto rhof = fv::interpolate(flow.rho(), fzi);
    return mag(vf) / sqrt(gamaf * pf / rhof);
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcFaceZoneTotalPressure(const flowModel &flow,
                                                               const integer fzi) {
    if (flow.mixtures().isSingular()) {
        const auto vf = fv::interpolate(flow.v(), fzi);
        const auto gamaf = fv::interpolate(flow.gama(), fzi);
        const auto pf = fv::interpolate(flow.p(), fzi);
        const auto rhof = fv::interpolate(flow.rho(), fzi);
        realArray ma(mag(vf) / sqrt(gamaf * pf / rhof));

        return pf * pow((real(1.0) + real(0.5) * (gamaf - real(1.0)) * sqr(ma)),
                        gamaf / (gamaf - real(1.0)));
    }
    const auto Ttf = calcFaceZoneTotalTemperature(flow, fzi);
    const auto Tf = fv::interpolate(flow.T(), fzi);
    const auto pf = fv::interpolate(flow.p(), fzi);
    auto yif = flow.mixtures().Yif(fzi);

    realArray ptf(pf.size(), Zero);
    for (integer i = 0; i < pf.size(); ++i) {
        ptf[i] = pf[i] * exp(flow.mixtures().thermalTable().inteCp0dT(Tf[i], Ttf[i], yif[i]) /
                             flow.mixtures().species().Rm(yif[i]));
    }

    return ptf;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcFaceZoneTotalPressureCCP(const flowModel &flow,
                                                                  const integer fzi) {
    const auto vf = fv::interpolate(flow.v(), fzi);
    const auto gamaf = fv::interpolate(flow.gama(), fzi);
    const auto pf = fv::interpolate(flow.p(), fzi);
    const auto rhof = fv::interpolate(flow.rho(), fzi);
    realArray ma(mag(vf) / sqrt(gamaf * pf / rhof));

    return pf * pow((real(1.0) + real(0.5) * (gamaf - real(1.0)) * sqr(ma)),
                    gamaf / (gamaf - real(1.0)));
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcFaceZoneTotalTemperature(const flowModel &flow,
                                                                  const integer fzi) {
    const auto Tf = fv::interpolate(flow.T(), fzi);
    const auto pf = fv::interpolate(flow.p(), fzi);
    const auto vf = fv::interpolate(flow.v(), fzi);

    if (flow.mixtures().isSingular()) {
        const auto gamaf = fv::interpolate(flow.gama(), fzi);
        const auto rhof = fv::interpolate(flow.rho(), fzi);
        realArray ma(mag(vf) / sqrt(gamaf * pf / rhof));

        return Tf * (real(1.0) + real(0.5) * (gamaf - real(1.0)) * sqr(ma));
    } else {
        realArray Ttf(pf.size(), Zero);
        auto yif = flow.mixtures().Yif(fzi);

        for (integer i = 0; i < pf.size(); ++i) {
            integer flag = 0;
            real pi = pf[i];
            real Ti = Tf[i];
            real h0 = const_cast<flowModel &>(flow).thermo().mixtures().thermalTable().ha_p(pi, Ti,
                                                                                            yif[i]);
            h0 += real(0.5) * vf[i].magSqr();
            Ttf[i] = const_cast<flowModel &>(flow).thermo().mixtures().thermalTable().THa_p(
                h0, pi, Ti, flag, yif[i]);
        }

        return Ttf;
    }
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcFaceZoneTotalTemperatureCCp(const flowModel &flow,
                                                                     const integer fzi) {
    const auto Tf = fv::interpolate(flow.T(), fzi);
    const auto pf = fv::interpolate(flow.p(), fzi);
    const auto vf = fv::interpolate(flow.v(), fzi);
    const auto gamaf = fv::interpolate(flow.gama(), fzi);
    const auto rhof = fv::interpolate(flow.rho(), fzi);
    realArray ma(mag(vf) / sqrt(gamaf * pf / rhof));

    return Tf * (real(1.0) + real(0.5) * (gamaf - real(1.0)) * sqr(ma));
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::wallFaceZoneHeatFlux(const flowModel &flow,
                                                          const integer zoneId) {
    const auto &fZL = flow.mesh().faceZones();
    const auto &faces = flow.mesh().faces();
    const auto &fA = flow.mesh().faceArea();
    const auto &fC = flow.mesh().faceCentre();
    const auto &cC = flow.mesh().cellCentre();
    if (!fZL[zoneId].isWall()) {
        checkWarning("Only for wall face zone");
        return realArray();
    }
    if (flow.kappal().size() == 0) {
        checkWarning("Cannot compute heat flux for not giving thermal conductivity");
        return realArray();
    }

    realArray heatFlux(fZL[zoneId].size());

    integer id = 0;
    for (integer facei = fZL[zoneId].firstIndex(); facei < fZL[zoneId].lastIndex() + 1; ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();
        const real Tw = 0.5 * (flow.T()[cl] + flow.T()[cr]);
        const real dd = (cC[cl] - fC[facei]) * fA[facei].normalized();
        const real kappaf = flow.kappal()(cl); // +cp()[cl] * mut()[cl] / prt();
        heatFlux[id] = -kappaf * (flow.T()[cl] - Tw) / dd;
        id++;
    }

    return heatFlux;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::wallFaceZoneFrictionCoefficient(const flowModel &flow,
                                                                     const integer zoneId) {
    const auto &fz = flow.mesh().faceZones()[zoneId];
    const auto &fA = flow.mesh().fA();
    const auto &fC = flow.mesh().faceCentre();
    const auto &cC = flow.mesh().cellCentre();
    const auto &faces = flow.mesh().faces();
    const auto &fS = flow.mesh().Iteration().refValues();
    if (!fz.isWall()) {
        checkWarning("Only for wall face zone");
        return realArray();
    }
    realArray cf(fz.size());
    integer id = 0;
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();
        const real muw = (flow.mul()[cl] + flow.mul()[cr]) * 0.5;
        const real dist = (cC[cl] - fC[facei]) * fA[facei].normalized();
        real uu = flow.v()[cl] * fA[facei].normalized();
        real ut = sqrt(flow.v()[cl].magSqr() - sqr(uu));
        real tauw = muw * ut / dist;
        cf[id] = tauw / (0.5 * fS.rho() * sqr(fS.vMag()));

        id++;
    }
    return cf;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::wallFaceZoneAbsPressCoefficient(const flowModel &flow,
                                                                     const integer zoneId) {
    const auto &fz = flow.mesh().faceZones()[zoneId];
    const auto &faces = flow.mesh().faces();
    const auto &fS = flow.mesh().Iteration().refValues();
    if (!fz.isWall()) {
        checkWarning("Only for wall face zone");
        return realArray();
    }
    realArray cpAbs(fz.size());
    integer id = 0;
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();
        real pw = (flow.p()[cl] + flow.p()[cr]) * 0.5;
        cpAbs[id] = pw / (0.5 * fS.rho() * sqr(fS.vMag()));

        id++;
    }
    return cpAbs;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::wallFaceZoneRelPressCoefficient(const flowModel &flow,
                                                                     const integer zoneId) {
    const auto &fz = flow.mesh().faceZones()[zoneId];
    const auto &faces = flow.mesh().faces();
    const auto &fS = flow.mesh().Iteration().refValues();
    if (!fz.isWall()) {
        checkWarning("Only for wall face zone");
        return realArray();
    }
    realArray cpRel(fz.size());
    integer id = 0;
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();
        real pw = (flow.p()[cl] + flow.p()[cr]) * 0.5;
        cpRel[id] = (pw - fS.p()) / (0.5 * fS.rho() * sqr(fS.vMag()));
        id++;
    }
    return cpRel;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::wallFaceZoneUPlus(const flowModel &flow,
                                                       const integer zoneId) {
    const auto &mesh = flow.mesh();
    const auto &fz = mesh.faceZones()[zoneId];
    const auto &fA = mesh.faceArea();
    const auto &faces = mesh.faces();

    realArray muw = flow.mul(zoneId);
    realArray upw(fz.size());

    const auto &fc = mesh.faceCentre();
    const auto &cc = mesh.cellCentre();
    integer id = 0;
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer &cl = faces[facei].leftCell();
        const integer &cr = faces[facei].rightCell();
        const auto nn = fA[facei].normalized();
        const auto dist = mag((cc[cl] - fc[facei]) * nn);
        real uu = flow.v()[cl] * nn;
        const real rhow = 0.5 * (flow.rho()[cl] + flow.rho()[cr]);
        real ut = sqrt(flow.v()[cl].magSqr() - sqr(uu));
        real tauw = muw[id] * ut / dist;
        real utau = sqrt(tauw / rhow);
        upw[id++] = ut / utau;
    }

    return upw;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::wallFaceZoneYPlus(const flowModel &flow,
                                                       const integer zoneId) {
    const auto &mesh = flow.mesh();
    const auto &fz = mesh.faceZones()[zoneId];
    const auto &fA = mesh.faceArea();
    const auto &faces = mesh.faces();
    realArray ypw(fz.size(), Zero);
    realArray muw = flow.mul(zoneId);

    const auto &fc = mesh.faceCentre();
    const auto &cc = mesh.cellCentre();
    integer id = 0;
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer &cl = faces[facei].leftCell();
        const auto &cr = faces[facei].rightCell();

        const auto nn = fA[facei].normalized();
        const auto dist = mag((cc[cl] - fc[facei]) * nn);
        real uu = flow.v()[cl] * nn;

        real rhow = real(0.5) * (flow.rho()[cl] + flow.rho()[cr]);
        real ut = sqrt(flow.v()[cl].magSqr() - sqr(uu));
        real tauw = muw[id] * ut / (dist);
        real utau = sqrt(tauw / rhow);
        ypw[id] = dist * utau / (muw[id] / rhow);
        id++;
    }

    return ypw;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::wallFaceZoneHeatCoefficient(const flowModel &flow,
                                                                 const integer zoneId) {
    const auto &mesh = flow.mesh();
    const auto &fz = mesh.faceZones()[zoneId];
    const auto &fA = mesh.faceArea();
    const auto &faces = mesh.faces();

    const auto &fS = flow.mesh().Iteration().refValues();

    realArray muw = flow.mul(zoneId);
    realArray ch(fz.size());

    const auto &fC = mesh.faceCentre();
    const auto &cC = mesh.cellCentre();
    integer id = 0;
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();
        const real Tw = 0.5 * (flow.T()[cl] + flow.T()[cr]);
        const real dd = (cC[cl] - fC[facei]) * fA[facei].normalized();

        ch[id] = muw[id] / (flow.Prl() * fS.rho() * fS.vMag() * (fS.Tt() - Tw)) *
                 (flow.T()[cl] - Tw) / dd;
        id++;
    }

    return ch;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcFaceZoneVorticity(const flowModel &flow,
                                                           const integer zoneId) {
    const auto &mesh = flow.mesh();
    const auto &fz = mesh.faceZones()[zoneId];
    const auto &faces = mesh.faces();
    const auto &fW = mesh.faceWeight();
    integer id = 0;
    realArray vor(fz.size(), Zero);
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();
        const auto vgf =
            fW[facei] * flow.v().grad()[cl] + (real(1.0) - fW[facei]) * flow.v().grad()[cr];
        vor[id] = skewMagnitude(vgf);
        id++;
    }
    return vor;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcViscousRatio(const flowModel &flow, const integer zoneId) {
    if (NullRefObj::isNullRef(flow.mut())) {
        LFatal("Can not compute turbulent viscous ratio due to a null "
               "array of turbulent viscocity");
    }
    const auto &mesh = flow.mesh();
    const auto &fz = mesh.faceZones()[zoneId];
    const auto &faces = mesh.faces();
    const auto &fW = mesh.faceWeight();
    integer id = 0;
    realArray mutf(fz.size(), Zero);
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();
        const auto mul = fW[facei] * flow.mul()[cl] + (real(1.0) - fW[facei]) * flow.mul()[cr];
        const auto mut = fW[facei] * flow.mut()[cl] + (real(1.0) - fW[facei]) * flow.mut()[cr];
        mutf[id] = mut / max(mul, veryTiny);
        id++;
    }
    return mutf;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcU0(const flowModel &flow, const integer zoneId) {
    const auto &mesh = flow.mesh();
    const auto &fz = mesh.faceZones()[zoneId];
    const auto &faces = mesh.faces();
    integer id = 0;
    realArray u0(fz.size(), Zero);
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        u0[id] = flow.v()[cl].x();
        id++;
    }
    return u0;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcV0(const flowModel &flow, const integer zoneId) {
    const auto &mesh = flow.mesh();
    const auto &fz = mesh.faceZones()[zoneId];
    const auto &faces = mesh.faces();
    integer id = 0;
    realArray v0(fz.size(), Zero);
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        v0[id] = flow.v()[cl].y();
        id++;
    }
    return v0;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFaceZoneVar::calcW0(const flowModel &flow, const integer zoneId) {
    const auto &mesh = flow.mesh();
    const auto &fz = mesh.faceZones()[zoneId];
    const auto &faces = mesh.faces();
    integer id = 0;
    realArray w0(fz.size(), Zero);
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        w0[id] = flow.v()[cl].z();
        id++;
    }
    return w0;
}