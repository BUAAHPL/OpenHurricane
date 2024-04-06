#include "calculateFieldVar.hpp"
/*!
 * \file writeInstantFieldVar.inl
 * \brief In-Line subroutines of the <i>writeInstantFieldVar.hpp</i> file.
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
#pragma once

inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcMachNumber(const vectorArray &v, const realArray &gama,
                                                 const realArray &p, const realArray &rho) {
    return mag(v) / sqrt(gama * p / rho);
}

inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcTotalPressure(const realArray &p, const realArray &gama,
                                                    const realArray &Ma, const flowModel &flow) {
    return p *
           pow((real(1.0) + real(0.5) * (gama - real(1.0)) * sqr(Ma)), gama / (gama - real(1.0)));
}

inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcTotalPressure(const realArray &p, const realArray &gama,
                                                    const vectorArray &v, const realArray &rho,
                                                    const flowModel &flow) {
    return p * pow(real(1.0) + real(0.5) * (gama - real(1.0)) * sqr(mag(v) / sqrt(gama * p / rho)),
                   gama / (gama - real(1.0)));
}

inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcTotalTemperature(const realArray &T, const realArray &gama,
                                                       const realArray &Ma, const flowModel &flow) {
    return T * (real(1.0) + real(0.5) * (gama - real(1.0)) * sqr(Ma));
}

inline OpenHurricane::realArray OpenHurricane::calculateFieldVar::calcTotalTemperature(
    const realArray &T, const realArray &gama, const vectorArray &v, const realArray &p,
    const realArray &rho, const flowModel &flow) {
    return T * (real(1.0) + real(0.5) * (gama - real(1.0)) * sqr(mag(v) / sqrt(gama * p / rho)));
}

hur_nodiscard inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::totalTemperature(const flowModel &flow, const realArray &p,
                                                   const realArray &T, const vectorArray &v,
                                                   const PtrList<cellRealArray> &yi) {
    realArray Tt(flow.mesh().nTotalCells());
    for (integer i = 0; i < flow.mesh().nCells(); i++) {
        integer flag = 0;
        real pi = p[i];
        real Ti = T[i];
        real h0 =
            const_cast<flowModel &>(flow).thermo().mixtures().thermalTable().ha_p(pi, Ti, yi, i);
        h0 += real(0.5) * v[i].magSqr();
        Tt[i] = const_cast<flowModel &>(flow).thermo().mixtures().thermalTable().THa_p(h0, pi, Ti,
                                                                                       flag, yi, i);
    }
    return Tt;
}

hur_nodiscard inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcViscousRatio(const flowModel &flow) {
    if (NullRefObj::isNullRef(flow.mut())) {
        LFatal("Can not compute turbulent viscous ratio due to a null "
               "array of turbulent viscocity");
    }
    return flow.mut() / flow.mul();
}

inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcHeatReleaseRate(const flowModel &flow,
                                                      const combustionModel *chemtryPtr) {
    if (!flow.mixtures().noReaction()) {
        return const_cast<combustionModel *>(chemtryPtr)->heatReleaseRate();
    } else {
        return realArray();
    }
}

inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcGFO(const flowModel &flow,
                                          const combustionModel *chemtryPtr) {
    if (!flow.mixtures().noReaction()) {
        if (const_cast<combustionModel *>(chemtryPtr)->fuelName().size() == 0) {
            checkWarning("The fuel is not given for compute flame index and would not compute");
            return realArray();
        }
        if (const_cast<combustionModel *>(chemtryPtr)->oxygenName().size() == 0) {
            checkWarning("The oxygen is not given for compute flame index and would not compute");
            return realArray();
        }
        return flow.GFO(const_cast<combustionModel *>(chemtryPtr)->fuelName(),
                        const_cast<combustionModel *>(chemtryPtr)->oxygenName());
    } else {
        return realArray();
    }
}

inline void OpenHurricane::calculateFieldVar::calcMoleSpecies(
    const flowModel &flow, std::map<std::string, object *> &outFieldVarMap, const string &type) {
    if (flow.mixtures().isSingular()) {
        return;
    }
    bool isOutXi = false;
    for (integer i = 0; i < flow.mixtures().species().size(); ++i) {
        auto iter = outFieldVarMap.find(flow.mixtures().species()[i].name() + type);
        if (iter != outFieldVarMap.end()) {
            isOutXi = true;
            break;
        }
    }
    if (!isOutXi) {
        return;
    }
    PtrList<cellRealArray> tmpXi;
    for (integer i = 0; i < flow.mixtures().species().size(); ++i) {
        tmpXi.append(new cellRealArray(
            object(flow.mixtures().species()[i].name() + type, flow.mesh(), object::WRITE_OUTPUT),
            flow.mesh()));
    }
    flow.mixtures().species().Yi2Xi(flow.mixtures().Yi(), tmpXi);
    for (integer i = 0; i < flow.mixtures().species().size(); ++i) {
        auto iter = outFieldVarMap.find(flow.mixtures().species()[i].name() + type);
        if (iter != outFieldVarMap.end()) {
            iter->second = tmpXi.set(i, nullptr);
        }
    }
}

inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcVorticity(const flowModel &flow) {
    auto vor = skewMagnitude(flow.v().grad());

    return vor;
}

inline OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcDamkohler(const flowModel &flow,
                                                const combustionModel *chemtryPtr) {
    if (!flow.mixtures().noReaction()) {
        return const_cast<combustionModel *>(chemtryPtr)->calcDamkohler();
    } else {
        return realArray();
    }
}
