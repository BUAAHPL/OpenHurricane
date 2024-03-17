/*!
 * \file calculateFieldVar.cpp
 * \brief Main subroutines for calculateFieldVar.
 * \author Chen Zhenyi
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
#include "calculateFieldVar.hpp"

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcTotalPressure(const flowModel &flow) {
    if (flow.mixtures().isSingular()) {
        return flow.p() *
               pow(real(1.0) + real(0.5) * (flow.gama() - real(1.0)) *
                                   sqr(mag(flow.v()) / sqrt(flow.gama() * flow.p() / flow.rho())),
                   flow.gama() / (flow.gama() - real(1.0)));
    }
    realArray pt(flow.mesh().nTotalCells());
    const auto Tt = totalTemperature(flow, flow.p(), flow.T(), flow.v(), flow.mixtures().Yi());
    for (integer n = 0; n < flow.mesh().nCells(); ++n) {
        pt[n] = flow.p()[n] * exp(flow.mixtures().thermalTable().inteCp0dT(
                                      flow.T()[n], Tt[n], flow.mixtures().Yi(), n) /
                                  flow.mixtures().species().Rm(flow.mixtures().Yi(), n));
    }
    return pt;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcTotalEnthalpy(const flowModel &flow) {
    realArray hat(flow.mesh().nTotalCells());
    for (integer i = 0; i < flow.mesh().nCells(); i++) {
        integer flag = 0;
        real pi = flow.p()[i];
        real Ti = flow.T()[i];
        hat[i] = const_cast<flowModel &>(flow).thermo().mixtures().thermalTable().ha_p(
            pi, Ti, flow.mixtures().Yi(), i);
        hat[i] += real(0.5) * flow.v()[i].magSqr();
    }
    return hat;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFieldVar::cellVolume(const flowModel &flow) {
    realArray cv(flow.mesh().nTotalCells());
    for (integer i = 0; i < flow.mesh().nCells(); i++) {
        cv[i] = flow.mesh().cellVol()[i];
    }
    return cv;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcQCriterion(const flowModel &flow) {
    realArray Q(flow.mesh().nTotalCells());
    const auto &v = flow.v();
    for (integer i = 0; i < flow.mesh().nCells(); i++) {
        auto S = symm(v.grad()[i]);
        auto R = skew(v.grad()[i]);
        auto P = -tr(S);

        Q[i] = 0.5 * (sqr(P) + R.magSqr() - S.magSqr());
    }
    return Q;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcOmegaCriterion(const flowModel &flow) {
    const real ep = 0.001 * max(calcQCriterion(flow));
    realArray Omega(flow.mesh().nTotalCells());
    const auto &v = flow.v();
    for (integer i = 0; i < flow.mesh().nCells(); i++) {
        const auto wwij = skew(v.grad()[i]).magSqr();
        const auto SSij = symm(v.grad()[i]).magSqr();
        Omega[i] = wwij / (wwij + SSij + ep);
    }
    return Omega;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::calculateFieldVar::calcDeltaCriterion(const flowModel &flow) {
    realArray delta(flow.mesh().nTotalCells());
    const auto &v = flow.v();
    for (integer i = 0; i < flow.mesh().nCells(); i++) {
        const auto Q = 0.5 * (skew(v.grad()[i]).magSqr() - symm(v.grad()[i]).magSqr());
        const auto P = -tr(v.grad()[i]);
        const auto R = det(v.grad()[i]);
        const real q = (Q - sqr(P)) / 3.0;
        const real r = R + 2.0 / 27.0 * pow3(P) - (Q * P) / 3.0;
        delta[i] = pow3(q) / 27.0 + sqr(r) / 4.0;
    }
    return delta;
}
