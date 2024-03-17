/*!
 * \file BS23.cpp
 * \brief Main subroutines for Bogacki-Shampine 2/3 ODEs embedded pair solver.
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

#include "BS23.hpp"

namespace OpenHurricane {
    createClassNameStr(BS23, "BS23");
    registerObjFty(ODEsSolver, BS23, controller);
} // namespace OpenHurricane

OpenHurricane::BS23::BS23(const integer nEqs)
    : ODEsSolver(nEqs), s2_(nEqs, Zero), s3_(nEqs, Zero), s4_(nEqs, Zero), error_(nEqs, Zero) {}

OpenHurricane::BS23::BS23(const integer nEqs, const real atol, const real rtol, const integer maxStep)
    : ODEsSolver(nEqs, atol, rtol, maxStep), s2_(nEqs, Zero), s3_(nEqs, Zero), s4_(nEqs, Zero),
      error_(nEqs, Zero) {}

OpenHurricane::BS23::BS23(const integer nEqs, const controller &cont)
    : ODEsSolver(nEqs, cont), s2_(nEqs, Zero), s3_(nEqs, Zero), s4_(nEqs, Zero),
      error_(nEqs, Zero) {}

bool OpenHurricane::BS23::reset(const integer nEqsNew) {
    if (ODEsSolver::reset(nEqsNew)) {
        resetArray(s2_);
        resetArray(s3_);
        resetArray(s4_);
        resetArray(error_);
        return true;
    } else {
        return false;
    }
}

OpenHurricane::real OpenHurricane::BS23::solve(const real t0, const realArray &y0, const realArray &dydt0,
                                       const real dt, realArray &y) {
    auto &s1 = dydt0;

    realArray yTmp;
    yTmp = y0 + real(0.5) * dt * s1;

    DyDt(t0 + real(0.5) * dt, yTmp, s2_);

    yTmp = y0 + real(0.75) * dt * s2_;
    DyDt(t0 + real(0.75) * dt, yTmp, s3_);

    auto &error = error_;

    y = y0 + dt / real(9.0) * (real(2.0) * s1 + real(4.0) * s3_ + real(3.0) * s2_);
    DyDt(t0 + dt, y, s4_);

    error = mag(dt / real(72) *
                (real(-5.0) * s1 + real(6.0) * s2_ + real(8.0) * s3_ + real(-9.0) * s4_));

    return maxError(y0, y, error);
}

void OpenHurricane::BS23::solveOnly(const real t0, const realArray &y0, const realArray &dydt0,
                                const real dt, realArray &y) {
    auto &s1 = dydt0;
    ;
    realArray yTmp;
    yTmp = y0 + real(0.5) * dt * s1;

    DyDt(t0 + real(0.5) * dt, yTmp, s2_);

    yTmp = y0 + real(0.75) * dt * s2_;
    DyDt(t0 + real(0.75) * dt, yTmp, s3_);

    y = y0 + dt / real(9.0) * (real(2.0) * s1 + real(4.0) * s3_ + real(3.0) * s2_);
}

OpenHurricane::real OpenHurricane::BS23::solve(const real t0, const realArray &y0, const realArray &dydt0,
                                       const real dt, const integerListList &Gml, const integer im,
                                       realArray &y) {
    if (im >= Gml.size()) {
        errorAbortStr(("The group index (im): " + toString(im) + " is larger than the group size " +
                       toString(Gml.size())));
    }
    LFatal("This functions cannot be used!");
    return real();
}
