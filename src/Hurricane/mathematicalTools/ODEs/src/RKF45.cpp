/*!
 * \file RKF45.cpp
 * \brief Main subroutines for 4/5th Order Runge-Kutta-Fehlberg ODE solver.
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
#include "RKF45.hpp"

namespace OpenHurricane {
    createClassNameStr(RKF45, "RKF45");
    registerObjFty(ODEsSolver, RKF45, controller);
} // namespace OpenHurricane

OpenHurricane::RKF45::RKF45(const integer nEqs)
    : ODEsSolver(nEqs), s2_(nEqs, Zero), s3_(nEqs, Zero), s4_(nEqs, Zero), s5_(nEqs, Zero),
      s6_(nEqs, Zero), error_(nEqs, Zero), yTmp_(nEqs, Zero) {}

OpenHurricane::RKF45::RKF45(const integer nEqs, const real atol, const real rtol, const integer maxStep)
    : ODEsSolver(nEqs, atol, rtol, maxStep), s2_(nEqs, Zero), s3_(nEqs, Zero), s4_(nEqs, Zero),
      s5_(nEqs, Zero), s6_(nEqs, Zero), error_(nEqs, Zero), yTmp_(nEqs, Zero) {}

OpenHurricane::RKF45::RKF45(const integer nEqs, const controller &cont)
    : ODEsSolver(nEqs, cont), s2_(nEqs, Zero), s3_(nEqs, Zero), s4_(nEqs, Zero), s5_(nEqs, Zero),
      s6_(nEqs, Zero), error_(nEqs, Zero), yTmp_(nEqs, Zero) {}

bool OpenHurricane::RKF45::reset(const integer nEqsNew) {
    if (ODEsSolver::reset(nEqsNew)) {
        resetArray(s2_);
        resetArray(s3_);
        resetArray(s4_);
        resetArray(s5_);
        resetArray(s6_);
        resetArray(error_);
        resetArray(yTmp_);
        return true;
    } else {
        return false;
    }
}

OpenHurricane::real OpenHurricane::RKF45::solve(const real t0, const realArray &y0, const realArray &dydt0,
                                        const real dt, realArray &y) {
    auto &s1 = dydt0;

    auto &yTmp = yTmp_;

    yTmp = y0 + a21 * dt * dydt0;
    DyDt(t0 + c2 * dt, yTmp, s2_);

    yTmp = y0 + dt * (a31 * s1 + a32 * s2_);
    DyDt(t0 + c3 * dt, yTmp, s3_);

    yTmp = y0 + dt * (a41 * s1 + a42 * s2_ + a43 * s3_);
    DyDt(t0 + c4 * dt, yTmp, s4_);

    yTmp = y0 + dt * (a51 * s1 + a52 * s2_ + a53 * s3_ + a54 * s4_);

    DyDt(t0 + c5 * dt, yTmp, s5_);

    yTmp = y0 + dt * (a61 * s1 + a62 * s2_ + a63 * s3_ + a64 * s4_ + a65 * s5_);

    DyDt(t0 + c6 * dt, yTmp, s6_);

    y = y0 + dt * (b1 * s1 + b3 * s3_ + b4 * s4_ + b5 * s5_ + b6 * s6_);
    error_ = dt * (e1 * s1 + e3 * s3_ + e4 * s4_ + e5 * s5_ + e6 * s6_);

    return maxError(y0, y, error_);
}

void OpenHurricane::RKF45::solveOnly(const real t0, const realArray &y0, const realArray &dydt0,
                                 const real dt, realArray &y) {
    auto &s1 = dydt0;

    auto &yTmp = yTmp_;

    yTmp = y0 + a21 * dt * dydt0;
    DyDt(t0 + c2 * dt, yTmp, s2_);

    yTmp = y0 + dt * (a31 * s1 + a32 * s2_);
    DyDt(t0 + c3 * dt, yTmp, s3_);

    yTmp = y0 + dt * (a41 * s1 + a42 * s2_ + a43 * s3_);
    DyDt(t0 + c4 * dt, yTmp, s4_);

    yTmp = y0 + dt * (a51 * s1 + a52 * s2_ + a53 * s3_ + a54 * s4_);

    DyDt(t0 + c5 * dt, yTmp, s5_);

    yTmp = y0 + dt * (a61 * s1 + a62 * s2_ + a63 * s3_ + a64 * s4_ + a65 * s5_);

    DyDt(t0 + c6 * dt, yTmp, s6_);

    y = y0 + dt * (b1 * s1 + b3 * s3_ + b4 * s4_ + b5 * s5_ + b6 * s6_);
}

OpenHurricane::real OpenHurricane::RKF45::solve(const real t0, const realArray &y0, const realArray &dydt0,
                                        const real dt, const integerListList &Gml, const integer im,
                                        realArray &y) {
    if (im >= Gml.size()) {
        errorAbortStr(("The group index (im): " + toString(im) + " is larger than the group size " +
                       toString(Gml.size())));
    }
    LFatal("This functions cannot be used!");
    return real();
}
