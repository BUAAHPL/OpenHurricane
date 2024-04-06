/*!
 * \file Euler.cpp
 * \brief Main subroutines for Euler ODE solver.
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

#include "Euler.hpp"

namespace OpenHurricane {
    createClassNameStr(Euler, "Euler");
    registerObjFty(ODEsSolver, Euler, controller);
} // namespace OpenHurricane

OpenHurricane::Euler::Euler(const integer nEqs)
    : ODEsSolver(nEqs), yTemph2_(nEqs, Zero), error_(nEqs, Zero), dydt2_(nEqs, Zero) {}

OpenHurricane::Euler::Euler(const integer nEqs, const real atol, const real rtol, const integer maxStep)
    : ODEsSolver(nEqs, atol, rtol, maxStep), yTemph2_(nEqs, Zero), error_(nEqs, Zero),
      dydt2_(nEqs, Zero) {}

OpenHurricane::Euler::Euler(const integer nEqs, const controller &cont)
    : ODEsSolver(nEqs, cont), yTemph2_(nEqs, Zero), error_(nEqs, Zero), dydt2_(nEqs, Zero) {}

bool OpenHurricane::Euler::reset(const integer nEqsNew) {
    if (ODEsSolver::reset(nEqsNew)) {
        resetArray(yTemph2_);
        resetArray(error_);
        resetArray(dydt2_);
        return true;
    } else {
        return false;
    }
}

OpenHurricane::real OpenHurricane::Euler::solve(const real t0, const realArray &y0, const realArray &dydt0,
                                        const real dt, realArray &y) {
    error_ = dt * dydt0;
    yTemph2_ = y0 + error_;

    y = y0 + 0.5 * error_;
    auto &dydth2 = dydt2_;
    DyDt(t0 + real(0.5) * dt, y, dydth2);
    y += 0.5 * dt * dydth2;

    error_ = 2.0 / (1.0 - 2.0) * (y - yTemph2_);

    return maxError(yTemph2_, y, error_);
}

void OpenHurricane::Euler::solveOnly(const real t0, const realArray &y0, const realArray &dydt0,
                                 const real dt, realArray &y) {
    y = y0 + dt * dydt0;
}

OpenHurricane::real OpenHurricane::Euler::solve(const real t0, const realArray &y0, const realArray &dydt0,
                                        const real dt, const integerListList &Gml, const integer im,
                                        realArray &y) {
    if (im >= Gml.size()) {
        errorAbortStr(("The group index (im): " + toString(im) + " is larger than the group size " +
                       toString(Gml.size())));
    }
    error_ = Zero;
    for (int m = 0; m <= im; ++m) {
        for (int i = 0; i < Gml[m].size(); ++i) {
            const int j = Gml[m][i];
            error_[j] = dt * dydt0[j];
            y[j] = y0[j] + error_[j];
        }
    }

    for (int m = Gml.size() - 1; m > im; --m) {
        for (int i = 0; i < Gml[m].size(); ++i) {
            const int j = Gml[m][i];
            y[j] = y0[j];
        }
    }

    return maxError(y0, y, error_, Gml, im);
}

OpenHurricane::real OpenHurricane::Euler::getError(const real t0, const real tbase, real &tReached,
                                           const realArray &dt, const integerListList &Gml,
                                           realArray &y) {
    real tn = tbase;
    real t05n = 0.5 * (tn - t0) + t0;
    real tR = t0;
    real t05R = t0;
    realArray y0 = y;
    realArray y05 = y;

    // Solve dt = h, t_initial = t0
    ODEsSolver::solveList(t0, tn, tR, dt, Gml, y0);

    // Solve dt = h/2, t_initial = t0
    ODEsSolver::solveList(t0, t05n, t05R, dt, Gml, y05);
    real tnh2 = tbase;
    // Solve dt = h/2, t_initial = t05R
    ODEsSolver::solveList(t05R, tnh2, t05R, dt, Gml, y05);

    realArray error = 2.0 * (1.0 - 2.0) * (y0 - y05);

    y = y0;
    tReached = tR;
    return maxError(y0, y05, error);
}
