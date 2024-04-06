/*!
 * \file ODEsSolver.cpp
 * \brief Main subroutines for the ODEs solver.
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

#include "ODEsSolver.hpp"

namespace OpenHurricane {
    createClassNameStr(ODEsSolver, "ODEsSolver");
    createObjFty(ODEsSolver, controller);
} // namespace OpenHurricane

OpenHurricane::real OpenHurricane::ODEsSolver::maxError(const realArray &y0, const realArray &y,
                                                const realArray &e) const {
    real maxE = Zero;
    for (integer i = 0; i < y0.size(); ++i) {
        real tolelence = ATOL_[i] + RTOL_[i] * max(fabs(y0[i]), fabs(y[i]));

        maxE = max(maxE, fabs(e[i]) / tolelence);
    }
    return maxE;
}

OpenHurricane::real OpenHurricane::ODEsSolver::maxError(const realArray &y0, const realArray &y,
                                                const realArray &e, const integerListList &Gml,
                                                const integer im) const {
    real maxE = Zero;
    for (integer j = 0; j < Gml[im].size(); ++j) {
        const integer i = Gml[im][j];

        real tolelence = ATOL_[i] + RTOL_[i] * max(fabs(y0[i]), fabs(y[i]));
        maxE = max(maxE, fabs(e[i]) / tolelence);
    }
    return maxE;
}

OpenHurricane::ODEsSolver::ODEsSolver(const integer nEqs)
    : nEqs_(nEqs), maxEqns_(nEqs), maxSteps_(10000), ATOL_(nEqs, tiny), RTOL_(nEqs, 1e-4),
      safeScale_(0.8), minScale_(0.25), maxScale_(2.0), alphaInc_(0.25),
#ifdef TEST_PROCESS_TIME
      countIter_(0),
#endif
      yTemp_(nEqs, Zero), dydt0_(nEqs, Zero) {
}

OpenHurricane::ODEsSolver::ODEsSolver(const integer nEqs, const real atol, const real rtol,
                                  const integer maxStep)
    : nEqs_(nEqs), maxEqns_(nEqs), maxSteps_(maxStep), ATOL_(nEqs, atol), RTOL_(nEqs, rtol),
      safeScale_(0.8), minScale_(0.25), maxScale_(2.0), alphaInc_(0.25),
#ifdef TEST_PROCESS_TIME
      countIter_(0),
#endif
      yTemp_(nEqs, Zero), dydt0_(nEqs, Zero) {
}

OpenHurricane::ODEsSolver::ODEsSolver(const integer nEqs, const controller &cont)
    : nEqs_(nEqs), maxEqns_(nEqs), maxSteps_(cont.findOrDefault<integer>("maxstep", 10000)),
      ATOL_(nEqs, cont.findOrDefault<real>("atol", tiny)),
      RTOL_(nEqs, cont.findOrDefault<real>("rtol", 1e-4)),
      safeScale_(cont.findOrDefault<real>("safescale", 0.8)),
      minScale_(cont.findOrDefault<real>("minscale", 0.25)),
      maxScale_(cont.findOrDefault<real>("maxScale", 2.0)),
      alphaInc_(cont.findOrDefault<real>("alphainc", 0.25)),
#ifdef TEST_PROCESS_TIME
      countIter_(0),
#endif
      yTemp_(nEqs, Zero), dydt0_(nEqs, Zero) {
}

OpenHurricane::uniquePtr<OpenHurricane::ODEsSolver> OpenHurricane::ODEsSolver::creator(const integer nEqs,
                                                                           const controller &cont) {
    string solverType = cont.findWord(ODEsSolver::className_);

    defineInObjCreator(ODEsSolver, static_cast<std::string>(solverType), controller, (nEqs, cont));
}

void OpenHurricane::ODEsSolver::solve(const real t0, const real tn, real &tReached, const real &dt,
                                  const integerListList &Gml, const integer im,
                                  const bool advancingEnd, realArray &y) {
    tReached = t0;
    real subDt = dt;
    bool last = false;
    integer count = 0;
    auto &dydt0 = dydt0_;
    while (true) {
        if ((tReached + subDt - tn) * (tReached + subDt - t0) > 0) {
            subDt = tn - tReached;
            last = true;
        }
        yTemp_ = y;
        DyDt(tReached, yTemp_, dydt0);
        real error = solve(tReached, yTemp_, dydt0, subDt, Gml, im, y);
        tReached += subDt;
        if (error < 1.0 && !advancingEnd) {
            break;
        }
        if (last) {
            break;
        }
        if (count++ > maxSteps_) {
            break;
        }
    }
}

void OpenHurricane::ODEsSolver::solveList(const real t0, const real tnn, real &tReached,
                                      const realArray &dt, const integerListList &Gml,
                                      realArray &y) {
    real tn = tnn;
    real dt0 = tnn - t0;
    if (dt.size() <= 1) {
        solve(t0, tn, dt0, y);
        tReached = tn;
        return;
    }

    tn = t0;
    for (integer im = Gml.size() - 1; im >= 0; --im) {
        solve(tn, tnn, tReached, min(dt[im], dt0), Gml, im, (im == 0), y);
        tn = tReached;
    }
}

void OpenHurricane::ODEsSolver::adaptiveSolve(real &t, real &dt0, realArray &y) {
    real dt = dt0;
    auto &dydt0 = dydt0_;

    DyDt(t, y, dydt0);

    real error = Zero;
    //realArray yTemp_(y.size());
    yTemp_ = y;

    do {

#ifdef TEST_PROCESS_TIME
        countIter_++;
#endif
        error = solve(t, y, dydt0, dt, yTemp_);
        if (!checkNewY(dt, y, yTemp_)) {
            error = 2;
            continue;
        }
        if (error > 1.0) {
            real scale = max(safeScale_ * pow(error, -alphaInc_), minScale_);
            dt *= scale;
            if (dt <= veryTiny) {
                errorAbortStr(("step size underflow: " + toString(dt)));
            }
        }
    } while (error > 1.0);

    t += dt;
    y = yTemp_;
    if (error > pow(maxScale_ / safeScale_, real(-1.0) / alphaInc_)) {
        dt0 = min(max(safeScale_ * pow(error, -alphaInc_), minScale_), maxScale_) * dt;
    } else {
        dt0 = max(safeScale_ * maxScale_ * dt, dt);
    }
}

void OpenHurricane::ODEsSolver::solve(const real t0, real &tn, real &dt0, realArray &y) {

#ifdef TEST_PROCESS_TIME
    countIter_ = 0;
#endif
    real t = t0;
    real dt = dt0;
    real dt00 = dt0;
    real dtReached = t0;
    bool last = false;
    for (integer n = 0; n < maxSteps_; ++n) {
        real dtTry0 = dt;
        if ((t + dt - tn) * (t + dt - t0) > 0) {
            dt = tn - t;
            last = true;
        } else if (t + dt == tn) {
            last = true;
        }
        real tt0 = t;
        adaptiveSolve(t, dt, y);
        dtReached = t;
        if ((t - tn) * (tn - t0) >= 0) {
            if (n > 0 && last) {
                dt0 = min(dtTry0, tn);
            }

            dt = dt0;

            return;
        }
    }
    tn = dtReached;
    return;
    errorAbortStr(("Up to the maximum steps: " + toString(maxSteps_) +
                   ". The total time step is: " + toString(tn) + "s, " +
                   "the time did is: " + toString(dtReached) + "s"));
}

void OpenHurricane::ODEsSolver::solve(const real t0, const real tn, real &tReached,
                                  const realArray &dtl, const integerListList &Gml, realArray &y) {
    if (dtl.size() <= 1) {
        real dt0 = tn - t0;
        real tnn = tn;
        solve(t0, tnn, dt0, y);
        tReached = tnn;
        return;
    }

    real t = t0;
    real dt = tn - t0;
    //real dt0 = tn - t0;
    real dtReached = t0;
    //bool last = false;
    realArray y0 = y;
    realArray y1(nEqs_);
    for (integer n = 0; n < maxSteps_; ++n) {
        if ((t + dt - tn) * (t + dt - t0) > 0) {
            dt = tn - t;
            //last = true;
        }
        real tt0 = t;
        //dt0 = dt;
        y1 = y0;
        real maxErr = getError(tt0, tt0 + dt, t, dtl, Gml, y0);
        dtReached = t;
        if ((t - tn) * (t - t0) >= 0) {
            if (maxErr < 1.0) {
                y = y0;
                return;
            } else {
                t = tt0;
                dt *= 0.5;
                y0 = y1;
            }
        } else {
            if (maxErr >= 1.0) {
                t = tt0;
                dt *= 0.5;
                y0 = y1;
            }
        }
    }
    errorAbortStr(("Up to the maximum steps: " + toString(maxSteps_) +
                   ". The time should be reached is: " + toString(tn) + "s, " +
                   "the time reached is: " + toString(dtReached) + "s"));
}
