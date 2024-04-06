/*!
 * \file ImplicitEuler.cpp
 * \brief Main subroutines for implicit Euler ODE solver.
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

#include "ImplicitEuler.hpp"
#include "LUDecompose.hpp"

#include "controllerSwitch.hpp"

namespace OpenHurricane {
    createClassNameStr(ImplicitEuler, "ImplicitEuler");
    registerObjFty(ODEsSolver, ImplicitEuler, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::ImplicitEuler::maxIncrement(const realArray &y0,
                                                                     const realArray &y,
                                                                     const realArray &e) const {
    real maxE = Zero;
    for (integer i = 0; i < y0.size(); ++i) {
        if (isnan(y[i])) {
            return 1.5;
        }
        real tolelence = ATolIncrement_ + RTolIncrement_ * max(fabs(y0[i]), fabs(y[i]));
        maxE = max(maxE, fabs(e[i]) / tolelence);
    }
    return maxE;
}

OpenHurricane::ImplicitEuler::ImplicitEuler(const integer nEqsNew)
    : ODEsSolver(nEqsNew), dfdt_(nEqsNew, Zero), dydt_(nEqsNew, Zero), dfdy_(nEqsNew, nEqsNew),
      pivotIndices_(nEqsNew, Zero), yk_(nEqsNew, Zero), yTemp_(nEqsNew, Zero), y05k_(nEqsNew, Zero),
      isAutoUpdateJacobian_(false), stepToUpdateJac_(4), ATolIncrement_(1e-7),
      RTolIncrement_(0.001), yklast_(2), countStep_(0), lastTimeStep_(2, Zero) {}

OpenHurricane::ImplicitEuler::ImplicitEuler(const integer nEqsNew, const real atol, const real rtol,
                                        const integer maxStep)
    : ODEsSolver(nEqsNew, atol, rtol, maxStep), dfdt_(nEqsNew, Zero), dydt_(nEqsNew, Zero),
      dfdy_(nEqsNew, nEqsNew), pivotIndices_(nEqsNew, Zero), yk_(nEqsNew, Zero),
      yTemp_(nEqsNew, Zero), y05k_(nEqsNew, Zero), isAutoUpdateJacobian_(false),
      stepToUpdateJac_(4), ATolIncrement_(1e-7), RTolIncrement_(0.001), yklast_(2), countStep_(0),
      lastTimeStep_(2, Zero) {}

OpenHurricane::ImplicitEuler::ImplicitEuler(const integer nEqsNew, const controller &cont)
    : ODEsSolver(nEqsNew, cont), dfdt_(nEqsNew, Zero), dydt_(nEqsNew, Zero),
      dfdy_(nEqsNew, nEqsNew), pivotIndices_(nEqsNew, Zero), yk_(nEqsNew, Zero),
      yTemp_(nEqsNew, Zero), y05k_(nEqsNew, Zero), isAutoUpdateJacobian_(false),
      stepToUpdateJac_(10), ATolIncrement_(1e-7), RTolIncrement_(0.001), yklast_(2), countStep_(0),
      lastTimeStep_(2, Zero) {
    if (cont.found("ImplicitEuler")) {
        const auto &impECont = cont.subController("ImplicitEuler");

        controllerSwitch myContS(impECont);

        isAutoUpdateJacobian_ = myContS("autoUpdateJacobian", isAutoUpdateJacobian_);
        stepToUpdateJac_ = impECont.findOrDefault<integer>("stepToUpdateJac", stepToUpdateJac_);
        ATolIncrement_ = impECont.findOrDefault<real>("ATolIncrement", ATolIncrement_);
        RTolIncrement_ = impECont.findOrDefault<real>("RTolIncrement", RTolIncrement_);
    }
}

bool OpenHurricane::ImplicitEuler::reset(const integer nEqsNew) {
    if (ODEsSolver::reset(nEqsNew)) {
        resetArray(dfdt_);
        resetArray(dydt_);
        dfdy_.shallowResize(ODEsSolver::nEqs(), ODEsSolver::nEqs());
        resetArray(pivotIndices_);
        resetArray(yk_);
        resetArray(y05k_);
        resetArray(yTemp_);
        return true;
    } else {
        return false;
    }
}

OpenHurricane::real OpenHurricane::ImplicitEuler::solve(const real t0, const realArray &y0,
                                                const realArray &dydt0, const real dt,
                                                realArray &y) {
    yk_ = y0;
    bool isConvergence = true;
    auto errorN = NewtonIteration(t0, y0, dydt0, dt, yk_, isConvergence);
    if (!isConvergence) {
        return errorN;
    }

    y = yk_;
    if (countStep_ == 0) {
        realArray err00(y0.size(), Zero);
        DyDt(t0 + dt, yk_, dydt_);
        for (integer i = 0; i < y0.size(); ++i) {
            err00[i] = yk_[i] - y0[i] - dt * dydt_[i];
        }

        return maxIncrement(y0, yk_, err00);
    } else {
        realArray err00(y0.size(), Zero);

        // Delta t_n
        const real dtn = lastTimeStep_[0];

        // Delta t_(n-1)
        const real dtn1 = lastTimeStep_[1];

        const real f1 = 2.0 / (dtn * (dtn + dtn1));
        const real f2 = 2.0 / (dtn1 * (dtn + dtn1));
        const real LTEF = sqr(dtn) * 0.5;

        for (integer i = 0; i < y0.size(); ++i) {
            err00[i] = f1 * (yk_[i] - y0[i]) - f2 * (y0[i] - yklast_[1][i]);
            err00[i] *= LTEF;
        }

        return maxIncrement(y0, yk_, err00);
    }
    return errorN;
}

void OpenHurricane::ImplicitEuler::solveOnly(const real t0, const realArray &y0, const realArray &dydt0,
                                         const real dt, realArray &y) {
    LFatal("This function could not be called in implicit Euler solver.");
}

OpenHurricane::real OpenHurricane::ImplicitEuler::solve(const real t0, const realArray &y0,
                                                const realArray &dydt0, const real dt,
                                                const integerListList &Gml, const integer im,
                                                realArray &y) {
    if (im >= Gml.size()) {
        errorAbortStr(("The group index (im): " + toString(im) + " is larger than the group size " +
                       toString(Gml.size())));
    }
    LFatal("This functions cannot be used!");
    return real();
}

OpenHurricane::real OpenHurricane::ImplicitEuler::NewtonIteration(const real t0, const realArray &y0,
                                                          const realArray &dydt0, const real dt,
                                                          realArray &ykp1, bool &isConvergence) {
    integer n = ykp1.size();
    integer count = 0;
    real maxErr = 0.0;

    // Convergence rate constant
    real lastErr = 0.0;
    do {
        if (isAutoUpdateJacobian_) {
            if (count == 0 || count % stepToUpdateJac_ == 0) {
                jacobian(t0, ykp1, dydt_, dfdy_);
                for (integer i = 0; i < n; ++i) {
                    for (integer j = 0; j < n; ++j) {
                        dfdy_(i, j) = -dt * dfdy_(i, j);
                    }

                    dfdy_(i, i) += 1.0;
                }
                LUDecomposeDPivoting(dfdy_, pivotIndices_);
            } else {
                DyDt(t0, ykp1, dydt_);
            }
        } else {
            jacobian(t0, ykp1, dydt_, dfdy_);
            for (integer i = 0; i < n; ++i) {
                for (integer j = 0; j < n; ++j) {
                    dfdy_(i, j) = -dt * dfdy_(i, j);
                }

                dfdy_(i, i) += 1.0;
            }
            LUDecomposeDPivoting(dfdy_, pivotIndices_);
        }

        auto error = -ykp1 + y0 + dt * dydt_;
        LUBacksubstituteDPivoting(dfdy_, error, pivotIndices_);

        for (integer i = 0; i < ykp1.size(); ++i) {
            if (isnan(error[i]) || isinf(error[i])) {
                isConvergence = false;
                return 1000.0;
            }
            ykp1[i] += error[i];
        }

        maxErr = maxError(ykp1, ykp1, error);

#ifdef TEST_PROCESS_TIME
        countIter_++;
#endif

        if (count > 0) {
            // Diverging test
            real diver = maxErr / max(lastErr, tiny);

            // The iteration is diverging if diver > 2
            if (isAutoUpdateJacobian_) {
                if (diver > 2 && count > max(3 * stepToUpdateJac_, 5)) {
                    isConvergence = false;
                    return 1000.0;
                }
            } else {
                if (diver > 2 && count > 5) {
                    isConvergence = false;
                    return 1000.0;
                }
            }
        }
        lastErr = maxErr;
        count++;
        if (isnan(maxErr)) {
            if (report) {
                LWarning("Maximum error of Newton iteration is not-a-number");
            }
            isConvergence = false;
            return 1000.0;
        }
        if (count > 1000) {
            if (report) {
                LWarning("Exceed the maximum steps of Newton iteration: %d", count);
            }
            isConvergence = false;
            return maxErr;
        }
    } while (maxErr > 1.0);

    isConvergence = true;
    return maxErr;
}

void OpenHurricane::ImplicitEuler::adaptiveSolve(real &t, real &dt0, realArray &y) {
    real dt = dt0;
    realArray dydt0(y.size(), Zero);

    real error = Zero;

    yTemp_ = y;

    lastTimeStep_[1] = lastTimeStep_[0];
    integer count = 0;
    do {
        lastTimeStep_[0] = dt;
        error = solve(t, y, dydt0, dt, yTemp_);

        if (!checkNewY(dt, y, yTemp_) && count < 100) {
            error = 2;
            count++;
            real scale = max(safeScale_ * pow(error, -alphaInc_), minScale_);
            dt *= scale;
            if (dt <= veryTiny) {
                errorAbortStr(("step size underflow: " + toString(dt)));
            }
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

void OpenHurricane::ImplicitEuler::solve(const real t0, real &tn, real &dt0, realArray &y) {
#ifdef TEST_PROCESS_TIME
    countIter_ = 0;
#endif
    real t = t0;
    real dt = 0.2 * dt0;
    real dtReached = t0;
    bool last = false;
    for (integer n = 0; n < maxSteps_; ++n) {
        countStep_ = n;
        if (n == 0) {
            if (yklast_[0].size() != y.size()) {
                yklast_[0].resize(y.size());
            }
            yklast_[0] = y;
        } else {
            if (yklast_[1].size() != y.size()) {
                yklast_[1].resize(y.size());
            }
            yklast_[1] = yklast_[0];
            yklast_[0] = y;
        }
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
