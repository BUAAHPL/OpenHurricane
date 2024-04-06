/*!
 * \file BDF23.cpp
 * \brief Main subroutines for BDF23 solver.
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

#include "BDF23.hpp"
#include "LUDecompose.hpp"

#include "controllerSwitch.hpp"

namespace OpenHurricane {
    createClassNameStr(BDF23, "BDF23");
    registerObjFty(ODEsSolver, BDF23, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real
OpenHurricane::BDF23::maxIncrement(const realArray &y0, const realArray &y, const realArray &e) const {
    real maxE = Zero;
    for (integer i = 0; i < y0.size(); ++i) {
        if (isnan(y[i])) {
            return 1.5;
        }
        real tolelence = ATolIncrement_ + RTolIncrement_ * max(fabs(y0[i]), fabs(y[i]));
        maxE = max(maxE, fabs(e[i]) / tolelence);
        //  maxE += sqr(fabs(e[i]) / tolelence);
    }
    return maxE;
    //return sqrt(maxE / y0.size());
}

void OpenHurricane::BDF23::setAlpha(const integer iord) {
    if (iord == 1) {
        alpha_.resize(2);
        alpha_[0] = 1.0;
        alpha_[1] = -1.0;
    } else if (iord == 2) {
        /**
         * The previous time step [s].
         * lastTimeStep_[0] = dt^n i.e. current time step size.
         * lastTimeStep_[1] = dt^(n-1) i.e. previous time step size.
         * lastTimeStep_[2] = dt^(n-2) i.e.  previous time step size.
         */
        const auto &dtn = lastTimeStep_;

        const real Dtn = dtn[0];
        const real Dtn1 = dtn[1];
        //const real Dtn2 = dtn[2];

        alpha_.resize(3);
        /*alpha_[0] = 3.0 / 2.0;
        alpha_[1] = -2.0;
        alpha_[2] = 0.5;*/

        //alpha_[0] = (2.0 * Dtn + Dtn1) / (Dtn * (Dtn + Dtn1));
        alpha_[1] = -(Dtn + Dtn1) / (Dtn * Dtn1);
        alpha_[2] = Dtn / (Dtn1 * (Dtn + Dtn1));
        alpha_[0] = -(alpha_[1] + alpha_[2]);
        alpha_ *= Dtn;
    } else if (iord == 3) {
        /**
         * The previous time step [s].
         * lastTimeStep_[0] = dt^n i.e. current time step size.
         * lastTimeStep_[1] = dt^(n-1) i.e. previous time step size.
         * lastTimeStep_[2] = dt^(n-2) i.e.  previous time step size.
         */
        const auto &dtn = lastTimeStep_;

        const real Dtn = dtn[0];
        const real Dtn1 = dtn[1];
        const real Dtn2 = dtn[2];

        alpha_.resize(4);
        /*alpha_[0] = 11.0 / 6.0;
        alpha_[1] = -3.0;
        alpha_[2] = 3.0 / 2.0;
        alpha_[3] = -1.0/3.0;*/

        alpha_[1] = -(Dtn + Dtn1) * (Dtn + Dtn1 + Dtn2) / (Dtn * Dtn1 * (Dtn1 + Dtn2));
        alpha_[2] = Dtn * (Dtn + Dtn1 + Dtn2) / (Dtn1 * Dtn2 * (Dtn + Dtn1));
        alpha_[3] = -Dtn * (Dtn + Dtn1) / (Dtn2 * (Dtn1 + Dtn2) * (Dtn + Dtn1 + Dtn2));
        alpha_[0] = -(alpha_[1] + alpha_[2] + alpha_[3]);
        alpha_ *= Dtn;
    } else {
        errorAbortStr(("Unsupported order:" + toString(iord) + " for BDF"));
    }
}

OpenHurricane::BDF23::BDF23(const integer nEqs)
    : ODEsSolver(nEqs), dfdt_(nEqs, Zero), dydt_(nEqs, Zero), dfdy_(nEqs, nEqs),
      pivotIndices_(nEqs, Zero), yk_(nEqs, Zero), yTemp_(nEqs, Zero), y05k_(nEqs, Zero), yklast_(4),
      isAutoUpdateJacobian_(false), stepToUpdateJac_(4), ATolIncrement_(1e-8),
      RTolIncrement_(0.001), alpha_(), iOrder_(2), iCorder_(1), lastTimeStep_(4), initialFct_(0.1),
      countJacStep_(0) {}

OpenHurricane::BDF23::BDF23(const integer nEqs, const real atol, const real rtol, const integer maxStep)
    : ODEsSolver(nEqs, atol, rtol, maxStep), dfdt_(nEqs, Zero), dydt_(nEqs, Zero),
      dfdy_(nEqs, nEqs), pivotIndices_(nEqs, Zero), yk_(nEqs, Zero), yTemp_(nEqs, Zero),
      y05k_(nEqs, Zero), yklast_(4), isAutoUpdateJacobian_(false), stepToUpdateJac_(4),
      ATolIncrement_(1e-8), RTolIncrement_(0.001), alpha_(), iOrder_(2), iCorder_(1),
      lastTimeStep_(4), initialFct_(0.1), countJacStep_(0) {}

OpenHurricane::BDF23::BDF23(const integer nEqs, const controller &cont)
    : ODEsSolver(nEqs, cont), dfdt_(nEqs, Zero), dydt_(nEqs, Zero), dfdy_(nEqs, nEqs),
      pivotIndices_(nEqs, Zero), yk_(nEqs, Zero), yTemp_(nEqs, Zero), y05k_(nEqs, Zero), yklast_(4),
      isAutoUpdateJacobian_(false), stepToUpdateJac_(4), ATolIncrement_(1e-8),
      RTolIncrement_(0.001), alpha_(), iOrder_(2), iCorder_(1), lastTimeStep_(4), initialFct_(0.1),
      countJacStep_(0) {
    if (cont.found("BDF23")) {
        const auto &impECont = cont.subController("BDF23");

        controllerSwitch myContS(impECont);

        isAutoUpdateJacobian_ = myContS("autoUpdateJacobian", isAutoUpdateJacobian_);
        stepToUpdateJac_ = impECont.findOrDefault<integer>("stepToUpdateJac", stepToUpdateJac_);
        iOrder_ = impECont.findOrDefault<integer>("Order", iOrder_);
        initialFct_ = impECont.findOrDefault<real>("initialFactor", initialFct_);
        if (iOrder_ != 2 && iOrder_ != 3) {
            checkWarningStr(("The order must be 2 or 3 but given: " + toString(iOrder_) + " in " +
                             impECont.name() + ". Default value: 2 will be used."));
        }

        ATolIncrement_ = impECont.findOrDefault<real>("ATolIncrement", ATolIncrement_);
        RTolIncrement_ = impECont.findOrDefault<real>("RTolIncrement", RTolIncrement_);
    }
}

bool OpenHurricane::BDF23::reset(const integer nEqsNew) {
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

OpenHurricane::real OpenHurricane::BDF23::solve(const real t0, const realArray &y0, const realArray &dydt0,
                                        const real dt, realArray &y) {
    y = 0;
    for (integer k = alpha_.size() - 1; k > 0; --k) {
        y += alpha_[k] * yklast_[k - 1];
    }

    yk_ = y0;
    bool isConvergence = true;
    auto errorN = NewtonIteration(t0, y, dydt0, dt, yk_, isConvergence);
    if (!isConvergence) {
        return errorN;
    }

    y = yk_;

    if (countJacStep_ <= iOrder_) {
        realArray err00(y0.size(), Zero);
        DyDt(t0 + dt, yk_, dydt_);

        for (integer i = 0; i < y0.size(); ++i) {
            err00[i] = yk_[i] - y0[i] - dt * dydt_[i];
        }

        return maxIncrement(y0, yk_, err00);
    } else {
        realArray err00(y0.size(), Zero);

        if (iOrder_ == 2) {
            // Delta t_n
            const real dtn = lastTimeStep_[0];

            // Delta t_(n-1)
            const real dtn1 = lastTimeStep_[1];

            // Delta t_(n-2)
            const real dtn2 = lastTimeStep_[2];

            const real f1 = inv(dtn);
            const real f2 = inv(dtn1) * (1.0 + dtn / dtn1);
            const real f3 = dtn / (dtn1 * dtn2);

            const real ltef = (dtn + dtn1) / 6.0;
            for (integer i = 0; i < y0.size(); ++i) {
                err00[i] = f1 * (yk_[i] - yklast_[0][i]) - f2 * (yklast_[0][i] - yklast_[1][i]) +
                           f3 * (yklast_[1][i] - yklast_[2][i]);

                err00[i] *= ltef;
            }

            /*real LTEF = sqr(dtn) * (dtn + dtn1) / 6.0;
            LTEF /= 6.0 * (dtn2 * dtn1 * (dtn2 + dtn1) * dtn * (dtn1 + dtn) * (dtn2 + dtn1 + dtn));

            real f1 = dtn2 * dtn1 * (dtn2 + dtn1);
            real f2 = dtn2 * (dtn1 + dtn) * (dtn2 + dtn1 + dtn);
            real f3 = dtn * (dtn2 + dtn1) * (dtn2 + dtn1 + dtn);
            real f4 = dtn * dtn1 * (dtn1 + dtn);

            f1 *= LTEF;
            f2 *= LTEF;
            f3 *= LTEF;
            f4 *= LTEF;
            for (integer i = 0; i < y0.size(); ++i)
            {
                err00[i] = f1 * yk_[i] - f2 * yklast_[0][i] + f3 * yklast_[1][i] - f3 * yklast_[2][i];
            }*/
        } else {
            // Delta t_n
            const real dtn = lastTimeStep_[0];

            // Delta t_(n-1)
            const real dtn1 = lastTimeStep_[1];

            // Delta t_(n-2)
            const real dtn2 = lastTimeStep_[2];

            // Delta t_(n-3)
            const real dtn3 = lastTimeStep_[3];

            real LTEF = alpha_[1] * pow4(dtn) + alpha_[2] * pow4(dtn + dtn1) +
                        alpha_[3] * pow4(dtn + dtn1 + dtn2);
            ;
            LTEF /= 24.0;

            real f1 = dtn * (dtn1 + dtn) * (dtn2 + dtn1 + dtn) * (dtn3 + dtn2 + dtn1 + dtn);
            real f2 = dtn * dtn1 * (dtn2 + dtn1) * (dtn3 + dtn2 + dtn1);
            real f3 = dtn2 * dtn1 * (dtn3 + dtn2) * (dtn1 + dtn);
            real f4 = dtn3 * dtn2 * (dtn1 + dtn2) * (dtn2 + dtn1 + dtn);
            real f5 = dtn3 * (dtn2 + dtn3) * (dtn1 + dtn2 + dtn3) * (dtn3 + dtn2 + dtn1 + dtn);
            f1 = inv(f1) * LTEF;
            f2 = inv(f2) * LTEF;
            f3 = inv(f3) * LTEF;
            f4 = inv(f4) * LTEF;
            f5 = inv(f5) * LTEF;

            for (integer i = 0; i < y0.size(); ++i) {
                err00[i] = f1 * yk_[i] - f2 * yklast_[0][i] + f3 * yklast_[1][i] -
                           f4 * yklast_[2][i] + f5 * yklast_[3][i];
            }
        }
        return maxIncrement(y0, yk_, err00);
    }
    return 0.55;
}

void OpenHurricane::BDF23::solveOnly(const real t0, const realArray &y0, const realArray &dydt0,
                                 const real dt, realArray &y) {
    LFatal("This function could not be called in implicit Euler solver.");
}

OpenHurricane::real OpenHurricane::BDF23::solve(const real t0, const realArray &y0, const realArray &dydt0,
                                        const real dt, const integerListList &Gml, const integer im,
                                        realArray &y) {
    LFatal("This functions cannot be used!");
    return real();
}

OpenHurricane::real OpenHurricane::BDF23::NewtonIteration(const real t0, const realArray &y0,
                                                  const realArray &dydt0, const real dt,
                                                  realArray &ykp1, bool &isConvergence) {
    integer n = ykp1.size();
    integer count = 0;
    real maxErr = 0.0;
    bool shouldUpdateJac = false;
    if (countJacStep_ != 0 && (lastTimeStep_[0] / lastTimeStep_[1] >= 3.0)) {
        shouldUpdateJac = true;
    }

    // Convergence rate constant
    real lastErr = 0.0;
    do {
        if (isAutoUpdateJacobian_) {
            if (countJacStep_ == 0 || countJacStep_ % stepToUpdateJac_ == 0 || shouldUpdateJac) {
                jacobian(t0, ykp1, dydt_, dfdy_);
                for (integer i = 0; i < n; ++i) {
                    for (integer j = 0; j < n; ++j) {
                        dfdy_(i, j) = -dt * dfdy_(i, j);
                    }

                    dfdy_(i, i) += alpha_[0];
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

                dfdy_(i, i) += alpha_[0];
            }
            LUDecomposeDPivoting(dfdy_, pivotIndices_);
        }

        auto error = -(alpha_[0] * ykp1 + y0 - dt * dydt_);

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
            if (diver > 2) {
                isConvergence = false;
                return 1000.0;
            }
        }
        lastErr = maxErr;

        count++;

        if (count > 15 && count % 5 == 0) {
            shouldUpdateJac = true;
        } else {
            shouldUpdateJac = false;
        }

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

void OpenHurricane::BDF23::adaptiveSolve(real &t, real &dt0, realArray &y) {
    real dt = dt0;
    realArray dydt0(y.size(), Zero);

    real error = Zero;

    // yTemp_ = y;
    lastTimeStep_[3] = lastTimeStep_[2];
    lastTimeStep_[2] = lastTimeStep_[1];
    lastTimeStep_[1] = lastTimeStep_[0];
    real lastError = large;
    integer count = 0;
    do {
        lastTimeStep_[0] = dt;
        setAlpha(iCorder_);

        error = solve(t, y, dydt0, dt, yTemp_);

        if (!checkNewY(dt, y, yTemp_) && count < 100) {
            count++;
            error = 2;
            continue;
        }
        if (lastError < error && (error < 5)) {
            break;
        }
        if (error > 1.0) {
            real scale = max(safeScale_ * pow(error, -alphaInc_), minScale_);
            dt *= scale;
            if (dt <= veryTiny) {
                errorAbortStr(("step size underflow: " + toString(dt)));
            }
        }
        lastError = error;
    } while (error > 1.0);

    t += dt;
    y = yTemp_;
    if (error > pow(maxScale_ / safeScale_, real(-1.0) / alphaInc_)) {
        dt0 = min(max(safeScale_ * pow(error, -alphaInc_), minScale_), maxScale_) * dt;
    } else {
        dt0 = max(safeScale_ * maxScale_ * dt, dt);
    }
}

void OpenHurricane::BDF23::solve(const real t0, real &tn, real &dt0, realArray &y) {
#ifdef TEST_PROCESS_TIME
    countIter_ = 0;
#endif
    real t = t0;
    real dt = dt0;
    real dt00 = dt0;
    real dtReached = t0;
    bool last = false;
    for (integer n = 0; n < maxSteps_; ++n) {
        countJacStep_ = n;
        if (n == 0) {
            dt = min(initialFct_ * tn, real(1e-8));
            iCorder_ = 1;
            if (yklast_[0].size() != y.size()) {
                yklast_[0].resize(y.size());
            }
            yklast_[0] = y;
        } else if (n == 1) {
            iCorder_ = 2;
            if (yklast_[1].size() != y.size()) {
                yklast_[1].resize(y.size());
            }
            yklast_[1] = yklast_[0];
            yklast_[0] = y;
        } else {
            iCorder_ = iOrder_;
            if (n <= 3) {
                if (yklast_[2].size() != y.size()) {
                    yklast_[2].resize(y.size());
                }
                yklast_[2] = yklast_[1];
                yklast_[1] = yklast_[0];
                yklast_[0] = y;
            } else {
                if (n <= 5) {
                    dt *= 2.0;
                }
                if (yklast_[3].size() != y.size()) {
                    yklast_[3].resize(y.size());
                }
                yklast_[3] = yklast_[2];
                yklast_[2] = yklast_[1];
                yklast_[1] = yklast_[0];
                yklast_[0] = y;
            }
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
