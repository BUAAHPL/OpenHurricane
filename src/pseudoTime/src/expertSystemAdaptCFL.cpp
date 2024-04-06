/*!
 * \file expertSystemAdaptCFL.cpp
 * \brief The subroutines and functions of expert-system adapt CFL number.
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

#include "expertSystemAdaptCFL.hpp"
#include "linearRegression.hpp"

namespace OpenHurricane {
    createClassName(expertSystemAdaptCFL);
    registerObjFty(CFL, expertSystemAdaptCFL, controller);
} // namespace OpenHurricane

OpenHurricane::expertSystemAdaptCFL::expertSystemAdaptCFL(const iteration &iter,
                                                          const runtimeMesh &mesh,
                                                          const controller &cont)
    : CFL(iter, mesh, cont),
      CFLFactor0_(
          cont.subController("expertSystemAdaptCFL").findOrDefault<real>("CFLFactor0", 0.01)),
      CFLFactor0Max_(5.0),
      breakdownFactor_(
          cont.subController("expertSystemAdaptCFL").findOrDefault<real>("breakdownFactor", 0.5)),
      divergenceFactor_(
          cont.subController("expertSystemAdaptCFL").findOrDefault<real>("divergenceFactor", 0.8)),
      residualJumpThreshold_(cont.subController("expertSystemAdaptCFL")
                                 .findOrDefault<real>("residualJumpThreshold", 0.5)),
      updateFactor_(
          cont.subController("expertSystemAdaptCFL").findOrDefault<real>("updateFactor", 2.0)),
      interval_(cont.subController("expertSystemAdaptCFL").findOrDefault<integer>("interval", 10)),
      countUnchange_(0), countResUp_(0), m2_(4, Zero), lastMaxRes_(1.0), lastAveRes_(1.0),
      lastEMaxRes_(1.0), lastEAveRes_(1.0), isBreakdown_(false) {
    if (iter_.totalStep() <= initialStage_) {
        cCFL_ = 0.1 * CFLMax_;
    } else {
        cCFL_ = CFLMax_;
    }
    cCFL_ = max(CFLMin_, cCFL_);
    if (!cont.found("stepForPrintCFL")) {
        stepForPrintCFL_ = 20;
    }
}

OpenHurricane::real OpenHurricane::expertSystemAdaptCFL::getCFL() const {
    cCFL_ = max(CFLMin_, cCFL_);
    cCFL_ = min(CFLMax_, cCFL_);
    const auto &rho = mesh().findObjectRef<cellRealArray>("rho");
    const auto &E = mesh().findObjectRef<cellRealArray>("E");
    rho.getAveAndMaxRHS();
    E.getAveAndMaxRHS();
    auto rhs0 = rho.rhsAve0();
    auto rhsAve = rho.rhsAve();
    auto rhsMax = rho.rhsMax();

    auto Erhs0 = E.rhsAve0();
    auto ErhsAve = E.rhsAve();
    auto ErhsMax = E.rhsMax();

    real mm1 = Zero;
    real mm3 = Zero;

    if (isBreakdown_) {
        countUnchange_ = 0;
        breakdowns();
        isBreakdown_ = false;
        cCFL_ = max(CFLMin_, cCFL_);
        countResUp_++;
    } else if ((log10(rhsMax / lastMaxRes_) > residualJumpThreshold_ ||
                log10(rhsAve / lastAveRes_) > residualJumpThreshold_) ||
               (log10(ErhsMax / lastEMaxRes_) > residualJumpThreshold_ ||
                log10(ErhsAve / lastEAveRes_) > residualJumpThreshold_)) {
        divergence();
        countUnchange_ = 0;
        cCFL_ = max(CFLMin_, cCFL_);
        countResUp_++;
    } else if (notChange_ && ((!iter().hasSubIteration() && iter_.totalStep() <= initialStage_) ||
                              iter().hasSubIteration())) {
        notChange_ = false;
        countUnchange_ = 0;
        //divergence();
        cCFL_ = max(CFLMin_, cCFL_);
    } else {
        if ((!iter().hasSubIteration() && iter_.totalStep() <= initialStage_)) {
            if ((rhsMax > lastMaxRes_ || rhsAve > lastAveRes_) ||
                (ErhsMax > lastEMaxRes_ || ErhsAve > lastEAveRes_)) {
                countUnchange_ = 0;
                countResUp_++;
                if (countResUp_ >= max(interval_, integer(100))) {
                    divergence();
                    countResUp_ = 0;
                }
                cCFL_ = max(CFLMin_, cCFL_);
            } else {
                countResUp_ = 0;
                countUnchange_++;
                //if (countUnchange_ >= max(integer(interval_ / 2), 1))
                if (activateAdaptCFL(countUnchange_)) {
                    slowConvergence();
                    mm1 = log10(rhs0 / rhsAve);
                    mm3 = slope(m2_);
                    if (mm1 > 4 && m2_[m2_.size() - 1] > 4 && mm3 > 0) {
                        CFLFactor0_ *= 2.0;
                    }
                    CFLFactor0_ = min(CFLFactor0Max_, CFLFactor0_);

                    cCFL_ *= (1.0 + CFLFactor0_);
                    cCFL_ = max(CFLMin_, cCFL_);
                    cCFL_ = min(CFLMax_, cCFL_);
                }
            }
        } else {
            if ((rhsMax > lastMaxRes_ || rhsAve > lastAveRes_) ||
                (ErhsMax > lastEMaxRes_ || ErhsAve > lastEAveRes_)) {
                countResUp_++;
            } else {
                countResUp_ = 0;
            }

            if (countResUp_ >= max(interval_, integer(100))) {
                divergence();
                countUnchange_ = 0;
                countResUp_ = 0;
                cCFL_ = max(CFLMin_, cCFL_);
            } else {
                countUnchange_++;
                //if (countUnchange_ >= max(integer(interval_ / 2), 1))
                if (activateAdaptCFL(countUnchange_)) {
                    slowConvergence();
                    mm1 = log10(rhs0 / rhsAve);
                    mm3 = slope(m2_);
                    if (mm1 > 4 && m2_[m2_.size() - 1] > 4 && mm3 > 0) {
                        CFLFactor0_ *= 2.0;
                    }
                    CFLFactor0_ = min(CFLFactor0Max_, CFLFactor0_);

                    cCFL_ *= (1.0 + CFLFactor0_);
                    cCFL_ = max(CFLMin_, cCFL_);
                    cCFL_ = min(CFLMax_, cCFL_);
                }
            }
        }
    }

    lastMaxRes_ = rhsMax;
    lastAveRes_ = rhsAve;
    lastEMaxRes_ = ErhsMax;
    lastEAveRes_ = ErhsAve;
    if (!iter().hasSubIteration()) {
        printCFL(cCFL_);
    }
    return cCFL_;
}

void OpenHurricane::expertSystemAdaptCFL::setRelativeCorrection(const real relCorre) const {
    for (integer i = 0; i < m2_.size() - 1; ++i) {
        m2_[i] = m2_[i + 1];
    }
    m2_[m2_.size() - 1] = -log10(relCorre);
}

void OpenHurricane::expertSystemAdaptCFL::setBreakdowns() const {
    isBreakdown_ = true;
}

void OpenHurricane::expertSystemAdaptCFL::breakdowns() const {
    CFLFactor0_ *= breakdownFactor_;
    cCFL_ *= breakdownFactor_;
}

void OpenHurricane::expertSystemAdaptCFL::divergence() const {
    cCFL_ *= divergenceFactor_;
}

void OpenHurricane::expertSystemAdaptCFL::slowConvergence() const {
    if (countUnchange_ >= interval_) {
        CFLFactor0_ *= updateFactor_;
    }
}

void OpenHurricane::expertSystemAdaptCFL::printCFL(const real cfl) const {
    if (iter_.totalStep() % stepForPrintCFL_ == 0) {
        real cflmin = cfl;
        real cflmax = cfl;
        HurMPI::reduce(cflmin, MPI_MIN);
        HurMPI::reduce(cflmax, MPI_MAX);
        Pout << std::endl
             << "    Adapting CLF Info: The current global CFL is " << cflmin << "." << std::endl
             << std::endl;
    }
}

bool OpenHurricane::expertSystemAdaptCFL::activateAdaptCFL(const integer icount) const {
    if (iter().hasSubIteration()) {
        return icount >= 2;
    } else {
        if (iter_.totalStep() <= initialStage_) {
            return icount >= max(integer(interval_ / 2), 5);
        } else if (iter_.totalStep() <= (initialStage_ + initialStage_ / 2)) {
            return icount >= 3;
        } else if (iter_.totalStep() <= (initialStage_ + initialStage_)) {
            return icount >= 2;
        }
    }
    return true;
}
