#include "ESDIRK.hpp"
/*!
 * \file ESDIRK.inl
 * \brief In-Line subroutines of the <i>ESDIRK.hpp</i> file.
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
#pragma once

template <class timeMethod>
inline OpenHurricane::ESDIRK<timeMethod>::ESDIRK(const controller &cont, const runtimeMesh &mesh,
                                                 const flowModel &flowM, solver &_solver,
                                                 cellVectorArray &v)
    : timeMethod(cont, mesh, flowM, _solver, v, true), gamma_(), b_(), a_(), stage_(), k_(0),
      type_(ESDIRK3), rhoOi_(-1), m0_(cont.subController("timeMethod")
                                          .subController("dualTimeStep")
                                          .findOrDefault<integer>("beginStep", 5)),
      iterMethod_(iterationMethod::pseudoTimeStepping),
      newExterpolate_(newTimeExterpolationMethod::lastTime) {
    setESDIRKType(cont.subController("timeMethod"));

    setABG();
    if (cont.subController("timeMethod").subController("dualTimeStep").found("iterationMethod")) {
        string iw = cont.subController("timeMethod")
                        .subController("dualTimeStep")
                        .findWord("iterationMethod");
        stringToUpperCase(iw);
        if (iw == "PSEUDOTIMESTEPPING") {
            iterMethod_ = iterationMethod::pseudoTimeStepping;
        } else if (iw == "NEWTON") {
            iterMethod_ = iterationMethod::NewtonIteration;
        } else {
            errorAbortStr(("Unknown iteration method: " + iw));
        }
    }
    if (cont.subController("timeMethod")
            .subController("dualTimeStep")
            .found("newTimeStepInitialize")) {
        string iw = cont.subController("timeMethod")
                        .subController("dualTimeStep")
                        .findWord("newTimeStepInitialize");
        stringToUpperCase(iw);
        if (iw == "LASTTIME") {
            newExterpolate_ = newTimeExterpolationMethod::lastTime;
        } else if (iw == "LAGRANGEEXTERPOLATION") {
            newExterpolate_ = newTimeExterpolationMethod::Lagrange;
        } else {
            errorAbortStr(("Unknown initialzing method: " + iw));
        }
    }
    setSolverWrite();
}

template <class timeMethod> inline OpenHurricane::ESDIRK<timeMethod>::~ESDIRK() noexcept {}

template <class timeMethod>
inline bool OpenHurricane::ESDIRK<timeMethod>::explicitSource() const noexcept {
    return timeMethod::explicitSource();
}

template <class timeMethod>
inline bool OpenHurricane::ESDIRK<timeMethod>::diagonalImpSource() const noexcept {
    return timeMethod::diagonalImpSource();
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::initializing() {
    timeMethod::initializing();
    for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
        object *ob = timeMethod::objectList_[oi];
        if (ob->name() == "rho") {
            rhoOi_ = oi;
            break;
        }
    }
    for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
        object *ob = timeMethod::objectList_[oi];
        if (ob->nElements() == 1) {
            cellRealArray *f = static_cast<cellRealArray *>(ob);
            (*f).setLastArrayArray(stage_, Zero);
        } else if (ob->nElements() == 3) {
            cellVectorArray *f = static_cast<cellVectorArray *>(ob);
            (*f).setLastArrayArray(stage_, Zero);
        }
    }
    if (newExterpolate_ == newTimeExterpolationMethod::Lagrange) {
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            object *ob = timeMethod::objectList_[oi];
            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f).lastArray().setLastArrayArray(2, Zero);
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f).lastArray().setLastArrayArray(2, Zero);
            }
        }
    }
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::updateOld() {
    // physical time step
    const real deltaT = getPhyTimeStep();

    timeMethod::alphaDT_ = 1.0 / (gamma_ * deltaT);
    cellRealArray *rho = static_cast<cellRealArray *>(timeMethod::objectList_[rhoOi_]);
    if (newExterpolate_ == newTimeExterpolationMethod::Lagrange) {
        for (integer cellI = 0; cellI < timeMethod::mesh().nCells(); ++cellI) {
            for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
                object *ob = timeMethod::objectList_[oi];
                if (oi == rhoOi_) {
                    cellRealArray *f = static_cast<cellRealArray *>(ob);
                    (*f).lastArray().lastArrayArray()[1][cellI] =
                        (*f).lastArray().lastArrayArray()[0][cellI];
                    (*f).lastArray().lastArrayArray()[0][cellI] = (*f).lastArray()[cellI];
                    continue;
                }
                if (ob->nElements() == 1) {
                    cellRealArray *f = static_cast<cellRealArray *>(ob);
                    (*f).lastArray().lastArrayArray()[1][cellI] =
                        (*f).lastArray().lastArrayArray()[0][cellI];
                    (*f).lastArray().lastArrayArray()[0][cellI] =
                        (*rho).lastArray()[cellI] * (*f).lastArray()[cellI];
                } else if (ob->nElements() == 3) {
                    cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                    (*f).lastArray().lastArrayArray()[1][cellI] =
                        (*f).lastArray().lastArrayArray()[0][cellI];
                    (*f).lastArray().lastArrayArray()[0][cellI] =
                        (*rho).lastArray()[cellI] * (*f).lastArray()[cellI];
                }
            }
        }
    }
    // store U^n/rho
    timeMethod::storeOldPrimitiveValue();
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::updateRhs() {
    // Cell volume list
    auto &cV = timeMethod::mesh().cellVolume();
    cellRealArray *rho = static_cast<cellRealArray *>(timeMethod::objectList_[rhoOi_]);
    for (integer cellI = 0; cellI < timeMethod::mesh().nCells(); ++cellI) {
        // physical time step
        const real deltaT = getPhyTimeStep();
        real temp = cV[cellI] / (gamma_ * deltaT);
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            object *ob = timeMethod::objectList_[oi];
            if (oi == rhoOi_) {
                continue;
            }
            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                real aR = Zero;

                for (integer j = 0; j <= k_; j++) {
                    aR += a_[k_][j] * (*f).lastArrayArray()[j][cellI];
                }
                (*f).rhs()[cellI] =
                    aR / gamma_ + temp * ((*rho).lastArray()[cellI] * (*f).lastArray()[cellI] -
                                          (*rho)[cellI] * (*f)[cellI]);
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                vector aR = Zero;

                for (integer j = 0; j <= k_; j++) {
                    aR += a_[k_][j] * (*f).lastArrayArray()[j][cellI];
                }
                (*f).rhs()[cellI] =
                    aR / gamma_ + temp * ((*rho).lastArray()[cellI] * (*f).lastArray()[cellI] -
                                          (*rho)[cellI] * (*f)[cellI]);
            }
        }
    }

    for (integer cellI = 0; cellI < timeMethod::mesh().nCells(); ++cellI) {
        // physical time step
        const real deltaT = getPhyTimeStep();
        real temp = cV[cellI] / (gamma_ * deltaT);
        real aR = Zero;
        for (integer j = 0; j <= k_; j++) {
            aR += a_[k_][j] * (*rho).lastArrayArray()[j][cellI];
        }
        (*rho).rhs()[cellI] = aR / gamma_ + temp * ((*rho).lastArray()[cellI] - (*rho)[cellI]);
    }
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::storeRhs() {
    for (integer cellI = 0; cellI < timeMethod::mesh().nCells(); ++cellI) {
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            object *ob = timeMethod::objectList_[oi];

            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f).lastArrayArray()[k_][cellI] = (*f).rhs()[cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f).lastArrayArray()[k_][cellI] = (*f).rhs()[cellI];
            }
        }
    }
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::marching() {
    stageUpdate();
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::stageUpdate() {
    iteration &iter = const_cast<iteration &>(timeMethod::iter());
    updateOld();
    initialSub();

    k_ = 0;
    storeRhs();

    for (k_ = 1; k_ < stage_; k_++) {
        Pout << " ESDIRK stage: " << k_ + 1 << " of iter: " << iter.cStep() << std::endl;

        while (iter.subLoop()) {
            timeMethod::timeStep();

            // Refresh
            timeMethod::solver_.iterRefresh();

            timeMethod::solver_.calculateFc();
            timeMethod::solver_.calculateSource();
            timeMethod::solver_.calculateFv();
            storeRhs();
            updateRhs();
            setPseudoTimeStepScaleFactor(iter.subIter().cSubStep());
            timeMethod::marching();
            timeMethod::solver_.updatePrimitives();

            iter.writeSubIterResiduals();
        }
        iter.resetSubiter();
    }
}

template <class timeMethod>
hur_nodiscard inline bool OpenHurricane::ESDIRK<timeMethod>::isESDIRKScheme() const noexcept {
    return true;
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::timeStep() {
    updateDyPhyTimeStep();
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::setABG() {
    if (type_ == ESDIRK3) {
        a_.resize(4);
        b_.resize(4);
        a_ = Zero;
        gamma_ = 1767732205903.0 / 4055673282236.0;
        a_[1][0] = 1767732205903.0 / 4055673282236.0;
        a_[1][1] = 1767732205903.0 / 4055673282236.0;
        a_[2][0] = 2746238789719.0 / 10658868560708.0;
        a_[2][1] = -640167445237.0 / 6845629431997.0;
        a_[2][2] = 1767732205903.0 / 4055673282236.0;
        a_[3][0] = 1471266399579.0 / 7840856788654.0;
        a_[3][1] = -4482444167858.0 / 7529755066697.0;
        a_[3][2] = 11266239266428.0 / 11593286722821.0;
        a_[3][3] = 1767732205903.0 / 4055673282236.0;
        b_[0] = a_[3][0];
        b_[1] = a_[3][1];
        b_[2] = a_[3][2];
        b_[3] = a_[3][3];
    } else if (type_ == ESDIRK4) {
        a_.resize(6);
        b_.resize(6);
        a_ = Zero;
        gamma_ = 1.0 / 4.0;
        a_[1][0] = 1.0 / 4.0;
        a_[1][1] = 1.0 / 4.0;
        a_[2][0] = 8611.0 / 62500.0;
        a_[2][1] = -1743.0 / 31250.0;
        a_[2][2] = 1.0 / 4.0;
        a_[3][0] = 5012029.0 / 34652500.0;
        a_[3][1] = -654441.0 / 2922500.0;
        a_[3][2] = 174375.0 / 388108.0;
        a_[3][3] = 1.0 / 4.0;
        a_[4][0] = 15267082809.0 / 155376265600.0;
        a_[4][1] = -171443401.0 / 120774400.0;
        a_[4][2] = 730878875.0 / 902184768.0;
        a_[4][3] = 2285395.0 / 8070912.0;
        a_[4][4] = 1.0 / 4.0;
        a_[5][0] = 82889.0 / 524892.0;
        a_[5][1] = 0.0;
        a_[5][2] = 15625.0 / 83664.0;
        a_[5][3] = 69875.0 / 102672.0;
        a_[5][4] = -2260.0 / 8211.0;
        a_[5][5] = 1.0 / 4.0;
        b_[0] = a_[5][0];
        b_[1] = a_[5][1];
        b_[2] = a_[5][2];
        b_[3] = a_[5][3];
        b_[4] = a_[5][4];
        b_[5] = a_[5][5];
    } else {
        LFatal("Unknown BDF type");
    }
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::initialSub() {
    if (newExterpolate_ == newTimeExterpolationMethod::lastTime) {
        return;
    }
    if (timeMethod::iter().cStep() <= 3) {
        return;
    }
    const auto &dtn = timeMethod::iter().pTStep().lastTimeStep();
    const real Lnp2 = (dtn[1] + dtn[0]) * dtn[0] / ((dtn[2] + dtn[1]) * dtn[2]);
    const real Lnp1 = ((dtn[2] + dtn[1] + dtn[0]) * dtn[0]) / (-dtn[2] * dtn[1]);
    const real Ln = ((dtn[2] + dtn[1] + dtn[0]) * (dtn[1] + dtn[0])) / ((dtn[2] + dtn[1]) * dtn[1]);
    for (integer cellI = 0; cellI < timeMethod::mesh().nCells(); ++cellI) {
        object *rhoOb = timeMethod::objectList_[rhoOi_];
        cellRealArray *rho = static_cast<cellRealArray *>(rhoOb);
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            object *ob = timeMethod::objectList_[oi];

            if (oi == rhoOi_) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f)[cellI] = Ln * (*f).lastArray()[cellI] +
                              Lnp1 * (*f).lastArray().lastArrayArray()[0][cellI] +
                              Lnp2 * (*f).lastArray().lastArrayArray()[1][cellI];

                continue;
            }
            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f)[cellI] = Ln * (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                              Lnp1 * (*f).lastArray().lastArrayArray()[0][cellI] +
                              Lnp2 * (*f).lastArray().lastArrayArray()[1][cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f)[cellI] = Ln * (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                              Lnp1 * (*f).lastArray().lastArrayArray()[0][cellI] +
                              Lnp2 * (*f).lastArray().lastArrayArray()[1][cellI];
            }
        }
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            if (oi == rhoOi_) {
                continue;
            }
            object *ob = timeMethod::objectList_[oi];

            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f)[cellI] /= (*rho)[cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f)[cellI] /= (*rho)[cellI];
            }
        }
    }
}

template <class timeMethod>
inline void
OpenHurricane::ESDIRK<timeMethod>::setPseudoTimeStepScaleFactor(const integer subSteps) {
    if (iterMethod_ == iterationMethod::NewtonIteration) {
        if (subSteps >= m0_) {
            timeMethod::scaleTimeStep_ *= pow(10.0, subSteps - m0_);
        } else {
            timeMethod::scaleTimeStep_ = 1.0;
        }
    }
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::setSolverWrite() {
    LFatal("Must be reproduced in specific class");
}

template <class timeMethod>
inline void OpenHurricane::ESDIRK<timeMethod>::setESDIRKType(const controller &timeMethodCont) {
    LFatal("Must be reproduced in specific class");
}

template <class timeMethod> inline void OpenHurricane::ESDIRK<timeMethod>::updateDyPhyTimeStep() {
    if (timeMethod::iter().pTStep().isDynamicTimeStep()) {
        timeMethod::calcFluxSpectRadius();

        const_cast<iteration &>(timeMethod::iter())
            .pTStep()
            .setTimeStep(
                timeMethod::pseudoTimes().getGlobalTimeStep(timeMethod::iter().pTStep().dyCFL()));
    }
}