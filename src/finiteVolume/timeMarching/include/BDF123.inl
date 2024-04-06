#include "BDF123.hpp"
/*!
 * \file BDF123.inl
 * \brief In-Line subroutines of the <i>BDF123.hpp</i> file.
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
#pragma once

template <class timeMethod>
inline OpenHurricane::BDF123<timeMethod>::BDF123(const controller &cont, const runtimeMesh &mesh,
                                                 const flowModel &flowM, solver &_solver,
                                                 cellVectorArray &v)
    : timeMethod(cont, mesh, flowM, _solver, v, true), type_(BDF2),
      iterMethod_(iterationMethod::pseudoTimeStepping),
      newExterpolate_(newTimeExterpolationMethod::lastTime),
      m0_(cont.subController("timeMethod")
              .subController("dualTimeStep")
              .findOrDefault<integer>("beginStep", 5)),
      storeSize_(1), alphaSize_(1),
      //dtn_(),
      rhoi_(-1), useLinearDtCFL_(false), dtCFLConst_(10), dtCFLLinear_(10), dtCFL0_(0),
      dtFinal_(0) {
    setBDFType(cont.subController("timeMethod"));

    setAlpha(type_);

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
    storeSize_ = alpha_.size();
    alphaSize_ = alpha_.size();
    storeSize_ = max(storeSize_, 4);
    setSolverWrite();

    if (cont.subController("timeMethod").found(className_)) {
        const auto &LUCont = cont.subController("timeMethod").subController(className_);
        controllerSwitch myCont(LUCont);
        useLinearDtCFL_ = myCont("useLinearDtCFL", useLinearDtCFL_);
        if (useLinearDtCFL_ && LUCont.found("linearDtCFL")) {
            Pout << "    Info: using linear time-step or CFL for initial "
                    "stage..."
                 << std::endl;
            const auto &ldcCont = LUCont.subController("linearDtCFL");
            dtCFLConst_ = ldcCont.findOrDefault<integer>("dtCFLConst", dtCFLConst_);
            dtCFLLinear_ = ldcCont.findOrDefault<integer>("dtCFLLinear", dtCFLLinear_);
        }
    }
}

template <class timeMethod> inline OpenHurricane::BDF123<timeMethod>::~BDF123() noexcept {}

template <class timeMethod>
inline bool OpenHurricane::BDF123<timeMethod>::explicitSource() const noexcept {
    return timeMethod::explicitSource();
}

template <class timeMethod>
inline bool OpenHurricane::BDF123<timeMethod>::diagonalImpSource() const noexcept {
    return timeMethod::diagonalImpSource();
}

template <class timeMethod> inline void OpenHurricane::BDF123<timeMethod>::initializing() {
    timeMethod::initializing();

    for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
        object *ob = timeMethod::objectList_[oi];
        if (ob->nElements() == 1) {
            cellRealArray *f = static_cast<cellRealArray *>(ob);
            (*f).setLastArrayArray(storeSize_, Zero);
        } else if (ob->nElements() == 3) {
            cellVectorArray *f = static_cast<cellVectorArray *>(ob);
            (*f).setLastArrayArray(storeSize_, Zero);
        }
    }
}

template <class timeMethod> inline void OpenHurricane::BDF123<timeMethod>::updateOld() {
    if (timeMethod::iter().cStep() <= alphaSize_ - 1) {
        if (timeMethod::iter().restartFromUnsteady()) {
            setAlpha(type_);
        } else if (type_ == BDF2) {
            if (timeMethod::iter().cStep() == 1) {
                setAlpha(BDF1);

            } else {
                setAlpha(type_);
            }
        } else if (type_ == BDF3) {
            if (timeMethod::iter().cStep() == 1) {
                setAlpha(BDF1);
            } else if (timeMethod::iter().cStep() == 2) {
                setAlpha(BDF2);
            } else {
                setAlpha(type_);
            }
        }
    } else {
        setAlpha(type_);
    }

    timeMethod::alphaDT_ = alpha_[0] / getPhyTimeStep();

    for (integer cellI = 0; cellI < timeMethod::mesh().nCells(); ++cellI) {

        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            object *ob = timeMethod::objectList_[oi];
            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f).rhs().rhs()[cellI] = Zero;
                for (integer i = storeSize_ - 1; i > 0; --i) {
                    (*f).lastArrayArray()[i][cellI] = (*f).lastArrayArray()[i - 1][cellI];
                }
                for (integer i = alpha_.size() - 1; i > 0; --i) {
                    (*f).rhs().rhs()[cellI] += alpha_[i] * (*f).lastArrayArray()[i][cellI];
                }
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f).rhs().rhs()[cellI] = Zero;
                for (integer i = storeSize_ - 1; i > 0; --i) {
                    (*f).lastArrayArray()[i][cellI] = (*f).lastArrayArray()[i - 1][cellI];
                }
                for (integer i = alpha_.size() - 1; i > 0; --i) {
                    (*f).rhs().rhs()[cellI] += alpha_[i] * (*f).lastArrayArray()[i][cellI];
                }
            }
        }
    }
}

template <class timeMethod> inline void OpenHurricane::BDF123<timeMethod>::updateRhs() {
    auto &cV = timeMethod::mesh().cellVolume();
    for (integer cellI = 0; cellI < timeMethod::mesh().nCells(); ++cellI) {
        real temp = cV[cellI] / getPhyTimeStep();
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            object *ob = timeMethod::objectList_[oi];
            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f).rhs()[cellI] -=
                    temp * ((*f).rhs().rhs()[cellI] + alpha_[0] * (*f).lastArrayArray()[0][cellI]);
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f).rhs()[cellI] -=
                    temp * ((*f).rhs().rhs()[cellI] + alpha_[0] * (*f).lastArrayArray()[0][cellI]);
            }
        }
    }
}

template <class timeMethod> inline void OpenHurricane::BDF123<timeMethod>::updateWStar() {
    if (rhoi_ == -1) {
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            object *ob = timeMethod::objectList_[oi];
            if (ob->name() == "rho") {
                rhoi_ = oi;
                break;
            }
        }
    }

    for (integer cellI = 0; cellI < timeMethod::mesh().nCells(); ++cellI) {
        object *rhoOb = timeMethod::objectList_[rhoi_];
        cellRealArray *rho = static_cast<cellRealArray *>(rhoOb);
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            if (oi == rhoi_) {
                (*rho).lastArrayArray()[0][cellI] = (*rho)[cellI];
                continue;
            }
            object *ob = timeMethod::objectList_[oi];

            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f).lastArrayArray()[0][cellI] = (*rho)[cellI] * (*f)[cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f).lastArrayArray()[0][cellI] = (*rho)[cellI] * (*f)[cellI];
            }
        }
    }
}

template <class timeMethod> inline void OpenHurricane::BDF123<timeMethod>::marching() {
    iteration &iter = const_cast<iteration &>(timeMethod::iter());

#ifdef TEST_PROCESS_TIME
    fileName resOut;
    if (iter.cont().found("testProcessTime")) {
        const auto &testCont = iter.cont().subController("testProcessTime");
        string testN = testCont.findWord("fileName");
        trim(testN);
        fileName outFile = testN;
        if (!outFile.isAbsolute()) {
            outFile = iter.outputPath() / outFile;
        }

        resOut = outFile;
    } else {
        const auto &cfn = iter.configName();
        const auto cfnn = cfn.name(true);
        fileName outFile = cfnn + "TestTime.dat";
        outFile = iter.outputPath() / outFile;
        resOut = outFile;
    }

    testProcessTime myTestPT(iter, resOut);
#endif // TEST_PROCESS_TIME

    updateWStar();
    updateDyPhyTimeStep();
    while (iter.iterating()) {
#ifdef TEST_PROCESS_TIME
        myTestPT.start(iter.cStep());
#endif // TEST_PROCESS_TIME
        updateOld();
        initialSub();
#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("initialSub");
#endif // TEST_PROCESS_TIME
        timeMethod::solver_.previousSource();
#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("previousSource");
        real timestepAndIterRefresh = 0;
        real calculateFc = 0;
        real calculateSource = 0;
        real calculateFv = 0;
        real marching = 0;
        real UpdatePrimitives = 0;
        real other = 0;
#endif // TEST_PROCESS_TIME
        while (iter.subLoop()) {
            timeMethod::timeStep();

            // Refresh
            timeMethod::solver_.iterRefresh();
#ifdef TEST_PROCESS_TIME
            timestepAndIterRefresh += myTestPT.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
            timeMethod::solver_.calculateFc();
#ifdef TEST_PROCESS_TIME
            calculateFc += myTestPT.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
            timeMethod::solver_.calculateSource();
#ifdef TEST_PROCESS_TIME
            calculateSource += myTestPT.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
            timeMethod::solver_.calculateFv();
#ifdef TEST_PROCESS_TIME
            calculateFv += myTestPT.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
            updateRhs();
            setPseudoTimeStepScaleFactor(iter.subIter().cSubStep());
#ifdef TEST_PROCESS_TIME
            other += myTestPT.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
            timeMethod::marching();
#ifdef TEST_PROCESS_TIME
            marching += myTestPT.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
            updateWStar();
#ifdef TEST_PROCESS_TIME
            other += myTestPT.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
            timeMethod::solver_.updatePrimitives();
#ifdef TEST_PROCESS_TIME
            UpdatePrimitives += myTestPT.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
            iter.writeSubIterResiduals();
#ifdef TEST_PROCESS_TIME
            other += myTestPT.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
        }
#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("Refresh", timestepAndIterRefresh);
        myTestPT.clockTimeIncrement("ConvectiveFlux", calculateFc);
        myTestPT.clockTimeIncrement("Source", calculateSource);
        myTestPT.clockTimeIncrement("ViscousFlux", calculateFv);
        myTestPT.clockTimeIncrement("Marching", marching);
        myTestPT.clockTimeIncrement("UpdatePrimitives", UpdatePrimitives);
        myTestPT.clockTimeIncrement("other", other);
#endif // TEST_PROCESS_TIME
        timeMethod::solver_.postSource();
#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("postSource");
#endif // TEST_PROCESS_TIME
        updateWStar();

#ifdef TEST_PROCESS_TIME
        myTestPT.stop();
#endif // TEST_PROCESS_TIME

        timeMethod::solver_.write();
        updateDyPhyTimeStep();
    }
}

template <class timeMethod>
hur_nodiscard inline bool OpenHurricane::BDF123<timeMethod>::isBDFScheme() const noexcept {
    return true;
}

template <class timeMethod>
inline void OpenHurricane::BDF123<timeMethod>::setAlpha(const BDFType typ) {
    if (typ == BDF1) {
        alpha_.resize(2);
        alpha_[0] = 1.0;
        alpha_[1] = -1.0;
    } else if (typ == BDF2) {
        /**
         * The previous time step [s].
         * lastTimeStep_[0] = dt^n i.e. current time step size.
         * lastTimeStep_[1] = dt^(n-1) i.e. previous time step size.
         * lastTimeStep_[2] = dt^(n-2) i.e.  previous time step size.
         */
        const auto &dtn = timeMethod::iter().pTStep().lastTimeStep();

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
    } else if (typ == BDF3) {
        /**
         * The previous time step [s].
         * lastTimeStep_[0] = dt^n i.e. current time step size.
         * lastTimeStep_[1] = dt^(n-1) i.e. previous time step size.
         * lastTimeStep_[2] = dt^(n-2) i.e.  previous time step size.
         */
        const auto &dtn = timeMethod::iter().pTStep().lastTimeStep();

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
        LFatal("Unknown BDF type");
    }
}

template <class timeMethod> inline void OpenHurricane::BDF123<timeMethod>::initialSub() {
    if (newExterpolate_ == newTimeExterpolationMethod::lastTime) {
        return;
    }
    if (timeMethod::iter().cStep() <= 3) {
        return;
    }
    if (rhoi_ == -1) {
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            object *ob = timeMethod::objectList_[oi];
            if (ob->name() == "rho") {
                rhoi_ = oi;
                break;
            }
        }
    }
    const auto &dtn = timeMethod::iter().pTStep().lastTimeStep();
    const real Lnp2 = (dtn[1] + dtn[0]) * dtn[0] / ((dtn[2] + dtn[1]) * dtn[2]);
    const real Lnp1 = ((dtn[2] + dtn[1] + dtn[0]) * dtn[0]) / (-dtn[2] * dtn[1]);
    const real Ln = ((dtn[2] + dtn[1] + dtn[0]) * (dtn[1] + dtn[0])) / ((dtn[2] + dtn[1]) * dtn[1]);

    for (integer cellI = 0; cellI < timeMethod::mesh().nCells(); ++cellI) {
        object *rhoOb = timeMethod::objectList_[rhoi_];
        cellRealArray *rho = static_cast<cellRealArray *>(rhoOb);
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {

            object *ob = timeMethod::objectList_[oi];

            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f).lastArrayArray()[0][cellI] = Ln * (*f).lastArrayArray()[1][cellI] +
                                                  Lnp1 * (*f).lastArrayArray()[2][cellI] +
                                                  Lnp2 * (*f).lastArrayArray()[3][cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f).lastArrayArray()[0][cellI] = Ln * (*f).lastArrayArray()[1][cellI] +
                                                  Lnp1 * (*f).lastArrayArray()[2][cellI] +
                                                  Lnp2 * (*f).lastArrayArray()[3][cellI];
            }
        }
        for (integer oi = 0; oi < timeMethod::objectList_.size(); ++oi) {
            if (oi == rhoi_) {
                if ((*rho).lastArrayArray()[0][cellI] > 0) {
                    (*rho)[cellI] = (*rho).lastArrayArray()[0][cellI];
                }
                continue;
            }
            object *ob = timeMethod::objectList_[oi];

            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f)[cellI] = (*f).lastArrayArray()[0][cellI] / (*rho).lastArrayArray()[0][cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f)[cellI] = (*f).lastArrayArray()[0][cellI] / (*rho).lastArrayArray()[0][cellI];
            }
        }
    }
    timeMethod::solver_.updatePrimitives(true);
}

template <class timeMethod>
inline void
OpenHurricane::BDF123<timeMethod>::setPseudoTimeStepScaleFactor(const integer subSteps) {
    if (iterMethod_ == iterationMethod::NewtonIteration) {
        if (subSteps >= m0_) {
            timeMethod::scaleTimeStep_ *= pow(10.0, subSteps - m0_);
        } else {
            timeMethod::scaleTimeStep_ = 1.0;
        }
    }
}

template <class timeMethod> inline void OpenHurricane::BDF123<timeMethod>::setSolverWrite() {
    LFatal("Must be reproduced in specific class");
}

template <class timeMethod>
inline void OpenHurricane::BDF123<timeMethod>::setBDFType(const controller &timeMethodCont) {
    LFatal("Must be reproduced in specific class");
}

template <class timeMethod> inline void OpenHurricane::BDF123<timeMethod>::updateDyPhyTimeStep() {
    if (timeMethod::iter().restartFromUnsteady()) {
        if (timeMethod::iter().pTStep().isDynamicTimeStep()) {
            timeMethod::calcFluxSpectRadius();

            const_cast<iteration &>(timeMethod::iter())
                .pTStep()
                .setTimeStep(timeMethod::pseudoTimes().getGlobalTimeStep(
                    timeMethod::iter().pTStep().dyCFL()));
        }
    } else {
        if (useLinearDtCFL_) {
            const auto cstep = timeMethod::iter().cStep();
            if (timeMethod::iter().pTStep().isDynamicTimeStep()) {
                timeMethod::calcFluxSpectRadius();

                if (cstep == 0) {
                    dtCFL0_ = 0.1 * timeMethod::iter().pTStep().dyCFL();
                }

                if (cstep <= dtCFLConst_) {
                    const_cast<iteration &>(timeMethod::iter())
                        .pTStep()
                        .setTimeStep(timeMethod::pseudoTimes().getGlobalTimeStep(dtCFL0_));
                } else if (cstep <= dtCFLLinear_) {
                    real CFL1 =
                        dtCFL0_ + 0.9 * timeMethod::iter().pTStep().dyCFL() *
                                      fabs(real(cstep) - real(dtCFLConst_)) /
                                      max(real(1.0), fabs(real(dtCFLLinear_) - real(dtCFLConst_)));
                    CFL1 = min(CFL1, timeMethod::iter().pTStep().dyCFL());
                    const_cast<iteration &>(timeMethod::iter())
                        .pTStep()
                        .setTimeStep(timeMethod::pseudoTimes().getGlobalTimeStep(CFL1));
                } else {
                    const_cast<iteration &>(timeMethod::iter())
                        .pTStep()
                        .setTimeStep(timeMethod::pseudoTimes().getGlobalTimeStep(
                            timeMethod::iter().pTStep().dyCFL()));
                }
            } else {
                if (cstep == 0) {
                    dtCFL0_ = 0.1 * timeMethod::iter().pTStep().pTimeStep();
                    dtFinal_ = timeMethod::iter().pTStep().pTimeStep();
                }

                if (cstep <= dtCFLConst_) {
                    const_cast<iteration &>(timeMethod::iter()).pTStep().setTimeStep(dtCFL0_);
                } else if (cstep <= dtCFLLinear_) {
                    real dt0 =
                        dtCFL0_ + 0.9 * dtFinal_ * fabs(real(cstep) - real(dtCFLConst_)) /
                                      max(real(1.0), fabs(real(dtCFLLinear_) - real(dtCFLConst_)));
                    dt0 = min(dt0, dtFinal_);
                    const_cast<iteration &>(timeMethod::iter()).pTStep().setTimeStep(dt0);
                } else if (cstep <= dtCFLLinear_ + 1) {
                    const_cast<iteration &>(timeMethod::iter()).pTStep().setTimeStep(dtFinal_);
                }
            }
        }
    }
}
