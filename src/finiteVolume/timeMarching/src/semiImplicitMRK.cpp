/*!
 * \file semiImplicitMRK.cpp
 * \brief Main subroutines for semiImplicitMRK schemes.
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
#include "semiImplicitMRK.hpp"
#include "solver.hpp"
namespace OpenHurricane {
    createClassNameStr(semiImplicitMRK, "semiImplicitMRK");
    registerObjFty(timeMarching, semiImplicitMRK, controller);
} // namespace OpenHurricane

OpenHurricane::semiImplicitMRK::semiImplicitMRK(const controller &cont, const runtimeMesh &mesh,
                                                const flowModel &flowM, solver &_solver,
                                                cellVectorArray &v)
    : timeMarching(cont, mesh, flowM, _solver, v), alpha_(), sigma_(), rhoOi_(-1),
      ep_(cont.subController("timeMethod")
              .subController("semiImplicitMRK")
              .findOrDefault<real>("ep", 0.5)),
      maxJacStep_(cont.subController("timeMethod")
                      .subController("semiImplicitMRK")
                      .findOrDefault<integer>("maxJacStep", 2)),
      cflRatioForImpResSmoo_(cont.subController("timeMethod")
                                 .subController("semiImplicitMRK")
                                 .findOrDefault<real>("cflRatioForImpResSmoo", 2.0)),
      impResSmooth_(false), dw_(), keyList_() {
    string stageCoefW = cont.subController("timeMethod")
                            .subController("semiImplicitMRK")
                            .findWordOrDefault("MRKStage", "f3");

    setStageCoeffs(stageCoefW);

    const auto &MRKCont = cont.subController("timeMethod").subController("semiImplicitMRK");
    impResSmooth_ = controllerSwitch(MRKCont)("impResSmooth", impResSmooth_);

    pseudoTimes().unsetIsStretchAc();

    solver_.setBDFUnsteadySolver();
    if (diagonalImpSource()) {
        LFatal("The Multi-stage Runge-Kutta dose not support diagonal "
               "source terms Jacobian treatment");
    }
}

OpenHurricane::semiImplicitMRK::~semiImplicitMRK() noexcept {}

void OpenHurricane::semiImplicitMRK::stageUpdate() {
    if (explicitSource()) {
        storeOldPrimitiveValue();
        auto &cV = mesh().cellVolume();
        cellRealArray *rho = static_cast<cellRealArray *>(objectList_[rhoOi_]);
        for (integer stageI = 0; stageI < alpha_.size(); ++stageI) {
            implicitResidualSmoothing();
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                real dtdV = dt_[cellI] / cV[cellI];

                for (integer oi = 0; oi < objectList_.size(); ++oi) {
                    if (oi == rhoOi_) {
                        (*rho)[cellI] =
                            (*rho).lastArray()[cellI] + alpha_[stageI] * dtdV * (*rho).rhs()[cellI];
                        continue;
                    }
                    object *ob = objectList_[oi];

                    if (ob->nElements() == 1) {
                        cellRealArray *f = static_cast<cellRealArray *>(ob);
                        (*f)[cellI] = (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                                      alpha_[stageI] * dtdV * (*f).rhs()[cellI];
                    } else if (ob->nElements() == 3) {
                        cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                        (*f)[cellI] = (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                                      alpha_[stageI] * dtdV * (*f).rhs()[cellI];
                    }
                }
            }
            updatePrimitives();
            if (stageI < alpha_.size() - 1) {
                solver_.updatePrimitives(true);
                solver_.calculateFc();
                solver_.calculateSource();
                solver_.calculateFv();
            }
        }
    } else {
        storeOldPrimitiveValue();
        auto &cV = mesh().cellVolume();
        cellRealArray *rho = static_cast<cellRealArray *>(objectList_[rhoOi_]);
        realSquareMatrix Jac;
        for (integer stageI = 0; stageI < alpha_.size(); ++stageI) {
            implicitResidualSmoothing();
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                real dtdV = dt_[cellI] / cV[cellI];

                if (fullJacobianSource()) {
                    Jac = -Jacobian()[cellI];
                    for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                        Jac(pi, pi) = inv(alpha_[stageI] * dtdV) + Jac(pi, pi);
                    }
                    Jac = inv(Jac);
                }
                dw_.setZero();
                integer count = 0;
                for (integer oi = 0; oi < objectList_.size(); ++oi) {
                    if (oi == rhoOi_) {
                        dw_[count] = ((*rho).rhs()[cellI] + (*rho).lastArray().lastArray()[cellI]);
                        count++;
                        continue;
                    }
                    object *ob = objectList_[oi];

                    if (ob->nElements() == 1) {
                        cellRealArray *f = static_cast<cellRealArray *>(ob);
                        dw_[count] = ((*f).rhs()[cellI] + (*f).lastArray().lastArray()[cellI]);
                        count++;
                    } else if (ob->nElements() == 3) {
                        cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                        dw_[count] =
                            ((*f).rhs()[cellI].x() + (*f).lastArray().lastArray()[cellI].x());
                        count++;
                        dw_[count] =
                            ((*f).rhs()[cellI].y() + (*f).lastArray().lastArray()[cellI].y());
                        count++;
                        dw_[count] =
                            ((*f).rhs()[cellI].z() + (*f).lastArray().lastArray()[cellI].z());
                        count++;
                    } else {
                        errorAbortStr(("Number of components " + toString(ob->nElements()) +
                                       " of parameter not supported yet"));
                    }
                }
                if (fullJacobianSource()) {
                    dw_ = Jac * dw_;
                } 

                count = 0;
                for (integer oi = 0; oi < objectList_.size(); ++oi) {
                    if (oi == rhoOi_) {
                        (*rho)[cellI] = (*rho).lastArray()[cellI] + dw_[count];
                        count++;
                        continue;
                    }
                    object *ob = objectList_[oi];

                    if (ob->nElements() == 1) {
                        cellRealArray *f = static_cast<cellRealArray *>(ob);
                        (*f)[cellI] =
                            (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] + dw_[count];
                        count++;
                    } else if (ob->nElements() == 3) {
                        cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                        (*f)[cellI].x() =
                            (*rho).lastArray()[cellI] * (*f).lastArray()[cellI].x() + dw_[count];
                        count++;
                        (*f)[cellI].y() =
                            (*rho).lastArray()[cellI] * (*f).lastArray()[cellI].y() + dw_[count];
                        count++;
                        (*f)[cellI].z() =
                            (*rho).lastArray()[cellI] * (*f).lastArray()[cellI].z() + dw_[count];
                        count++;
                    } else {
                        errorAbortStr(("Number of components " + toString(ob->nElements()) +
                                       " of parameter not supported yet"));
                    }
                }
            }
            updatePrimitives();
            if (stageI < alpha_.size() - 1) {
                solver_.updatePrimitives(true);
                solver_.calculateFc();
                //solver_.calculateSource();
                solver_.calculateFv();
            }
        }
        if (fullJacobianSource()) {
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                Jacobian()[cellI].setZero();
            }
        }
    }
}

void OpenHurricane::semiImplicitMRK::updatePrimitives() {
    cellRealArray *rho = static_cast<cellRealArray *>(objectList_[rhoOi_]);
    for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
        for (integer oi = 0; oi < objectList_.size(); ++oi) {
            object *ob = objectList_[oi];
            if (oi == rhoOi_) {
                continue;
            }

            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f)[cellI] = (*f)[cellI] / (*rho)[cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f)[cellI] = (*f)[cellI] / (*rho)[cellI];
            } else {
                errorAbortStr(("Number of components " + toString(ob->nElements()) +
                               " of parameter not supported yet"));
            }
        }
    }
}

void OpenHurricane::semiImplicitMRK::restoreConserves() {
    cellRealArray *rho = static_cast<cellRealArray *>(objectList_[rhoOi_]);
    for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
        for (integer oi = 0; oi < objectList_.size(); ++oi) {
            object *ob = objectList_[oi];
            if (oi == rhoOi_) {
                continue;
            }

            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f)[cellI] = (*f)[cellI] * (*rho)[cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f)[cellI] = (*f)[cellI] * (*rho)[cellI];
            }
        }
    }
}

void OpenHurricane::semiImplicitMRK::setStageCoeffs(const string &w) {
    if (w == "f3") {
        alpha_.resize(3);
        alpha_[0] = 0.1481;
        alpha_[1] = 0.400;
        alpha_[2] = 1.0;
        sigma_ = 1.5;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "f4") {
        alpha_.resize(4);
        alpha_[0] = 0.0833;
        alpha_[1] = 0.2069;
        alpha_[2] = 0.4265;
        alpha_[3] = 1.0;
        sigma_ = 2.0;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "f5") {
        alpha_.resize(5);
        alpha_[0] = 0.0533;
        alpha_[1] = 0.1263;
        alpha_[2] = 0.2375;
        alpha_[3] = 0.4414;
        alpha_[4] = 1.0;
        sigma_ = 2.5;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "s3") {
        alpha_.resize(3);
        alpha_[0] = 0.1918;
        alpha_[1] = 0.4929;
        alpha_[2] = 1.0;
        sigma_ = 0.69;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "s4") {
        alpha_.resize(4);
        alpha_[0] = 0.1084;
        alpha_[1] = 0.2602;
        alpha_[2] = 0.5052;
        alpha_[3] = 1.0;
        sigma_ = 0.92;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "s5") {
        alpha_.resize(5);
        alpha_[0] = 0.0695;
        alpha_[1] = 0.1602;
        alpha_[2] = 0.2898;
        alpha_[3] = 0.5060;
        alpha_[4] = 1.0;
        sigma_ = 1.15;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else {
        alpha_.resize(3);
        alpha_[0] = 0.1918;
        alpha_[1] = 0.4929;
        alpha_[2] = 1.0;
        sigma_ = 0.69;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    }
    if (explicitSource()) {
        pseudoTimes().cfl().setCFLMax(min(sigma_, pseudoTimes().cfl().CFLMax()));
    } else {
        pseudoTimes().cfl().setCFLMax(pseudoTimes().cfl().CFLMax());
    }
}

void OpenHurricane::semiImplicitMRK::timeStep() {
    timeMarching::timeStep();
    if (impResSmooth_) {
        for (integer n = 0; n < mesh().nCells(); ++n) {
            dt_[n] *= cflRatioForImpResSmoo_;
        }
    }
}

void OpenHurricane::semiImplicitMRK::initializing() {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];
        if (ob->name() == "rho") {
            rhoOi_ = oi;
            break;
        }
    }
    if (explicitSource()) {
        return;
    }
    timeMarching::initializing();
    dq_.clear();
    dw_.resize(countParam_);

    if (fullJacobianSourceTable()) {
        keyList_.resize(mesh().nCells());
    }
}

void OpenHurricane::semiImplicitMRK::marching() {
    iteration &iters = const_cast<iteration &>(iter());

#ifdef TEST_PROCESS_TIME
    fileName resOut;
    if (iter().cont().found("testProcessTime")) {
        const auto &testCont = iter().cont().subController("testProcessTime");
        string testN = testCont.findWord("fileName");
        trim(testN);
        fileName outFile = testN;
        if (!outFile.isAbsolute()) {
            outFile = iter().outputPath() / outFile;
        }

        resOut = outFile;
    } else {
        const auto &cfn = iter().configName();
        const auto cfnn = cfn.name(true);
        fileName outFile = cfnn + "TestTime.dat";
        outFile = iter().outputPath() / outFile;
        resOut = outFile;
    }

    testProcessTime myTestPT(iter(), resOut);
#endif // TEST_PROCESS_TIME

    timeStep();

    Pout << "    Begin iteration..." << std::endl;

    // Third step: Begin the iteration.
    while (iters.iterating()) {
#ifdef TEST_PROCESS_TIME
        myTestPT.start(iters.cStep());
#endif // TEST_PROCESS_TIME

        // Refresh
        solver_.iterRefresh();
#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("Refresh");
#endif // TEST_PROCESS_TIME

        solver_.previousSource();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("previousSource");
#endif // TEST_PROCESS_TIME

        solver_.calculateFc(); // Firstly, compute the convective flux Fc

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("ConvectiveFlux");
#endif // TEST_PROCESS_TIME

        if (explicitSource()) {
            solver_.calculateSource(); // Secondly, evaluate the source terms S
        } else {
            initializeSource();
            solver_.calculateSource(); // Secondly, evaluate the source terms S
            getSource();
        }
#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("Source");
#endif // TEST_PROCESS_TIME

        solver_.calculateFv(); // Finally, calculate the viscous flux Fv

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("ViscousFlux");
#endif // TEST_PROCESS_TIME

        stageUpdate();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("Marching");
#endif // TEST_PROCESS_TIME

        // Update flow field
        solver_.updatePrimitives(!updatedTemperature());

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("UpdatePrimitives");
#endif // TEST_PROCESS_TIME

        solver_.postSource();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("postSource");
#endif // TEST_PROCESS_TIME

        solver_.calculateOtherTerms();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("OtherTerms");
#endif // TEST_PROCESS_TIME

        if (!explicitSource()) {
            restoreSource();
        }
        // Write solution
        solver_.write();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("Write");
#endif // TEST_PROCESS_TIME

        timeStep();
#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("TimeStep");
#endif // TEST_PROCESS_TIME

#ifdef TEST_PROCESS_TIME
        myTestPT.stop();
#endif // TEST_PROCESS_TIME
    }
}

void OpenHurricane::semiImplicitMRK::implicitResidualSmoothing(cellRealArray &Q) const {
    realArray Rm(Q.mesh().nTotalCells(), Zero); // R_m
    // R_(m-1)
    cellRealArray Rmm1(object(Q.name() + "_Rmm1", Q.mesh(), object::NOT_WRITE, object::TEMPORARY),
                       Q.mesh(), Zero);

    for (integer i = 0; i < mesh().nCells(); ++i) {
        Rmm1[i] = Q.rhs()[i];
    }

    const auto &faces = mesh().faces();
    const auto &cells = mesh().cells();
    for (integer m = 0; m < maxJacStep_; ++m) {
        //fv::transfer(Rmm1, true);
        realTransfer myTransfer(mesh(), Rmm1, false, true);
        myTransfer.transferInit();
        myTransfer.transferring();
        for (integer i = 0; i < mesh().nCells(); ++i) {
            const auto &fl = cells[i].facesList();

            real Rjm = Zero;
            integer Nf = fl.size();
            for (integer fli = 0; fli < fl.size(); ++fli) {
                const integer fi = fl[fli];
                const auto &cl = faces[fi].leftCell();
                const auto &cr = faces[fi].rightCell();

                const integer j = cl + cr - i;

                if (j < mesh().nCells() || (j >= mesh().nCells() && mag(Rmm1[j]) > tiny)) {
                    Nf--;
                    Rjm += Rmm1[j];
                }
            }

            Rm[i] = (Q.rhs()[i] + ep_ * Rjm) / (real(1.0) + ep_ * real(Nf));
        }

        Rmm1 = Rm;
    }

    for (integer i = 0; i < Rm.size(); ++i) {
        Q.rhs()[i] = Rm[i];
    }
}

void OpenHurricane::semiImplicitMRK::implicitResidualSmoothing(cellVectorArray &Q) const {
    vectorArray Rm(mesh().nCells()); // R_m

    // R_(m-1)
    cellVectorArray Rmm1(object(Q.name() + "_Rmm1", Q.mesh(), object::NOT_WRITE, object::TEMPORARY),
                         Q.mesh(), Zero);
    for (integer i = 0; i < mesh().nCells(); ++i) {
        Rmm1[i] = Q.rhs()[i];
    }

    const auto &faces = mesh().faces();
    const auto &cells = mesh().cells();
    for (integer m = 0; m < maxJacStep_; ++m) {
        //fv::transfer(Rmm1, true);
        vectorTransfer myTransfer(mesh(), Rmm1, false, true);
        myTransfer.transferInit();
        myTransfer.transferring();
        for (integer i = 0; i < mesh().nCells(); ++i) {
            const auto &fl = cells[i].facesList();

            vector Rjm = Zero;
            integer Nf = fl.size();
            for (integer fli = 0; fli < fl.size(); ++fli) {
                const integer fi = fl[fli];
                const auto &cl = faces[fi].leftCell();
                const auto &cr = faces[fi].rightCell();

                const integer j = cl + cr - i;

                if (j < mesh().nCells() || (j >= mesh().nCells() && mag(Rmm1[j]) > tiny)) {
                    Nf--;
                    Rjm += Rmm1[j];
                }
            }

            Rm[i] = (Q.rhs()[i] + ep_ * Rjm) / (real(1.0) + ep_ * real(Nf));
        }

        Rmm1 = Rm;
    }

    for (integer i = 0; i < Rm.size(); ++i) {
        Q.rhs()[i] = Rm[i];
    }
}

void OpenHurricane::semiImplicitMRK::implicitResidualSmoothing() {
    if (!impResSmooth_) {
        return;
    }

    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];

        if (ob->nElements() == 1) {
            cellRealArray *Qptr = static_cast<cellRealArray *>(ob);
            implicitResidualSmoothing(*Qptr);
        } else if (ob->nElements() == 3) {
            cellVectorArray *Qptr = static_cast<cellVectorArray *>(ob);
            implicitResidualSmoothing(*Qptr);
        }
    }
}

void OpenHurricane::semiImplicitMRK::initializeSource() {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];

        if (ob->nElements() == 1) {
            cellRealArray *f = static_cast<cellRealArray *>(ob);
            auto &Q = f->lastArray().lastArray();
            auto &R = f->rhs();
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                Q[cellI] = R[cellI];
            }
        } else if (ob->nElements() == 3) {
            cellVectorArray *f = static_cast<cellVectorArray *>(ob);
            auto &Q = f->lastArray().lastArray();
            auto &R = f->rhs();
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                Q[cellI] = R[cellI];
            }
        } else {
            errorAbortStr(("Number of components " + toString(ob->nElements()) +
                           " of parameter not supported yet"));
        }
    }
}

void OpenHurricane::semiImplicitMRK::getSource() {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];

        if (ob->nElements() == 1) {
            cellRealArray *f = static_cast<cellRealArray *>(ob);
            auto &Q = f->lastArray().lastArray();
            auto &R = f->rhs();
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                auto tmp = Q[cellI];
                Q[cellI] = R[cellI] - Q[cellI];
                R[cellI] = tmp;
            }
        } else if (ob->nElements() == 3) {
            cellVectorArray *f = static_cast<cellVectorArray *>(ob);
            auto &Q = f->lastArray().lastArray();
            auto &R = f->rhs();
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                vector tmp = Q[cellI];
                Q[cellI] = R[cellI] - Q[cellI];
                R[cellI] = tmp;
            }
        } else {
            errorAbortStr(("Number of components " + toString(ob->nElements()) +
                           " of parameter not supported yet"));
        }
    }
}

void OpenHurricane::semiImplicitMRK::restoreSource() {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];

        if (ob->nElements() == 1) {
            cellRealArray *f = static_cast<cellRealArray *>(ob);
            auto &Q = f->lastArray().lastArray();
            auto &R = f->rhs();
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                R[cellI] += Q[cellI];
            }
        } else if (ob->nElements() == 3) {
            cellVectorArray *f = static_cast<cellVectorArray *>(ob);
            auto &Q = f->lastArray().lastArray();
            auto &R = f->rhs();
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                R[cellI] += Q[cellI];
            }
        } else {
            errorAbortStr(("Number of components " + toString(ob->nElements()) +
                           " of parameter not supported yet"));
        }
    }
}
