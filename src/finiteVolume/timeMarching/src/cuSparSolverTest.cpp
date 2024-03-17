/*!
 * \file cuSparSolverTest.cpp
 * \brief Main subroutines for cuSparSolverTest.
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
#include "cuSparSolverTest.hpp"
#include "cudaStreams.hpp"
#include "solver.hpp"

#ifdef CUDA_PARALLEL
namespace OpenHurricane {
    createClassNameStr(cuSparSolverTest, "cuSparSolverTest");
}
namespace OpenHurricane {
    registerObjFty(timeMarching, cuSparSolverTest, controller);
}

OpenHurricane::cuSparSolverTest::cuSparSolverTest(const controller &cont, const runtimeMesh &mesh,
                                                  const flowModel &flowM, solver &_solver,
                                                  cellVectorArray &v,
                                                  const bool modifiedDiagnalTime)
    : timeMarching(cont, mesh, flowM, _solver, v),
      sparSolCont_(cont.subController("timeMethod")
                       .subController("cuSparSolverTest")
                       .subController("sparSolvers")),
      omega_(cont.subController("timeMethod")
                 .subController("cuSparSolverTest")
                 .findOrDefault<real>("omega", 1.05)),
      betaT_(cont.subController("timeMethod")
                 .subController("cuSparSolverTest")
                 .findOrDefault<real>("betaT", 0.2)),
      APtr_(nullptr), dqPtr_(nullptr), sparSolver_(nullptr),
      modifiedDiagnalTime_(modifiedDiagnalTime), alphaDT_(0.0), scaleTimeStep_(1.0),
      nIter_(cont.subController("timeMethod")
                 .subController("cuSparSolverTest")
                 .findOrDefault<integer>("nIter", 5)),
      isReserveGPUTmpMem_(false) {
    const auto &sscont = cont.subController("timeMethod").subController("cuSparSolverTest");
    controllerSwitch ssconts(sscont);
    isReserveGPUTmpMem_ = ssconts("isReserveGPUTmpMem", isReserveGPUTmpMem_);
}

OpenHurricane::cuSparSolverTest::~cuSparSolverTest() noexcept {
    sparSolver_.clear();
}

void OpenHurricane::cuSparSolverTest::update() {
    const real coef1 = -0.20;
    const real coef2 = 2.0;
    realArray qold(countParam());
    auto &rho = const_cast<flowModel &>(flow_).rho();
    auto &p = const_cast<flowModel &>(flow_).p();
    auto &T = const_cast<flowModel &>(flow_).T();
    auto &E = const_cast<flowModel &>(flow_).E();
    auto &iflag = const_cast<flowModel &>(flow_).temperatureFlag();
    real sumRelCorrec = Zero;

    integer countNan = 0;

    for (integer n = 0; n < mesh().nCells(); ++n) {
        qold[0] = rho[n];
        qold[1] = v_[n][0];
        qold[2] = v_[n][1];
        qold[3] = v_[n][2];
        qold[4] = E[n];
        for (integer j = 5; j < countParam(); j++) {
            qold[j] = vistParam(j, n);
        }
        real sigmaT = 1.0;
        integer countMT = 0;

        while (true) {
            real rhon = qold[0];
            real rhoOld = rhon;
            real rhou = rhon * qold[1];
            real rhov = rhon * qold[2];
            real rhow = rhon * qold[3];
            real rhoE = rhon * qold[4];
            /**\brief New value of conservative variables.*/
            rhon += sigmaT * (*dqPtr_)[n](0, 0);
            rhou += sigmaT * (*dqPtr_)[n](0, 1);
            rhov += sigmaT * (*dqPtr_)[n](0, 2);
            rhow += sigmaT * (*dqPtr_)[n](0, 3);
            rhoE += sigmaT * (*dqPtr_)[n](0, 4);
            real de = (*dqPtr_)[n](0, 4);
            if (rhon < 0.0) {
                pseudoTimes().cfl().setBreakdowns();
                break;
            }
            rho[n] = rhon;
            v_[n][0] = rhou / rhon;
            v_[n][1] = rhov / rhon;
            v_[n][2] = rhow / rhon;
            E[n] = rhoE / rhon;

            for (integer j = 5; j < countParam(); j++) {
                real qOld = qold[j];
                real rhoq = rhoOld * qOld;
                rhoq += sigmaT * (*dqPtr_)[n](0, j);
                real qnew = rhoq / rhon;
                if (isnan(qnew)) {
                    qnew = qOld;
                    countNan++;
                }
                //qnew = max(min(qnew, real((1.0 + 0.05)) * qold[j]), real((1.0 - 0.05)) * qold[j]);
                pushParam(j, n, qnew);
            }

            // New temperature
            real newT =
                flow_.mixtures().thermalTable().TEa_rho((E[n] - real(0.5) * v_[n].magSqr()), rho[n],
                                                        T[n], iflag[n], flow_.mixtures().Yi(), n);
            real Td = T[n];
            real deltaT = mag(Td - newT);
            if (deltaT > betaT_ * Td && countMT < 2) {
                countMT++;
                sigmaT = (betaT_ * Td / deltaT) * sigmaT;
            } else {
                T[n] = newT;
                sumRelCorrec += sqr(sigmaT * (*dqPtr_)[n](0, 0) / rhoOld);
                break;
            }
        }

        if (iflag[n] == 0 || iflag[n] == 2) {
            E[n] = (flow_.mixtures().thermalTable().ea_rho(rho[n], T[n], flow_.mixtures().Yi(), n) +
                    real(0.5) * v_[n].magSqr());
        }
    }
    HurMPI::allReduce(sumRelCorrec, MPI_SUM);
    pseudoTimes().cfl().setRelativeCorrection(
        sqrt(sumRelCorrec / real(mesh().allMeshCellNumber())));

    HurMPI::allReduce(countNan, MPI_SUM);
    if (countNan != 0) {
        if (countNan > 50) {
            errorAbortStr(("Too many NaN: " + toString(countNan)));
        } else {
            checkWarningStr(("NaN occurs in results: " + toString(countNan)));
        }
    }
}

void OpenHurricane::cuSparSolverTest::initializing() {
    timeMarching::initializing();
    dq_.clear();
    if (isReserveGPUTmpMem_) {
        auto &csradd = mesh().CRSAdrr();
        APtr_.reset(new cuCRSBlockMatrix<real>(csradd.nRows(), csradd.nColumns(), csradd.NNZ(),
                                               csradd.rowPtr().data(), csradd.col().data(),
                                               countParam()));
        sparSolver_ = sparseSolver::cuRealBlockSparSolver::creator(*APtr_, sparSolCont_);
        dqPtr_.reset(new cuBlockArray1<cu_real, cu_integer>(csradd.nRows(), 1, countParam()));
        APtr_->buildAddressingInDevice();
    }
}

void OpenHurricane::cuSparSolverTest::marching() {
    int deviceId = HurGPU::GPUDev().first();
    HurGPU::setDevice(deviceId);
    cudaStreams st;
    st.create();
    if (!isReserveGPUTmpMem_) {
        auto &csradd = mesh().CRSAdrr();
        APtr_.reset(new cuCRSBlockMatrix<real>(csradd.nRows(), csradd.nColumns(), csradd.NNZ(),
                                               csradd.rowPtr().data(), csradd.col().data(),
                                               countParam()));
        sparSolver_ = sparseSolver::cuRealBlockSparSolver::creator(*APtr_, sparSolCont_);
        dqPtr_.reset(new cuBlockArray1<cu_real, cu_integer>(csradd.nRows(), 1, countParam()));
        APtr_->buildAddressingInDeviceAsync(st());
    }
    constructOperator();

    hrClock myClock;

    APtr_->copyFromHostAsync(st());

    const integer nCells = mesh().nCells();
    cuBlockArray1<cu_real, cu_integer> b(nCells, 1, countParam());

    for (integer n = 0; n < nCells; ++n) {
        for (integer pi = 0; pi < paramMap_.size(); ++pi) {
            b[n](0, pi) = vistRhs(pi, n);
            //(*dqPtr_)[n](0, pi) = 0;
        }
    }
    b.copyFromHostAsync(st());
    //(*dqPtr_).copyFromHostAsync(st());
    sparSolver_->solve(*dqPtr_, b, nIter_, deviceId, st());
    st.destroy();

    real eleT = myClock.elapsedClockTime();
    HurMPIBase::reduce(eleT, MPI_MAX);
    Pout.setReal();
    Pout << " gpusparSolver_Time = " << eleT << std::endl;
    Pout.unsetReal();

    update();

    // Free GPU memory
    destroyCuArray(b);
    if (!isReserveGPUTmpMem_) {
        destroyCuArray((*APtr_));
        destroyCuArray((*dqPtr_));
        sparSolver_.clear();
        APtr_.clear();
        dqPtr_.clear();
    }
}

void OpenHurricane::cuSparSolverTest::constructOperator() {
    auto &vv = APtr_->value();
    auto &diagId = APtr_->diagIndex();
    const integer nCells = mesh().nCells();
    const auto &cV = mesh().cellVolume();
    const auto &fA = mesh().faceArea();
    const auto &cC = mesh().cellCentre();
    const auto &faceL = mesh().faces();
    const auto &cellL = mesh().cells();

    for (integer n = 0; n < nCells; ++n) {
        real dm;
        dm = cV[n] / (scaleTimeStep_ * dt_[n]);
        if (modifiedDiagnalTime_) {
            dm += alphaDT_ * cV[n];
        }
        auto id = diagId[n];
        auto Aii = vv[id];
        for (integer fi = 0; fi < cellL[n].faceSize(); fi++) {
            const integer faceI = cellL[n].facei(fi);
            const auto &cl = faceL[faceI].leftCell();
            const auto &cr = faceL[faceI].rightCell();

            //const integer m = faceL[faceI].oppositeCell(n);
            const integer m = (n == cl) ? cr : cl;

            if (n == cl) {
                dm += 0.5 * (omega_ * rai_[faceI][0]) + rav_[faceI][0];
            } else {
                dm += 0.5 * (omega_ * rai_[faceI][1]) + rav_[faceI][1];
            }
        }

        if (explicitSource()) {
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                Aii(pi, pi) = dm;
            }
        } else if (diagonalImpSource()) {
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                Aii(pi, pi) = (dm - vistDiagSource(pi, n));
            }
        } else if (fullJacobianSource()) {
            LFatal("Unsupported source term treatment");
        } else if (fullJacobianSourceTable()) {
            LFatal("Unsupported source term treatment");
        } else {
            LFatal("Unsupported source term treatment");
        }
    }

    const auto &leftCellIndex = mesh().CRSAdrr().leftCellIndex();
    const auto &rightCellIndex = mesh().CRSAdrr().rightCellIndex();
    realSquareMatrix ACji(paramMap_.size());
    realSquareMatrix ACij(paramMap_.size());
    for (integer fi = 0; fi < mesh().nInteriorFaces(); ++fi) {

        // index of right cell
        const auto i = mesh().CRSAdrr().col()[leftCellIndex[fi]];

        // index of left cell
        const auto j = mesh().CRSAdrr().col()[rightCellIndex[fi]];

        auto aji = vv[leftCellIndex[fi]];
        auto aij = vv[rightCellIndex[fi]];

        real rami;
        real ramj;
        ACji.setZero();
        ACij.setZero();
        solver_.Ac(i, -fA[fi], ACji);
        solver_.Ac(j, fA[fi], ACij);
        rami = 0.5 * rai_[fi][1] + rav_[fi][1];
        ramj = 0.5 * rai_[fi][0] + rav_[fi][0];

        ACji *= 0.5;
        ACij *= 0.5;
        for (integer pi = 0; pi < paramMap_.size(); ++pi) {
            for (integer pj = 0; pj < paramMap_.size(); ++pj) {
                aji(pi, pj) = ACji(pi, pj);
                aij(pi, pj) = ACij(pi, pj);
            }
        }
        for (integer pi = 0; pi < paramMap_.size(); ++pi) {
            aji(pi, pi) -= rami;
            aij(pi, pi) -= ramj;
        }
    }
}
#endif