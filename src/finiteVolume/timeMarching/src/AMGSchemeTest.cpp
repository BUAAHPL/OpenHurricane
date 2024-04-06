/*!
 * \file AMGSchemeTest.cpp
 * \brief Main subroutines for AMG AMGSchemeTest.
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
#include "AMGSchemeTest.hpp"
#include "solver.hpp"
namespace OpenHurricane {
    createClassNameStr(AMGSchemeTest, "AMGSchemeTest");
    registerObjFty(timeMarching, AMGSchemeTest, controller);
} // namespace OpenHurricane

OpenHurricane::AMGSchemeTest::AMGSchemeTest(const controller &cont, const runtimeMesh &mesh,
                                            const flowModel &flowM, solver &_solver,
                                            cellVectorArray &v, const bool modifiedDiagnalTime)
    : timeMarching(cont, mesh, flowM, _solver, v), omega_(cont.subController("timeMethod")
                                                              .subController("AMGSchemeTest")
                                                              .findOrDefault<real>("omega", 1.05)),
      betaT_(cont.subController("timeMethod")
                 .subController("AMGSchemeTest")
                 .findOrDefault<real>("betaT", 0.2)),
      A_(mesh.CRSAdrr()), sparSolver_(nullptr), modifiedDiagnalTime_(modifiedDiagnalTime),
      alphaDT_(0.0), scaleTimeStep_(1.0), nIter_(cont.subController("timeMethod")
                                                     .subController("AMGSchemeTest")
                                                     .findOrDefault<integer>("nIter", 5)) {
    const auto &AMGTeCont = cont.subController("timeMethod").subController("AMGSchemeTest");
    const auto &sparSolCont = AMGTeCont.subController("sparSolvers");
    sparSolver_ = sparseSolver::realBlockSparSolver::creator(A_, sparSolCont);
}

OpenHurricane::AMGSchemeTest::~AMGSchemeTest() noexcept {
    sparSolver_.clear();
}

void OpenHurricane::AMGSchemeTest::update() {
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
            rhon += sigmaT * dq_[n][0];
            rhou += sigmaT * dq_[n][1];
            rhov += sigmaT * dq_[n][2];
            rhow += sigmaT * dq_[n][3];
            rhoE += sigmaT * dq_[n][4];
            real de = dq_[n][4];
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
                rhoq += sigmaT * dq_[n][j];
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
                sumRelCorrec += sqr(sigmaT * dq_[n][0] / rhoOld);
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

void OpenHurricane::AMGSchemeTest::initializing() {
    timeMarching::initializing();
    for (integer i = 0; i < A_.NNZ(); ++i) {
        A_.data()[i].resize(countParam_);
        A_.data()[i].setZero();
    }

    for (auto &e : A_.interfaceA()) {
        e.resize(countParam_);
        e.setZero();
    }
}

void OpenHurricane::AMGSchemeTest::marching() {
    constructOperator();

    const integer nTotalCells = mesh().nTotalCells();
    for (integer ci = 0; ci < nTotalCells; ++ci) {
        dq_[ci] = Zero;
    }
    const integer nCells = mesh().nCells();
    realArrayArray b(nCells);

    for (integer n = 0; n < nCells; ++n) {
        b[n].resize(countParam_, Zero);
        for (integer pi = 0; pi < paramMap_.size(); ++pi) {
            b[n][pi] = vistRhs(pi, n);
        }
    }
    //AMGSolver_.coarsening();
    //AMGSolver_.firstSolveFinetsLevel(dq_, b);

    hrClock myClock;
    sparSolver_->solve(dq_, b, nIter_);
    real eleT = myClock.elapsedClockTime();
    HurMPIBase::reduce(eleT, MPI_MAX);
    Pout.setReal();
    Pout << " sparSolver_Time = " << eleT << std::endl;
    Pout.unsetReal();
    update();
}

void OpenHurricane::AMGSchemeTest::constructOperator() {
    auto vv = A_.data();
    const integer nCells = mesh().nCells();
    const auto &cV = mesh().cellVolume();
    const auto &fA = mesh().faceArea();
    const auto &cC = mesh().cellCentre();
    const auto &faceL = mesh().faces();
    const auto &cellL = mesh().cells();
    if (fullJacobianSourceTable()) {
        A_.setInvAii();
    }
    for (integer n = 0; n < nCells; ++n) {
        real dm;
        dm = cV[n] / (scaleTimeStep_ * dt_[n]);
        if (modifiedDiagnalTime_) {
            dm += alphaDT_ * cV[n];
        }
        auto &Aii = A_.Aii(n);
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
            Aii.setZero();
            Aii.setOnlyDiagFillNonZero();
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                Aii(pi, pi) = dm;
            }
        } else if (diagonalImpSource()) {
            Aii.setZero();
            Aii.setOnlyDiagFillNonZero();
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                Aii(pi, pi) = (dm - vistDiagSource(pi, n));
            }
        } else if (fullJacobianSource()) {
            Aii = -Jacobian()[n];
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                Aii(pi, pi) += dm;
            }
        } else {
            LFatal("Unsupported source term treatment");
        }
    }

    const auto &leftCellIndex = A_.leftCellIndex();
    const auto &rightCellIndex = A_.rightCellIndex();
    for (integer fi = 0; fi < mesh().nInteriorFaces(); ++fi) {

        // index of right cell
        const auto i = A_.col()[leftCellIndex[fi]];

        // index of left cell
        const auto j = A_.col()[rightCellIndex[fi]];

        auto &aji = vv[leftCellIndex[fi]];
        auto &aij = vv[rightCellIndex[fi]];
        aji.setZero();
        aij.setZero();

        real rami;
        real ramj;

        solver_.Ac(i, -fA[fi], aji);
        solver_.Ac(j, fA[fi], aij);
        rami = 0.5 * rai_[fi][1] + rav_[fi][1];
        ramj = 0.5 * rai_[fi][0] + rav_[fi][0];

        aji *= 0.5;
        aij *= 0.5;
        for (integer pi = 0; pi < paramMap_.size(); ++pi) {
            aji(pi, pi) -= rami;
            aij(pi, pi) -= ramj;
        }
    }

    const auto &fzl = mesh().faceZones();
    integer cnt = 0;
    auto &intA = A_.interfaceA();
    for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
        if (fzl[fzi].isCutFace() || fzl[fzi].isPeriodic() || fzl[fzi].isPeriodicShadow()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
                const auto cl = faceL(fi).leftCell();
                const auto cr = faceL(fi).rightCell();

                auto &aij = intA[cnt++];
                real ramj;
                solver_.Ac(cr, -fA[fi], aij);
                ramj = 0.5 * rai_[fi][1] + rav_[fi][1];
                aij *= 0.5;
                for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                    aij(pi, pi) -= ramj;
                }
            }
        }
    }

    if (fullJacobianSource()) {
        for (integer n = 0; n < nCells; ++n) {
            Jacobian()[n].setZero();
        }
    }
}
