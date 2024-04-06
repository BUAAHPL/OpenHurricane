/*!
 * \file CUDAChemistrySource.cpp
 * \brief Main subroutines for chemistry source in CUDA platform.
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

#ifdef CUDA_PARALLEL
#include "CUDAChemistrySource.hpp"
#include "sourceTermCUDAKernel.hpp"

namespace OpenHurricane {
    createClassNameStr(CUDAChemistrySource, "CUDAChemistrySource");
}

void OpenHurricane::CUDAChemistrySource::makeCUDAChemODEs() const {
    if (CUDAChemODEsPtr_ == nullptr) {
        short itpv = 0;
        if (this->isConstPressure()) {
            itpv = 0;
        } else if (this->isConstVolume()) {
            itpv = 1;
        }

        createReactionTable();

        CUDAChemODEsPtr_ = new CUDAChemistrySourceODEs(*reactionsCUDAPtr_, itpv);
    }
}

void OpenHurricane::CUDAChemistrySource::makeEulerSolver() const {
    if (EulerSolverPtr_ == nullptr) {
        const auto &topCont = this->mesh().Iteration().cont();
        const auto &reactCont =
            topCont.subController("flow").subController("mixture").subController("reactions");
        const auto &chemSorCont = reactCont.subController("chemistrySource");

        real atol = 1e-11;
        real rtol = 1e-4;
        integer maxStep = 10000;
        if (chemSorCont.found("solver")) {
            const auto &solCont = chemSorCont.subController("solver");
            atol = solCont.findOrDefault<real>("atol", atol);
            rtol = solCont.findOrDefault<real>("rtol", rtol);
            maxStep = solCont.findOrDefault<integer>("maxstep", maxStep);
        }
        EulerSolverPtr_ = new EulerCUDA(atol, rtol, 10000);
    }
}

void OpenHurricane::CUDAChemistrySource::getYiRhoT(real *hur_restrict hostYiRhoTPtr_) const {
    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp(); ++i) {
            hostYiRhoTPtr_[n * (nsp() + 2) + i] = yi()[i][n];
        }
        hostYiRhoTPtr_[n * (nsp() + 2) + nsp()] = flows_.rho()[n];
        hostYiRhoTPtr_[n * (nsp() + 2) + nsp() + 1] = flows_.T()[n];
    }
}

void OpenHurricane::CUDAChemistrySource::getYiRhoTSlice(real *hur_restrict hostYiRhoTPtr_,
                                                    const integer pos, const integer offset) const {
    for (integer n = 0; n < offset; ++n) {
        for (integer i = 0; i < nsp(); ++i) {
            hostYiRhoTPtr_[n * (nsp() + 2) + i] = yi()[i][n + pos];
        }
        hostYiRhoTPtr_[n * (nsp() + 2) + nsp()] = flows_.rho()[n + pos];
        hostYiRhoTPtr_[n * (nsp() + 2) + nsp() + 1] = flows_.T()[n + pos];
    }
}

OpenHurricane::CUDAChemistrySource::CUDAChemistrySource(flowModel &flows, const controller &cont)
    : chemistrySource(flows, cont), reactionsCUDAPtr_(nullptr), CUDAChemODEsPtr_(nullptr),
      EulerSolverPtr_(nullptr), ODEsSolverType_(ODEsSolverType::Euler),
      nTB_(nsp(), mesh().nCells()), subDt_() {
    if (cont.found("CUDAChemistrySource")) {
        const auto &CUDACont = cont.subController("CUDAChemistrySource");
        if (CUDACont.found("ODEsSolverType")) {
            auto ew = CUDACont.findWord("ODEsSolverType");
            if (ew == "Euler") {
                ODEsSolverType_ = ODEsSolverType::Euler;
            } else if (ew == "MTSEuler") {
                ODEsSolverType_ = ODEsSolverType::MTSEuler;
            } else {
                errorAbortStr(("Unknown string: " + ew + " in " + CUDACont.name()));
            }
        }
    }
    /*reactionsCUDAPtr_ = new CUDAReactions::reactionTable
    (
            mixtures().reactionCUDA(),
            mixtures().speciesCUDA()
    );*/
}

OpenHurricane::CUDAChemistrySource::~CUDAChemistrySource() noexcept {
    HurDelete(CUDAChemODEsPtr_);
    HurDelete(EulerSolverPtr_);
    if (reactionsCUDAPtr_ != nullptr) {
        reactionsCUDAPtr_->destroy();
    }
    HurDelete(reactionsCUDAPtr_);
}

void OpenHurricane::CUDAChemistrySource::calculateSourceTerms(const bool withoutLastSpec) {
    createReactionTable();
    real *hostYiRhoTPtr_ = nullptr;
    size_t sizeOfH_YiRhoT =
        (static_cast<unsigned long long>(nsp()) + 2) * mesh().nCells() * sizeof(real);
    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostYiRhoTPtr_, sizeOfH_YiRhoT, cudaHostAllocDefault));
    getYiRhoT(hostYiRhoTPtr_);
    real *RRi = nullptr;
    size_t sizeOfRRi = static_cast<unsigned long long>(nsp()) * mesh().nCells() * sizeof(real);
    checkCUDAGPUError(cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

    onlyChemistrySource::calcOnlyChemSourceTerms(hostYiRhoTPtr_, *reactionsCUDAPtr_, nsp(), nrc(),
                                                 mesh().nCells(), RRi);

    for (integer i = 0; i < nsp(); ++i) {
        for (integer n = 0; n < mesh().nCells(); ++n) {
            Ri()[i][n] = RRi[n * nsp() + i];
        }
    }

    checkCUDAGPUError(cudaFreeHost(hostYiRhoTPtr_));
    checkCUDAGPUError(cudaFreeHost(RRi));
    destroyReactionTable();
}

void OpenHurricane::CUDAChemistrySource::calculateSourceTermsImp(const bool withoutLastSpec) {
    createReactionTable();
    real *hostYiRhoTPtr_ = nullptr;
    size_t sizeOfH_YiRhoT =
        (static_cast<unsigned long long>(nsp()) + 2) * mesh().nCells() * sizeof(real);
    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostYiRhoTPtr_, sizeOfH_YiRhoT, cudaHostAllocDefault));

    getYiRhoT(hostYiRhoTPtr_);

    real *RRi = nullptr;
    size_t sizeOfRRi = static_cast<unsigned long long>(nsp()) * mesh().nCells() * sizeof(real);
    checkCUDAGPUError(cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

    real *dRRidrhoyi = nullptr;
    size_t sizeOfdRRi =
        (static_cast<unsigned long long>(nsp()) - 1) * mesh().nCells() * sizeof(real);
    checkCUDAGPUError(cudaHostAlloc((void **)&dRRidrhoyi, sizeOfdRRi, cudaHostAllocDefault));

    onlyChemistrySource::calcChemSourceTermsKimDiagonal(hostYiRhoTPtr_, *reactionsCUDAPtr_, nsp(),
                                                        nrc(), mesh().nCells(), RRi, dRRidrhoyi);

    for (integer i = 0; i < nsp(); ++i) {
        for (integer n = 0; n < mesh().nCells(); ++n) {
            Ri()[i][n] = RRi[n * nsp() + i];
            if (i < nsp() - 1) {
                yi_[i].diagSource()[n] = dRRidrhoyi[n * (nsp() - 1) + i];
            }
        }
    }
    checkCUDAGPUError(cudaFreeHost(hostYiRhoTPtr_));
    checkCUDAGPUError(cudaFreeHost(RRi));
    checkCUDAGPUError(cudaFreeHost(dRRidrhoyi));
    destroyReactionTable();
}

void OpenHurricane::CUDAChemistrySource::calculateSourceTermsAsync(real *hur_restrict hostYiRhoTPtr_,
                                                               cu2DArray<cu_real> &dYiRhoT,
                                                               real *hur_restrict RRi,
                                                               cu2DArray<cu_real> &dRi,
                                                               const cudaStreams &streams) {
    getYiRhoT(hostYiRhoTPtr_);
    onlyChemistrySource::calcOnlyChemSourceTermsAsync(hostYiRhoTPtr_, dYiRhoT, *reactionsCUDAPtr_,
                                                      nsp(), nrc(), mesh().nCells(), nTB_, RRi, dRi,
                                                      streams);
}

void OpenHurricane::CUDAChemistrySource::calculateSourceTermsAsyncSlice(
    real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT, real *hur_restrict RRi,
    cu2DArray<cu_real> &dRi, const cudaStreams &streams, const integer pos,
    const integer offset) {
    nTB_.setBlockAndGridSize(nsp(), offset);
    nTB_.setSharedMemSize1<cu_real>(3 * nsp() * nTB_.blockSize().y);

    getYiRhoTSlice(hostYiRhoTPtr_, pos, offset);
    onlyChemistrySource::calcOnlyChemSourceTermsAsync(
        hostYiRhoTPtr_, dYiRhoT, *reactionsCUDAPtr_, nsp(), nrc(), offset, nTB_, RRi, dRi, streams);
}

void OpenHurricane::CUDAChemistrySource::calculateSourceTermsImpAsync(
    real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT, real *hur_restrict RRi,
    cu2DArray<cu_real> &dRi, real *hur_restrict dRidrhoyi, cu2DArray<cu_real> &ddRidrhoyi,
    const cudaStreams &streams) {
    getYiRhoT(hostYiRhoTPtr_);
    onlyChemistrySource::calcChemSourceTermsKimDiagonalAsync(
        hostYiRhoTPtr_, dYiRhoT, *reactionsCUDAPtr_, nsp(), nrc(), mesh().nCells(), nTB_, RRi,
        dRidrhoyi, dRi, ddRidrhoyi, streams);
}

void OpenHurricane::CUDAChemistrySource::calculateSourceTermsImpAsyncSlice(
    real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT, real *hur_restrict RRi,
    cu2DArray<cu_real> &dRi, real *hur_restrict dRidrhoyi, cu2DArray<cu_real> &ddRidrhoyi,
    const cudaStreams &streams, const integer pos, const integer offset) {
    nTB_.setBlockAndGridSize(nsp(), offset);
    nTB_.setSharedMemSize1<cu_real>(3 * nsp() * nTB_.blockSize().y);

    getYiRhoTSlice(hostYiRhoTPtr_, pos, offset);
    onlyChemistrySource::calcChemSourceTermsKimDiagonalAsync(
        hostYiRhoTPtr_, dYiRhoT, *reactionsCUDAPtr_, nsp(), nrc(), offset, nTB_, RRi, dRidrhoyi,
        dRi, ddRidrhoyi, streams);
}

void OpenHurricane::CUDAChemistrySource::calculateSourceTermsImpAsyncHybrid(
    real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT, real *hur_restrict RRi,
    cu2DArray<cu_real> &dRi, cu_float *hur_restrict dRidrhoyi,
    cu2DArray<cu_float> &ddRidrhoyi, const cudaStreams &streams) {

    getYiRhoT(hostYiRhoTPtr_);
    onlyChemistrySource::calcChemSourceTermsKimDiagonalAsyncHybrid(
        hostYiRhoTPtr_, dYiRhoT, *reactionsCUDAPtr_, nsp(), nrc(), mesh().nCells(), nTB_, RRi,
        dRidrhoyi, dRi, ddRidrhoyi, streams);
}

void OpenHurricane::CUDAChemistrySource::calculateSourceTermsImpAsyncHybridSlice(
    real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT, real *hur_restrict RRi,
    cu2DArray<cu_real> &dRi, cu_float *hur_restrict dRidrhoyi,
    cu2DArray<cu_float> &ddRidrhoyi, const cudaStreams &streams, const integer pos,
    const integer offset) {
    nTB_.setBlockAndGridSize(nsp(), offset);
    nTB_.setSharedMemSize1<cu_real>(3 * nsp() * nTB_.blockSize().y);

    getYiRhoTSlice(hostYiRhoTPtr_, pos, offset);
    onlyChemistrySource::calcChemSourceTermsKimDiagonalAsyncHybrid(
        hostYiRhoTPtr_, dYiRhoT, *reactionsCUDAPtr_, nsp(), nrc(), offset, nTB_, RRi, dRidrhoyi,
        dRi, ddRidrhoyi, streams);
}

void OpenHurricane::CUDAChemistrySource::createReactionTable() const {
    if (reactionsCUDAPtr_ == nullptr) {
        reactionsCUDAPtr_ =
            new cuChem::reactionTable(mixtures().reactionCUDA(), mixtures().speciesCUDA());
    }
}

void OpenHurricane::CUDAChemistrySource::createReactionTableAsync(const cudaStreams &streams) const {
    if (reactionsCUDAPtr_ == nullptr) {
        reactionsCUDAPtr_ = new cuChem::reactionTable(mixtures().reactionCUDA(),
                                                      mixtures().speciesCUDA(streams), streams);
    }
}

void OpenHurricane::CUDAChemistrySource::destroyReactionTable() const {
    if (reactionsCUDAPtr_ != nullptr) {
        reactionsCUDAPtr_->destroy();
    }
    HurDelete(reactionsCUDAPtr_);
    mixtures().destroySpeciesCUDA();
}

OpenHurricane::real OpenHurricane::CUDAChemistrySource::solveTest(const real t, const real dT, real &_p,
                                                          real &_T, real &_rho, realArray &yi,
                                                          fileOsstream &fos) {
    p0_ = _p;
    rho0_ = _rho;
    p_ = p0_;
    rho_ = rho0_;
    T_ = _T;
    real leftTime = t;
    real subDt = dT;
    real *Tci = nullptr;
    checkCUDAError(cudaHostAlloc((void **)&Tci, (nsp_ + 1) * sizeof(real), cudaHostAllocDefault));

    for (integer i = 0; i < nsp_; ++i) {
        Tci[i + 1] = rho0_ * yi[i] / therm()[i].Wi();
    }
    Tci[0] = T_;
    /*HurGPU::setSharedMenBankSize();
    HurGPU::setCachePreferShared();*/
    while (leftTime > 0.0) {
        real dt = leftTime;
        if (ODEsSolverType_ == ODEsSolverType::Euler) {
            calcChemODEs(nsp_, CUDAChemODEs(), EulerSolver(), 0, dt, subDt, Tci);
        } else if (ODEsSolverType_ == ODEsSolverType::MTSEuler) {
            calcChemODEsMTS(nsp_, CUDAChemODEs(), EulerSolver(), 0, dt, subDt, Tci);
        }
        leftTime -= dt;
    }

    for (integer i = 0; i < nsp_; ++i) {
        Tci[i + 1] = max(real(0.0), Tci[i + 1]);
    }
    T_ = Tci[0];

    if (reactionFlowType_ == ReactionFlowType::ConstantPressure) {
        rho_ = therm().rho(p_, T_, yi);
    } else {
        p_ = therm().p(rho_, T_, yi);
    }
    real yiSum = 0;
    for (integer i = 0; i < nsp_; ++i) {
        yi[i] = Tci[i + 1] * species().W(i) / max(rho_, tiny);
        yiSum += yi[i];
    }
    if (yiSum != 1) {
        yi /= yiSum;
    }
    _rho = rho_;
    _p = p_;
    _T = T_;

    checkCUDAError(cudaFreeHost(Tci));
    return subDt;
}

OpenHurricane::real OpenHurricane::CUDAChemistrySource::solve(const realArray &timeStep,
                                                      const real dtFactor) {
    if (!isReacting_) {
        return large;
    }
    const auto &T = flows_.T();
    const auto &p = flows_.p();
    const auto &rho = flows_.rho();
    nThreadsAndBlocks &nTBa = nTB_;

    const auto nTotal = nTBa.gridSize().x * nTBa.blockSize().y;

    if (subDt_.size() == 0) {
        subDt_.resize(nTotal, Zero);

        for (integer n = 0; n < mesh().nCells(); ++n) {
            subDt_[n] = dtFactor * timeStep[n];
        }
    }

    real *hostSubDtPtr_ = nullptr;
    size_t sizeOfSubDt = nTotal * sizeof(real);

    checkCUDAGPUError(cudaHostAlloc((void **)&hostSubDtPtr_, sizeOfSubDt, cudaHostAllocDefault));

    /*checkCUDAGPUError
    (
            cudaHostRegister(subDt_.data(), subDt_.size(), cudaHostAllocDefault)
    );*/

    real *hostTYiRhoPtr_ = nullptr;
    size_t sizeOfH_TYiRho = (static_cast<unsigned long long>(nsp()) + 2) * nTotal * sizeof(real);

    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostTYiRhoPtr_, sizeOfH_TYiRho, cudaHostAllocDefault));

    real *hostDtPtr_ = nullptr;
    size_t sizeOfH_dt = nTotal * sizeof(real);
    checkCUDAGPUError(cudaHostAlloc((void **)&hostDtPtr_, sizeOfH_dt, cudaHostAllocDefault));

    for (integer n = 0; n < mesh().nCells(); ++n) {
        hostSubDtPtr_[n] = subDt_[n];
        hostDtPtr_[n] = dtFactor * timeStep[n];
        hostTYiRhoPtr_[n * (nsp() + 2)] = T[n];
        for (integer i = 0; i < nsp(); ++i) {
            hostTYiRhoPtr_[n * (nsp() + 2) + 1 + i] = yi_[i][n];
        }
        hostTYiRhoPtr_[n * (nsp() + 2) + nsp() + 1] = rho[n];
    }
    if (ODEsSolverType_ == ODEsSolverType::Euler) {
        calcChemODEsArray(hostTYiRhoPtr_, nsp(), mesh().nCells(), CUDAChemODEs(), EulerSolver(),
                          nTBa, hostDtPtr_, hostSubDtPtr_, true);
    } else if (ODEsSolverType_ == ODEsSolverType::MTSEuler) {
        calcChemODEsArrayMTS(hostTYiRhoPtr_, nsp(), mesh().nCells(), CUDAChemODEs(), EulerSolver(),
                             nTBa, hostDtPtr_, hostSubDtPtr_, true);
    }

    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp(); ++i) {
            real ci = hostTYiRhoPtr_[n * (nsp() + 2) + 1 + i];
            real ci0 = rho[n] * yi_[i][n] / species_.W(i);
            Ri()[i][n] = (ci - ci0) / (timeStep[n] * dtFactor) * species_.W(i);
        }
        subDt_[n] = hostSubDtPtr_[n];
    }

    checkCUDAError(cudaFreeHost(hostSubDtPtr_));
    checkCUDAError(cudaFreeHost(hostTYiRhoPtr_));
    checkCUDAError(cudaFreeHost(hostDtPtr_));
    //cudaHostUnregister(subDt_.data());
    return large;
}

void OpenHurricane::CUDAChemistrySource::solveUpdate(const realArray &timeStep, const real dtFactor) {
    if (!isReacting_) {
        return;
    }
    auto &T = flows_.T();
    auto &p = flows_.p();
    auto &rho = flows_.rho();
    nThreadsAndBlocks &nTBa = nTB_;

    const auto nTotal = nTBa.gridSize().x * nTBa.blockSize().y;

    if (subDt_.size() == 0) {
        subDt_.resize(nTotal, Zero);

        for (integer n = 0; n < mesh().nCells(); ++n) {
            subDt_[n] = dtFactor * timeStep[n];
        }
    }

    cudaHostRegister(subDt_.data(), subDt_.size(), cudaHostRegisterPortable);
    real *hostTYiRhoPtr_ = nullptr;
    size_t sizeOfH_TYiRho = (static_cast<unsigned long long>(nsp()) + 2) * nTotal * sizeof(real);

    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostTYiRhoPtr_, sizeOfH_TYiRho, cudaHostAllocDefault));

    real *hostDtPtr_ = nullptr;
    size_t sizeOfH_dt = nTotal * sizeof(real);
    checkCUDAGPUError(cudaHostAlloc((void **)&hostDtPtr_, sizeOfH_dt, cudaHostAllocDefault));

    for (integer n = 0; n < mesh().nCells(); ++n) {
        hostDtPtr_[n] = dtFactor * timeStep[n];
        hostTYiRhoPtr_[n * (nsp() + 2)] = T[n];
        for (integer i = 0; i < nsp(); ++i) {
            hostTYiRhoPtr_[n * (nsp() + 2) + 1 + i] = yi_[i][n];
        }
        hostTYiRhoPtr_[n * (nsp() + 2) + nsp() + 1] = rho[n];
    }
    if (ODEsSolverType_ == ODEsSolverType::Euler) {
        calcChemODEsArray(hostTYiRhoPtr_, nsp(), mesh().nCells(), CUDAChemODEs(), EulerSolver(),
                          nTBa, hostDtPtr_, subDt_.data(), false);
    } else if (ODEsSolverType_ == ODEsSolverType::MTSEuler) {
        calcChemODEsArrayMTS(hostTYiRhoPtr_, nsp(), mesh().nCells(), CUDAChemODEs(), EulerSolver(),
                             nTBa, hostDtPtr_, subDt_.data(), false);
    }

    realArray yi0(nsp());
    for (integer n = 0; n < mesh().nCells(); ++n) {
        T[n] = hostTYiRhoPtr_[n * (nsp() + 2)];
        for (integer i = 0; i < nsp(); ++i) {
            yi_[i][n] = hostTYiRhoPtr_[n * (nsp() + 2) + 1 + i];
            yi0[i] = yi_[i][n];
        }
        if (isConstPressure()) {
            rho[n] = therm().rho(p[n], T[n], yi0);
        } else {
            p[n] = therm().p(rho[n], T[n], yi0);
        }
    }

    checkCUDAError(cudaFreeHost(hostTYiRhoPtr_));
    checkCUDAError(cudaFreeHost(hostDtPtr_));
    cudaHostUnregister(subDt_.data());
}

void OpenHurricane::CUDAChemistrySource::solve(cellRealArray &rhoi, cellRealArray &pi,
                                           cellRealArray &Ti, PtrList<cellRealArray> &yii,
                                           const realArray &dt, realArray &subDt,
                                           const real dtFactor) {
    if (!isReacting_) {
        return;
    }

    nThreadsAndBlocks nTBa(nsp(), mesh().nCells());

    cudaHostRegister(subDt.data(), subDt.size(), cudaHostRegisterPortable);
    real *hostTYiRhoPtr_ = nullptr;
    size_t sizeOfH_TYiRho =
        (static_cast<unsigned long long>(nsp()) + 2) * mesh().nCells() * sizeof(real);
    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostTYiRhoPtr_, sizeOfH_TYiRho, cudaHostAllocDefault));

    real *hostDtPtr_ = nullptr;
    size_t sizeOfH_dt = mesh().nCells() * sizeof(real);
    checkCUDAGPUError(cudaHostAlloc((void **)&hostDtPtr_, sizeOfH_dt, cudaHostAllocDefault));

    for (integer n = 0; n < mesh().nCells(); ++n) {
        hostDtPtr_[n] = dtFactor * dt[n];
        hostTYiRhoPtr_[n * (nsp() + 2)] = Ti[n];
        for (integer i = 0; i < nsp(); ++i) {
            hostTYiRhoPtr_[n * (nsp() + 2) + 1 + i] = yii[i][n];
        }
        hostTYiRhoPtr_[n * (nsp() + 2) + nsp() + 1] = rhoi[n];
    }

    calcChemODEsArray(hostTYiRhoPtr_, nsp(), mesh().nCells(), CUDAChemODEs(), EulerSolver(), nTB_,
                      hostDtPtr_, subDt.data());

    realArray yi0(nsp());
    for (integer n = 0; n < mesh().nCells(); ++n) {
        Ti[n] = hostTYiRhoPtr_[n * (nsp() + 2)];
        for (integer i = 0; i < nsp(); ++i) {
            yii[i][n] = hostTYiRhoPtr_[n * (nsp() + 2) + 1 + i];
            yi0[i] = yii[i][n];
        }
        if (isConstPressure()) {
            rhoi[n] = therm().rho(pi[n], Ti[n], yi0);
        } else {
            pi[n] = therm().p(rhoi[n], Ti[n], yi0);
        }
    }

    checkCUDAError(cudaFreeHost(hostTYiRhoPtr_));
    checkCUDAError(cudaFreeHost(hostDtPtr_));
    cudaHostUnregister(subDt.data());
}

void OpenHurricane::CUDAChemistrySource::solve(cellRealArray &rhoi, cellRealArray &pi,
                                           cellRealArray &Ti, PtrList<cellRealArray> &yii,
                                           const realArray &dt, realArray &subDt,
                                           const real dtFactor, const realArray &odeFactor) {
    if (!isReacting_) {
        return;
    }
    HurGPU::setSharedMenBankSize();
    HurGPU::setCachePreferShared();
    cudaHostRegister(subDt.data(), subDt.size(), cudaHostRegisterPortable);
    real *hostTYiRhoPtr_ = nullptr;
    size_t sizeOfH_TYiRho =
        (static_cast<unsigned long long>(nsp()) + 2) * mesh().nCells() * sizeof(real);
    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostTYiRhoPtr_, sizeOfH_TYiRho, cudaHostAllocDefault));

    real *hostDtPtr_ = nullptr;
    size_t sizeOfH_dt = mesh().nCells() * sizeof(real);
    checkCUDAGPUError(cudaHostAlloc((void **)&hostDtPtr_, sizeOfH_dt, cudaHostAllocDefault));

    real *hostOdeFactorPtr_ = nullptr;
    size_t sizeOfH_odeFactor = mesh().nCells() * sizeof(real);
    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostOdeFactorPtr_, sizeOfH_odeFactor, cudaHostAllocDefault));

    for (integer n = 0; n < mesh().nCells(); ++n) {
        hostDtPtr_[n] = dt[n];
        hostOdeFactorPtr_[n] = odeFactor[n];
        hostTYiRhoPtr_[n * (nsp() + 2)] = Ti[n];
        for (integer i = 0; i < nsp(); ++i) {
            hostTYiRhoPtr_[n * (nsp() + 2) + 1 + i] = yii[i][n];
        }
        hostTYiRhoPtr_[n * (nsp() + 2) + nsp() + 1] = rhoi[n];
    }

    calcChemODEsArray(hostTYiRhoPtr_, nsp(), mesh().nCells(), CUDAChemODEs(), EulerSolver(),
                      hostDtPtr_, subDt.data(), hostOdeFactorPtr_);

    realArray yi0(nsp());
    for (integer n = 0; n < mesh().nCells(); ++n) {
        Ti[n] = hostTYiRhoPtr_[n * (nsp() + 2)];
        for (integer i = 0; i < nsp(); ++i) {
            yii[i][n] = hostTYiRhoPtr_[n * (nsp() + 2) + 1 + i];
            yi0[i] = yii[i][n];
        }
        if (isConstPressure()) {
            rhoi[n] = therm().rho(pi[n], Ti[n], yi0);
        } else {
            pi[n] = therm().p(rhoi[n], Ti[n], yi0);
        }
    }

    checkCUDAError(cudaFreeHost(hostTYiRhoPtr_));
    checkCUDAError(cudaFreeHost(hostDtPtr_));
    checkCUDAError(cudaFreeHost(hostOdeFactorPtr_));
    cudaHostUnregister(subDt.data());
    HurDelete(CUDAChemODEsPtr_);
    HurDelete(EulerSolverPtr_);
}

#endif // CUDA_PARALLEL
