/*!
 * \file laminarSpeciesSolver.cpp
 * \brief Main subroutines of laminar species transport Solver.
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

#include "laminarSpeciesSolver.hpp"
#include "calculateFieldVar.hpp"
#include "laminar.hpp"
#include "laminarFlow.hpp"
#include "solutionWrite.hpp"
#include "viscousFlux.hpp"

namespace OpenHurricane {
    createClassName(laminarSpeciesSolver);
         registerObjFty(solver, laminarSpeciesSolver, controller);
} // namespace OpenHurricane

OpenHurricane::laminarSpeciesSolver::laminarSpeciesSolver(iteration &iter, const runtimeMesh &mesh)
    : solver(iter, mesh), chemtryPtr_(nullptr), rhoId_(0), rhouId_(1), rhoEId_(4), rhoYi0Id_(5)
#ifdef CUDA_PARALLEL
      ,
      fzBoundStartPtr_(nullptr)
#endif //CUDA_PARALLEL
{

    turbPtr_.reset(new laminar(iter.cont(), *flowPtr_));
    if (iter.cont().found("flow")) {
        const auto &flowCont = iter.cont().subController("flow");
        if (!mixtures().noReaction()) {
            if (flowCont.found("mixture")) {
                const auto &reacCont = flowCont.subController("mixture").subController("reactions");
                if (reacCont.subController("combustionModel").findWord("type") == "finiteRate") {
                    chemtryPtr_ = combustionModel::creator(*flowPtr_, reacCont, *turbPtr_);
                } else {
                    const auto w = reacCont.subController("combustionModel").findWord("type");
                    errorAbortStr(("Only can use finiteRate in laminar reacting "
                                   "flows solver in current version but not " +
                                   w));
                }
            } else {
                LFatal("Cannot find mixture setting section in flow section");
            }
        }
    }
}

void OpenHurricane::laminarSpeciesSolver::solving() {
    const integer nsp = mixtures().species().size();

    // First step: Register the primitives parameter in the time marching method.
    marching().addObject(rho());
    marching().addObject(v());
    marching().addObject(E());

    for (integer i = 0; i < nsp - 1; ++i) {
        marching().addObject(yi()[i]);
    }
    rhoId_ = 0;
    rhouId_ = 1;
    rhoEId_ = 4;
    rhoYi0Id_ = 5;

    solver::solving();

    //marching().initializing();

    //// Second step: Initialize the whole flow field
    //initialize();

    //marching().timeStep();
    //// Third step: Begin the iteration.
    //while (iter_.iterating())
    //{
    //    //marching().timeStep();
    //    // 20210322 ��˼�� ����iterRefresh()����
    //    // Refresh
    //    iterRefresh();

    //    // Calculate the right hand side of the governing equations using the following three steps:
    //    calculateFc();// Firstly, compute the convective flux Fc
    //    calculateSource();// Finally, evaluate the source terms S
    //    calculateFv();// Secondly, calculate the viscous flux Fv

    //    // Time marching
    //    marching().marching();

    //    // Update flow field
    //    updatePrimitives();

    //    calculateOtherTerms();

    //    // Write solution
    //    write();
    //    marching().timeStep();
    //}
}

void OpenHurricane::laminarSpeciesSolver::BDFSolve() {
    const integer nsp = mixtures().species().size();

    // First step: Register the primitives parameter in the time marching method.
    marching().addObject(rho());
    marching().addObject(v());
    marching().addObject(E());
    for (integer i = 0; i < nsp - 1; ++i) {
        marching().addObject(yi()[i]);
    }

    rhoId_ = 0;
    rhouId_ = 1;
    rhoEId_ = 4;
    rhoYi0Id_ = 5;
    marching().initializing();

    // Second step: Initialize the whole flow field
    initialize();

    solver::BDFSolve();
    //marching().marching();
}

void OpenHurricane::laminarSpeciesSolver::clear() noexcept {
    flowPtr_.clear();
    invFluxPtr_.clear();
    chemtryPtr_.clear();
}

void OpenHurricane::laminarSpeciesSolver::bc() {
    const integer nsp = mixtures().species().size();
    rho().updateBoundary();
    v().updateBoundary();
    p().updateBoundary();

    for (integer isp = 0; isp < nsp; ++isp) {
        yi()[isp].updateBoundary();
    }

    T().updateBoundary();
    //E().updateBoundary();
    // 20210322 ��˼�� ȡ������solver::bc()����
    // solver::bc()�����еģ�RiemannValue::updated_ = false;������������solver::iterRefresh()
    //solver::bc();
}

void OpenHurricane::laminarSpeciesSolver::timeStep(realArray &dt) {}

void OpenHurricane::laminarSpeciesSolver::updateProperties() {
    mixtures().E(p(), T(), v(), E(), false, true);
    mixtures().gammaBoundary(p(), T(), gama(), true);
    mixtures().gasProperties(rho(), p(), T(), mu(), kappal(), cp(), Diff(), hi(), false, true);
    /*mixtures().muKappa(rho(), p(), T(), mu(), kappal(), false, true);
    flowPtr_->mixtures().Diff(p(), T(), Diff(), false, true);
    flowPtr_->mixtures().updateHai(p(), T(), false, true);*/
}

void OpenHurricane::laminarSpeciesSolver::initialize() {
    solver::initialize();
    const integer nsp = mixtures().species().size();
    if (!iter().restart()) {
        real rhoi = rho().initialize();
        vector vi = v().initialize();
        real pi = p().initialize();
        real Ti = T().initialize();

        real Ei = E().initialize();
        for (integer isp = 0; isp < nsp; ++isp) {
            yi()[isp].initialize();
        }

        for (integer celli = 0; celli < mesh().nCells(); ++celli) {
            rho()[celli] = mixtures().rho(p()[celli], T()[celli], celli);
        }
    } else {
        for (integer isp = 0; isp < nsp; ++isp) {
            yi()[isp] = Zero;
        }
        const_cast<iteration &>(iter_).readRelay();
    }
    if (initPatching()) {
        for (integer celli = 0; celli < mesh().nCells(); ++celli) {
            rho()[celli] = mixtures().rho(p()[celli], T()[celli], celli);
        }
    }

    mixtures().gamma(p(), T(), gama(), true);

    mixtures().E(p(), T(), v(), E(), true, false);
    bc();
    mixtures().gammaBoundary(p(), T(), gama(), true);

    updateProperties();
}

void OpenHurricane::laminarSpeciesSolver::calculateFc() {
    // The size of species.
    const integer nsp = mixtures().species().size();

    // Caculate the inviscous fluxes for the continuity, momentum and energy equations
    invFluxPtr_->basicFlux();

    // Calculate the convective fluxes for the species equations.
    invFluxPtr_->invFluxSpecies(yi(), false);

    // Calculate the gradient of last species.
    invFluxPtr_->grad(yi()[nsp - 1]);
}

void OpenHurricane::laminarSpeciesSolver::calculateFv() {
    const auto &spT = mixtures().species();

    faceTensorArray deltafV(fv::gradf(v()));

    faceVectorArray Vf(fv::interpolate(v()));
    faceRealArray rhof(fv::interpolate(rho()));

    faceSymmTensorArray tau(
        object("tau", mesh(), object::NOT_WRITE, object::TEMPORARY), mesh(),
        //fv::interpolate(mu()) * (twoSymm(deltafV) - real(2.0 / 3.0) * (diag(deltafV)))
        fv::interpolate(mu()) *
            (twoSymm(deltafV) - (real(2.0 / 3.0) * (div(diagToVector(deltafV))) * I)));

    //v().rhs() += visFlux(tau);
    visFlux(tau, v());

    const faceRealArray kappaf(fv::interpolate(kappal()));

    Vf = tau * Vf;
    Vf += (kappaf * fv::gradf(T()));
    faceRealArray Dimf(object("interpolation(Dimf)", mesh(), object::NOT_WRITE, object::TEMPORARY),
                       mesh());
    /* cellRealArray hi
     (
         object
         (
             "hitmp",
             mesh(),
             object::NOT_WRITE,
             object::TEMPORARY
         ),
         mesh()
     );*/

    faceVectorArray Ji(object("Ji", mesh(), object::NOT_WRITE, object::TEMPORARY), mesh());

    
    for (integer i = 0; i < spT.size(); ++i) {
        //const faceRealArray Dimf(fv::interpolate(Diff()[i]));
        fv::interpolate(Diff()[i], Dimf);
        Dimf *= rhof;
        fv::gradf(yi()[i], Ji);
        Ji *= Dimf;
        //Ji = (rhof * Dimf) * fv::gradf(yi()[i]);
        if (i < spT.size() - 1) {
            visFlux<real>(Ji, yi()[i]);
        }
        auto &hif = Dimf;
        fv::interpolate(hi()[i], hif);
        Ji *= hif;
        Vf += Ji;
    }
    //E().rhs() += visFlux<real>(Vf);
    visFlux<real>(Vf, E());
}

void OpenHurricane::laminarSpeciesSolver::calculateSource() {
    
    // To get chemical source
    if (!mixtures().noReaction()) {
        chemtryPtr_->evaluateSource(*timeMarcingPtr_, rhoId_, rhouId_, rhoEId_, rhoYi0Id_);
    }
    // To get other source terms
    if (sorcTermPtr_) {
        const integer nsp = mixtures().species().size();
        sorcTerm().addSourceTerms(rho());
        sorcTerm().addSourceTerms(rho(), v());
        sorcTerm().addSourceTerms(rho(), E());
        for (integer i = 0; i < nsp - 1; ++i) {
            sorcTerm().addSourceTerms(rho(), yi()[i]);
        }
    }
}

void OpenHurricane::laminarSpeciesSolver::updatePrimitives(const bool shouldUpdateTemp) {
    mixtures().lastSpeAndNormalized();

    solver::updatePrimitives(shouldUpdateTemp);

    limits();
#ifdef CUDA_PARALLEL
    if (HurGPU::useGPU()) {
        updatePropertiesCUDA();
    } else {
        mixtures().gamma(p(), T(), gama(), true);
        bc();
        updateProperties();
    }
#else
    mixtures().gamma(p(), T(), gama(), true);
    bc();
    updateProperties();

#endif //CUDA_PARALLEL
}

void OpenHurricane::laminarSpeciesSolver::updateFlowOld() {}

void OpenHurricane::laminarSpeciesSolver::write() {

    if (iter().solWrite().writeNow()) {
        flowPtr_->mixtures().cp(p(), T(), cp(), true);
        if (!mixtures().noReaction()) {
            addSolOutput(iter(), "heatReleaseRate", chemtryPtr_->heatReleaseRate);
            addSolOutput(iter(), "tcFRR", chemtryPtr_->tcFRR);
            addSolOutput(iter(), "tcJacDT", chemtryPtr_->tcJacDT);
            addSolOutput(iter(), "tcSFR", chemtryPtr_->tcSFR);
            addSolOutput(iter(), "tcGSPR", chemtryPtr_->tcGSPR);
            if (iter().cont().subController("iteration").found("cellLoadWeight")) {
                const auto &clwCont =
                    iter().cont().subController("iteration").subController("cellLoadWeight");

                if (clwCont.found("weightTypes")) {
                    const auto wtw = clwCont.findWord("weightTypes");

                    if (wtw == "DAC") {
                        chemtryPtr_->chemistry().setCellLoadWeights();
                    }
                }
            }

            if (iter().solWrite().found("mixtureFractionZ")) {
                chemtryPtr_->calcMixtureFraction();
            }
            if (iter().solWrite().found("mixtureDisspationRate")) {
                iter().solWrite().setOutputField("mixtureDisspationRate",
                                                 chemtryPtr_->mixtureDisspationRate(*invFluxPtr_));
            }
        }
    }
    iter_.write();
}

#ifdef CUDA_PARALLEL
#include "thermoListCUDA.hpp"
#include "transportListCUDA.hpp"

void OpenHurricane::laminarSpeciesSolver::updatePropertiesCUDA() {
    nThreadsAndBlocks nTBForTran_, nTBForHaiCp_;
    nTBForTran_.setBlockAndGridSize(specTable().size(), mesh().nCells());
    nTBForHaiCp_.setBlockAndGridSize(specTable().size(), mesh().nCells());

    nTBForTran_.setSharedMemSize1<real>(2 * specTable().size() * nTBForTran_.blockSize().y);
    nTBForHaiCp_.setSharedMemSize1<real>(specTable().size() * nTBForTran_.blockSize().y);

    // First step: Create CUDA stream
    cudaStreams streamTran;
    streamTran.create();

    // Second step: Create transport table in memory of the device (GPU)
    mixtures().transportCUDA(streamTran);

    // Third step: (1) Alloc memory in the host (CPU) for array "hostYiTPPtr", "hostDimmMuKap" and "hostHai".
    //             The pointer "hostYiTPPtr" is used for storing yi (species mass fraction), T (temperature) and p (pressure).
    //             The pointer "hostDimmMuKap" is used for storing Dim, mu and kappa.
    //             The pointer "hostHai" is used for storing hai.
    //             (2) Alloc memory in the device (GPU) for array "dYiTP", "dDimmMuKap" and "dHai".
    real *hostYiTPPtr = nullptr;
    const integer nsp = flowPtr_->mixtures().species().size();
    size_t sizeOfH_YiTP =
        (static_cast<unsigned long long>(nsp) + 2) * mesh().nCells() * sizeof(real);

    checkCUDAGPUError(cudaHostAlloc((void **)&hostYiTPPtr, sizeOfH_YiTP, cudaHostAllocDefault));

    getYiTP(hostYiTPPtr);
    real *hostDimmMuKap = nullptr;
    size_t sizeOfH_DimmMuKap =
        (static_cast<unsigned long long>(nsp) + 2) * mesh().nCells() * sizeof(real);

    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostDimmMuKap, sizeOfH_DimmMuKap, cudaHostAllocDefault));
    cu2DArray<cu_real> dYiTP(nsp + 2, mesh().nCells());
    cu2DArray<cu_real> dDimmMuKap(nsp + 2, mesh().nCells());

    real *hostHaiCP = nullptr;
    size_t sizeOfH_HaiCP =
        (static_cast<unsigned long long>(nsp + 1)) * mesh().nCells() * sizeof(real);

    checkCUDAGPUError(cudaHostAlloc((void **)&hostHaiCP, sizeOfH_HaiCP, cudaHostAllocDefault));
    cu2DArray<cu_real> dHaiCp(nsp + 1, mesh().nCells());

    // Fourth step: Compute the hai in the device (GPU)
    CUDAThermo::calchaicpAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA().species(),
                               nTBForHaiCp_, nsp, mesh().nCells(), hostHaiCP, dHaiCp, streamTran);

    // Fifth step: Compute the Dim, mu and kappa in the device (GPU)
    CUDATransport::calcTransportPropertiesAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA(),
                                                nTBForTran_, nsp, mesh().nCells(), hostDimmMuKap,
                                                dDimmMuKap, streamTran, false);

    // Sixth step: Calculate boundary conditions in the host (CPU)
    mixtures().gamma(p(), T(), gama(), true);
    bc();
    mixtures().E(p(), T(), v(), E(), false, true);
    mixtures().gammaBoundary(p(), T(), gama(), true);

    // Eighth step: Destroy memory pointer, get result from arrays.
    streamTran.synchronize();

    mixtures().destroySpeciesCUDA();
    mixtures().destroyTransportCUDA();
    dYiTP.clear();
    dHaiCp.clear();
    dDimmMuKap.clear();
    streamTran.destroy();

    checkCUDAGPUError(cudaFreeHost(hostYiTPPtr));
    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp; ++i) {
            mixtures().hi(i)[n] = hostHaiCP[n * (nsp + 1) + i];
        }
        cp()[n] = hostHaiCP[n * (nsp + 1) + nsp];
    }
    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp; ++i) {
            Diff()[i][n] = hostDimmMuKap[n * (nsp + 2) + i];
        }
        mu()[n] = hostDimmMuKap[n * (nsp + 2) + nsp];
        kappal()[n] = hostDimmMuKap[n * (nsp + 2) + nsp + 1];
    }
    checkCUDAGPUError(cudaFreeHost(hostHaiCP));
    checkCUDAGPUError(cudaFreeHost(hostDimmMuKap));

    muKappaDiffBoundary();
}

void OpenHurricane::laminarSpeciesSolver::getYiTP(real *hur_restrict hostYiTPPtr_) const {
    const auto nsp = specTable().size();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp; ++i) {
            hostYiTPPtr_[n * (nsp + 2) + i] = yi()[i][n];
        }
        hostYiTPPtr_[n * (nsp + 2) + nsp] = T()[n];
        hostYiTPPtr_[n * (nsp + 2) + nsp + 1] = p()[n];
    }
}

void OpenHurricane::laminarSpeciesSolver::muKappaDiffBoundary() {
    nThreadsAndBlocks nTBForTran_;
    nTBForTran_.setBlockAndGridSize(specTable().size(), mesh().nCells());
    nTBForTran_.setSharedMemSize1<real>(2 * specTable().size() * nTBForTran_.blockSize().y);

    const integer nsp = flowPtr_->mixtures().species().size();
    nThreadsAndBlocks nTBForTranBnd(nsp, fzBoundStart().last());
    nThreadsAndBlocks nTBForHaiCpBnd(nsp, fzBoundStart().last());
    nTBForTranBnd.setSharedMemSize1<real>(2 * nsp * nTBForTranBnd.blockSize().y);
    nTBForHaiCpBnd.setSharedMemSize1<real>(specTable().size() * nTBForTran_.blockSize().y);

    cudaStreams streamTran;
    streamTran.create();

    mixtures().transportCUDA(streamTran);

    real *hostYiTPPtr = nullptr;
    size_t sizeOfH_YiTP =
        (static_cast<unsigned long long>(nsp) + 2) * fzBoundStart().last() * sizeof(real);

    checkCUDAGPUError(cudaHostAlloc((void **)&hostYiTPPtr, sizeOfH_YiTP, cudaHostAllocDefault));

    getYiTPBoundary(hostYiTPPtr);
    real *hostDimmMuKap = nullptr;
    size_t sizeOfH_DimmMuKap =
        (static_cast<unsigned long long>(nsp) + 2) * fzBoundStart().last() * sizeof(real);

    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostDimmMuKap, sizeOfH_DimmMuKap, cudaHostAllocDefault));
    cu2DArray<cu_real> dYiTP(nsp + 2, fzBoundStart().last());
    cu2DArray<cu_real> dDimmMuKap(nsp + 2, fzBoundStart().last());

    real *hostHaiCp = nullptr;
    size_t sizeOfH_HaiCp =
        (static_cast<unsigned long long>(nsp + 1)) * fzBoundStart().last() * sizeof(real);

    checkCUDAGPUError(cudaHostAlloc((void **)&hostHaiCp, sizeOfH_HaiCp, cudaHostAllocDefault));
    cu2DArray<cu_real> dHaiCp(nsp + 1, fzBoundStart().last());

    CUDAThermo::calchaicpAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA().species(),
                               nTBForHaiCpBnd, nsp, fzBoundStart().last(), hostHaiCp, dHaiCp,
                               streamTran);
    CUDATransport::calcTransportPropertiesAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA(),
                                                nTBForTranBnd, nsp, fzBoundStart().last(),
                                                hostDimmMuKap, dDimmMuKap, streamTran, false);

    transferTranP();

    streamTran.synchronize();
    mixtures().destroySpeciesCUDA();
    mixtures().destroyTransportCUDA();
    dYiTP.clear();
    dHaiCp.clear();
    dDimmMuKap.clear();
    streamTran.destroy();
    checkCUDAGPUError(cudaFreeHost(hostYiTPPtr));

    setExtrapolate(hostDimmMuKap, hostHaiCp);

    checkCUDAGPUError(cudaFreeHost(hostHaiCp));
    checkCUDAGPUError(cudaFreeHost(hostDimmMuKap));
}

void OpenHurricane::laminarSpeciesSolver::makeFzBoundStart() const {
    const auto &fZL = mesh().faceZones();
    integer countFZB = Zero;
    for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
        if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
            !fZL[fzi].isPeriodicShadow()) {
            ++countFZB;
        }
    }
    fzBoundStartPtr_ = new integerList(countFZB + 1);

    countFZB = Zero;
    (*fzBoundStartPtr_)[countFZB] = 0;
    for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
        if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
            !fZL[fzi].isPeriodicShadow()) {
            (*fzBoundStartPtr_)[countFZB + 1] = (*fzBoundStartPtr_)[countFZB] + fZL[fzi].size();
            ++countFZB;
        }
    }
}

void OpenHurricane::laminarSpeciesSolver::getYiTPBoundary(real *hur_restrict hostYiTPPtr_) const {
    const auto &fZL = mesh().faceZones();
    const auto nsp = specTable().size();
    integer countFZB = Zero;

    for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
        if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
            !fZL[fzi].isPeriodicShadow()) {
            const auto tf = fv::interpolate(T(), fzi);
            const auto pf = fv::interpolate(p(), fzi);
            const auto yif = mixtures().Yif(fzi);
            integer cntz = 0;
            for (integer n = fzBoundStart()[countFZB]; n < fzBoundStart()[countFZB + 1]; ++n) {
                for (integer i = 0; i < nsp; ++i) {
                    hostYiTPPtr_[n * (nsp + 2) + i] = yif[cntz][i];
                }
                hostYiTPPtr_[n * (nsp + 2) + nsp] = tf[cntz];
                hostYiTPPtr_[n * (nsp + 2) + nsp + 1] = pf[cntz];
                ++cntz;
            }
            ++countFZB;
        }
    }
}

void OpenHurricane::laminarSpeciesSolver::setExtrapolate(const real *hur_restrict hostDimmMuKap,
                                                     const real *hur_restrict hostHaiCp) {
    const auto nsp = specTable().size();
    const auto &fZL = mesh().faceZones();
    const auto &fs = mesh().faces();
    integer countFZB = Zero;
    for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
        if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
            !fZL[fzi].isPeriodicShadow()) {
            integer cntz = fzBoundStart()[countFZB];
            for (integer fi = fZL[fzi].firstIndex(); fi <= fZL[fzi].lastIndex(); ++fi) {
                const auto &cl = fs[fi].leftCell();
                const auto &cr = fs[fi].rightCell();
                for (integer isp = 0; isp < nsp; ++isp) {
                    mixtures().hi(isp)[cr] =
                        real(2) * hostHaiCp[cntz * (nsp + 1) + isp] - mixtures().hi(isp)[cl];
                    Diff()[isp][cr] =
                        real(2.0) * hostDimmMuKap[cntz * (nsp + 2) + isp] - Diff()[isp][cl];
                }

                cp()[cr] = real(2) * hostHaiCp[cntz * (nsp + 1) + nsp] - cp()[cl];
                mu()[cr] = real(2.0) * hostDimmMuKap[cntz * (nsp + 2) + nsp] - mu()[cl];
                kappal()[cr] = real(2.0) * hostDimmMuKap[cntz * (nsp + 2) + nsp + 1] - kappal()[cl];
                ++cntz;
            }

            ++countFZB;
        }
    }
}

void OpenHurricane::laminarSpeciesSolver::transferTranP() {
    const auto nsp = specTable().size();
    realTransfer myTransfer1(mesh(), mu(), false, true);
    realTransfer myTransfer2(mesh(), kappal(), false, true);
    realTransfer myTransfer3(mesh(), cp(), false, true);

    myTransfer1.transferInit();
    myTransfer2.transferInit();
    myTransfer3.transferInit();
    for (integer i = 0; i < nsp; ++i) {
        realTransfer myTransfer(mesh(), mixtures().hi(i), false, true);
        realTransfer myTransferDiff(mesh(), Diff()[i], false, true);
        myTransfer.transferInit();
        myTransferDiff.transferInit();

        myTransfer.transferring();
        myTransferDiff.transferring();
    }
    myTransfer1.transferring();
    myTransfer2.transferring();
    myTransfer3.transferring();
}
#endif //CUDA_PARALLEL