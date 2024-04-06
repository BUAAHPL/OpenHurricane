/*!
 * \file laminarSpeciesSolverCUDA.cpp
 * \brief Main subroutines of laminar species transport Solver using CUDA accelerating chemical source terms calculation.
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

#include "laminarSpeciesSolverCUDA.hpp"
#ifdef CUDA_PARALLEL

#include "calculateFieldVar.hpp"
#include "laminar.hpp"
#include "laminarFlow.hpp"
#include "solutionWrite.hpp"
#include "viscousFlux.hpp"

namespace OpenHurricane {
    createClassName(laminarSpeciesSolverCUDA);
    registerObjFty(solver, laminarSpeciesSolverCUDA, controller);
} // namespace OpenHurricane

OpenHurricane::laminarSpeciesSolverCUDA::laminarSpeciesSolverCUDA(iteration &iter,
                                                              const runtimeMesh &mesh)
    : solver(iter, mesh), chemtryPtr_(nullptr), rhoId_(0), rhouId_(1), rhoEId_(4), rhoYi0Id_(5) {

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

    if (HurGPU::GPUDev().size() != 1) {
        LFatal("Only support single GPU mode in this class");
    }
    HurGPU::setDevice(HurGPU::GPUDev()[0]);
}

void OpenHurricane::laminarSpeciesSolverCUDA::solving() {
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

void OpenHurricane::laminarSpeciesSolverCUDA::BDFSolve() {
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

void OpenHurricane::laminarSpeciesSolverCUDA::clear() noexcept {
    invFluxPtr_.clear();
    chemtryPtr_.clear();
}

void OpenHurricane::laminarSpeciesSolverCUDA::bc() {
    const integer nsp = mixtures().species().size();

    rho().updateBoundary();
    v().updateBoundary();
    p().updateBoundary();

    //20210409 ��˼�� �������������yi�߽���µ������¶�T֮ǰ
    for (integer isp = 0; isp < nsp; ++isp) {
        yi()[isp].updateBoundary();
    }
    T().updateBoundary();
    //E().updateBoundary();
    // 20210322 ��˼�� ȡ������solver::bc()����
    // solver::bc()�����еģ�RiemannValue::updated_ = false;������������solver::iterRefresh()
    //solver::bc();
}

void OpenHurricane::laminarSpeciesSolverCUDA::timeStep(realArray &dt) {}

void OpenHurricane::laminarSpeciesSolverCUDA::updateProperties() {
    mixtures().E(p(), T(), v(), E(), false, true);
    mixtures().gammaBoundary(p(), T(), gama(), true);
    mixtures().gasProperties(rho(), p(), T(), mu(), kappal(), cp(), Diff(), hi(), false, true);
    /*mixtures().muKappa(rho(), p(), T(), mu(), kappal(), false, true);
    flowPtr_->mixtures().Diff(p(), T(), Diff(), false, true);
    flowPtr_->mixtures().updateHai(p(), T(), false, true);*/
}

void OpenHurricane::laminarSpeciesSolverCUDA::initialize() {
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

void OpenHurricane::laminarSpeciesSolverCUDA::calculateFc() {
    // The size of species.
    const integer nsp = mixtures().species().size();

    // Caculate the inviscous fluxes for the continuity, momentum and energy equations
    invFluxPtr_->basicFlux();

    // Calculate the convective fluxes for the species equations.
    invFluxPtr_->invFluxSpecies(yi(), false);

    // Calculate the gradient of last species.
    invFluxPtr_->grad(yi()[nsp - 1]);
}

void OpenHurricane::laminarSpeciesSolverCUDA::calculateFv() {}

void OpenHurricane::laminarSpeciesSolverCUDA::actualCalFv() {
    const auto &spT = mixtures().species();

    faceTensorArray deltafV(fv::gradf(v()));

    faceVectorArray Vf(fv::interpolate(v()));
    faceRealArray rhof(fv::interpolate(rho()));

    faceSymmTensorArray tau(
        object("tau", mesh(), object::NOT_WRITE, object::TEMPORARY), mesh(),
        fv::interpolate(mu()) *
            (twoSymm(deltafV) - (real(2.0 / 3.0) * (div(diagToVector(deltafV))) * I)));

    //v().rhs() += visFlux(tau);
    visFlux(tau, v());

    const faceRealArray kappaf(fv::interpolate(kappal()));

    Vf = tau * Vf;
    Vf += (kappaf * fv::gradf(T()));
    faceRealArray Dimf(object("interpolation(Dimf)", mesh(), object::NOT_WRITE, object::TEMPORARY),
                       mesh());

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

void OpenHurricane::laminarSpeciesSolverCUDA::calculateSource() {
    
    const integer nsp = flowPtr_->mixtures().species().size();
    if (timeMarcingPtr_->explicitSource()) {
        if (!mixtures().noReaction()) {
            HurGPU::setSharedMenBankSize();
            HurGPU::setCachePreferShared();
            // First step: Create CUDA stream
            cudaStreams stream;
            stream.create();

            // Second step: Create reaction table in memory of the device (GPU).
            chemtryPtr_->chemistry().createReactionTableAsync(stream);

            // Third step: Alloc memory in the host (CPU) for array "hostYiRhoTPtr" and "RRi".
            //             The pointer "hostYiRhoTPtr" is used for storing yi (species mass fraction), rho (density) and T (temperature).
            //             The pointer "RRi" is used for storing Ri (chemical source terms).
            real *hostYiRhoTPtr = nullptr;
            size_t sizeOfH_YiRhoT =
                (static_cast<unsigned long long>(nsp) + 2) * mesh().nCells() * sizeof(real);
            checkCUDAGPUError(
                cudaHostAlloc((void **)&hostYiRhoTPtr, sizeOfH_YiRhoT, cudaHostAllocDefault));
            real *RRi = nullptr;
            size_t sizeOfRRi =
                static_cast<unsigned long long>(nsp) * mesh().nCells() * sizeof(real);
            checkCUDAGPUError(cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

            // Fourth step: Alloc memory in the device (GPU) for array "dYiRhoT" and "dRi".
            //             The pointer "dYiRhoT" is used for storing yi (species mass fraction), rho (density) and T (temperature).
            //             The pointer "dRi" is used for storing Ri (chemical source terms).
            cu2DArray<cu_real> dYiRhoT(nsp + 2, mesh().nCells());
            cu2DArray<cu_real> dRi(nsp, mesh().nCells());

            // Fifth step: Compute the chemical source terms in the device (GPU)
            chemtryPtr_->chemistry().calculateSourceTermsAsync(hostYiRhoTPtr, dYiRhoT, RRi, dRi,
                                                               stream);

            // Sixth step: Calculate the viscous flux in the host (CPU)
            actualCalFv();

            // Seventh step: Synchronize the computational tasks of host (CPU) and device (GPU)
            stream.synchronize();

            // Eighth step: Destroy memory pointer and get chemical source terms from array RRi.
            chemtryPtr_->chemistry().destroyReactionTable();
            dYiRhoT.clear();
            dRi.clear();
            stream.destroy();
            checkCUDAGPUError(cudaFreeHost(hostYiRhoTPtr));
            for (integer i = 0; i < nsp; ++i) {
                for (integer n = 0; n < mesh().nCells(); ++n) {
                    chemtryPtr_->chemistry().Ri()[i][n] = RRi[n * nsp + i];
                }
            }
            checkCUDAGPUError(cudaFreeHost(RRi));
            chemtryPtr_->getChemistrySource();
        }
    } else if (timeMarcingPtr_->diagonalImpSource()) {
        if (!mixtures().noReaction()) {
            HurGPU::setSharedMenBankSize();
            HurGPU::setCachePreferShared();
            // First step: Create CUDA stream
            cudaStreams stream;
            stream.create();

            // Second step: Create reaction table in memory of the device (GPU).
            chemtryPtr_->chemistry().createReactionTableAsync(stream);

            // Third step: Alloc memory in the host (CPU) for array "hostYiRhoTPtr", "RRi" and "dRRidrhoyi".
            //             The pointer "hostYiRhoTPtr" is used for storing yi (species mass fraction), rho (density) and T (temperature).
            //             The pointer "RRi" is used for storing Ri (chemical source terms).
            //             The pointer "dRRidrhoyi" is used for storing the diagonal elements of the Jacobian matrix of the chemical source terms.
            real *hostYiRhoTPtr = nullptr;
            size_t sizeOfH_YiRhoT =
                (static_cast<unsigned long long>(nsp) + 2) * mesh().nCells() * sizeof(real);
            checkCUDAGPUError(
                cudaHostAlloc((void **)&hostYiRhoTPtr, sizeOfH_YiRhoT, cudaHostAllocDefault));

            real *RRi = nullptr;
            size_t sizeOfRRi =
                static_cast<unsigned long long>(nsp) * mesh().nCells() * sizeof(real);
            checkCUDAGPUError(cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

            real *dRRidrhoyi = nullptr;
            size_t sizeOfdRRi =
                (static_cast<unsigned long long>(nsp) - 1) * mesh().nCells() * sizeof(real);
            checkCUDAGPUError(
                cudaHostAlloc((void **)&dRRidrhoyi, sizeOfdRRi, cudaHostAllocDefault));

            // Fourth step: Alloc memory in the device (GPU) for array "dYiRhoT", "dRi" and "device_dRidrhoyi".
            //             The pointer "dYiRhoT" is used for storing yi (species mass fraction), rho (density) and T (temperature).
            //             The pointer "dRi" is used for storing Ri (chemical source terms).
            //             The pointer "dRi" is used for storing Ri the diagonal elements of the Jacobian matrix of the chemical source terms.
            cu2DArray<cu_real> dYiRhoT(nsp + 2, mesh().nCells());
            cu2DArray<cu_real> dRi(nsp, mesh().nCells());
            cu2DArray<cu_real> device_dRidrhoyi(nsp - 1, mesh().nCells());

            // Fifth step: Compute the chemical source terms in the device (GPU)
            chemtryPtr_->chemistry().calculateSourceTermsImpAsync(
                hostYiRhoTPtr, dYiRhoT, RRi, dRi, dRRidrhoyi, device_dRidrhoyi, stream);

            // Sixth step: Calculate the viscous flux in the host (CPU)
            actualCalFv();

            // Seventh step: Synchronize the computational tasks of host (CPU) and device (GPU)
            stream.synchronize();

            // Eighth step: Destroy memory pointer, get chemical source terms from array RRi and the diagonal elements of the chemical Jacobian matrix from dRRidrhoyi.
            chemtryPtr_->chemistry().destroyReactionTable();
            dYiRhoT.clear();
            dRi.clear();
            device_dRidrhoyi.clear();
            stream.destroy();
            checkCUDAGPUError(cudaFreeHost(hostYiRhoTPtr));
            for (integer i = 0; i < nsp; ++i) {
                for (integer n = 0; n < mesh().nCells(); ++n) {
                    chemtryPtr_->chemistry().Ri()[i][n] = RRi[n * nsp + i];
                    if (i < nsp - 1) {
                        yi()[i].diagSource()[n] = dRRidrhoyi[n * (nsp - 1) + i];
                    }
                }
            }
            checkCUDAGPUError(cudaFreeHost(RRi));
            checkCUDAGPUError(cudaFreeHost(dRRidrhoyi));
            chemtryPtr_->getImpChemistrySource(dt_, !iter().isGlobalTimeStep());
        }
    } else if (timeMarcingPtr_->fullJacobianSource()) {
        LFatal("This solver does not support full Jacobian "
                   "point-implicit operation");
    }

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

void OpenHurricane::laminarSpeciesSolverCUDA::updatePrimitives(const bool shouldUpdateTemp) {
    mixtures().lastSpeAndNormalized();

    solver::updatePrimitives(shouldUpdateTemp);

    limits();
    mixtures().gamma(p(), T(), gama(), true);
    bc();

    updateProperties();
}

void OpenHurricane::laminarSpeciesSolverCUDA::updateFlowOld() {}

void OpenHurricane::laminarSpeciesSolverCUDA::write() {
    if (iter().solWrite().writeNow()) {
        flowPtr_->mixtures().cp(p(), T(), cp(), true);
        if (!mixtures().noReaction()) {
            addSolOutput(iter(), "heatReleaseRate", chemtryPtr_->heatReleaseRate);
            addSolOutput(iter(), "tcFRR", chemtryPtr_->tcFRR);
            addSolOutput(iter(), "tcJacDT", chemtryPtr_->tcJacDT);
            addSolOutput(iter(), "tcSFR", chemtryPtr_->tcSFR);
            addSolOutput(iter(), "tcGSPR", chemtryPtr_->tcGSPR);
        }
    }
    iter_.write();
}

#endif // CUDA_PARALLEL