/*!
 * \file laminarSpeciesSolverCUDATran.cpp
 * \brief Main subroutines of laminar species transport Solver using CUDA accelerating chemical source terms calculation.
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
#ifdef CUDA_PARALLEL
#include "laminarSpeciesSolverCUDATran.hpp"
#include "calculateFieldVar.hpp"
#include "laminar.hpp"
#include "laminarFlow.hpp"
#include "solutionWrite.hpp"
#include "viscousFlux.hpp"

#ifdef TEST_PROCESS_TIME
#include "cudaEvents.hpp"
#endif // TEST_PROCESS_TIME

namespace OpenHurricane {
    createClassName(laminarSpeciesSolverCUDATran);
    registerObjFty(solver, laminarSpeciesSolverCUDATran, controller);
} // namespace OpenHurricane

OpenHurricane::laminarSpeciesSolverCUDATran::laminarSpeciesSolverCUDATran(iteration &iter,
                                                                      const runtimeMesh &mesh)
    : solver(iter, mesh), chemtryPtr_(nullptr), rhoId_(0), rhouId_(1), rhoEId_(4), rhoYi0Id_(5),
      nTBForTran_(), nTBForHaiCp_(), fzBoundStartPtr_(nullptr) {

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
    nTBForTran_.setBlockAndGridSize(specTable().size(), mesh.nCells());
    nTBForHaiCp_.setBlockAndGridSize(specTable().size(), mesh.nCells());

    nTBForTran_.setSharedMemSize1<real>(2 * specTable().size() * nTBForTran_.blockSize().y);

    if (HurGPU::GPUDev().size() != 1) {
        LFatal("Only support single GPU mode in this class");
    }
    HurGPU::setDevice(HurGPU::GPUDev()[0]);
}

void OpenHurricane::laminarSpeciesSolverCUDATran::solving() {
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

void OpenHurricane::laminarSpeciesSolverCUDATran::BDFSolve() {
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

void OpenHurricane::laminarSpeciesSolverCUDATran::clear() noexcept {
    invFluxPtr_.clear();
    chemtryPtr_.clear();
    fzBoundStartPtr_.clear();
}

void OpenHurricane::laminarSpeciesSolverCUDATran::bc() {
    const integer nsp = mixtures().species().size();
    rho().updateBoundary();
    v().updateBoundary();
    p().updateBoundary();

    //20210409 ��˼�� �������������yi�߽���µ������¶�T֮ǰ
    for (integer isp = 0; isp < nsp; ++isp) {
        yi()[isp].updateBoundary();
    }

    T().updateBoundary();
}

void OpenHurricane::laminarSpeciesSolverCUDATran::timeStep(realArray &dt) {}

void OpenHurricane::laminarSpeciesSolverCUDATran::updateProperties() {
    /*HurGPU::setSharedMenBankSize();
    HurGPU::setCachePreferShared();*/
    // First step: Create CUDA stream
    cudaStreams streamTran;
    streamTran.create();

#ifdef TEST_PROCESS_TIME
    HurMPI::barrier();
    cudaEvents myStart;
    cudaEvents myEnd;
    myStart.create();
    myEnd.create();
    myStart.record(streamTran);
    hrClock recordGPU;
    integer ava, tot;
    HurGPU::getMemInfo(ava, tot);
    Pout("    Info: available GPU memory: %d MB, total GPU memory: %d MB\n", ava, tot);
#endif // TEST_PROCESS_TIME
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

    real *hostHai = nullptr;
    size_t sizeOfH_Hai = (static_cast<unsigned long long>(nsp)) * mesh().nCells() * sizeof(real);

    checkCUDAGPUError(cudaHostAlloc((void **)&hostHai, sizeOfH_Hai, cudaHostAllocDefault));
#ifdef TEST_PROCESS_TIME
    integer allUseMem =
        static_cast<integer>((sizeOfH_YiTP + sizeOfH_DimmMuKap + sizeOfH_Hai) / 1024 / 1024);
    HurMPIBase::reduce(allUseMem, MPI_SUM, HurGPU::sameGPUCommtor().masterNoInSub(),
                       HurGPU::sameGPUCommtor());
    if (HurGPU::sameGPUCommtor().isMasterInSub()) {
        printf("    Info: At least used GPU memory: %d MB for gas properties\n", allUseMem);
    }
#endif // TEST_PROCESS_TIME

    cu2DArray<cu_real> dYiTP(nsp + 2, mesh().nCells());
    cu2DArray<cu_real> dDimmMuKap(nsp + 2, mesh().nCells());
    cu2DArray<cu_real> dHai(nsp, mesh().nCells());

#ifdef TEST_PROCESS_TIME
    auto recordGPUBeforeCalc = recordGPU.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

    // Fourth step: Compute the hai in the device (GPU)
    CUDAThermo::calchaiAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA().species(), nTBForHaiCp_,
                             nsp, mesh().nCells(), hostHai, dHai, streamTran);

    // Fifth step: Compute the Dim, mu and kappa in the device (GPU)
    CUDATransport::calcTransportPropertiesAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA(),
                                                nTBForTran_, nsp, mesh().nCells(), hostDimmMuKap,
                                                dDimmMuKap, streamTran, false);

#ifdef TEST_PROCESS_TIME
    myEnd.record(streamTran);
    hrClock cpuStart;
#endif // TEST_PROCESS_TIME

    // Sixth step: Calculate boundary conditions in the host (CPU)
    mixtures().gamma(p(), T(), gama(), true);
    bc();
    mixtures().E(p(), T(), v(), E(), false, true);
    mixtures().gammaBoundary(p(), T(), gama(), true);

#ifdef TEST_PROCESS_TIME
    auto cpuTime1 = cpuStart.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

    // Eighth step: Destroy memory pointer, get result from arrays.
    streamTran.synchronize();

#ifdef TEST_PROCESS_TIME
    auto recordGPUOnlyCalc = recordGPU.clockTimeIncrement();
    auto recordGPUCalc = recordGPU.elapsedClockTime();
    auto cpuTime2 = cpuStart.clockTimeIncrement();
    float gpuTime1;
    checkCUDAGPUError(cudaEventSynchronize(myEnd()));
    checkCUDAGPUError(cudaEventElapsedTime(&gpuTime1, myStart(), myEnd()));
#endif // TEST_PROCESS_TIME

    mixtures().destroySpeciesCUDA();
    mixtures().destroyTransportCUDA();
    dYiTP.clear();
    dHai.clear();
    dDimmMuKap.clear();
    streamTran.destroy();

#ifdef TEST_PROCESS_TIME
    myStart.destroy();
    myEnd.destroy();
#endif // TEST_PROCESS_TIME

    checkCUDAGPUError(cudaFreeHost(hostYiTPPtr));
    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp; ++i) {
            mixtures().hi(i)[n] = hostHai[n * nsp + i];
        }
    }
    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp; ++i) {
            Diff()[i][n] = hostDimmMuKap[n * (nsp + 2) + i];
        }
        mu()[n] = hostDimmMuKap[n * (nsp + 2) + nsp];
        kappal()[n] = hostDimmMuKap[n * (nsp + 2) + nsp + 1];
    }
    checkCUDAGPUError(cudaFreeHost(hostHai));
    checkCUDAGPUError(cudaFreeHost(hostDimmMuKap));

#ifdef TEST_PROCESS_TIME
    auto cpuTime3 = cpuStart.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

    muKappaDiffBoundary();

#ifdef TEST_PROCESS_TIME
    auto cpuTime4 = cpuStart.clockTimeIncrement();
    auto recordTotalCalc = recordGPU.elapsedClockTime();
    auto cpuTime12 = cpuTime1 + cpuTime2;

    HurMPI::allReduce(cpuTime1, MPI_MAX);
    HurMPI::allReduce(cpuTime2, MPI_MAX);
    HurMPI::allReduce(cpuTime3, MPI_MAX);
    HurMPI::allReduce(cpuTime4, MPI_MAX);
    HurMPI::allReduce(cpuTime12, MPI_MAX);

    HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
    HurMPI::allReduce(recordTotalCalc, MPI_MAX);
    HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
    HurMPI::allReduce(recordGPUCalc, MPI_MAX);

    real gputime11 = (real)gpuTime1 / 1000.0;
    HurMPI::allReduce(gputime11, MPI_MAX);
    if (HurMPI::master()) {
        Pout << "    Info: testing processing time in \"updateProperties\"" << std::endl;
        Pout << "          CPU_bc_gamma_E  : " << cpuTime1 << "[s]" << std::endl;
        Pout << "          CPU_sync        : " << cpuTime2 << "[s]" << std::endl;
        Pout << "          CPU_convert data: " << cpuTime3 << "[s]" << std::endl;
        Pout << "          CPU_bc_HaiMuKDim: " << cpuTime4 << "[s]" << std::endl;
        Pout << "          GPU_HaiMuKDim   : " << gputime11 << "[s]" << std::endl;
        Pout << "          Total_calc_Propeeer  : " << recordTotalCalc << "[s]" << std::endl;
        Pout << "          Total_before_calc    : " << recordGPUBeforeCalc << "[s]" << std::endl;
        Pout << "          Total_GPUCal_to_sync : " << recordGPUOnlyCalc << "[s]" << std::endl;
        Pout << "          Total_begin_to_sync  : " << recordGPUCalc << "[s]" << std::endl;
        Pout << "          Total_Fvturb_to_sync : " << cpuTime12 << "[s]" << std::endl;
    }

#endif // TEST_PROCESS_TIME
}

void OpenHurricane::laminarSpeciesSolverCUDATran::initialize() {
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
    //20210409 ��˼�� �޸�gama�����˳�򡪡�������ڳ�����>����ԭʼ�����ı߽�ֵ����>����gama�ı߽�ֵ
    //mixtures().gamma(p(), T(), gama(), true);
    //20210422 ��˼�� ����E��ʼ��
    mixtures().E(p(), T(), v(), E(), true, false);
    //bc();
    //mixtures().gammaBoundary(p(), T(), gama(), true);

    updateProperties();
}

void OpenHurricane::laminarSpeciesSolverCUDATran::calculateFc() {
    // The size of species.
    const integer nsp = mixtures().species().size();

    // Caculate the inviscous fluxes for the continuity, momentum and energy equations
    invFluxPtr_->basicFlux();

    // Calculate the convective fluxes for the species equations.
    invFluxPtr_->invFluxSpecies(yi(), false);

    // Calculate the gradient of last species.
    invFluxPtr_->grad(yi()[nsp - 1]);
}

void OpenHurricane::laminarSpeciesSolverCUDATran::calculateFv() {}

void OpenHurricane::laminarSpeciesSolverCUDATran::actualCalFv() {
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

void OpenHurricane::laminarSpeciesSolverCUDATran::calculateSource() {    
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
            for (integer n = 0; n < mesh().nCells(); ++n) {
                for (integer i = 0; i < nsp; ++i) {
                    chemtryPtr_->chemistry().Ri()[i][n] = RRi[n * nsp + i];
                }
            }
            checkCUDAGPUError(cudaFreeHost(RRi));
            chemtryPtr_->getChemistrySource();
        }
    } else if (timeMarcingPtr_->diagonalImpSource()) {
        if (!mixtures().noReaction()) {
            bool hybridPrecision = false;
            if (!hybridPrecision) {
                Pout("    Using double precision\n");
                /*HurGPU::setSharedMenBankSize();
                HurGPU::setCachePreferShared();*/
                // First step: Create CUDA stream
                cudaStreams stream;
                stream.create();
#ifdef TEST_PROCESS_TIME
                HurMPI::barrier();
                cudaEvents myStart;
                cudaEvents myEnd;
                myStart.create();
                myEnd.create();
                myStart.record(stream);
                hrClock recordGPU;
#endif // TEST_PROCESS_TIME

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

#ifdef TEST_PROCESS_TIME
                integer allUseMem =
                    static_cast<integer>((sizeOfH_YiRhoT + sizeOfRRi + sizeOfdRRi) / 1024 / 1024);
                HurMPIBase::reduce(allUseMem, MPI_SUM, HurGPU::sameGPUCommtor().masterNoInSub(),
                                   HurGPU::sameGPUCommtor());
                if (HurGPU::sameGPUCommtor().isMasterInSub()) {
                    printf("    Info: At least used GPU memory: %d MB\n", allUseMem);
                }
#endif // TEST_PROCESS_TIME

                // Fourth step: Alloc memory in the device (GPU) for array "dYiRhoT", "dRi" and "device_dRidrhoyi".
                //             The pointer "dYiRhoT" is used for storing yi (species mass fraction), rho (density) and T (temperature).
                //             The pointer "dRi" is used for storing Ri (chemical source terms).
                //             The pointer "dRi" is used for storing Ri the diagonal elements of the Jacobian matrix of the chemical source terms.
                cu2DArray<cu_real> dYiRhoT(nsp + 2, mesh().nCells());
                cu2DArray<cu_real> dRi(nsp, mesh().nCells());
                cu2DArray<cu_real> device_dRidrhoyi(nsp - 1, mesh().nCells());

#ifdef TEST_PROCESS_TIME
                auto recordGPUBeforeCalc = recordGPU.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

                // Fifth step: Compute the chemical source terms in the device (GPU)
                chemtryPtr_->chemistry().calculateSourceTermsImpAsync(
                    hostYiRhoTPtr, dYiRhoT, RRi, dRi, dRRidrhoyi, device_dRidrhoyi, stream);

#ifdef TEST_PROCESS_TIME
                myEnd.record(stream);
                hrClock cpuStart;
#endif // TEST_PROCESS_TIME

                // Sixth step: Calculate the viscous flux in the host (CPU)
                actualCalFv();

#ifdef TEST_PROCESS_TIME
                auto cpuTime1 = cpuStart.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
                // Seventh step: Synchronize the computational tasks of host (CPU) and device (GPU)
                stream.synchronize();

#ifdef TEST_PROCESS_TIME
                auto recordGPUOnlyCalc = recordGPU.clockTimeIncrement();
                auto recordGPUCalc = recordGPU.elapsedClockTime();
                auto cpuTime2 = cpuStart.clockTimeIncrement();
                float gpuTime1;
                checkCUDAGPUError(cudaEventSynchronize(myEnd()));
                checkCUDAGPUError(cudaEventElapsedTime(&gpuTime1, myStart(), myEnd()));
#endif // TEST_PROCESS_TIME

                // Eighth step: Destroy memory pointer, get chemical source terms from array RRi and the diagonal elements of the chemical Jacobian matrix from dRRidrhoyi.
                chemtryPtr_->chemistry().destroyReactionTable();
                dYiRhoT.clear();
                dRi.clear();
                device_dRidrhoyi.clear();
                stream.destroy();
#ifdef TEST_PROCESS_TIME
                myStart.destroy();
                myEnd.destroy();
#endif // TEST_PROCESS_TIME

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

#ifdef TEST_PROCESS_TIME
                auto cpuTime3 = cpuStart.clockTimeIncrement();
                auto recordTotalCalc = recordGPU.elapsedClockTime();
                auto cpuTime12 = cpuTime1 + cpuTime2;
                HurMPI::allReduce(cpuTime1, MPI_MAX);
                HurMPI::allReduce(cpuTime2, MPI_MAX);
                HurMPI::allReduce(cpuTime3, MPI_MAX);
                HurMPI::allReduce(cpuTime12, MPI_MAX);

                HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
                HurMPI::allReduce(recordTotalCalc, MPI_MAX);
                HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
                HurMPI::allReduce(recordGPUCalc, MPI_MAX);

                real gputime11 = (real)gpuTime1 / 1000.0;
                HurMPI::allReduce(gputime11, MPI_MAX);
                if (HurMPI::master()) {
                    Pout << "    Info: testing processing time in "
                            "\"calculateSource\""
                         << std::endl;
                    Pout << "          CPU_FV               : " << cpuTime1 << "[s]" << std::endl;
                    Pout << "          CPU_sync             : " << cpuTime2 << "[s]" << std::endl;
                    Pout << "          CPU_convert data     : " << cpuTime3 << "[s]" << std::endl;
                    Pout << "          GPU_ChemicalSource   : " << gputime11 << "[s]" << std::endl;
                    Pout << "          Total_calc_source    : " << recordTotalCalc << "[s]"
                         << std::endl;
                    Pout << "          Total_before_calc    : " << recordGPUBeforeCalc << "[s]"
                         << std::endl;
                    Pout << "          Total_GPUCal_to_sync : " << recordGPUOnlyCalc << "[s]"
                         << std::endl;
                    Pout << "          Total_begin_to_sync  : " << recordGPUCalc << "[s]"
                         << std::endl;
                    Pout << "          Total_Fv_to_sync     : " << cpuTime12 << "[s]" << std::endl;
                }
#endif // TEST_PROCESS_TIME
            } else {
                Pout("    Using hybrid precision\n");
                /*HurGPU::setSharedMenBankSize();
                HurGPU::setCachePreferShared();*/
                // First step: Create CUDA stream
                cudaStreams stream;
                stream.create();
#ifdef TEST_PROCESS_TIME
                HurMPI::barrier();
                cudaEvents myStart;
                cudaEvents myEnd;
                myStart.create();
                myEnd.create();
                myStart.record(stream);
                hrClock recordGPU;
#endif // TEST_PROCESS_TIME

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

                cu_float *dRRidrhoyi = nullptr;
                size_t sizeOfdRRi = (static_cast<unsigned long long>(nsp) - 1) * mesh().nCells() *
                                    sizeof(cu_float);
                checkCUDAGPUError(
                    cudaHostAlloc((void **)&dRRidrhoyi, sizeOfdRRi, cudaHostAllocDefault));

#ifdef TEST_PROCESS_TIME
                integer allUseMem =
                    static_cast<integer>((sizeOfH_YiRhoT + sizeOfRRi + sizeOfdRRi) / 1024 / 1024);
                HurMPIBase::reduce(allUseMem, MPI_SUM, HurGPU::sameGPUCommtor().masterNoInSub(),
                                   HurGPU::sameGPUCommtor());
                if (HurGPU::sameGPUCommtor().isMasterInSub()) {
                    printf("    Info: At least used GPU memory: %d MB\n", allUseMem);
                }
#endif // TEST_PROCESS_TIME
                // Fourth step: Alloc memory in the device (GPU) for array "dYiRhoT", "dRi" and "device_dRidrhoyi".
                //             The pointer "dYiRhoT" is used for storing yi (species mass fraction), rho (density) and T (temperature).
                //             The pointer "dRi" is used for storing Ri (chemical source terms).
                //             The pointer "dRi" is used for storing Ri the diagonal elements of the Jacobian matrix of the chemical source terms.
                cu2DArray<cu_real> dYiRhoT(nsp + 2, mesh().nCells());
                cu2DArray<cu_real> dRi(nsp, mesh().nCells());
                cu2DArray<cu_float> device_dRidrhoyi(nsp - 1, mesh().nCells());

#ifdef TEST_PROCESS_TIME
                auto recordGPUBeforeCalc = recordGPU.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

                // Fifth step: Compute the chemical source terms in the device (GPU)
                chemtryPtr_->chemistry().calculateSourceTermsImpAsyncHybrid(
                    hostYiRhoTPtr, dYiRhoT, RRi, dRi, dRRidrhoyi, device_dRidrhoyi, stream);

#ifdef TEST_PROCESS_TIME
                myEnd.record(stream);
                hrClock cpuStart;
#endif // TEST_PROCESS_TIME

                // Sixth step: Calculate the viscous flux in the host (CPU)
                actualCalFv();

#ifdef TEST_PROCESS_TIME
                auto cpuTime1 = cpuStart.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
                // Seventh step: Synchronize the computational tasks of host (CPU) and device (GPU)
                stream.synchronize();

#ifdef TEST_PROCESS_TIME
                auto recordGPUOnlyCalc = recordGPU.clockTimeIncrement();
                auto recordGPUCalc = recordGPU.elapsedClockTime();
                auto cpuTime2 = cpuStart.clockTimeIncrement();
                float gpuTime1;
                checkCUDAGPUError(cudaEventSynchronize(myEnd()));
                checkCUDAGPUError(cudaEventElapsedTime(&gpuTime1, myStart(), myEnd()));
#endif // TEST_PROCESS_TIME

                // Eighth step: Destroy memory pointer, get chemical source terms from array RRi and the diagonal elements of the chemical Jacobian matrix from dRRidrhoyi.
                chemtryPtr_->chemistry().destroyReactionTable();
                dYiRhoT.clear();
                dRi.clear();
                device_dRidrhoyi.clear();
                stream.destroy();
#ifdef TEST_PROCESS_TIME
                myStart.destroy();
                myEnd.destroy();
#endif // TEST_PROCESS_TIME

                checkCUDAGPUError(cudaFreeHost(hostYiRhoTPtr));

                for (integer i = 0; i < nsp; ++i) {
                    for (integer n = 0; n < mesh().nCells(); ++n) {
                        chemtryPtr_->chemistry().Ri()[i][n] = RRi[n * nsp + i];
                        if (i < nsp - 1) {
                            yi()[i].diagSource()[n] =
                                static_cast<real>(dRRidrhoyi[n * (nsp - 1) + i]);
                        }
                    }
                }
                checkCUDAGPUError(cudaFreeHost(RRi));
                checkCUDAGPUError(cudaFreeHost(dRRidrhoyi));
                chemtryPtr_->getImpChemistrySource(dt_, !iter().isGlobalTimeStep());

#ifdef TEST_PROCESS_TIME
                auto cpuTime3 = cpuStart.clockTimeIncrement();
                auto recordTotalCalc = recordGPU.elapsedClockTime();
                auto cpuTime12 = cpuTime1 + cpuTime2;
                HurMPI::allReduce(cpuTime1, MPI_MAX);
                HurMPI::allReduce(cpuTime2, MPI_MAX);
                HurMPI::allReduce(cpuTime3, MPI_MAX);
                HurMPI::allReduce(cpuTime12, MPI_MAX);

                HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
                HurMPI::allReduce(recordTotalCalc, MPI_MAX);
                HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
                HurMPI::allReduce(recordGPUCalc, MPI_MAX);

                real gputime11 = (real)gpuTime1 / 1000.0;
                HurMPI::allReduce(gputime11, MPI_MAX);
                if (HurMPI::master()) {
                    Pout << "    Info: testing processing time in "
                            "\"calculateSource\""
                         << std::endl;
                    Pout << "          CPU_FV               : " << cpuTime1 << "[s]" << std::endl;
                    Pout << "          CPU_sync             : " << cpuTime2 << "[s]" << std::endl;
                    Pout << "          CPU_convert data     : " << cpuTime3 << "[s]" << std::endl;
                    Pout << "          GPU_ChemicalSource   : " << gputime11 << "[s]" << std::endl;
                    Pout << "          Total_calc_source    : " << recordTotalCalc << "[s]"
                         << std::endl;
                    Pout << "          Total_before_calc    : " << recordGPUBeforeCalc << "[s]"
                         << std::endl;
                    Pout << "          Total_GPUCal_to_sync : " << recordGPUOnlyCalc << "[s]"
                         << std::endl;
                    Pout << "          Total_begin_to_sync  : " << recordGPUCalc << "[s]"
                         << std::endl;
                    Pout << "          Total_Fv_to_sync     : " << cpuTime12 << "[s]" << std::endl;
                }
#endif // TEST_PROCESS_TIME
            }
        }
    } else if (timeMarcingPtr_->fullJacobianSource()) {
        LFatal("This solver does not support full Jacobian point-implicit operation");
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

void OpenHurricane::laminarSpeciesSolverCUDATran::updatePrimitives(const bool shouldUpdateTemp) {
    mixtures().lastSpeAndNormalized();

    solver::updatePrimitives(shouldUpdateTemp);

    limits();
    updateProperties();
}

void OpenHurricane::laminarSpeciesSolverCUDATran::updateFlowOld() {}

void OpenHurricane::laminarSpeciesSolverCUDATran::write() {
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

void OpenHurricane::laminarSpeciesSolverCUDATran::getYiTP(real *hur_restrict hostYiTPPtr_) const {
    const auto nsp = specTable().size();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp; ++i) {
            hostYiTPPtr_[n * (nsp + 2) + i] = yi()[i][n];
        }
        hostYiTPPtr_[n * (nsp + 2) + nsp] = T()[n];
        hostYiTPPtr_[n * (nsp + 2) + nsp + 1] = p()[n];
    }
}

void OpenHurricane::laminarSpeciesSolverCUDATran::makeFzBoundStart() const {
    const auto &fZL = mesh().faceZones();
    integer countFZB = Zero;
    for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
        if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
            !fZL[fzi].isPeriodicShadow()) {
            ++countFZB;
        }
    }
    fzBoundStartPtr_.reset(new integerList(countFZB + 1));

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

void OpenHurricane::laminarSpeciesSolverCUDATran::muKappaDiffBoundary() {
    /*HurGPU::setSharedMenBankSize();
    HurGPU::setCachePreferShared();*/
    const integer nsp = flowPtr_->mixtures().species().size();
    nThreadsAndBlocks nTBForTranBnd(nsp, fzBoundStart().last());
    nThreadsAndBlocks nTBForHaiCpBnd(nsp, fzBoundStart().last());
    nTBForTranBnd.setSharedMemSize1<real>(2 * nsp * nTBForTranBnd.blockSize().y);

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

    real *hostHai = nullptr;
    size_t sizeOfH_Hai =
        (static_cast<unsigned long long>(nsp)) * fzBoundStart().last() * sizeof(real);

    checkCUDAGPUError(cudaHostAlloc((void **)&hostHai, sizeOfH_Hai, cudaHostAllocDefault));
    cu2DArray<cu_real> dHai(nsp, fzBoundStart().last());

    CUDAThermo::calchaiAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA().species(),
                             nTBForHaiCpBnd, nsp, fzBoundStart().last(), hostHai, dHai, streamTran);
    CUDATransport::calcTransportPropertiesAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA(),
                                                nTBForTranBnd, nsp, fzBoundStart().last(),
                                                hostDimmMuKap, dDimmMuKap, streamTran, false);

    transferTranP();

    streamTran.synchronize();
    mixtures().destroySpeciesCUDA();
    mixtures().destroyTransportCUDA();
    dYiTP.clear();
    dHai.clear();
    dDimmMuKap.clear();
    streamTran.destroy();
    checkCUDAGPUError(cudaFreeHost(hostYiTPPtr));

    setExtrapolate(hostDimmMuKap, hostHai);

    checkCUDAGPUError(cudaFreeHost(hostHai));
    checkCUDAGPUError(cudaFreeHost(hostDimmMuKap));
}

void OpenHurricane::laminarSpeciesSolverCUDATran::getYiTPBoundary(
    real *hur_restrict hostYiTPPtr_) const {
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

void OpenHurricane::laminarSpeciesSolverCUDATran::setExtrapolate(const real *hur_restrict hostDimmMuKap,
                                                             const real *hur_restrict hostHai) {
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
                        real(2) * hostHai[cntz * nsp + isp] - mixtures().hi(isp)[cl];
                    Diff()[isp][cr] =
                        real(2.0) * hostDimmMuKap[cntz * (nsp + 2) + isp] - Diff()[isp][cl];
                }
                mu()[cr] = real(2.0) * hostDimmMuKap[cntz * (nsp + 2) + nsp] - mu()[cl];
                kappal()[cr] = real(2.0) * hostDimmMuKap[cntz * (nsp + 2) + nsp + 1] - kappal()[cl];
                ++cntz;
            }

            ++countFZB;
        }
    }
}

void OpenHurricane::laminarSpeciesSolverCUDATran::transferTranP() {
    const auto nsp = specTable().size();
    realTransfer myTransfer1(mesh(), mu(), false, true);
    realTransfer myTransfer2(mesh(), kappal(), false, true);

    myTransfer1.transferInit();
    myTransfer2.transferInit();
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
}

#endif // CUDA_PARALLEL