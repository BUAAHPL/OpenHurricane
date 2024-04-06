/*!
 * \file turbulentSpeciesSolverCUDATran.cpp
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

#ifdef CUDA_PARALLEL
#include "turbulentSpeciesSolverCUDATran.hpp"
#include "calculateFieldVar.hpp"
#include "laminar.hpp"
#include "laminarFlow.hpp"
#include "solutionWrite.hpp"
#include "viscousFlux.hpp"

#ifdef TEST_PROCESS_TIME
#include "cudaEvents.hpp"
#endif // TEST_PROCESS_TIME

namespace OpenHurricane {
    createClassName(turbulentSpeciesSolverCUDATran);
    registerObjFty(solver, turbulentSpeciesSolverCUDATran, controller);
} // namespace OpenHurricane

OpenHurricane::turbulentSpeciesSolverCUDATran::turbulentSpeciesSolverCUDATran(iteration &iter,
                                                                          const runtimeMesh &mesh)
    : solver(iter, mesh), chemtryPtr_(nullptr), rhoId_(0), rhouId_(1), rhoEId_(4), rhoYi0Id_(6),
      rhoTurb0Id_(5), nTBForTran_(), nTBForHaiCp_(), fzBoundStartPtr_(nullptr)
#ifdef TEST_PROCESS_TIME
      ,
      gasProFOS_(IOsstream::streamFormat::ASCII_FORMAT, IOsstream::openOptions::ONLY_MASTER,
                 std::ios_base::out),
      gasProProcessTime_(), chemSorFOS_(IOsstream::streamFormat::ASCII_FORMAT,
                                        IOsstream::openOptions::ONLY_MASTER, std::ios_base::out),
      chemSorProcessTime_()
#endif // TEST_PROCESS_TIME
{

    turbPtr_.reset(new laminar(iter.cont(), *flowPtr_));
    if (iter.cont().found("flow")) {
        const auto &flowCont = iter.cont().subController("flow");
        if (flowCont.found("turbulence")) {
            auto &turbCont = flowCont.subController("turbulence");
            turbPtr_ = turbulenceModel::creator(turbCont, *flowPtr_);
        } else {
            LFatal("Cannot find turbulence setting section in flow section");
        }
        if (!mixtures().noReaction()) {
            if (flowCont.found("mixture")) {
                const auto &reacCont = flowCont.subController("mixture").subController("reactions");
                if (reacCont.subController("combustionModel").findWord("type") == "finiteRate" ||
                    reacCont.subController("combustionModel").findWord("type") == "PaSR") {
                    chemtryPtr_ = combustionModel::creator(*flowPtr_, reacCont, *turbPtr_);
                } else {
                    const auto w = reacCont.subController("combustionModel").findWord("type");
                    errorAbortStr(("Only can use finiteRate or PaSR in reacting "
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
    nTBForHaiCp_.setSharedMemSize1<real>(specTable().size() * nTBForHaiCp_.blockSize().y);

#ifdef TEST_PROCESS_TIME
    if (!mixtures().noReaction()) {
        fileName outFile = iter.configName().name(true) + "" + "TestGPUChemTime" + ".dat";
        outFile = iter.outputPath() / outFile;
        chemSorFOS_.changeFileName(outFile);

        chemSorFOS_.open();
    } else {
        chemSorFOS_.close();
    }
    fileName outFile = iter.configName().name(true) + "" + "TestGPUGasProTime" + ".dat";
    outFile = iter.outputPath() / outFile;
    gasProFOS_.changeFileName(outFile);

    gasProFOS_.open();
#endif // TEST_PROCESS_TIME

    if (HurGPU::GPUDev().size() != 1) {
        LFatal("Only support single GPU mode in this class");
    }
    HurGPU::setDevice(HurGPU::GPUDev()[0]);
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::solving() {
    const integer nsp = mixtures().species().size();

    // First step: Register the primitives parameter in the time marching method.
    marching().addObject(rho());
    marching().addObject(v());
    marching().addObject(E());
    if (turbPtr_->isCoupled()) {
        for (integer i = 0; i < turbPtr_->nEq(); ++i) {
            marching().addObject(turbPtr_->var(i));
        }
    }
    for (integer i = 0; i < nsp - 1; ++i) {
        marching().addObject(yi()[i]);
    }
    rhoId_ = 0;
    rhouId_ = 1;
    rhoEId_ = 4;
    if (turbPtr_->isCoupled()) {
        rhoTurb0Id_ = 5;
        rhoYi0Id_ = 5 + turbPtr_->nEq();
    } else {
        rhoYi0Id_ = 5;
    }
    solver::solving();
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::BDFSolve() {
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
    rhoTurb0Id_ = 5;
    rhoYi0Id_ = 5 + turbPtr_->nEq();
    marching().initializing();

    // Second step: Initialize the whole flow field
    initialize();

    solver::BDFSolve();
    //marching().marching();
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::clear() noexcept {
    invFluxPtr_.clear();
    chemtryPtr_.clear();
    fzBoundStartPtr_.clear();
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::bc() {
    const integer nsp = mixtures().species().size();
    rho().updateBoundary();
    v().updateBoundary();
    p().updateBoundary();

    for (integer isp = 0; isp < nsp; ++isp) {
        yi()[isp].updateBoundary();
    }

    T().updateBoundary();
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::timeStep(realArray &dt) {}

void OpenHurricane::turbulentSpeciesSolverCUDATran::updateProperties() {
    const integer nsp = flowPtr_->mixtures().species().size();
    size_t sizeOfH_YiTP =
        (static_cast<unsigned long long>(nsp) + 2) * mesh().nCells() * sizeof(real);
    size_t sizeOfH_DimmMuKap =
        (static_cast<unsigned long long>(nsp) + 2) * mesh().nCells() * sizeof(real);
    size_t sizeOfH_HaiCP =
        (static_cast<unsigned long long>(nsp + 1)) * mesh().nCells() * sizeof(real);

    integer ava, tot;
    HurGPU::getMemInfo(ava, tot);
    integer allUseMem =
        static_cast<integer>((sizeOfH_YiTP + sizeOfH_DimmMuKap + sizeOfH_HaiCP) / 1024 / 1024);
    integer allUseMems;
    HurMPIBase::allReduce(&allUseMem, &allUseMems, 1, feature<integer>::MPIType, MPI_SUM,
                          HurGPU::sameGPUComm());

#ifdef TEST_PROCESS_TIME
    Pout("    Info: available GPU memory: %d MB, total GPU memory: %d MB\n", ava, tot);
    if (HurGPU::sameGPUCommtor().isMasterInSub()) {
        printf("    Info: At least used GPU memory: %d MB at GPU: %d for gas properties\n",
               allUseMems, HurGPU::GPUDev().first());
    }
#endif // TEST_PROCESS_TIME
    integer nSlice = allUseMems / ava + 1;
    if (nSlice > 1) {
#ifdef TEST_PROCESS_TIME
        if (HurGPU::sameGPUCommtor().isMasterInSub()) {
            printf("    Info: Using %d slices at GPU: %d for gas properties\n", nSlice,
                   HurGPU::GPUDev().first());
        }
#endif
        updatePropertiesBySlice(nSlice);
        return;
    }

    //HurGPU::setSharedMenBankSize();
    //HurGPU::setCachePreferShared();
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
#endif // TEST_PROCESS_TIME

    // Second step: Create transport table in memory of the device (GPU)
    mixtures().transportCUDA(streamTran);

    // Third step: (1) Alloc memory in the host (CPU) for array "hostYiTPPtr", "hostDimmMuKap" and "hostHai".
    //             The pointer "hostYiTPPtr" is used for storing yi (species mass fraction), T (temperature) and p (pressure).
    //             The pointer "hostDimmMuKap" is used for storing Dim, mu and kappa.
    //             The pointer "hostHai" is used for storing hai.
    //             (2) Alloc memory in the device (GPU) for array "dYiTP", "dDimmMuKap" and "dHai".
    real *hostYiTPPtr = nullptr;
    checkCUDAGPUError(cudaHostAlloc((void **)&hostYiTPPtr, sizeOfH_YiTP, cudaHostAllocDefault));
    getYiTP(hostYiTPPtr);

    real *hostDimmMuKap = nullptr;
    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostDimmMuKap, sizeOfH_DimmMuKap, cudaHostAllocDefault));

    real *hostHaiCP = nullptr;
    checkCUDAGPUError(cudaHostAlloc((void **)&hostHaiCP, sizeOfH_HaiCP, cudaHostAllocDefault));

    cu2DArray<cu_real> dYiTP(nsp + 2, mesh().nCells());
    cu2DArray<cu_real> dDimmMuKap(nsp + 2, mesh().nCells());
    cu2DArray<cu_real> dHaiCp(nsp + 1, mesh().nCells());

#ifdef TEST_PROCESS_TIME
    auto recordGPUBeforeCalc = recordGPU.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

    // Fourth step: Compute the hai in the device (GPU)
    CUDAThermo::calchaicpAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA().species(),
                               nTBForHaiCp_, nsp, mesh().nCells(), hostHaiCP, dHaiCp, streamTran);

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
    dHaiCp.clear();
    dDimmMuKap.clear();
    streamTran.destroy();

#ifdef TEST_PROCESS_TIME
    myStart.destroy();
    myEnd.destroy();
#endif // TEST_PROCESS_TIME

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

#ifdef TEST_PROCESS_TIME
    auto cpuTime3 = cpuStart.clockTimeIncrement();
#endif // TEST_PROCESS_TIME
    muKappaDiffBoundary();

#ifdef TEST_PROCESS_TIME
    auto cpuTime4 = cpuStart.clockTimeIncrement();
    auto recordTotalCalc = recordGPU.elapsedClockTime();
    auto cpuTime12 = cpuTime1 + cpuTime2;
    gasProProcessTime_.CPU_bc_gamma_E_ = cpuTime1;
    gasProProcessTime_.CPU_sync_ = cpuTime2;
    gasProProcessTime_.CPU_convert_data_ = cpuTime3;
    gasProProcessTime_.CPU_bc_HaiMuKDim_ = cpuTime4;
    gasProProcessTime_.Total_Fvturb_to_sync_ = cpuTime12;
    gasProProcessTime_.Total_before_calc_ = recordGPUBeforeCalc;
    gasProProcessTime_.Total_calc_Propeeer_ = recordTotalCalc;
    gasProProcessTime_.Total_GPUCal_to_sync_ = recordGPUOnlyCalc;
    gasProProcessTime_.Total_begin_to_sync_ = recordGPUCalc;
    //HurMPI::allReduce(cpuTime1, MPI_MAX);
    //HurMPI::allReduce(cpuTime2, MPI_MAX);
    //HurMPI::allReduce(cpuTime3, MPI_MAX);
    //HurMPI::allReduce(cpuTime4, MPI_MAX);
    //HurMPI::allReduce(cpuTime12, MPI_MAX);

    //HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
    //HurMPI::allReduce(recordTotalCalc, MPI_MAX);
    //HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
    //HurMPI::allReduce(recordGPUCalc, MPI_MAX);
    gasProProcessTime_.GPU_HaiMuKDim_ = (real)gpuTime1 / 1000.0;
    /*real gputime11 = (real)gpuTime1 / 1000.0;
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
    }*/

#endif // TEST_PROCESS_TIME
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::initialize() {
    solver::initialize();
    const integer nsp = mixtures().species().size();
    PtrList<cellRealArray> tmpCellArray;
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

    mixtures().E(p(), T(), v(), E(), true, false);

    updateProperties();

    if (!iter().restart()) {
        turbPtr_->initialize();
    } else {        
        turbPtr_->initializeRestart();
    }

#ifdef TEST_PROCESS_TIME
    if (!mixtures().noReaction()) {
        chemSorProcessTime_.writeTecplotHead(chemSorFOS_);
    }
    gasProProcessTime_.writeTecplotHead(gasProFOS_);
#endif // TEST_PROCESS_TIME
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::calculateFc() {
    // The size of species.
    const integer nsp = mixtures().species().size();

    // Caculate the inviscous fluxes for the continuity, momentum and energy equations
    invFluxPtr_->basicFlux();
    // Calculate the convective fluxes for the turbulent equations for turbulent coupling solver.
    if (turbPtr_->isCoupled()) {
        for (integer i = 0; i < turbPtr_->nEq(); ++i) {
            invFluxPtr_->invFlux(turbPtr_->var(i),
                                 realBounded::lowerBound(real(0), realBounded::BOUNDED_FLAG));
        }
    }
    // Calculate the convective fluxes for the species equations.
    invFluxPtr_->invFluxSpecies(yi(), false);

    // Calculate the gradient of turbulent variables for splitting solver.
    turbPtr_->calcGrad(*invFluxPtr_);
    // Calculate the gradient of last species.
    invFluxPtr_->grad(yi()[nsp - 1]);
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::calculateFv() {}

void OpenHurricane::turbulentSpeciesSolverCUDATran::actualCalFv() {
    const auto &spT = flowPtr_->mixtures().species();
    const faceTensorArray deltafV(fv::gradf(v()));

    faceVectorArray Vf(fv::interpolate(v()));
    const faceRealArray rhof(fv::interpolate(rho()));
    const faceRealArray muf(fv::interpolate(mu()));
    faceRealArray mutf(fv::interpolate(mut()));
    const faceRealArray cpf(fv::interpolate(cp()));

    const faceSymmTensorArray tau(turbPtr_->tauEff(rhof, muf, mutf, deltafV));

    //v().rhs() += visFlux(tau);
    visFlux(tau, v());
    turbPtr_->visFlux(rhof, muf, mutf, mu(), mut(), rho());

    const faceRealArray kappaf(fv::interpolate(kappal()));
    const realArray kappatf(mutf * cpf / prt());
    Vf = tau * Vf;
    Vf += ((kappaf + kappatf) * fv::gradf(T()));

    mutf /= sct();

    faceRealArray Dimf(object("interpolation(Dimf)", mesh(), object::NOT_WRITE, object::TEMPORARY),
                       mesh());

    faceVectorArray Ji(object("Ji", mesh(), object::NOT_WRITE, object::TEMPORARY), mesh());

    for (integer i = 0; i < spT.size(); ++i) {
        //const faceRealArray Dimf(fv::interpolate(Diff()[i]));
        fv::interpolate(Diff()[i], Dimf);
        Dimf *= rhof;
        Dimf += mutf;
        fv::gradf(yi()[i], Ji);
        Ji *= Dimf;
        //Ji = (rhof * Dimf + mutf) * fv::gradf(yi()[i]);
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

    turbPtr_->correctEnergyEqVisFlux(E());
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::calculateSource() {   
    const integer nsp = flowPtr_->mixtures().species().size();
    if (timeMarcingPtr_->explicitSource()) {
        if (!mixtures().noReaction()) {
            size_t sizeOfH_YiRhoT =
                (static_cast<unsigned long long>(nsp) + 2) * mesh().nCells() * sizeof(real);
            size_t sizeOfRRi =
                static_cast<unsigned long long>(nsp) * mesh().nCells() * sizeof(real);

            integer ava, tot;
            HurGPU::getMemInfo(ava, tot);
            integer allUseMem = static_cast<integer>((sizeOfH_YiRhoT + sizeOfRRi) / 1024 / 1024);
            integer allUseMems;
            HurMPIBase::allReduce(&allUseMem, &allUseMems, 1, feature<integer>::MPIType, MPI_SUM,
                                  HurGPU::sameGPUComm());

#ifdef TEST_PROCESS_TIME
            Pout("    Info: available GPU memory: %d MB, total GPU memory: %d MB\n", ava, tot);
            if (HurGPU::sameGPUCommtor().isMasterInSub()) {
                printf("    Info: At least used GPU memory: %d MB at GPU: %d for chemical source\n",
                       allUseMems, HurGPU::GPUDev().first());
            }
#endif // TEST_PROCESS_TIME
            integer nSlice = allUseMems / ava + 1;
            if (nSlice > 1) {
#ifdef TEST_PROCESS_TIME
                if (HurGPU::sameGPUCommtor().isMasterInSub()) {
                    printf("    Info: Using %d slices at GPU: %d for chemical source\n", nSlice,
                           HurGPU::GPUDev().first());
                }
#endif
                updateChemSourceBySlice(nSlice);
            } else {
                /* HurGPU::setSharedMenBankSize();
                 HurGPU::setCachePreferShared();*/
                // First step: Create CUDA stream
                cudaStreams stream;
                stream.create();

                // Second step: Create reaction table in memory of the device (GPU).
                chemtryPtr_->chemistry().createReactionTableAsync(stream);

                // Third step: Alloc memory in the host (CPU) for array "hostYiRhoTPtr" and "RRi".
                //             The pointer "hostYiRhoTPtr" is used for storing yi (species mass fraction), rho (density) and T (temperature).
                //             The pointer "RRi" is used for storing Ri (chemical source terms).
                real *hostYiRhoTPtr = nullptr;
                checkCUDAGPUError(
                    cudaHostAlloc((void **)&hostYiRhoTPtr, sizeOfH_YiRhoT, cudaHostAllocDefault));
                real *RRi = nullptr;
                checkCUDAGPUError(cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

                // Fourth step: Alloc memory in the device (GPU) for array "dYiRhoT" and "dRi".
                //             The pointer "dYiRhoT" is used for storing yi (species mass fraction), rho (density) and T (temperature).
                //             The pointer "dRi" is used for storing Ri (chemical source terms).
                cu2DArray<cu_real> dYiRhoT(nsp + 2, mesh().nCells());
                cu2DArray<cu_real> dRi(nsp, mesh().nCells());

                // Fifth step: Compute the chemical source terms in the device (GPU)
                chemtryPtr_->chemistry().calculateSourceTermsAsync(hostYiRhoTPtr, dYiRhoT, RRi, dRi,
                                                                   stream);

                // Sixth step: Calculate the turbulence source and viscous flux in the host (CPU)
                turbPtr_->expSource();
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
        } else {
            turbPtr_->expSource();
        }
    } else if (timeMarcingPtr_->diagonalImpSource()) {
        if (!mixtures().noReaction()) {
            bool hybridPrecision = false;
            if (!hybridPrecision) {
                size_t sizeOfH_YiRhoT =
                    (static_cast<unsigned long long>(nsp) + 2) * mesh().nCells() * sizeof(real);
                size_t sizeOfRRi =
                    static_cast<unsigned long long>(nsp) * mesh().nCells() * sizeof(real);
                size_t sizeOfdRRi =
                    (static_cast<unsigned long long>(nsp) - 1) * mesh().nCells() * sizeof(real);
                integer ava, tot;
                HurGPU::getMemInfo(ava, tot);
                integer allUseMem =
                    static_cast<integer>((sizeOfH_YiRhoT + sizeOfRRi + sizeOfdRRi) / 1024 / 1024);
                integer allUseMems;
                HurMPIBase::allReduce(&allUseMem, &allUseMems, 1, feature<integer>::MPIType,
                                      MPI_SUM, HurGPU::sameGPUComm());

#ifdef TEST_PROCESS_TIME
                Pout("    Info: available GPU memory: %d MB, total GPU memory: %d MB\n", ava, tot);
                if (HurGPU::sameGPUCommtor().isMasterInSub()) {
                    printf("    Info: At least used GPU memory: %d MB at GPU: %d for chemical "
                           "source (Imp-real)\n",
                           allUseMems, HurGPU::GPUDev().first());
                }
#endif // TEST_PROCESS_TIME
                integer nSlice = allUseMems / ava + 1;
                if (nSlice > 1) {
#ifdef TEST_PROCESS_TIME
                    if (HurGPU::sameGPUCommtor().isMasterInSub()) {
                        printf("    Info: Using %d slices at GPU: %d for chemical source\n", nSlice,
                               HurGPU::GPUDev().first());
                    }
#endif
                    updateChemSourceImpBySlice(nSlice, hybridPrecision);
                } else {
                    //HurGPU::setSharedMenBankSize();
                    //HurGPU::setCachePreferShared();
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
                    checkCUDAGPUError(cudaHostAlloc((void **)&hostYiRhoTPtr, sizeOfH_YiRhoT,
                                                    cudaHostAllocDefault));
                    real *RRi = nullptr;
                    checkCUDAGPUError(
                        cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

                    real *dRRidrhoyi = nullptr;
                    checkCUDAGPUError(
                        cudaHostAlloc((void **)&dRRidrhoyi, sizeOfdRRi, cudaHostAllocDefault));

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

                    // Sixth step: Calculate the turbulence source and viscous flux in the host (CPU)
                    turbPtr_->impSource();
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
                    chemSorProcessTime_.CPU_FV_turbSource_ = cpuTime1;
                    chemSorProcessTime_.CPU_sync_ = cpuTime2;
                    chemSorProcessTime_.CPU_convert_data_ = cpuTime3;
                    chemSorProcessTime_.Total_Fvturb_to_sync_ = cpuTime12;
                    chemSorProcessTime_.Total_before_calc_ = recordGPUBeforeCalc;
                    chemSorProcessTime_.Total_calc_source_ = recordTotalCalc;
                    chemSorProcessTime_.Total_GPUCal_to_sync_ = recordGPUOnlyCalc;
                    chemSorProcessTime_.Total_begin_to_sync_ = recordGPUCalc;
                    //HurMPI::allReduce(cpuTime1, MPI_MAX);
                    //HurMPI::allReduce(cpuTime2, MPI_MAX);
                    //HurMPI::allReduce(cpuTime3, MPI_MAX);
                    //HurMPI::allReduce(cpuTime12, MPI_MAX);

                    //HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
                    //HurMPI::allReduce(recordTotalCalc, MPI_MAX);
                    //HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
                    //HurMPI::allReduce(recordGPUCalc, MPI_MAX);
                    chemSorProcessTime_.GPU_ChemicalSource_ = (real)gpuTime1 / 1000.0;
                    /*real gputime11 = (real)gpuTime1 / 1000.0;
                    HurMPI::allReduce(gputime11, MPI_MAX);*/
                    /*if (HurMPI::master()) {
                        Pout << "    Info: testing processing time in "
                                "\"calculateSource\""
                             << std::endl;
                        Pout << "          CPU_FV_turbSource    : " << cpuTime1 << "[s]"
                             << std::endl;
                        Pout << "          CPU_sync             : " << cpuTime2 << "[s]"
                             << std::endl;
                        Pout << "          CPU_convert data     : " << cpuTime3 << "[s]"
                             << std::endl;
                        Pout << "          GPU_ChemicalSource   : " << gputime11 << "[s]"
                             << std::endl;
                        Pout << "          Total_calc_source    : " << recordTotalCalc << "[s]"
                             << std::endl;
                        Pout << "          Total_before_calc    : " << recordGPUBeforeCalc << "[s]"
                             << std::endl;
                        Pout << "          Total_GPUCal_to_sync : " << recordGPUOnlyCalc << "[s]"
                             << std::endl;
                        Pout << "          Total_begin_to_sync  : " << recordGPUCalc << "[s]"
                             << std::endl;
                        Pout << "          Total_Fvturb_to_sync : " << cpuTime12 << "[s]"
                             << std::endl;
                    }*/
#endif // TEST_PROCESS_TIME
                }
            } else {
                size_t sizeOfH_YiRhoT =
                    (static_cast<unsigned long long>(nsp) + 2) * mesh().nCells() * sizeof(real);
                size_t sizeOfRRi =
                    static_cast<unsigned long long>(nsp) * mesh().nCells() * sizeof(real);
                size_t sizeOfdRRi = (static_cast<unsigned long long>(nsp) - 1) * mesh().nCells() *
                                    sizeof(cu_float);
                integer ava, tot;
                HurGPU::getMemInfo(ava, tot);
                integer allUseMem =
                    static_cast<integer>((sizeOfH_YiRhoT + sizeOfRRi + sizeOfdRRi) / 1024 / 1024);
                integer allUseMems;
                HurMPIBase::allReduce(&allUseMem, &allUseMems, 1, feature<integer>::MPIType,
                                      MPI_SUM, HurGPU::sameGPUComm());

#ifdef TEST_PROCESS_TIME
                Pout("    Info: available GPU memory: %d MB, total GPU memory: %d MB\n", ava, tot);
                if (HurGPU::sameGPUCommtor().isMasterInSub()) {
                    printf("    Info: At least used GPU memory: %d MB at GPU: %d for chemical "
                           "source (imp-hybrid) \n",
                           allUseMems, HurGPU::GPUDev().first());
                }
#endif // TEST_PROCESS_TIME
                integer nSlice = allUseMems / ava + 1;
                if (nSlice > 1) {
#ifdef TEST_PROCESS_TIME
                    if (HurGPU::sameGPUCommtor().isMasterInSub()) {
                        printf("    Info: Using %d slices at GPU: %d for chemical source\n", nSlice,
                               HurGPU::GPUDev().first());
                    }
#endif
                    updateChemSourceImpBySlice(nSlice, hybridPrecision);
                } else {
                    //HurGPU::setSharedMenBankSize();
                    //HurGPU::setCachePreferShared();
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
                    checkCUDAGPUError(cudaHostAlloc((void **)&hostYiRhoTPtr, sizeOfH_YiRhoT,
                                                    cudaHostAllocDefault));

                    real *RRi = nullptr;
                    checkCUDAGPUError(
                        cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

                    cu_float *dRRidrhoyi = nullptr;
                    checkCUDAGPUError(
                        cudaHostAlloc((void **)&dRRidrhoyi, sizeOfdRRi, cudaHostAllocDefault));

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

                    // Sixth step: Calculate the turbulence source and viscous flux in the host (CPU)
                    turbPtr_->impSource();
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
                    chemSorProcessTime_.CPU_FV_turbSource_ = cpuTime1;
                    chemSorProcessTime_.CPU_sync_ = cpuTime2;
                    chemSorProcessTime_.CPU_convert_data_ = cpuTime3;
                    chemSorProcessTime_.Total_Fvturb_to_sync_ = cpuTime12;
                    chemSorProcessTime_.Total_before_calc_ = recordGPUBeforeCalc;
                    chemSorProcessTime_.Total_calc_source_ = recordTotalCalc;
                    chemSorProcessTime_.Total_GPUCal_to_sync_ = recordGPUOnlyCalc;
                    chemSorProcessTime_.Total_begin_to_sync_ = recordGPUCalc;
                    //HurMPI::allReduce(cpuTime1, MPI_MAX);
                    //HurMPI::allReduce(cpuTime2, MPI_MAX);
                    //HurMPI::allReduce(cpuTime3, MPI_MAX);
                    //HurMPI::allReduce(cpuTime12, MPI_MAX);

                    //HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
                    //HurMPI::allReduce(recordTotalCalc, MPI_MAX);
                    //HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
                    //HurMPI::allReduce(recordGPUCalc, MPI_MAX);
                    chemSorProcessTime_.GPU_ChemicalSource_ = (real)gpuTime1 / 1000.0;
                    /*real gputime11 = (real)gpuTime1 / 1000.0;
                    HurMPI::allReduce(gputime11, MPI_MAX);
                    if (HurMPI::master()) {
                        Pout << "    Info: testing processing time in "
                                "\"calculateSource\""
                             << std::endl;
                        Pout << "          CPU_FV_turbSource    : " << cpuTime1 << "[s]"
                             << std::endl;
                        Pout << "          CPU_sync             : " << cpuTime2 << "[s]"
                             << std::endl;
                        Pout << "          CPU_convert data     : " << cpuTime3 << "[s]"
                             << std::endl;
                        Pout << "          GPU_ChemicalSource   : " << gputime11 << "[s]"
                             << std::endl;
                        Pout << "          Total_calc_source    : " << recordTotalCalc << "[s]"
                             << std::endl;
                        Pout << "          Total_before_calc    : " << recordGPUBeforeCalc << "[s]"
                             << std::endl;
                        Pout << "          Total_GPUCal_to_sync : " << recordGPUOnlyCalc << "[s]"
                             << std::endl;
                        Pout << "          Total_begin_to_sync  : " << recordGPUCalc << "[s]"
                             << std::endl;
                        Pout << "          Total_Fvturb_to_sync : " << cpuTime12 << "[s]"
                             << std::endl;
                    }*/
#endif // TEST_PROCESS_TIME
                }
            }
        } else {
            turbPtr_->impSource();
        }
    } else if (timeMarcingPtr_->fullJacobianSource()) {
        LFatal("This solver does not support full Jacobian point-implicit operation");
    }

    if (sorcTermPtr_) {
        const integer nsp = mixtures().species().size();
        sorcTerm().addSourceTerms(rho());
        sorcTerm().addSourceTerms(rho(), v());
        sorcTerm().addSourceTerms(rho(), E());
        if (turbPtr_->isCoupled()) {
            for (integer i = 0; i < turbPtr_->nEq(); ++i) {
                sorcTerm().addSourceTerms(rho(), turbPtr_->var(i));
            }
        }
        for (integer i = 0; i < nsp - 1; ++i) {
            sorcTerm().addSourceTerms(rho(), yi()[i]);
        }
    }
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::updatePrimitives(const bool shouldUpdateTemp) {
    mixtures().lastSpeAndNormalized();

    solver::updatePrimitives(shouldUpdateTemp);

    limits();
    updateProperties();
    if (turbPtr_->isCoupled()) {
        turbPtr_->limitAndUpdateBoundary();
        turbPtr_->update();
    }
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::calculateOtherTerms() {
    while (turbPtr_->loop()) {
        turbPtr_->solving(dt_);
        turbPtr_->limitAndUpdateBoundary();
        turbPtr_->update();
    }
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::updateFlowOld() {}

void OpenHurricane::turbulentSpeciesSolverCUDATran::write() {
#ifdef TEST_PROCESS_TIME
    if (!mixtures().noReaction()) {
        if (iter().cStep() != 0) {
            chemSorProcessTime_.writeTecplot(iter().cStep(), chemSorFOS_);
        }
    }
    if (iter().cStep() != 0) {
        gasProProcessTime_.writeTecplot(iter().cStep(), gasProFOS_);
    }
#endif // TEST_PROCESS_TIME
    if (iter().solWrite().writeNow()) {
        if (!mixtures().noReaction()) {
            addSolOutput(iter(), "Damkohler", chemtryPtr_->calcDamkohler);

            if (iter().solWrite().found("GFO")) {
                if (chemtryPtr_->fuelName().size() == 0) {
                    checkWarning("The fuel is not given for compute flame "
                                 "index, and would not compute");
                } else if (chemtryPtr_->oxygenName().size() == 0) {
                    checkWarning("The oxygen is not given for compute flame "
                                 "index, and would not compute");
                } else {
                    iter().solWrite().setOutputField(
                        "GFO", flowPtr_->GFO(chemtryPtr_->fuelName(), chemtryPtr_->oxygenName()));
                }
            }

            if (iter().solWrite().found("nGFO")) {
                if (chemtryPtr_->fuelName().size() == 0) {
                    checkWarning("The fuel is not given for compute flame "
                                 "index, and would not compute");
                } else if (chemtryPtr_->oxygenName().size() == 0) {
                    checkWarning("The oxygen is not given for compute flame "
                                 "index, and would not compute");
                } else {
                    iter().solWrite().setOutputField(
                        "nGFO", flowPtr_->nGFO(chemtryPtr_->fuelName(), chemtryPtr_->oxygenName()));
                }
            }
            addSolOutput(iter(), "flameIndex", chemtryPtr_->flameIndex);
            addSolOutput(iter(), "heatReleaseRate", chemtryPtr_->heatReleaseRate);
            addSolOutput(iter(), "tcFRR", chemtryPtr_->tcFRR);
            addSolOutput(iter(), "tcSFR", chemtryPtr_->tcSFR);
            addSolOutput(iter(), "tcGSPR", chemtryPtr_->tcGSPR);
            addSolOutput(iter(), "tcJacDT", chemtryPtr_->tcJacDT);
        }

        addSolOutput(iter(), "Ret", turbPtr_->Ret);
        addSolOutput(iter(), "KolmogorovLengthScale", turbPtr_->KolmogorovLengthScale);
        addSolOutput(iter(), "KolmogorovTimeScale", turbPtr_->KolmogorovTimeScale);
        addSolOutput(iter(), "KolmogorovVelocityScale", turbPtr_->KolmogorovVelocityScale);
        addSolOutput(iter(), "integralLengthScale", turbPtr_->integralLengthScale);
        addSolOutput(iter(), "integralTimeScale", turbPtr_->integralTimeScale);
        addSolOutput(iter(), "integralVelocityScale", turbPtr_->integralVelocityScale);

        addSolOutput(iter(), "turbulentKineticEnergy", turbPtr_->k);
        addSolOutput(iter(), "turbulentDissipationRate", turbPtr_->epsilon);

        addSolOutputWithOutName(
            iter(), "ReynoldsStressTensor",
            string("\"ReynoldsStress_uu\",\"ReynoldsStress_uv\",\"ReynoldsStress_"
                   "uw\","
                   "\"ReynoldsStress_vv\",\"ReynoldsStress_vw\","
                   "\"ReynoldsStress_ww\""),
            turbPtr_->ReynoldsStressTensor);
    }
    iter_.write();
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::getYiTP(real *hur_restrict hostYiTPPtr_) const {
    const auto nsp = specTable().size();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp; ++i) {
            hostYiTPPtr_[n * (nsp + 2) + i] = yi()[i][n];
        }
        hostYiTPPtr_[n * (nsp + 2) + nsp] = T()[n];
        hostYiTPPtr_[n * (nsp + 2) + nsp + 1] = p()[n];
    }
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::getYiTPSlice(real *hur_restrict hostYiTPPtr_,
                                                             const integer pos,
                                                             const integer offset) const {
    const auto nsp = specTable().size();
    for (integer n = 0; n < offset; ++n) {
        for (integer i = 0; i < nsp; ++i) {
            hostYiTPPtr_[n * (nsp + 2) + i] = yi()[i][n + pos];
        }
        hostYiTPPtr_[n * (nsp + 2) + nsp] = T()[n + pos];
        hostYiTPPtr_[n * (nsp + 2) + nsp + 1] = p()[n + pos];
    }
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::makeFzBoundStart() const {
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

void OpenHurricane::turbulentSpeciesSolverCUDATran::muKappaDiffBoundary() {
    /* HurGPU::setSharedMenBankSize();
     HurGPU::setCachePreferShared();*/
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

void OpenHurricane::turbulentSpeciesSolverCUDATran::getYiTPBoundary(
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

void OpenHurricane::turbulentSpeciesSolverCUDATran::setExtrapolate(
    const real *hur_restrict hostDimmMuKap, const real *hur_restrict hostHaiCp) {
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

void OpenHurricane::turbulentSpeciesSolverCUDATran::transferTranP() {
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

void OpenHurricane::turbulentSpeciesSolverCUDATran::gasPropertiesBySlice(const integer pos,
                                                                     const integer offset,
                                                                     const bool firstCall) {
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
#endif // TEST_PROCESS_TIME

    nTBForTran_.setBlockAndGridSize(specTable().size(), offset);
    nTBForHaiCp_.setBlockAndGridSize(specTable().size(), offset);

    nTBForTran_.setSharedMemSize1<real>(2 * specTable().size() * nTBForTran_.blockSize().y);
    nTBForHaiCp_.setSharedMemSize1<real>(specTable().size() * nTBForHaiCp_.blockSize().y);

    // Third step: (1) Alloc memory in the host (CPU) for array "hostYiTPPtr", "hostDimmMuKap" and "hostHai".
    //             The pointer "hostYiTPPtr" is used for storing yi (species mass fraction), T (temperature) and p (pressure).
    //             The pointer "hostDimmMuKap" is used for storing Dim, mu and kappa.
    //             The pointer "hostHai" is used for storing hai.
    //             (2) Alloc memory in the device (GPU) for array "dYiTP", "dDimmMuKap" and "dHai".
    real *hostYiTPPtr = nullptr;
    const integer nsp = flowPtr_->mixtures().species().size();
    size_t sizeOfH_YiTP = (static_cast<unsigned long long>(nsp) + 2) * offset * sizeof(real);

    checkCUDAGPUError(cudaHostAlloc((void **)&hostYiTPPtr, sizeOfH_YiTP, cudaHostAllocDefault));

    getYiTPSlice(hostYiTPPtr, pos, offset);
    real *hostDimmMuKap = nullptr;
    size_t sizeOfH_DimmMuKap = (static_cast<unsigned long long>(nsp) + 2) * offset * sizeof(real);

    checkCUDAGPUError(
        cudaHostAlloc((void **)&hostDimmMuKap, sizeOfH_DimmMuKap, cudaHostAllocDefault));

    real *hostHaiCP = nullptr;
    size_t sizeOfH_HaiCP = (static_cast<unsigned long long>(nsp + 1)) * offset * sizeof(real);

    checkCUDAGPUError(cudaHostAlloc((void **)&hostHaiCP, sizeOfH_HaiCP, cudaHostAllocDefault));

    cu2DArray<cu_real> dYiTP(nsp + 2, offset);
    cu2DArray<cu_real> dDimmMuKap(nsp + 2, offset);
    cu2DArray<cu_real> dHaiCp(nsp + 1, offset);

#ifdef TEST_PROCESS_TIME
    auto recordGPUBeforeCalc = recordGPU.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

    // Fourth step: Compute the hai in the device (GPU)
    CUDAThermo::calchaicpAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA().species(),
                               nTBForHaiCp_, nsp, offset, hostHaiCP, dHaiCp, streamTran);

    // Fifth step: Compute the Dim, mu and kappa in the device (GPU)
    CUDATransport::calcTransportPropertiesAsync(hostYiTPPtr, dYiTP, mixtures().transportCUDA(),
                                                nTBForTran_, nsp, offset, hostDimmMuKap, dDimmMuKap,
                                                streamTran, false);

#ifdef TEST_PROCESS_TIME
    myEnd.record(streamTran);
    hrClock cpuStart;
#endif // TEST_PROCESS_TIME
    if (firstCall) {
        // Sixth step: Calculate boundary conditions in the host (CPU)
        mixtures().gamma(p(), T(), gama(), true);
        bc();
        mixtures().E(p(), T(), v(), E(), false, true);
        mixtures().gammaBoundary(p(), T(), gama(), true);
    }
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

    dYiTP.clear();
    dHaiCp.clear();
    dDimmMuKap.clear();
    streamTran.destroy();

#ifdef TEST_PROCESS_TIME
    myStart.destroy();
    myEnd.destroy();
#endif // TEST_PROCESS_TIME

    checkCUDAGPUError(cudaFreeHost(hostYiTPPtr));
    for (integer n = 0; n < offset; ++n) {
        for (integer i = 0; i < nsp; ++i) {
            mixtures().hi(i)[n + pos] = hostHaiCP[n * (nsp + 1) + i];
        }
        cp()[n + pos] = hostHaiCP[n * (nsp + 1) + nsp];
    }
    for (integer n = 0; n < offset; ++n) {
        for (integer i = 0; i < nsp; ++i) {
            Diff()[i][n + pos] = hostDimmMuKap[n * (nsp + 2) + i];
        }
        mu()[n + pos] = hostDimmMuKap[n * (nsp + 2) + nsp];
        kappal()[n + pos] = hostDimmMuKap[n * (nsp + 2) + nsp + 1];
    }
    checkCUDAGPUError(cudaFreeHost(hostHaiCP));
    checkCUDAGPUError(cudaFreeHost(hostDimmMuKap));

#ifdef TEST_PROCESS_TIME
    auto cpuTime3 = cpuStart.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

#ifdef TEST_PROCESS_TIME
    //auto cpuTime4 = cpuStart.clockTimeIncrement();
    auto recordTotalCalc = recordGPU.elapsedClockTime();
    auto cpuTime12 = cpuTime1 + cpuTime2;
    gasProProcessTime_.CPU_bc_gamma_E_ += cpuTime1;
    gasProProcessTime_.CPU_sync_ += cpuTime2;
    gasProProcessTime_.CPU_convert_data_ += cpuTime3;
    //gasProProcessTime_.CPU_bc_HaiMuKDim_ += cpuTime4;
    gasProProcessTime_.Total_Fvturb_to_sync_ += cpuTime12;
    gasProProcessTime_.Total_before_calc_ += recordGPUBeforeCalc;
    gasProProcessTime_.Total_calc_Propeeer_ += recordTotalCalc;
    gasProProcessTime_.Total_GPUCal_to_sync_ += recordGPUOnlyCalc;
    gasProProcessTime_.Total_begin_to_sync_ += recordGPUCalc;
    //HurMPI::allReduce(cpuTime1, MPI_MAX);
    //HurMPI::allReduce(cpuTime2, MPI_MAX);
    //HurMPI::allReduce(cpuTime3, MPI_MAX);
    //HurMPI::allReduce(cpuTime4, MPI_MAX);
    //HurMPI::allReduce(cpuTime12, MPI_MAX);

    //HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
    //HurMPI::allReduce(recordTotalCalc, MPI_MAX);
    //HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
    //HurMPI::allReduce(recordGPUCalc, MPI_MAX);
    gasProProcessTime_.GPU_HaiMuKDim_ += (real)gpuTime1 / 1000.0;

    /*real gputime11 = (real)gpuTime1 / 1000.0;
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
    }*/

#endif // TEST_PROCESS_TIME
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::updatePropertiesBySlice(const integer nSlice) {
#ifdef TEST_PROCESS_TIME
    gasProProcessTime_.clear();
#endif // TEST_PROCESS_TIME
    mixtures().transportCUDA();

    integer pos = 0;
    integer offset = mesh().nCells() / nSlice;
    for (integer i = 0; i < nSlice; ++i) {
        gasPropertiesBySlice(pos, offset, i == 0);
        pos += offset;
        if (i + 1 == nSlice - 1) {
            offset = mesh().nCells() - pos;
        }
    }

    mixtures().destroySpeciesCUDA();
    mixtures().destroyTransportCUDA();

#ifdef TEST_PROCESS_TIME
    hrClock CPUStart;
#endif // TEST_PROCESS_TIME

    muKappaDiffBoundary();

#ifdef TEST_PROCESS_TIME
    gasProProcessTime_.CPU_bc_HaiMuKDim_ = CPUStart.elapsedClockTime();
#endif // TEST_PROCESS_TIME
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::chemicalSourceBySlice(const integer pos,
                                                                      const integer offset,
                                                                      const bool firstCall) {
    const integer nsp = flowPtr_->mixtures().species().size();
    //HurGPU::setSharedMenBankSize();
    //HurGPU::setCachePreferShared();
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

    // Third step: Alloc memory in the host (CPU) for array "hostYiRhoTPtr", "RRi" and "dRRidrhoyi".
    //             The pointer "hostYiRhoTPtr" is used for storing yi (species mass fraction), rho (density) and T (temperature).
    //             The pointer "RRi" is used for storing Ri (chemical source terms).
    //             The pointer "dRRidrhoyi" is used for storing the diagonal elements of the Jacobian matrix of the chemical source terms.
    real *hostYiRhoTPtr = nullptr;
    size_t sizeOfH_YiRhoT = (static_cast<unsigned long long>(nsp) + 2) * offset * sizeof(real);
    checkCUDAGPUError(cudaHostAlloc((void **)&hostYiRhoTPtr, sizeOfH_YiRhoT, cudaHostAllocDefault));

    real *RRi = nullptr;
    size_t sizeOfRRi = static_cast<unsigned long long>(nsp) * offset * sizeof(real);
    checkCUDAGPUError(cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

    // Fourth step: Alloc memory in the device (GPU) for array "dYiRhoT", "dRi" and "device_dRidrhoyi".
    //             The pointer "dYiRhoT" is used for storing yi (species mass fraction), rho (density) and T (temperature).
    //             The pointer "dRi" is used for storing Ri (chemical source terms).
    //             The pointer "dRi" is used for storing Ri the diagonal elements of the Jacobian matrix of the chemical source terms.
    cu2DArray<cu_real> dYiRhoT(nsp + 2, offset);
    cu2DArray<cu_real> dRi(nsp, offset);

#ifdef TEST_PROCESS_TIME
    auto recordGPUBeforeCalc = recordGPU.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

    // Fifth step: Compute the chemical source terms in the device (GPU)
    chemtryPtr_->chemistry().calculateSourceTermsAsyncSlice(hostYiRhoTPtr, dYiRhoT, RRi, dRi,
                                                            stream, pos, offset);

#ifdef TEST_PROCESS_TIME
    myEnd.record(stream);
    hrClock cpuStart;
#endif // TEST_PROCESS_TIME

    if (firstCall) {
        // Sixth step: Calculate the turbulence source and viscous flux in the host (CPU)
        turbPtr_->impSource();
        actualCalFv();
    }

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

    dYiRhoT.clear();
    dRi.clear();
    stream.destroy();

#ifdef TEST_PROCESS_TIME
    myStart.destroy();
    myEnd.destroy();
#endif // TEST_PROCESS_TIME

    checkCUDAGPUError(cudaFreeHost(hostYiRhoTPtr));

    for (integer i = 0; i < nsp; ++i) {
        for (integer n = 0; n < offset; ++n) {
            chemtryPtr_->chemistry().Ri()[i][n + pos] = RRi[n * nsp + i];
        }
    }
    checkCUDAGPUError(cudaFreeHost(RRi));

#ifdef TEST_PROCESS_TIME
    auto cpuTime3 = cpuStart.clockTimeIncrement();
    auto recordTotalCalc = recordGPU.elapsedClockTime();
    auto cpuTime12 = cpuTime1 + cpuTime2;
    chemSorProcessTime_.CPU_FV_turbSource_ += cpuTime1;
    chemSorProcessTime_.CPU_sync_ += cpuTime2;
    chemSorProcessTime_.CPU_convert_data_ += cpuTime3;
    chemSorProcessTime_.Total_Fvturb_to_sync_ += cpuTime12;
    chemSorProcessTime_.Total_before_calc_ += recordGPUBeforeCalc;
    chemSorProcessTime_.Total_calc_source_ += recordTotalCalc;
    chemSorProcessTime_.Total_GPUCal_to_sync_ += recordGPUOnlyCalc;
    chemSorProcessTime_.Total_begin_to_sync_ += recordGPUCalc;
    //HurMPI::allReduce(cpuTime1, MPI_MAX);
    //HurMPI::allReduce(cpuTime2, MPI_MAX);
    //HurMPI::allReduce(cpuTime3, MPI_MAX);
    //HurMPI::allReduce(cpuTime12, MPI_MAX);

    //HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
    //HurMPI::allReduce(recordTotalCalc, MPI_MAX);
    //HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
    //HurMPI::allReduce(recordGPUCalc, MPI_MAX);
    chemSorProcessTime_.GPU_ChemicalSource_ += (real)gpuTime1 / 1000.0;
    /*real gputime11 = (real)gpuTime1 / 1000.0;
    HurMPI::allReduce(gputime11, MPI_MAX);
    if (HurMPI::master()) {
        Pout << "    Info: testing processing time in "
                "\"calculateSource\""
             << std::endl;
        Pout << "          CPU_FV_turbSource    : " << cpuTime1 << "[s]" << std::endl;
        Pout << "          CPU_sync             : " << cpuTime2 << "[s]" << std::endl;
        Pout << "          CPU_convert data     : " << cpuTime3 << "[s]" << std::endl;
        Pout << "          GPU_ChemicalSource   : " << gputime11 << "[s]" << std::endl;
        Pout << "          Total_calc_source    : " << recordTotalCalc << "[s]" << std::endl;
        Pout << "          Total_before_calc    : " << recordGPUBeforeCalc << "[s]" << std::endl;
        Pout << "          Total_GPUCal_to_sync : " << recordGPUOnlyCalc << "[s]" << std::endl;
        Pout << "          Total_begin_to_sync  : " << recordGPUCalc << "[s]" << std::endl;
        Pout << "          Total_Fvturb_to_sync : " << cpuTime12 << "[s]" << std::endl;
    }*/
#endif // TEST_PROCESS_TIME
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::chemicalSourceImpBySlice(
    const integer pos, const integer offset, const bool firstCall, const bool isHybridPrecision) {
    const integer nsp = flowPtr_->mixtures().species().size();
    if (!isHybridPrecision) {
        //HurGPU::setSharedMenBankSize();
        //HurGPU::setCachePreferShared();
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

        // Third step: Alloc memory in the host (CPU) for array "hostYiRhoTPtr", "RRi" and "dRRidrhoyi".
        //             The pointer "hostYiRhoTPtr" is used for storing yi (species mass fraction), rho (density) and T (temperature).
        //             The pointer "RRi" is used for storing Ri (chemical source terms).
        //             The pointer "dRRidrhoyi" is used for storing the diagonal elements of the Jacobian matrix of the chemical source terms.
        real *hostYiRhoTPtr = nullptr;
        size_t sizeOfH_YiRhoT = (static_cast<unsigned long long>(nsp) + 2) * offset * sizeof(real);
        checkCUDAGPUError(
            cudaHostAlloc((void **)&hostYiRhoTPtr, sizeOfH_YiRhoT, cudaHostAllocDefault));

        real *RRi = nullptr;
        size_t sizeOfRRi = static_cast<unsigned long long>(nsp) * offset * sizeof(real);
        checkCUDAGPUError(cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

        real *dRRidrhoyi = nullptr;
        size_t sizeOfdRRi = (static_cast<unsigned long long>(nsp) - 1) * offset * sizeof(real);
        checkCUDAGPUError(cudaHostAlloc((void **)&dRRidrhoyi, sizeOfdRRi, cudaHostAllocDefault));

        // Fourth step: Alloc memory in the device (GPU) for array "dYiRhoT", "dRi" and "device_dRidrhoyi".
        //             The pointer "dYiRhoT" is used for storing yi (species mass fraction), rho (density) and T (temperature).
        //             The pointer "dRi" is used for storing Ri (chemical source terms).
        //             The pointer "dRi" is used for storing Ri the diagonal elements of the Jacobian matrix of the chemical source terms.
        cu2DArray<cu_real> dYiRhoT(nsp + 2, offset);
        cu2DArray<cu_real> dRi(nsp, offset);
        cu2DArray<cu_real> device_dRidrhoyi(nsp - 1, offset);

#ifdef TEST_PROCESS_TIME
        auto recordGPUBeforeCalc = recordGPU.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

        // Fifth step: Compute the chemical source terms in the device (GPU)
        chemtryPtr_->chemistry().calculateSourceTermsImpAsyncSlice(
            hostYiRhoTPtr, dYiRhoT, RRi, dRi, dRRidrhoyi, device_dRidrhoyi, stream, pos, offset);

#ifdef TEST_PROCESS_TIME
        myEnd.record(stream);
        hrClock cpuStart;
#endif // TEST_PROCESS_TIME

        if (firstCall) {
            // Sixth step: Calculate the turbulence source and viscous flux in the host (CPU)
            turbPtr_->impSource();
            actualCalFv();
        }

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
            for (integer n = 0; n < offset; ++n) {
                chemtryPtr_->chemistry().Ri()[i][n + pos] = RRi[n * nsp + i];
                if (i < nsp - 1) {
                    yi()[i].diagSource()[n + pos] = dRRidrhoyi[n * (nsp - 1) + i];
                }
            }
        }
        checkCUDAGPUError(cudaFreeHost(RRi));
        checkCUDAGPUError(cudaFreeHost(dRRidrhoyi));

#ifdef TEST_PROCESS_TIME
        auto cpuTime3 = cpuStart.clockTimeIncrement();
        auto recordTotalCalc = recordGPU.elapsedClockTime();
        auto cpuTime12 = cpuTime1 + cpuTime2;
        chemSorProcessTime_.CPU_FV_turbSource_ += cpuTime1;
        chemSorProcessTime_.CPU_sync_ += cpuTime2;
        chemSorProcessTime_.CPU_convert_data_ += cpuTime3;
        chemSorProcessTime_.Total_Fvturb_to_sync_ += cpuTime12;
        chemSorProcessTime_.Total_before_calc_ += recordGPUBeforeCalc;
        chemSorProcessTime_.Total_calc_source_ += recordTotalCalc;
        chemSorProcessTime_.Total_GPUCal_to_sync_ += recordGPUOnlyCalc;
        chemSorProcessTime_.Total_begin_to_sync_ += recordGPUCalc;
        //HurMPI::allReduce(cpuTime1, MPI_MAX);
        //HurMPI::allReduce(cpuTime2, MPI_MAX);
        //HurMPI::allReduce(cpuTime3, MPI_MAX);
        //HurMPI::allReduce(cpuTime12, MPI_MAX);

        //HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
        //HurMPI::allReduce(recordTotalCalc, MPI_MAX);
        //HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
        //HurMPI::allReduce(recordGPUCalc, MPI_MAX);
        chemSorProcessTime_.GPU_ChemicalSource_ += (real)gpuTime1 / 1000.0;

        /*real gputime11 = (real)gpuTime1 / 1000.0;
        HurMPI::allReduce(gputime11, MPI_MAX);
        if (HurMPI::master()) {
            Pout << "    Info: testing processing time in "
                    "\"calculateSource\""
                 << std::endl;
            Pout << "          CPU_FV_turbSource    : " << cpuTime1 << "[s]" << std::endl;
            Pout << "          CPU_sync             : " << cpuTime2 << "[s]" << std::endl;
            Pout << "          CPU_convert data     : " << cpuTime3 << "[s]" << std::endl;
            Pout << "          GPU_ChemicalSource   : " << gputime11 << "[s]" << std::endl;
            Pout << "          Total_calc_source    : " << recordTotalCalc << "[s]" << std::endl;
            Pout << "          Total_before_calc    : " << recordGPUBeforeCalc << "[s]"
                 << std::endl;
            Pout << "          Total_GPUCal_to_sync : " << recordGPUOnlyCalc << "[s]" << std::endl;
            Pout << "          Total_begin_to_sync  : " << recordGPUCalc << "[s]" << std::endl;
            Pout << "          Total_Fvturb_to_sync : " << cpuTime12 << "[s]" << std::endl;
        }*/
#endif // TEST_PROCESS_TIME
    } else {
        //HurGPU::setSharedMenBankSize();
        //HurGPU::setCachePreferShared();
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

        // Third step: Alloc memory in the host (CPU) for array "hostYiRhoTPtr", "RRi" and "dRRidrhoyi".
        //             The pointer "hostYiRhoTPtr" is used for storing yi (species mass fraction), rho (density) and T (temperature).
        //             The pointer "RRi" is used for storing Ri (chemical source terms).
        //             The pointer "dRRidrhoyi" is used for storing the diagonal elements of the Jacobian matrix of the chemical source terms.
        real *hostYiRhoTPtr = nullptr;
        size_t sizeOfH_YiRhoT = (static_cast<unsigned long long>(nsp) + 2) * offset * sizeof(real);
        checkCUDAGPUError(
            cudaHostAlloc((void **)&hostYiRhoTPtr, sizeOfH_YiRhoT, cudaHostAllocDefault));

        real *RRi = nullptr;
        size_t sizeOfRRi = static_cast<unsigned long long>(nsp) * offset * sizeof(real);
        checkCUDAGPUError(cudaHostAlloc((void **)&RRi, sizeOfRRi, cudaHostAllocDefault));

        cu_float *dRRidrhoyi = nullptr;
        size_t sizeOfdRRi =
            (static_cast<unsigned long long>(nsp) - 1) * offset * sizeof(cu_float);
        checkCUDAGPUError(cudaHostAlloc((void **)&dRRidrhoyi, sizeOfdRRi, cudaHostAllocDefault));

        // Fourth step: Alloc memory in the device (GPU) for array "dYiRhoT", "dRi" and "device_dRidrhoyi".
        //             The pointer "dYiRhoT" is used for storing yi (species mass fraction), rho (density) and T (temperature).
        //             The pointer "dRi" is used for storing Ri (chemical source terms).
        //             The pointer "dRi" is used for storing Ri the diagonal elements of the Jacobian matrix of the chemical source terms.
        cu2DArray<cu_real> dYiRhoT(nsp + 2, offset);
        cu2DArray<cu_real> dRi(nsp, offset);
        cu2DArray<cu_float> device_dRidrhoyi(nsp - 1, offset);

#ifdef TEST_PROCESS_TIME
        auto recordGPUBeforeCalc = recordGPU.clockTimeIncrement();
#endif // TEST_PROCESS_TIME

        // Fifth step: Compute the chemical source terms in the device (GPU)
        chemtryPtr_->chemistry().calculateSourceTermsImpAsyncHybridSlice(
            hostYiRhoTPtr, dYiRhoT, RRi, dRi, dRRidrhoyi, device_dRidrhoyi, stream, pos, offset);

#ifdef TEST_PROCESS_TIME
        myEnd.record(stream);
        hrClock cpuStart;
#endif // TEST_PROCESS_TIME

        if (firstCall) {
            // Sixth step: Calculate the turbulence source and viscous flux in the host (CPU)
            turbPtr_->impSource();
            actualCalFv();
        }

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
            for (integer n = 0; n < offset; ++n) {
                chemtryPtr_->chemistry().Ri()[i][n + pos] = RRi[n * nsp + i];
                if (i < nsp - 1) {
                    yi()[i].diagSource()[n + pos] =
                        static_cast<real>(dRRidrhoyi[n * (nsp - 1) + i]);
                }
            }
        }
        checkCUDAGPUError(cudaFreeHost(RRi));
        checkCUDAGPUError(cudaFreeHost(dRRidrhoyi));

#ifdef TEST_PROCESS_TIME
        auto cpuTime3 = cpuStart.clockTimeIncrement();
        auto recordTotalCalc = recordGPU.elapsedClockTime();
        auto cpuTime12 = cpuTime1 + cpuTime2;
        chemSorProcessTime_.CPU_FV_turbSource_ += cpuTime1;
        chemSorProcessTime_.CPU_sync_ += cpuTime2;
        chemSorProcessTime_.CPU_convert_data_ += cpuTime3;
        chemSorProcessTime_.Total_Fvturb_to_sync_ += cpuTime12;
        chemSorProcessTime_.Total_before_calc_ += recordGPUBeforeCalc;
        chemSorProcessTime_.Total_calc_source_ += recordTotalCalc;
        chemSorProcessTime_.Total_GPUCal_to_sync_ += recordGPUOnlyCalc;
        chemSorProcessTime_.Total_begin_to_sync_ += recordGPUCalc;
        //HurMPI::allReduce(cpuTime1, MPI_MAX);
        //HurMPI::allReduce(cpuTime2, MPI_MAX);
        //HurMPI::allReduce(cpuTime3, MPI_MAX);
        //HurMPI::allReduce(cpuTime12, MPI_MAX);

        //HurMPI::allReduce(recordGPUBeforeCalc, MPI_MAX);
        //HurMPI::allReduce(recordTotalCalc, MPI_MAX);
        //HurMPI::allReduce(recordGPUOnlyCalc, MPI_MAX);
        //HurMPI::allReduce(recordGPUCalc, MPI_MAX);
        chemSorProcessTime_.GPU_ChemicalSource_ += (real)gpuTime1 / 1000.0;

        /*real gputime11 = (real)gpuTime1 / 1000.0;
        HurMPI::allReduce(gputime11, MPI_MAX);
        if (HurMPI::master()) {
            Pout << "    Info: testing processing time in "
                    "\"calculateSource\""
                 << std::endl;
            Pout << "          CPU_FV_turbSource    : " << cpuTime1 << "[s]" << std::endl;
            Pout << "          CPU_sync             : " << cpuTime2 << "[s]" << std::endl;
            Pout << "          CPU_convert data     : " << cpuTime3 << "[s]" << std::endl;
            Pout << "          GPU_ChemicalSource   : " << gputime11 << "[s]" << std::endl;
            Pout << "          Total_calc_source    : " << recordTotalCalc << "[s]" << std::endl;
            Pout << "          Total_before_calc    : " << recordGPUBeforeCalc << "[s]"
                 << std::endl;
            Pout << "          Total_GPUCal_to_sync : " << recordGPUOnlyCalc << "[s]" << std::endl;
            Pout << "          Total_begin_to_sync  : " << recordGPUCalc << "[s]" << std::endl;
            Pout << "          Total_Fvturb_to_sync : " << cpuTime12 << "[s]" << std::endl;
        }*/
#endif // TEST_PROCESS_TIME
    }
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::updateChemSourceBySlice(const integer nSlice) {
#ifdef TEST_PROCESS_TIME
    chemSorProcessTime_.clear();
#endif // TEST_PROCESS_TIME

    chemtryPtr_->chemistry().createReactionTable();
    integer pos = 0;
    integer offset = mesh().nCells() / nSlice;
    for (integer i = 0; i < nSlice; ++i) {
        chemicalSourceBySlice(pos, offset, i == 0);
        pos += offset;
        if (i + 1 == nSlice - 1) {
            offset = mesh().nCells() - pos;
        }
    }
    chemtryPtr_->getChemistrySource();
    chemtryPtr_->chemistry().destroyReactionTable();
}

void OpenHurricane::turbulentSpeciesSolverCUDATran::updateChemSourceImpBySlice(
    const integer nSlice, const bool isHybridPrecision) {
#ifdef TEST_PROCESS_TIME
    chemSorProcessTime_.clear();
#endif // TEST_PROCESS_TIME

    chemtryPtr_->chemistry().createReactionTable();
    integer pos = 0;
    integer offset = mesh().nCells() / nSlice;
    for (integer i = 0; i < nSlice; ++i) {
        chemicalSourceImpBySlice(pos, offset, i == 0, isHybridPrecision);
        pos += offset;
        if (i + 1 == nSlice - 1) {
            offset = mesh().nCells() - pos;
        }
    }
    chemtryPtr_->getImpChemistrySource(dt_, !iter().isGlobalTimeStep());
    chemtryPtr_->chemistry().destroyReactionTable();
}

#endif // CUDA_PARALLEL