/*!
 * \file turbulentSpeciesSolver.cpp
 * \brief Main subroutines of turbulent Species Solver Solver.
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

#include "turbulentSpeciesSolver.hpp"
#include "calculateFieldVar.hpp"
#include "solutionWrite.hpp"
#include "viscousFlux.hpp"

namespace OpenHurricane {
    createClassName(turbulentSpeciesSolver);
    registerObjFty(solver, turbulentSpeciesSolver, controller);
} // namespace OpenHurricane

OpenHurricane::turbulentSpeciesSolver::turbulentSpeciesSolver(iteration &iter, const runtimeMesh &mesh)
    : solver(iter, mesh), chemtryPtr_(nullptr), rhoId_(0), rhouId_(1), rhoEId_(4), rhoYi0Id_(6),
      rhoTurb0Id_(5)
#ifdef TEST_PROCESS_TIME
      ,
      fosTestChemTime_(IOsstream::streamFormat::ASCII_FORMAT, IOsstream::openOptions::ONLY_MASTER,
                       std::ios_base::out)
#endif
      ,
      sorTime_(0), otherTime_(0), calcLoadWeight_(false)
#ifdef CUDA_PARALLEL
      ,
      fzBoundStartPtr_(nullptr)
#endif //CUDA_PARALLEL
{
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
                chemtryPtr_ = combustionModel::creator(*flowPtr_, reacCont, *turbPtr_);
            } else {
                LFatal("Cannot find mixture setting section in flow section");
            }
        }
    }
#ifdef TEST_PROCESS_TIME
    if (!mixtures().noReaction()) {
        fileName outFile = iter.configName().name(true) + "" + "TestChemTime" + ".dat";
        outFile = iter.outputPath() / outFile;
        fosTestChemTime_.changeFileName(outFile);

        fosTestChemTime_.open();
    } else {
        fosTestChemTime_.close();
    }
#endif // TEST_PROCESS_TIME

    if (iter.cont().subController("iteration").found("cellLoadWeight")) {
        const auto &clwCont =
            iter.cont().subController("iteration").subController("cellLoadWeight");

        if (clwCont.found("weightTypes")) {
            const auto wtw = clwCont.findWord("weightTypes");

            if (wtw == "DAC") {
                calcLoadWeight_ = true;
            }
        }
    }
#ifdef CUDA_PARALLEL
    if (HurGPU::useGPU()) {
        if (HurGPU::GPUDev().size() != 1) {
            LFatal("Only support single GPU mode in this class");
        }
        HurGPU::setDevice(HurGPU::GPUDev()[0]);
    }
#endif //CUDA_PARALLEL
}

void OpenHurricane::turbulentSpeciesSolver::solving() {
    const integer nsp = flowPtr_->mixtures().species().size();
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

    HurMPI::barrier();

    hrClock myClock;

    solver::solving();
    HurMPI::barrier();

    Pout << "    Finish solving time: " << myClock.elapsedClockTime() << std::endl;

    //  marching().initializing();

    //  // Second step: Initialize the whole flow field
    //  initialize();
    //  marching().timeStep();
    //  Pout << "    Begin iteration..." << std::endl;

    //  // Third step: Begin the iteration.
    //  while (iter_.iterating())
    //  {
    //      // Refresh
    //      iterRefresh();

    //      previousSource();

    //// Firstly, compute the convective flux Fc
    //      calculateFc();

    //      // Secondly, evaluate the source terms S
    //      calculateSource();

    //      // Finally, calculate the viscous flux Fv
    //      calculateFv();

    //      marching().marching();

    //      // Update flow field
    //      updatePrimitives(!timeMarcingPtr_->updatedTemperature());

    //      postSource();

    //      calculateOtherTerms();

    //      // Write solution
    //      write();

    //      marching().timeStep();
    //  }
}

void OpenHurricane::turbulentSpeciesSolver::BDFSolve() {
    const integer nsp = flowPtr_->mixtures().species().size();

    // First step: Register the primitives parameter in the time marching method.
    marching().addObject(rho());
    marching().addObject(v());
    marching().addObject(E());
    for (integer i = 0; i < turbPtr_->nEq(); ++i) {
        marching().addObject(turbPtr_->var(i));
    }

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
#ifdef HUR_DEBUG
    write();
#endif // HUR_DEBUG

    solver::BDFSolve();
    //marching().marching();
}

void OpenHurricane::turbulentSpeciesSolver::bc() {
    const integer nsp = flowPtr_->mixtures().species().size();

    rho().updateBoundary();
    v().updateBoundary();
    p().updateBoundary();

    for (integer isp = 0; isp < nsp; ++isp) {
        yi()[isp].updateBoundary();
    }

    T().updateBoundary();
}

void OpenHurricane::turbulentSpeciesSolver::updateProperties() {
    mixtures().E(p(), T(), v(), E(), false, true);
    mixtures().gammaBoundary(p(), T(), gama(), true);
    mixtures().gasProperties(rho(), p(), T(), mu(), kappal(), cp(), Diff(), hi(), false, true);
    /*flowPtr_->mixtures().muKappaCp(rho(), p(), T(), mu(), kappal(), cp(), false, true);
    flowPtr_->mixtures().Diff(p(), T(), Diff(), false, true);
    flowPtr_->mixtures().updateHai(p(), T(), false, true);*/
}

void OpenHurricane::turbulentSpeciesSolver::initialize() {
    solver::initialize();
    const integer nsp = flowPtr_->mixtures().species().size();
    PtrList<cellRealArray> tmpCellArray;
    if (!iter().restart()) {
        real rhoi = rho().initialize();
        vector vi = v().initialize();
        real pi = p().initialize();
        real Ti = T().initialize();
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

    mixtures().E(p(), T(), v(), E(), false, true);
    flowPtr_->mixtures().muKappaCp(rho(), p(), T(), mu(), kappal(), cp(), false, true);

    flowPtr_->mixtures().Diff(p(), T(), Diff(), false, true);

    flowPtr_->mixtures().updateHai(p(), T(), false, true);
    if (!iter().restart()) {
        turbPtr_->initialize();
    } else {
        turbPtr_->initializeRestart();
    }

#ifdef TEST_PROCESS_TIME
    if (!mixtures().noReaction()) {
        chemtryPtr_->chemistry().writeTimeTitle(fosTestChemTime_);
    }
#endif // TEST_PROCESS_TIME
}

void OpenHurricane::turbulentSpeciesSolver::timeStep(realArray &dt) {}

void OpenHurricane::turbulentSpeciesSolver::calculateFc() {
    hrClock myClocks;
    // The size of species.
    const integer nsp = flowPtr_->mixtures().species().size();

    // Caculate the inviscous fluxes for the continuity, momentum and energy equations
    invFluxPtr_->basicFlux();

    // Calculate the convective fluxes for the turbulent equations for turbulent coupling solver.
    if (turbPtr_->isCoupled()) {
        for (integer i = 0; i < turbPtr_->nEq(); ++i) {
            invFluxPtr_->invFlux(turbPtr_->var(i),
                                 realBounded::lowerBound(real(0), realBounded::BOUNDED_FLAG));
        }
    }

    invFluxPtr_->invFluxSpecies(yi(), false);

    // Calculate the gradient of turbulent variables for splitting solver.
    turbPtr_->calcGrad(*invFluxPtr_);

    // Calculate the gradient of last species.
    invFluxPtr_->grad(yi()[nsp - 1]);

    if (isUpdatingCalcTime()) {
        HurMPI::barrier();
        otherTime_ = myClocks.elapsedClockTime();
    }
}

void OpenHurricane::turbulentSpeciesSolver::calculateFv() {
    hrClock myClocks;
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

    if (isUpdatingCalcTime()) {
        HurMPI::barrier();
        otherTime_ += myClocks.elapsedClockTime();
    }
}

void OpenHurricane::turbulentSpeciesSolver::calculateSource() {  

    // To get turbulent source
    if (timeMarcingPtr_->explicitSource()) {
        turbPtr_->expSource();
        /*if (!mixtures().noReaction())
        {
            chemtryPtr_->expChemistrySource(dt_);
        }*/
    } else if (timeMarcingPtr_->diagonalImpSource()) {
        turbPtr_->impSource();
        /*if (!mixtures().noReaction())
        {
            chemtryPtr_->impChemistrySource(dt_, !iter().isGlobalTimeStep());
        }*/
    } else if (timeMarcingPtr_->fullJacobianSource() ||
               timeMarcingPtr_->fullJacobianSourceTable()) {
        turbPtr_->fullImpSource(timeMarcingPtr_->Jacobian(), rhoId_, rhoTurb0Id_);
       
    }

    // To get chemical source
    if (!mixtures().noReaction()) {
        hrClock myClocks;
        chemtryPtr_->evaluateSource(*timeMarcingPtr_, rhoId_, rhouId_, rhoEId_, rhoYi0Id_);
        if (isUpdatingCalcTime()) {
            HurMPI::barrier();
            sorTime_ += myClocks.elapsedClockTime();
        }
    }

    // To get other source
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

void OpenHurricane::turbulentSpeciesSolver::updatePrimitives(const bool shouldUpdateTemp) {
    hrClock myClocks;
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
    if (turbPtr_->isCoupled()) {
        turbPtr_->limitAndUpdateBoundary();
        turbPtr_->update();
    }

    if (isUpdatingCalcTime()) {
        HurMPI::barrier();
        otherTime_ = myClocks.elapsedClockTime();
    }
}

void OpenHurricane::turbulentSpeciesSolver::calculateOtherTerms() {
    while (turbPtr_->loop()) {
        turbPtr_->solving(dt_);
        turbPtr_->limitAndUpdateBoundary();
        turbPtr_->update();
    }
}

void OpenHurricane::turbulentSpeciesSolver::updateFlowOld() {}

void OpenHurricane::turbulentSpeciesSolver::write() {

#ifdef TEST_PROCESS_TIME
    if (!mixtures().noReaction()) {
        if (iter().cStep() != 0) {
            chemtryPtr_->chemistry().writeSolveODETime(fosTestChemTime_);
        }
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
                    checkWarning(
                        "The fuel is not given for compute flame index, and would not compute");
                } else if (chemtryPtr_->oxygenName().size() == 0) {
                    checkWarning(
                        "The oxygen is not given for compute flame index, and would not compute");
                } else {
                    iter().solWrite().setOutputField(
                        "nGFO", flowPtr_->nGFO(chemtryPtr_->fuelName(), chemtryPtr_->oxygenName()));
                }
            }

            if (iter().solWrite().found("mixtureFractionZ")) {
                chemtryPtr_->calcMixtureFraction();
            }
            if (iter().solWrite().found("mixtureDisspationRate")) {
                iter().solWrite().setOutputField("mixtureDisspationRate",
                                                 chemtryPtr_->mixtureDisspationRate(*invFluxPtr_));
            }

            addSolOutput(iter(), "flameIndex", chemtryPtr_->flameIndex);
            addSolOutput(iter(), "heatReleaseRate", chemtryPtr_->heatReleaseRate);
            addSolOutput(iter(), "tcFRR", chemtryPtr_->tcFRR);
            addSolOutput(iter(), "tcSFR", chemtryPtr_->tcSFR);
            addSolOutput(iter(), "tcGSPR", chemtryPtr_->tcGSPR);
            addSolOutput(iter(), "tcJacDT", chemtryPtr_->tcJacDT);

            setCellLoadWeight(sorTime_, otherTime_);
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

hur_deprecated void OpenHurricane::turbulentSpeciesSolver::testChemical() {}

void OpenHurricane::turbulentSpeciesSolver::setCellLoadWeight(const real sorTime,
                                                          const real otherTime) {
    if (calcLoadWeight_ && iter().solWrite().writeNow()) {
        if (!mixtures().noReaction()) {
            integer weightF1 = 2;
            integer weightF2 = 1;
            if (iter().cont().subController("iteration").found("cellLoadWeight")) {
                const auto &clwCont =
                    iter().cont().subController("iteration").subController("cellLoadWeight");
                weightF1 = clwCont.findOrDefault<integer>("weightFactor1", weightF1);
                weightF2 = clwCont.findOrDefault<integer>("weightFactor2", weightF2);
            }
            real ot = otherTime;
            HurMPI::allReduce(ot, MPI_SUM);
            integer nc = mesh().nCells();
            HurMPI::allReduce(nc, MPI_SUM);

            real st = sorTime;
            HurMPI::allReduce(st, MPI_SUM);

            ot /= nc;
            st /= nc;

            const real nff = st / ot;
            const auto &csct = chemtryPtr_->chemistry().cellSourceCalTime();

            auto &cellLWgt = const_cast<runtimeMesh &>(mesh()).cellLoadWeights();

            integer maxWgt = 0;
            for (integer n = 0; n < cellLWgt.size(); ++n) {
                cellLWgt[n] = integer(weightF1 + weightF2 * nff * integer(csct[n] / ot));
                maxWgt = max(maxWgt, cellLWgt[n]);
            }
            HurMPI::allReduce(maxWgt, MPI_MAX);
            if (maxWgt > 500) {
                real fct = 500.0 / ((real)maxWgt);
                for (integer n = 0; n < cellLWgt.size(); ++n) {
                    cellLWgt[n] = max(integer(1), (integer(fct * cellLWgt[n])));
                }
            }

            if (iter().solWrite().found("cellLoadWeight")) {
                realArray cwgt(mesh().nTotalCells(), Zero);
                for (integer n = 0; n < mesh().nCells(); ++n) {
                    cwgt[n] = cellLWgt[n];
                }
                iter().solWrite().setOutputField("cellLoadWeight", cwgt);
            }
        }
    }
}

hur_nodiscard bool OpenHurricane::turbulentSpeciesSolver::isUpdatingCalcTime() const noexcept {
    return calcLoadWeight_ && iter().solWrite().writeNow();
}

#ifdef CUDA_PARALLEL
#include "thermoListCUDA.hpp"
#include "transportListCUDA.hpp"

#ifdef TEST_PROCESS_TIME
#include "cudaEvents.hpp"
#endif // TEST_PROCESS_TIME

void OpenHurricane::turbulentSpeciesSolver::updatePropertiesCUDA() {
    nThreadsAndBlocks nTBForTran_, nTBForHaiCp_;
    nTBForTran_.setBlockAndGridSize(specTable().size(), mesh().nCells());
    nTBForHaiCp_.setBlockAndGridSize(specTable().size(), mesh().nCells());

    nTBForTran_.setSharedMemSize1<real>(2 * specTable().size() * nTBForTran_.blockSize().y);
    nTBForHaiCp_.setSharedMemSize1<real>(specTable().size() * nTBForTran_.blockSize().y);
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

void OpenHurricane::turbulentSpeciesSolver::getYiTP(real *hur_restrict hostYiTPPtr_) const {
    const auto nsp = specTable().size();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        for (integer i = 0; i < nsp; ++i) {
            hostYiTPPtr_[n * (nsp + 2) + i] = yi()[i][n];
        }
        hostYiTPPtr_[n * (nsp + 2) + nsp] = T()[n];
        hostYiTPPtr_[n * (nsp + 2) + nsp + 1] = p()[n];
    }
}

void OpenHurricane::turbulentSpeciesSolver::muKappaDiffBoundary() {
    nThreadsAndBlocks nTBForTran_;
    nTBForTran_.setBlockAndGridSize(specTable().size(), mesh().nCells());
    nTBForTran_.setSharedMemSize1<real>(2 * specTable().size() * nTBForTran_.blockSize().y);

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

void OpenHurricane::turbulentSpeciesSolver::makeFzBoundStart() const {
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

void OpenHurricane::turbulentSpeciesSolver::getYiTPBoundary(real *hur_restrict hostYiTPPtr_) const {
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

void OpenHurricane::turbulentSpeciesSolver::setExtrapolate(const real *hur_restrict hostDimmMuKap,
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

void OpenHurricane::turbulentSpeciesSolver::transferTranP() {
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
#endif // CUDA_PARALLEL