/*!
 * \file turbulentSpeciesSolverCUDA.cpp
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
#ifdef CUDA_PARALLEL
#include "turbulentSpeciesSolverCUDA.hpp"
#include "calculateFieldVar.hpp"
#include "cudaEvents.hpp"
#include "solutionWrite.hpp"
#include "viscousFlux.hpp"

namespace OpenHurricane {
    createClassName(turbulentSpeciesSolverCUDA);
         registerObjFty(solver, turbulentSpeciesSolverCUDA, controller);
} // namespace OpenHurricane

OpenHurricane::turbulentSpeciesSolverCUDA::turbulentSpeciesSolverCUDA(iteration &iter,
                                                                  const runtimeMesh &mesh)
    : solver(iter, mesh), chemtryPtr_(nullptr), rhoId_(0), rhouId_(1), rhoEId_(4), rhoYi0Id_(6),
      rhoTurb0Id_(5) {
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

    if (HurGPU::GPUDev().size() != 1) {
        LFatal("Only support single GPU mode in this class");
    }
    HurGPU::setDevice(HurGPU::GPUDev()[0]);
}

void OpenHurricane::turbulentSpeciesSolverCUDA::solving() {
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
    solver::solving();
    //
    //    marching().initializing();
    //
    //#ifdef TEST_PROCESS_TIME
    //    testProcessTime myTestPT
    //    (
    //        iter(),
    //        iter().cont().subController("testProcessTime")
    //    );
    //#endif // TEST_PROCESS_TIME
    //    // Second step: Initialize the whole flow field
    //    initialize();
    //
    //    marching().timeStep();
    //    // Third step: Begin the iteration.
    //    while (iter_.iterating())
    //    {
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.start(iter_.cStep());
    //#endif // TEST_PROCESS_TIME
    //
    //        // Refresh
    //        iterRefresh();
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.clockTimeIncrement("Refresh");
    //#endif // TEST_PROCESS_TIME
    //
    //        calculateFc();// Firstly, compute the convective flux Fc
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.clockTimeIncrement("ConvectiveFlux");
    //#endif // TEST_PROCESS_TIME
    //
    //        calculateSource();// Secondly, evaluate the source terms S
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.clockTimeIncrement("Source");
    //#endif // TEST_PROCESS_TIME
    //
    //        calculateFv();// Finally, calculate the viscous flux Fv
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.clockTimeIncrement("ViscousFlux");
    //#endif // TEST_PROCESS_TIME
    //
    //        marching().marching();
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.clockTimeIncrement("Marching");
    //#endif // TEST_PROCESS_TIME
    //
    //        // Update flow field
    //        updatePrimitives();
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.clockTimeIncrement("UpdatePrimitives");
    //#endif // TEST_PROCESS_TIME
    //
    //        calculateOtherTerms();
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.clockTimeIncrement("OtherTerms");
    //#endif // TEST_PROCESS_TIME
    //
    //        // Write solution
    //        write();
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.clockTimeIncrement("Write");
    //#endif // TEST_PROCESS_TIME
    //
    //        marching().timeStep();
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.clockTimeIncrement("TimeStep");
    //#endif // TEST_PROCESS_TIME
    //
    //#ifdef TEST_PROCESS_TIME
    //        myTestPT.stop();
    //#endif // TEST_PROCESS_TIME
    //    }
}

void OpenHurricane::turbulentSpeciesSolverCUDA::BDFSolve() {
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

    solver::BDFSolve();
    //marching().marching();
}

void OpenHurricane::turbulentSpeciesSolverCUDA::bc() {
    const integer nsp = flowPtr_->mixtures().species().size();

    rho().updateBoundary();
    v().updateBoundary();
    p().updateBoundary();

    for (integer isp = 0; isp < nsp; ++isp) {
        yi()[isp].updateBoundary();
    }

    T().updateBoundary();
}

void OpenHurricane::turbulentSpeciesSolverCUDA::updateProperties() {
    mixtures().E(p(), T(), v(), E(), false, true);
    mixtures().gammaBoundary(p(), T(), gama(), true);
    mixtures().gasProperties(rho(), p(), T(), mu(), kappal(), cp(), Diff(), hi(), false, true);
    /*flowPtr_->mixtures().muKappaCp(rho(), p(), T(), mu(), kappal(), cp(), false, true);
    flowPtr_->mixtures().Diff(p(), T(), Diff(), false, true);
    flowPtr_->mixtures().updateHai(p(), T(), false, true);*/
}

void OpenHurricane::turbulentSpeciesSolverCUDA::initialize() {
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
    flowPtr_->mixtures().muKappaCp(rho(), p(), T(), mu(), kappal(), cp());
    //mixtures().gamma(p(), T(), gama());
    flowPtr_->mixtures().Diff(p(), T(), Diff());
    //updateProperties();
    flowPtr_->mixtures().updateHai(p(), T(), false, true);
    if (!iter().restart()) {
        turbPtr_->initialize();
    } else {       
        turbPtr_->initializeRestart();
    }
}

void OpenHurricane::turbulentSpeciesSolverCUDA::timeStep(realArray &dt) {}

void OpenHurricane::turbulentSpeciesSolverCUDA::calculateFc() {
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

    // Calculate the convective fluxes for the species equations.
    invFluxPtr_->invFluxSpecies(yi(), false);

    // Calculate the gradient of turbulent variables for splitting solver.
    turbPtr_->calcGrad(*invFluxPtr_);

    // Calculate the gradient of last species.
    invFluxPtr_->grad(yi()[nsp - 1]);
}

void OpenHurricane::turbulentSpeciesSolverCUDA::calculateFv() {}

void OpenHurricane::turbulentSpeciesSolverCUDA::actualCalFv() {
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

void OpenHurricane::turbulentSpeciesSolverCUDA::calculateSource() {
    
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
        } else {
            turbPtr_->expSource();
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

            // Sixth step: Calculate the turbulence source and viscous flux in the host (CPU)
            turbPtr_->impSource();
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
        } else {
            turbPtr_->impSource();
        }
    } else if (timeMarcingPtr_->fullJacobianSource()) {
        LFatal("This solver does not support full Jacobian point-implicit operation");
    }

    if (sorcTermPtr_ ) {
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

void OpenHurricane::turbulentSpeciesSolverCUDA::updatePrimitives(const bool shouldUpdateTemp) {
    mixtures().lastSpeAndNormalized();
    solver::updatePrimitives(shouldUpdateTemp);
    limits();
    mixtures().gamma(p(), T(), gama(), true);
    bc();
    updateProperties();

    if (turbPtr_->isCoupled()) {
        turbPtr_->limitAndUpdateBoundary();
        turbPtr_->update();
    }
}

void OpenHurricane::turbulentSpeciesSolverCUDA::calculateOtherTerms() {
    while (turbPtr_->loop()) {
        turbPtr_->solving(dt_);
        turbPtr_->limitAndUpdateBoundary();
        turbPtr_->update();
    }
}

void OpenHurricane::turbulentSpeciesSolverCUDA::updateFlowOld() {}

void OpenHurricane::turbulentSpeciesSolverCUDA::write() {
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

        addSolOutputWithOutName(iter(), "ReynoldsStressTensor",
                                string("\"ReynoldsStress_uu\",\"ReynoldsStress_uv\",\"ReynoldsStress_"
                                     "uw\","
                                     "\"ReynoldsStress_vv\",\"ReynoldsStress_vw\","
                                     "\"ReynoldsStress_ww\""),
                                turbPtr_->ReynoldsStressTensor);
    }
    iter_.write();
}

hur_deprecated void OpenHurricane::turbulentSpeciesSolverCUDA::testChemical() {}
#endif // CUDA_PARALLEL