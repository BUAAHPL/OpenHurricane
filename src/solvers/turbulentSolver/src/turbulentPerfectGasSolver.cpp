/*!
 * \file turbulentPerfectGasSolver.cpp
 * \brief Main subroutines of turbulent Perfect Gas Solver.
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

#include "turbulentPerfectGasSolver.hpp"
#include "SpalartAllmaras.hpp"
#include "calculateFieldVar.hpp"
#include "solutionWrite.hpp"
#include "viscousFlux.hpp"

namespace OpenHurricane {
    createClassName(turbulentPerfectGasSolver);
         registerObjFty(solver, turbulentPerfectGasSolver, controller);
} // namespace OpenHurricane

OpenHurricane::turbulentPerfectGasSolver::turbulentPerfectGasSolver(iteration &iter,
                                                                const runtimeMesh &mesh)
    : solver(iter, mesh), SAModelPtr_(nullptr), SSTModel_(nullptr) {
    if (iter.cont().found("flow")) {
        const auto &flowCont = iter.cont().subController("flow");
        if (flowCont.found("turbulence")) {
            auto &turbCont = flowCont.subController("turbulence");
            turbPtr_ = turbulenceModel::creator(turbCont, *flowPtr_);
        } else {
            LFatal("Cannot find turbulence setting section in flow section");
        }
    }
    //kappal().setNoWrite();
}

void OpenHurricane::turbulentPerfectGasSolver::solving() {
    // First step: Register the primitives parameter in the time marching method.
    marching().addObject(rho());
    marching().addObject(v());
    marching().addObject(E());
    if (turbPtr_->isCoupled()) {
        for (integer i = 0; i < turbPtr_->nEq(); ++i) {
            marching().addObject(turbPtr_->var(i));
        }
    }
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
    //
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

void OpenHurricane::turbulentPerfectGasSolver::BDFSolve() {
    // First step: Register the primitives parameter in the time marching method.
    marching().addObject(rho());
    marching().addObject(v());
    marching().addObject(E());
    for (integer i = 0; i < turbPtr_->nEq(); ++i) {
        marching().addObject(turbPtr_->var(i));
    }
    marching().initializing();

    // Second step: Initialize the whole flow field
    initialize();

    solver::BDFSolve();
    //marching().marching();
}

void OpenHurricane::turbulentPerfectGasSolver::clear() noexcept {
    flowPtr_.clear();
    /*if (flowPtr_ != nullptr)
    {
        delete flowPtr_;
        flowPtr_ = nullptr;
    }*/
    invFluxPtr_.clear();
    /*if (invFluxPtr_ != nullptr)
    {
        delete invFluxPtr_;
        invFluxPtr_ = nullptr;
    }*/
    HurDelete(SAModelPtr_);
    /*if (SAModelPtr_ != nullptr)
    {
        delete SAModelPtr_;
        SAModelPtr_ = nullptr;
    }*/
    HurDelete(SSTModel_);
    /*if (SSTModel_ != nullptr)
    {
        delete SSTModel_;
        SSTModel_ = nullptr;
    }*/
}

void OpenHurricane::turbulentPerfectGasSolver::bc() {
    rho().updateBoundary();
    v().updateBoundary();
    p().updateBoundary();
    T().updateBoundary();
    // E().updateBoundary();

    // 20210318 ���� updatePrimitives()
    if (turbPtr_->isCoupled()) {
        turbPtr_->updateBoundary();
    }

    // 20210322 ��˼�� ȡ������solver::bc()����
    // solver::bc()�����еģ�RiemannValue::updated_ = false;������������solver::iterRefresh()
    //solver::bc();
}

void OpenHurricane::turbulentPerfectGasSolver::updateProperties() {
    mixtures().E(p(), T(), v(), E(), false, true);
    mixtures().gammaBoundary(p(), T(), gama(), true);

    //mixtures().mu(p(), T(), mu());
    mixtures().muKappaCp(rho(), p(), T(), mu(), kappal(), cp(), false, true);

    //mixtures().cp(rho(), T(), cp());

    // 20210318 ���� updatePrimitives()
    /*if (turbPtr_->isCoupled())
    {
        turbPtr_->update();
    }*/
}

void OpenHurricane::turbulentPerfectGasSolver::initialize() {
    solver::initialize();
    PtrList<cellRealArray> tmpCellArray;
    if (!iter().restart()) {
        real rhoi = rho().initialize();
        vector vi = v().initialize();
        real pi = p().initialize();
        real Ti = T().initialize();

        for (integer celli = 0; celli < mesh().nCells(); ++celli) {
            rho()[celli] = mixtures().rho(p()[celli], T()[celli], celli);
        }
    } else {
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
    mixtures().gammaBoundary(p(), T(), gama());

    mixtures().E(p(), T(), v(), E(), false, true);
    //mixtures().gamma(p(), T(), gama());

    //mixtures().mu(p(), T(), mu());
    mixtures().muKappaCp(rho(), p(), T(), mu(), kappal(), cp());
    /*updateProperties();*/

    if (!iter().restart()) {
        turbPtr_->initialize();
    } else {
        turbPtr_->initializeRestart();
    }
}

void OpenHurricane::turbulentPerfectGasSolver::timeStep(realArray &dt) {}

void OpenHurricane::turbulentPerfectGasSolver::calculateFc() {
    // Caculate the inviscous fluxes for the continuity, momentum and energy equations
    invFluxPtr_->basicFlux();

    // Calculate the convective fluxes for the turbulent equations for turbulent coupling solver.
    if (turbPtr_->isCoupled()) {
        for (integer i = 0; i < turbPtr_->nEq(); ++i) {
            invFluxPtr_->invFlux(turbPtr_->var(i),
                                 realBounded::lowerBound(real(0), realBounded::BOUNDED_FLAG));
        }
    }

    // Calculate the gradient of turbulent variables for splitting solver.
    turbPtr_->calcGrad(*invFluxPtr_);
}

void OpenHurricane::turbulentPerfectGasSolver::calculateFv() {
    const faceTensorArray deltafV(fv::gradf(v()));

    const faceRealArray cpf(fv::interpolate(cp()));

    const faceRealArray muf(fv::interpolate(mu()));
    const faceRealArray mutf(fv::interpolate(mut()));
    faceVectorArray Vf(fv::interpolate(v()));
    const faceRealArray rhof(fv::interpolate(rho()));

    //(muf + mutf) * (twoSymm(deltafV) - (real(2.0 / 3.0) * (div(diagToVector(deltafV))) * I))
    const faceSymmTensorArray tau(turbPtr_->tauEff(rhof, muf, mutf, deltafV));

    //v().rhs() += visFlux(tau);
    visFlux(tau, v());

    turbPtr_->visFlux(rhof, muf, mutf, mu(), mut(), rho());

    const faceRealArray kappaf(fv::interpolate(kappal()));
    const faceRealArray kappatf(object("kappatf", mesh(), object::NOT_WRITE, object::TEMPORARY),
                                mesh(), mutf * cpf / prt());

    Vf = tau * Vf;
    Vf += ((kappaf + kappatf) * fv::gradf(T()));

    //E().rhs() += visFlux<real>(Vf);
    visFlux<real>(Vf, E());

    turbPtr_->correctEnergyEqVisFlux(E());
}

void OpenHurricane::turbulentPerfectGasSolver::calculateSource() {
   
    if (timeMarcingPtr_->explicitSource()) {
        turbPtr_->expSource();
    } else if (timeMarcingPtr_->diagonalImpSource()) {
        turbPtr_->impSource();
    }

    if (sorcTermPtr_ ) {
        sorcTerm().addSourceTerms(rho());
        sorcTerm().addSourceTerms(rho(), v());
        sorcTerm().addSourceTerms(rho(), E());
        if (turbPtr_->isCoupled()) {
            for (integer i = 0; i < turbPtr_->nEq(); ++i) {
                sorcTerm().addSourceTerms(rho(), turbPtr_->var(i));
            }
        }
    }
}

void OpenHurricane::turbulentPerfectGasSolver::updatePrimitives(const bool shouldUpdateTemp) {
    solver::updatePrimitives(shouldUpdateTemp);
    limits();
    mixtures().gamma(p(), T(), gama(), true);
    bc();
    updateProperties();

    // 20210318 ��֤�������������ٸ�����������
    if (turbPtr_->isCoupled()) {
        turbPtr_->limitAndUpdateBoundary();
        turbPtr_->update();
    }
}

void OpenHurricane::turbulentPerfectGasSolver::calculateOtherTerms() {
    while (turbPtr_->loop()) {
        turbPtr_->solving(dt_);
        turbPtr_->limitAndUpdateBoundary();
        turbPtr_->update();
    }
}

void OpenHurricane::turbulentPerfectGasSolver::updateFlowOld() {}

//void OpenHurricane::turbulentPerfectGasSolver::calcOutput() const
//{
//    OpenHurricane::solver::calcOutput();
//    calcViscousRatio();
//    calcVorticity();
//    calcOtherVorticity();
//}

void OpenHurricane::turbulentPerfectGasSolver::write() {
    if (iter().solWrite().writeNow()) {
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
