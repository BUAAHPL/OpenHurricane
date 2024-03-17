/*!
 * \file laminarPerfectGasSolver.cpp
 * \brief Main subroutines of laminar Perfect Gas Solver.
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

#include "laminarPerfectGasSolver.hpp"
#include "laminar.hpp"
#include "laminarFlow.hpp"
#include "viscousFlux.hpp"

#include "calculateFieldVar.hpp"
#include "writeFieldVars.hpp"
namespace OpenHurricane {
    createClassName(laminarPerfectGasSolver);
         registerObjFty(solver, laminarPerfectGasSolver, controller);
} // namespace OpenHurricane

OpenHurricane::laminarPerfectGasSolver::laminarPerfectGasSolver(iteration &iter,
                                                            const runtimeMesh &mesh)
    : solver(iter, mesh) {

    //kappal().setNoWrite();
}

void OpenHurricane::laminarPerfectGasSolver::solving() {
    // First step: Register the primitives parameter in the time marching method.
    marching().addObject(rho());
    marching().addObject(v());
    marching().addObject(E());

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

void OpenHurricane::laminarPerfectGasSolver::BDFSolve() {
    // First step: Register the primitives parameter in the time marching method.
    marching().addObject(rho());
    marching().addObject(v());
    marching().addObject(E());
    marching().initializing();

    // Second step: Initialize the whole flow field
    initialize();

    solver::BDFSolve();
    //marching().marching();
}

void OpenHurricane::laminarPerfectGasSolver::clear() noexcept {
    flowPtr_.clear();
    invFluxPtr_.clear();
    /*if (flowPtr_ != nullptr)
    {
        delete flowPtr_;
        flowPtr_ = nullptr;
    }
    if (invFluxPtr_ != nullptr)
    {
        delete invFluxPtr_;
        invFluxPtr_ = nullptr;
    }*/
}

void OpenHurricane::laminarPerfectGasSolver::bc() {
    rho().updateBoundary();
    v().updateBoundary();
    p().updateBoundary();
    T().updateBoundary();
    //E().updateBoundary();

    // 20210322 ��˼�� ȡ������solver::bc()����
    // solver::bc()�����еģ�RiemannValue::updated_ = false;������������solver::iterRefresh()
    //solver::bc();
}

void OpenHurricane::laminarPerfectGasSolver::updateProperties() {
    mixtures().E(p(), T(), v(), E(), false, true);
    mixtures().gammaBoundary(p(), T(), gama(), true);

    /*mixtures().mu(p(), T(), mu());

    mixtures().cp(rho(), T(), cp());*/

    mixtures().muKappaCp(rho(), p(), T(), mu(), kappal(), cp(), false, true);
}

void OpenHurricane::laminarPerfectGasSolver::initialize() {
    solver::initialize();
    if (!iter().restart()) {
        real rhoi = rho().initialize();
        vector vi = v().initialize();
        real pi = p().initialize();
        real Ti = T().initialize();

        real Ei = E().initialize();


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
    mixtures().gammaBoundary(p(), T(), gama(), true);

    mixtures().E(p(), T(), v(), E(), false, true);
    updateProperties();
}

void OpenHurricane::laminarPerfectGasSolver::timeStep(realArray &dt) {}

void OpenHurricane::laminarPerfectGasSolver::calculateFc() {
    // Caculate the inviscous fluxes for the continuity, momentum and energy equations.
    invFluxPtr_->basicFlux();
}

void OpenHurricane::laminarPerfectGasSolver::calculateFv() {
    faceTensorArray deltafV(fv::gradf(v()));

    const faceRealArray cpf(fv::interpolate(cp()));
    const faceRealArray muf(fv::interpolate(mu()));
    faceVectorArray Vf(fv::interpolate(v()));

    faceSymmTensorArray tau(
        object("tau", mesh(), object::NOT_WRITE, object::TEMPORARY), mesh(),
        muf * (twoSymm(deltafV) - (real(2.0 / 3.0) * (div(diagToVector(deltafV))) * I)));

    v().rhs() += visFlux(tau);

    const faceRealArray kappaf(fv::interpolate(kappal()));

    Vf = tau * Vf;
    Vf += (kappaf * fv::gradf(T()));

    E().rhs() += visFlux<real>(Vf);
}

void OpenHurricane::laminarPerfectGasSolver::calculateSource() {
    if (sorcTermPtr_ ) {
        sorcTerm().addSourceTerms(rho());
        sorcTerm().addSourceTerms(rho(), v());
        sorcTerm().addSourceTerms(rho(), E());
    }
}

void OpenHurricane::laminarPerfectGasSolver::updatePrimitives(const bool shouldUpdateTemp) {
    solver::updatePrimitives(shouldUpdateTemp);
    limits();
    mixtures().gamma(p(), T(), gama(), true);
    bc();
    updateProperties();
}

void OpenHurricane::laminarPerfectGasSolver::updateFlowOld() {}

void OpenHurricane::laminarPerfectGasSolver::write() {
    iter_.write();
}
