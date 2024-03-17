/*!
 * \file EulerSpeciesSolver.cpp
 * \brief Main subroutines of Euler species reacting Solver.
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

#include "EulerSpeciesSolver.hpp"
#include "laminar.hpp"
#include "laminarFlow.hpp"
#include "viscousFlux.hpp"

#include "calculateFieldVar.hpp"
#include "solutionWrite.hpp"

namespace OpenHurricane {
    createClassName(EulerSpeciesSolver);
    registerObjFty(solver, EulerSpeciesSolver, controller);
} // namespace OpenHurricane

OpenHurricane::EulerSpeciesSolver::EulerSpeciesSolver(iteration &iter, const runtimeMesh &mesh)
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
}

void OpenHurricane::EulerSpeciesSolver::solving() {
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

void OpenHurricane::EulerSpeciesSolver::BDFSolve() {
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

void OpenHurricane::EulerSpeciesSolver::clear() noexcept {
    flowPtr_.clear();
    invFluxPtr_.clear();
    chemtryPtr_.clear();
    
}

void OpenHurricane::EulerSpeciesSolver::bc() {
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

void OpenHurricane::EulerSpeciesSolver::timeStep(realArray &dt) {}

void OpenHurricane::EulerSpeciesSolver::updateProperties() {
    mixtures().E(p(), T(), v(), E(), false, true);
    //mixtures().muKappa(rho(), p(), T(), mu(), kappal());
    //mixtures().gamma(p(), T(), gama());
    mixtures().gammaBoundary(p(), T(), gama(), true);

    //flowPtr_->mixtures().Diff(p(), T(), Diff());
    flowPtr_->mixtures().updateHai(p(), T(), false, true);
}

void OpenHurricane::EulerSpeciesSolver::initialize() {
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

void OpenHurricane::EulerSpeciesSolver::calculateFc() {
    // The size of species.
    const integer nsp = mixtures().species().size();

    // Caculate the inviscous fluxes for the continuity, momentum and energy equations
    invFluxPtr_->basicFlux();

    // Calculate the convective fluxes for the species equations.    
    invFluxPtr_->invFluxSpecies(yi(), false);

    // Calculate the gradient of last species.
    invFluxPtr_->grad(yi()[nsp - 1]);
}

void OpenHurricane::EulerSpeciesSolver::calculateFv() {}

void OpenHurricane::EulerSpeciesSolver::calculateSource() {

    // To get chemical source
    if (!mixtures().noReaction()) {
        chemtryPtr_->evaluateSource(*timeMarcingPtr_, rhoId_, rhouId_, rhoEId_, rhoYi0Id_);
    }
       
    // To get other source
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

void OpenHurricane::EulerSpeciesSolver::updatePrimitives(const bool shouldUpdateTemp) {
    mixtures().lastSpeAndNormalized();

    solver::updatePrimitives(shouldUpdateTemp);

    limits();
    mixtures().gamma(p(), T(), gama(), true);
    bc();

    updateProperties();
}

void OpenHurricane::EulerSpeciesSolver::updateFlowOld() {}

void OpenHurricane::EulerSpeciesSolver::write() {

    if (iter().solWrite().writeNow()) {
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
