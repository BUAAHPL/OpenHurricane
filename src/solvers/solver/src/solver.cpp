/*!
 * \file solver.cpp
 * \brief Main subroutines of base class of solver.
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

#include "solver.hpp"
#include "AUSM.hpp"
#include "HLL.hpp"
#include "HLLC.hpp"
#include "LUSGS.hpp"
#include "SpalartAllmaras.hpp"
#include "viscousFlux.hpp"
#include "viscousThermo.hpp"
#include "writeFaceZone.hpp"
#include "writeFieldVars.hpp"

namespace OpenHurricane {
    createClassName(solver);
    createObjFty(solver, controller);
} // namespace OpenHurricane

void OpenHurricane::solver::makeSorcTerm(const controller &sorcCont) {
    if (!sorcTermPtr_) {
        sorcTermPtr_.reset(new sourceTerms(*flowPtr_, iter_, sorcCont));
    }
}

OpenHurricane::solver::solver(iteration &iter, const runtimeMesh &mesh)
    : reqAutoPatch_(false), autoPatch_(), iter_(iter), flowPtr_(nullptr), invFluxPtr_(nullptr),
      timeMarcingPtr_(nullptr), turbPtr_(nullptr),
      dt_(object("timeStep", mesh, object::WRITE_OUTPUT), mesh, Zero), nLowTemp_(0), nHighTemp_(0),
      nLowPress_(0), nHighPress_(0), hasTPLmtCells_(false), timeType_(timeTypes::UsualSolver),
      sorcTermPtr_(nullptr), useLowMachPrecon_(false), a4Ptr_(nullptr), a5Ptr_(nullptr) {
    if (iter.cont().found("flow")) {
        flowPtr_ = flowModel::creator(mesh, iter.cont().subController("flow"));
    }

    if (iter.cont().found("spatialScheme")) {
        if (!flowPtr_.isNull()) {
            invFluxPtr_ =
                spatialScheme::creator(iter.cont().subController("spatialScheme"), mesh, *flowPtr_);
        } else {
            LFatal("The flowPtr_ is null pointer");
        }
    }
    if (iter.cont().found("iteration")) {
        const controller &iterCont = iter.cont().subController("iteration");
        timeMarcingPtr_ = timeMarching::creator(iterCont, mesh, *flowPtr_, *this, flowPtr_->v());
    }
    setAutoPatch();

    if (iter.cont().found("sourceTerms")) {
        const auto &sorcCont = iter.cont().subController("sourceTerms");
        makeSorcTerm(sorcCont);
    }
}

void OpenHurricane::solver::setAutoPatch() {
    if (!iter_.cont().found("initialization")) {
        return;
    }

    const auto &autoCont = iter_.cont().subController("initialization");
    if (!autoCont.found("autoPatch")) {
        return;
    }
    std::string rnl = autoCont.findText("autoPatch");
    replaceAllMarks(rnl, "\n", " ");
    if (!rnl.empty()) {
        size_t pos = 0;
        stdStringList rll;
        split(rnl, rll, ",");
        for (integer i = 0; i < rll.size(); ++i) {
            trim(rll[i]);
            if (!autoCont.found(rll[i])) {
                checkWarningStr(("Cannot find region: " + rll[i]));
                continue;
            }
            const auto &regionCont = autoCont.subController(rll[i]);
            autoPatch_.append(new patching(mesh(), regionCont));
        }
    }
}

void OpenHurricane::solver::solving() {
    marching().initializing();

#ifdef TEST_PROCESS_TIME
    fileName resOut;
    if (iter().cont().found("testProcessTime")) {
        const auto &testCont = iter().cont().subController("testProcessTime");
        string testN = testCont.findWord("fileName");
        trim(testN);
        fileName outFile = testN;
        if (!outFile.isAbsolute()) {
            outFile = iter().outputPath() / outFile;
        }

        resOut = outFile;
    } else {
        const auto &cfn = iter().configName();
        const auto cfnn = cfn.name(true);
        fileName outFile = cfnn + "TestTime.dat";
        outFile = iter().outputPath() / outFile;
        resOut = outFile;
    }

    testProcessTime myTestPT(iter(), resOut);
#endif // TEST_PROCESS_TIME

    // Second step: Initialize the whole flow field
    initialize();

    marching().timeStep();

    Pout << "    Begin iteration..." << std::endl;

    // Third step: Begin the iteration.
    while (iter_.iterating()) {
#ifdef TEST_PROCESS_TIME
        myTestPT.start(iter_.cStep());
#endif // TEST_PROCESS_TIME

        // Refresh
        iterRefresh();
#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("Refresh");
#endif // TEST_PROCESS_TIME

        previousSource();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("previousSource");
#endif // TEST_PROCESS_TIME

        calculateFc(); // Firstly, compute the convective flux Fc

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("ConvectiveFlux");
#endif // TEST_PROCESS_TIME

        calculateSource(); // Secondly, evaluate the source terms S

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("Source");
#endif // TEST_PROCESS_TIME

        calculateFv(); // Finally, calculate the viscous flux Fv

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("ViscousFlux");
#endif // TEST_PROCESS_TIME

        marching().marching();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("Marching");
#endif // TEST_PROCESS_TIME

        // Update flow field
        updatePrimitives(!timeMarcingPtr_->updatedTemperature());

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("UpdatePrimitives");
#endif // TEST_PROCESS_TIME

        postSource();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("postSource");
#endif // TEST_PROCESS_TIME

        calculateOtherTerms();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("OtherTerms");
#endif // TEST_PROCESS_TIME

        // Write solution
        write();

#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("Write");
#endif // TEST_PROCESS_TIME

        marching().timeStep();
#ifdef TEST_PROCESS_TIME
        myTestPT.clockTimeIncrement("TimeStep");
#endif // TEST_PROCESS_TIME

#ifdef TEST_PROCESS_TIME
        myTestPT.stop();
#endif // TEST_PROCESS_TIME
    }
}

void OpenHurricane::solver::BDFSolve() {
    Pout << "    Begin iteration..." << std::endl;
    marching().marching();
}

void OpenHurricane::solver::clear() noexcept {
    flowPtr_.clear();
    invFluxPtr_.clear();
    timeMarcingPtr_.clear();
    turbPtr_.clear();
    HurDelete(a4Ptr_);
    HurDelete(a5Ptr_);
}

void OpenHurricane::solver::bc() {}

bool OpenHurricane::solver::initPatching() {    

    if (!iter_.cont().found("initialization")) {
        return false;
    }

    const auto &initCont = iter_.cont().subController("initialization");

    if (!initCont.found("initPatch")) {
        return false;
    }
    std::string rnl = initCont.findText("initPatch");
    replaceAllMarks(rnl, "\n", " ");
    if (!rnl.empty()) {
        size_t pos = 0;
        stdStringList rll;
        split(rnl, rll, ",");
        for (integer i = 0; i < rll.size(); ++i) {
            trim(rll[i]);
            if (!initCont.found(rll[i])) {
                checkWarningStr(("Cannot find region: " + rll[i]));
                continue;
            }
            const auto &regionCont = initCont.subController(rll[i]);
            patching initPatch(mesh(), regionCont);
            Pout << "    Info: patching region: " << initPatch.regionName() << std::endl;
            initPatch.patchRegion();
        }
    }
    return true;
}

void OpenHurricane::solver::autoPatching(const integer istep) {
    if (!reqAutoPatch_) {
        return;
    }
    for (integer i = 0; i < autoPatch_.size(); ++i) {
        autoPatch_[i].patchRegion(istep);
    }
}

void OpenHurricane::solver::timeStep(realArray &dt) {}

void OpenHurricane::solver::calculateFc() {
    invFluxPtr_->basicFlux();
}

void OpenHurricane::solver::updatePrimitives(const bool shouldUpdateTemp) {
    if (shouldUpdateTemp) { // Update the temperature
        flowPtr_->thermo().TEa(E(), v(), rho(), T(), temperatureFlag());
    }

    // Update pressure
    flowPtr_->thermo().p(rho(), T(), p());

    // Only for free flame
    //integer nCell = T().mesh().nCells();
    /*real minDx = veryLarge;
        integer minCellI = 0;
    for (integer cellI = 0; cellI < nCell; ++cellI)
    {
        real dx = mag(mesh().cellCentre()[cellI].x() - 0.0095);
        if (minDx > dx)
        {
            minDx = dx;
            minCellI = cellI;
        }
    }*/
    /*real rho0 = 0;
    vector u0 = 0;*/
    /*if (iter().cStep() > 100)
    {
        T()[minCellI] = 700;
        rho0 = rho()[minCellI];
        rho()[minCellI] = mixtures().thermalTable().eos().rhom
        (
            p()[minCellI],
            T()[minCellI],
            mixtures().Yi(),
            minCellI
        );
        u0 = v()[minCellI];
        if (iter().cStep() < 150)
        {
            v()[minCellI].x() = rho0 * u0.x() / rho()[minCellI];
        }
        const auto& fzl = mesh().faceZones();
        for (integer fzi = 0; fzi < fzl.size(); ++fzi)
        {
            if (fzl[fzi].name() == "INLET")
            {
                const auto fi = fzl[fzi].firstIndex();
                const auto cl = mesh().faces()[fi].leftCell();
                const auto cr = mesh().faces()[fi].rightCell();
                real rhoi0 = 0.5 * (rho()[cl] + rho()[cr]);
                real ui0 = rho0 * u0.x() / rhoi0;
                vector vi(ui0, 0, 0);
                v()[cr] = 2.0 * vi - v()[cl];
            }
        }
    }
    else
    {
        const auto& fzl = mesh().faceZones();
        for (integer fzi = 0; fzi < fzl.size(); ++fzi)
        {
            if (fzl[fzi].name() == "INLET")
            {
                const auto fi = fzl[fzi].firstIndex();
                const auto cl = mesh().faces()[fi].leftCell();
                const auto cr = mesh().faces()[fi].rightCell();
                rho0 = 0.5 * (rho()[cl] + rho()[cr]);
                u0 = 0.5 * (v()[cl] + v()[cr]);
            }
        }
    }*/
    /*const auto& fzl = mesh().faceZones();
    for (integer fzi = 0; fzi < fzl.size(); ++fzi)
    {
        if (fzl[fzi].name() == "INLET")
        {
            const auto fi = fzl[fzi].firstIndex();
            const auto cl = mesh().faces()[fi].leftCell();
            const auto cr = mesh().faces()[fi].rightCell();
            rho0 = 0.5 * (rho()[cl] + rho()[cr]);
            u0 = 0.5 * (v()[cl] + v()[cr]);
        }
    }
    for (integer cellI = 0; cellI < nCell; ++cellI)
    {
        rho()[cellI] = mixtures().thermalTable().eos().rhom
        (
            p()[cellI],
            T()[cellI],
            mixtures().Yi(),
            cellI
        );
        E()[cellI] -= 0.5 * v()[cellI].magSqr();
        v()[cellI] = rho0 / rho()[cellI] * u0;
        E()[cellI] += 0.5 * v()[cellI].magSqr();
    }*/

    //autoPatching(iter().cStep());
}

void OpenHurricane::solver::updateFlowOld() {
    flowPtr_->thermo().mu(p(), T(), mu());
}

void OpenHurricane::solver::limits() {
    nLowTemp_ = 0;
    nHighTemp_ = 0;
    nLowPress_ = 0;
    nHighPress_ = 0;

    const integer nCell = mesh().nCells();
    const auto &cells = mesh().cells();
    const auto &faces = mesh().faces();
    const auto &cellCtr = mesh().cellCentre();
    hasTPLmtCells_ = false;
    for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
        if (temperatureFlag()[cellI] == 0) {
            nLowTemp_++;
        } else if (temperatureFlag()[cellI] == 2) {
            nHighTemp_++;
        }
        pressureFlag()[cellI] = 1;

        if (p()[cellI] > flowPtr_->pHigh()) {
            nHighPress_++;
            pressureFlag()[cellI] = 2;
            E()[cellI] -= p()[cellI] / rho()[cellI];
            p()[cellI] = flowPtr_->pHigh();
            //flowPtr_->thermo().correctLimits(p(), cellI, flowPtr_->pLow() + SMALL, flowPtr_->pHigh(), pressureFlag_[cellI]);
            rho()[cellI] = flowPtr_->mixtures().rho(p()[cellI], T()[cellI], cellI);

            E()[cellI] += p()[cellI] / rho()[cellI];
        } else if (p()[cellI] < flowPtr_->pLow()) {
            nLowPress_++;
            pressureFlag()[cellI] = 0;
            E()[cellI] -= p()[cellI] / rho()[cellI];
            p()[cellI] = flowPtr_->pLow() + tiny;
            //flowPtr_->thermo().correctLimits(p(), cellI, flowPtr_->pLow() + SMALL, flowPtr_->pHigh(), pressureFlag_[cellI]);
            rho()[cellI] = flowPtr_->mixtures().rho(p()[cellI], T()[cellI], cellI);
            E()[cellI] += p()[cellI] / rho()[cellI];
        }
    }
    HurMPI::allReduce(nLowTemp_, MPI_SUM);
    HurMPI::allReduce(nHighTemp_, MPI_SUM);
    HurMPI::allReduce(nLowPress_, MPI_SUM);
    HurMPI::allReduce(nHighPress_, MPI_SUM);

    if (nLowTemp_ != 0) {
        hasTPLmtCells_ = true;
        if (HurMPI::master()) {
            std::cerr << "    Temperature limited to " << std::scientific << flowPtr_->TLow()
                      << " in " << nLowTemp_ << " cells" << std::endl;
        }
    }

    if (nHighTemp_ != 0) {
        hasTPLmtCells_ = true;
        if (HurMPI::master()) {
            std::cerr << "    Temperature limited to " << flowPtr_->THigh() << " in " << nHighTemp_
                      << " cells" << std::endl;
        }
    }
    if (nLowPress_ != 0) {
        hasTPLmtCells_ = true;
        if (HurMPI::master()) {
            std::cerr << "    Pressure limited to " << flowPtr_->pLow() << " in " << nLowPress_
                      << " cells" << std::endl;
        }
    }

    if (nHighPress_ != 0) {
        hasTPLmtCells_ = true;
        if (HurMPI::master()) {
            std::cerr << "    Pressure limited to " << flowPtr_->pHigh() << " in " << nHighPress_
                      << " cells" << std::endl;
        }
    }

    HurMPI::barrier(HurMPI::getComm());
}
