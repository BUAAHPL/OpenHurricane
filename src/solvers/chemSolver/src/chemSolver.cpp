/*!
 * \file chemSolver.cpp
 * \brief Main subroutines of solver for chemical ODEs.
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

#include "chemSolver.hpp"

namespace OpenHurricane {
    createClassName(chemSolver);
    registerObjFty(solver, chemSolver, controller);
} // namespace OpenHurricane

void OpenHurricane::chemSolver::calcPhiAndPhiProgress(const real rho, const realArray &yi, real &phi,
                                                  real &phip, real &phil) const {
    const auto &spcs = chemistryPtr_->species();

    realArray ci(spcs.size(), Zero);

    for (integer isp = 0; isp < spcs.size(); ++isp) {
        ci[isp] = rho * yi[isp] / spcs[isp].W();
    }

    real NCo = 0;
    real NHo = 0;
    real NOo = 0;

    real NC = 0;
    real NH = 0;
    real NO = 0;

    real NCL = 0;
    real NHL = 0;
    real NOL = 0;
    for (integer isp = 0; isp < spcs.size(); ++isp) {
        NCo += NCs_[isp] * ci[isp];
        NHo += NHs_[isp] * ci[isp];
        NOo += NOs_[isp] * ci[isp];
        if (spcs[isp].name() != "CO2" && spcs[isp].name() != "H2O") {
            NC += NCs_[isp] * ci[isp];
            NH += NHs_[isp] * ci[isp];
            NO += NOs_[isp] * ci[isp];
        }

        if (spcs[isp].name() == "O2") {
            NCL += NCs_[isp] * ci[isp];
            NHL += NHs_[isp] * ci[isp];
            NOL += NOs_[isp] * ci[isp];
        }
    }
    for (integer i = 0; i < fuelSpecId_.size(); ++i) {
        const integer isp = fuelSpecId_[i];
        NCL += NCs_[isp] * ci[isp];
        NHL += NHs_[isp] * ci[isp];
        NOL += NOs_[isp] * ci[isp];
    }

    phi = (2.0 * NCo + NHo / 2.0 - zp_ * NCo) / (NOo - zp_ * NCo);
    phip = (2.0 * NC + NH / 2.0 - zp_ * NC) / (NO - zp_ * NC);
    //phip = (2.0 * NC + NH / 2.0 - zp_ * NC) / (NOo - zp_ * NCo);
    //phil = (2.0 * NCL + NHL / 2.0) / NOL;
    phil = (2.0 * NCL + NHL / 2.0) / NOo;
}

OpenHurricane::chemSolver::chemSolver(iteration &iter, const runtimeMesh &mesh)
    : solver(iter, mesh), chemistryPtr_(nullptr), yyi_(), p_(Zero), rho_(Zero), T_(Zero),
      timeStep_(Zero), endTime_(Zero), writeMolarFraction_(false), printStep_(10), writeStep_(1),
      fuelSpecId_(), NCs_(), NHs_(), NOs_(), NNs_(), zp_(), writePhiProgress_(false),
      writeHeatReleaseRate_(false), writeSpeciesTimeScale_(false), ignitTemThreshold_(400) {
    const auto &flowCont = iter.cont().subController("flow");
    const auto &mixtureCont = flowCont.subController("mixture");
    const auto &reactionsCont = mixtureCont.subController("reactions");
    if (reactionsCont.found("chemistrySource")) {
        chemistryPtr_ =
            chemistrySource::creator(*flowPtr_, reactionsCont.subController("chemistrySource"));
    } else {
        errorAbortStr(("Chemistry source must be set in " + reactionsCont.name()));
    }
    if (iter.cont().found("chemSolvers")) {
        auto &chemCont = iter.cont().subController("chemSolvers");
        if (chemCont.found("type")) {
            auto w = chemCont.findWord("type");
            trim(w);
            if (w == "constPressure") {
                chemistryPtr_->setConstantPressure();
            } else if (w == "constVolume") {
                chemistryPtr_->setConstantVolume();
            } else {
                errorAbortStr(("Unknown type " + w + "  in " + chemCont.name()));
            }
        }
        auto yi = getBoundariesFromController::getSpeciesMassFractions(chemCont,
                                                                       flowPtr_->mixtures(), false);
        yyi_.resize(yi.size());
        yyi_ = yi;
        p_ = chemCont.findOrDefault<real>("p", 1.01325e5);
        T_ = chemCont.findOrDefault<real>("T", 1200.0);
        rho_ = flowPtr_->mixtures().rho(p_, T_, yyi_);
        timeStep_ = chemCont.findOrDefault<real>("timeStep", 1e-6);
        endTime_ = chemCont.findOrDefault<real>("endTime", 0.0001);
        if (chemistryPtr_) {
            if (chemistryPtr_->isConstPressure()) {
                chemistryPtr_->hea0() = mixtures().thermalTable().ha_p(p_, T_, yyi_);
            } else if (chemistryPtr_->isConstVolume()) {
                chemistryPtr_->hea0() = mixtures().thermalTable().ea_p(p_, T_, yyi_);
            }
        }
        ignitTemThreshold_ = chemCont.findOrDefault<real>("ignitTemThreshold", ignitTemThreshold_);
    } else {
        errorAbortStr(("chemSolvers must be set in " + iter.cont().name()));
    }

    if (iter.cont().found("writeControl")) {
        const auto &writeCont = iter.cont().subController("writeControl");
        if (writeCont.found("chemSolverOut")) {
            const auto &chemSolverOutCont = writeCont.subController("chemSolverOut");
            if (chemSolverOutCont.found("writeSpecies")) {
                auto w = chemSolverOutCont.findWord("writeSpecies");
                trim(w);
                if (w == "massFraction") {
                    writeMolarFraction_ = false;
                } else if (w == "molarFraction") {
                    writeMolarFraction_ = true;
                } else {
                    errorAbortStr(("Unknown type " + w + "  in " + chemSolverOutCont.name()));
                }
            }
        }
        printStep_ = writeCont.findOrDefault<integer>("printStep", printStep_);
        printStep_ = max(integer(1), printStep_);

        writeStep_ = writeCont.findOrDefault<integer>("writeStep", writeStep_);
        writeStep_ = max(integer(1), writeStep_);

        controllerSwitch myCont(writeCont);
        writePhiProgress_ = myCont("writePhiProgress", writePhiProgress_);
        writeHeatReleaseRate_ = myCont("writeHeatReleaseRate", writeHeatReleaseRate_);
        writeSpeciesTimeScale_ = myCont("writeSpeciesTimeScale", writeSpeciesTimeScale_);
    }

    if (writePhiProgress_) {
        if (!chemistryPtr_) {
            LFatal("The chemistry pointer is null");
        }
        const auto &spcs = chemistryPtr_->species();
        NCs_.resize(spcs.size(), Zero);
        NHs_.resize(spcs.size(), Zero);
        NOs_.resize(spcs.size(), Zero);
        NNs_.resize(spcs.size(), Zero);
        for (integer isp = 0; isp < spcs.size(); ++isp) {
            const auto &ele = spcs[isp].elementList();

            for (integer iele = 0; iele < ele.size(); ++iele) {
                if (ele[iele].name() == "C") {
                    NCs_[isp] += ele[iele].nAtoms();
                } else if (ele[iele].name() == "H") {
                    NHs_[isp] += ele[iele].nAtoms();
                } else if (ele[iele].name() == "O") {
                    NOs_[isp] += ele[iele].nAtoms();
                } else if (ele[iele].name() == "N") {
                    NNs_[isp] += ele[iele].nAtoms();
                }
            }
        }

        bool success = false;
        realArray fuelMassFraction;
        if (reactionsCont.found("defineFuelAndOxygen")) {
            const auto &defCont = reactionsCont.subController("defineFuelAndOxygen");

            if (defCont.found("fuel")) {
                success = true;
                auto fuelName_ = defCont.findTextStr("fuel");

                fuelSpecId_.resize(fuelName_.size());
                fuelMassFraction.resize(fuelName_.size(), Zero);
                for (integer isp = 0; isp < spcs.size(); ++isp) {
                    for (integer j = 0; j < fuelName_.size(); ++j) {
                        if (spcs[isp].name() == fuelName_[j]) {
                            fuelSpecId_[j] = isp;
                        }
                    }
                }
                if (defCont.found("fuelMassFraction")) {
                    const auto &fmfcont = defCont.subController("fuelMassFraction");
                    for (integer j = 0; j < fuelName_.size(); ++j) {
                        fuelMassFraction[j] = fmfcont.findOrDefault<real>(fuelName_[j], 0);
                    }
                } else {
                    errorAbortStr(("fuelMassFraction must be defined in " + defCont.name()));
                }
            }
        }

        if (!success) {
            errorAbortStr(("Fule species must be defined in " + reactionsCont.name()));
        }

        real Wf = 0;
        for (integer j = 0; j < fuelSpecId_.size(); ++j) {
            const auto isp = fuelSpecId_[j];
            Wf += fuelMassFraction[j] / spcs[isp].W();
        }
        Wf = inv(Wf);

        real WC = 0;
        real WO = 0;
        for (integer j = 0; j < fuelSpecId_.size(); ++j) {
            const auto isp = fuelSpecId_[j];
            WC += fuelMassFraction[j] * Wf / spcs[isp].W() * NCs_[isp];
            WO += fuelMassFraction[j] * Wf / spcs[isp].W() * NOs_[isp];
        }
        zp_ = WO / WC;
    }
}

void OpenHurricane::chemSolver::solve() {
    fileName outN = iter().outputName();
    auto pathOut = outN.parentPath();
    auto fname = outN.name(true);
    string fext;
    fext += outN.ext();

    outN = fname + fext;
    outN = pathOut / outN;

    fileOsstream fos(outN);
    fos.os() << "variables = "
                "\"time\",\"rho\",\"p\",\"T\",\"nsp\",\"nrc\",\"elapsedTime\"";

#ifdef TEST_PROCESS_TIME
    fileOsstream fosTestChemTime_(IOsstream::streamFormat::ASCII_FORMAT,
                                  IOsstream::openOptions::ONLY_MASTER, std::ios_base::out);
    fileName outFile = iter().configName().name(true) + "" + "TestChemTime" + ".dat";
    outFile = iter().outputPath() / outFile;
    fosTestChemTime_.changeFileName(outFile);

    fosTestChemTime_.open();

    if (!fosTestChemTime_.isOpen()) {
        integer ccount = 0;
        integer kij = 1;
        while (true) {
            fileName outFile1 =
                iter().configName().name(true) + "" + "TestChemTime-" + toString(kij) + ".dat";
            outFile1 = iter().outputPath() / outFile1;
            fosTestChemTime_.changeFileName(outFile1);

            fosTestChemTime_.open();
            if (fosTestChemTime_.isOpen()) {
                PLInfo("    Info: cannot open file: %s.", outFile.name().c_str());
                PLInfo("    Info: open file: %s instead.\n", outFile1.name().c_str());
                break;
            }
            if (ccount++ > 20) {
                LFatal("Cannot open file: %s\n", outFile.name().c_str());
            }
        }
    }

    chemistryPtr_->writeTimeTitleSingle(fosTestChemTime_);

    fos.os() << ",\"foundInTable\"";
#endif

    if (writePhiProgress_) {
        fos.os() << ",\"phi\",\"phiProgress\",\"phil\"";
    }
    if (writeHeatReleaseRate_) {
        fos.os() << ",\"heatReleaseRate[j/m^3-sec]\"";
    }

    for (integer i = 0; i < specTable().size(); ++i) {
        fos.os() << ",\"" << specTable()[i].name().c_str() << "\"";
    }

#ifdef TEST_PROCESS_TIME
    fos.os() << ",\"nSpcGroup\",\"nFastSpc\",\"odeCountIter\"";
#endif

    if (writeSpeciesTimeScale_) {
        for (integer i = 0; i < specTable().size(); ++i) {
            fos.os() << ",\"" << specTable()[i].name().c_str() << "_timescale[s]"
                     << "\"";
        }
    }

    fos.os() << std::endl;
    fos.os() << "zone T = " << iter().configName().name(true) << std::endl;
    fos.setRealPrecision();

    real timedid = Zero;
    Pout.setScientic();
    integer count = 0;
    real ddT = timeStep_;
    real deltaT = timeStep_;
    real T0 = T_;
    real Tig = T0 + ignitTemThreshold_;
    real lastT = T0;
    real lastTime = 0;
    bool ignite = false;
    real tign = 0;
    hrClock myCl;
    bool last = false;
    while (true) {
        if (deltaT + timedid >= endTime_) {
            deltaT = endTime_ - timedid;
            last = true;
        }
        hrClock myCl2;

        ddT = chemistryPtr_->solveTest(deltaT, ddT, p_, T_, rho_, yyi_, fos);
        real elapsedTime = myCl2.elapsedClockTime();
        timedid += deltaT;
        if (!ignite) {
            if (T_ >= Tig) {
                ignite = true;
                if (T_ == Tig) {
                    tign = timedid;
                } else {
                    tign = lastTime + (Tig - lastT) * (timedid - lastTime) / (T_ - lastT);
                }
            }
            lastT = T_;
            lastTime = timedid;
        }
        if (count == 0 || count % (printStep_ * 10) == 0) {
            Pout << std::endl
                 << "step      time did(s)         rho(kg/m^3)         p(Pa)   "
                    "            T(K)                remaining time(s)"
                 << std::endl;
        }

        if (count % printStep_ == 0 || last) {
            std::cout << std::left << std::setfill(' ') << std::setw(10) << count;
            std::cout.setf(std::ios::showpoint);
            std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(8)
                      << timedid;
            std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(8)
                      << rho_;
            std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(8)
                      << p_;
            std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(8)
                      << T_;
            std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(8)
                      << (endTime_ - timedid);

            std::cout.unsetf(std::ios::showpoint);
            std::cout << std::right;
            Pout << std::endl;
        }
        if (count % writeStep_ == 0 || last) {
            fos.os() << timedid << " " << rho_ << " " << p_ << " " << T_ << " "
                     << chemistryPtr_->nSpc() << " " << chemistryPtr_->nRact() << " "
                     << elapsedTime;

#ifdef TEST_PROCESS_TIME
            fos.os() << " " << chemistryPtr_->findInTable();
#endif
            if (writePhiProgress_) {
                real phi, phip, phil;
                calcPhiAndPhiProgress(rho_, yyi_, phi, phip, phil);
                fos.os() << " " << phi;
                fos.os() << " " << phip;
                fos.os() << " " << phil;
            }
            if (writeHeatReleaseRate_) {
                realArray ci(specTable().size());
                for (integer i = 0; i < specTable().size(); ++i) {
                    ci[i] = rho_ * yyi_[i] / specTable()[i].W();
                }
                fos.os() << " " << chemistryPtr_->heatReleaseRate(p_, T_, ci);
            }

            if (writeMolarFraction_) {
                const auto xi = mixtures().species().Yi2Xi(yyi_);
                for (integer i = 0; i < specTable().size(); ++i) {
                    fos.os() << " " << xi[i];
                }
            } else {
                for (integer i = 0; i < specTable().size(); ++i) {
                    fos.os() << " " << yyi_[i];
                }
            }
#ifdef TEST_PROCESS_TIME
            fos.os() << " " << chemistryPtr_->nSpcGroup() << " "
                     << chemistryPtr_->nEachSpcGroup().first() << " "
                     << chemistryPtr_->odeCountIter();

            chemistryPtr_->writeSolveODETimeSingle(fosTestChemTime_, count);

#endif // TEST_PROCESS_TIME
            if (writeSpeciesTimeScale_) {
                realArray ci(specTable().size());
                realArray ts(specTable().size());
                for (integer i = 0; i < specTable().size(); ++i) {
                    ci[i] = rho_ * yyi_[i] / specTable()[i].W();
                }
                chemistryPtr_->timescale2(p_, T_, ci, ts);
                //chemistryPtr_->timescaleReacMode(p_, T_, ci, ts);
                for (integer i = 0; i < specTable().size(); ++i) {
                    fos.os() << " " << ts[i];
                }
            }
            fos.os() << std::endl;
        }
        if (last) {
            break;
        }
        count++;
    }
    const auto calcTime = myCl.elapsedClockTime();
    fos.close();
    Pout.setReal();
    Pout << "    The autoignition delay time: " << tign << " [s]..." << std::endl;
    Pout << "    The calculation elapsed time: " << calcTime << " [s]..." << std::endl;
    Pout.unsetReal();
#ifdef TEST_PROCESS_TIME
    fosTestChemTime_.close();
    chemistryPtr_->printTestTime();
#endif // TEST_PROCESS_TIME
}

void OpenHurricane::chemSolver::solving() {}

void OpenHurricane::chemSolver::BDFSolve() {}

void OpenHurricane::chemSolver::clear() noexcept {
    flowPtr_.clear();
    invFluxPtr_.clear();
    chemistryPtr_.clear();
}

void OpenHurricane::chemSolver::bc() {}

void OpenHurricane::chemSolver::updateProperties() {}

void OpenHurricane::chemSolver::initialize() {}

void OpenHurricane::chemSolver::calculateFc() {}

void OpenHurricane::chemSolver::updatePrimitives(const bool shouldUpdateTemp) {}

void OpenHurricane::chemSolver::write() {}
