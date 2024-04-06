/*!
 * \file combustionModel.cpp
 * \brief Main subroutines for combustion model.
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

#include "combustionModel.hpp"
#include "Lapacians.hpp"
#include "makeChemistryODESolvers.hpp"
#include "spatialScheme.hpp"
#include "timeMarching.hpp"

std::map<std::string, OpenHurricane::combustionModel::chemistryOption>
    OpenHurricane::combustionModel::validOptions;

namespace OpenHurricane {
    createClassNameStr(combustionModel, "combustionModel");
    createObjFty(combustionModel, controller);
} // namespace OpenHurricane

OpenHurricane::combustionModel::initValidTables::initValidTables() {
    combustionModel::addValidOptions("COUPLED",
                                     OpenHurricane::combustionModel::chemistryOption::coupled);
    combustionModel::addValidOptions(
        "STRANGSPLITTED", OpenHurricane::combustionModel::chemistryOption::strangSplitted);
    combustionModel::addValidOptions("INTEGRATED",
                                     OpenHurricane::combustionModel::chemistryOption::integrated);
}

OpenHurricane::combustionModel::initValidTables dummyChemInitValidTables;

void OpenHurricane::combustionModel::addValidOptions(const std::string &opt,
                                                     const chemistryOption param) {
    validOptions.emplace(opt, param);
}

OpenHurricane::combustionModel::combustionModel(flowModel &flows, const controller &cont,
                                                turbulenceModel &turb)
    : flows_(flows), chemistryPtr_(nullptr), turbulence_(turb),
      maxSubStep_(cont.findOrDefault<integer>("maxSubStep", 3, true)), fuelName_(), oxygenName_(),
      productionName_(), minDtFactor_(0.1), fuelOxygenProduIdPtr_(nullptr),
      chemOptions_(chemistryOption::coupled), tauIntPtr_(nullptr), mixZPtr_(nullptr),
      fuelMassFraction_(), oxygenMassFraction_(), fuelId_(), oxygenId_(), stoichiometricRatio_(0) {
    if (cont.found("chemistrySource")) {
        if (cont.isController("chemistrySource")) {
            chemistryPtr_ = chemistrySource::creator(flows, cont.subController("chemistrySource"));
        } else {
            string ty = cont.findWord("chemistrySource");
            if (ty != OpenHurricane::chemistryNoSolver::className_) {
                errorAbortStr(("Wrong setting chemistry: " + ty));
            }
            chemistryPtr_.reset(new chemistryNoSolver(flows, cont));
        }
    }
    if (cont.found("defineFuelAndOxygen")) {
        const auto &defCont = cont.subController("defineFuelAndOxygen");

        if (defCont.found("fuel")) {
            fuelName_ = defCont.findTextStr("fuel");
            fuelId_.resize(fuelName_.size(), -1);
            for (integer isp = 0; isp < species().size(); ++isp) {
                for (integer j = 0; j < fuelName_.size(); ++j) {
                    if (fuelName_[j] == species()[isp].name()) {
                        fuelId_[j] = isp;
                    }
                }
            }

            if (defCont.found("fuelMassFraction")) {
                fuelMassFraction_.resize(fuelName_.size(), -1);
                const auto &fmfcont = defCont.subController("fuelMassFraction");
                for (integer j = 0; j < fuelName_.size(); ++j) {
                    fuelMassFraction_[j] = fmfcont.findOrDefault<real>(fuelName_[j], 0);
                }
            }
        }

        if (defCont.found("oxygen")) {
            oxygenName_ = defCont.findTextStr("oxygen");
            oxygenId_.resize(oxygenName_.size(), -1);
            for (integer isp = 0; isp < species().size(); ++isp) {
                for (integer j = 0; j < oxygenName_.size(); ++j) {
                    if (oxygenName_[j] == species()[isp].name()) {
                        oxygenId_[j] = isp;
                    }
                }
            }
            if (defCont.found("oxygenMassFraction")) {
                oxygenMassFraction_.resize(oxygenName_.size(), -1);
                const auto &fmfcont = defCont.subController("oxygenMassFraction");
                for (integer j = 0; j < oxygenName_.size(); ++j) {
                    oxygenMassFraction_[j] = fmfcont.findOrDefault<real>(oxygenName_[j], 0);
                    ;
                }
            }
        }
        stoichiometricRatio_ =
            defCont.findOrDefault<real>("stoichiometricRatio", stoichiometricRatio_);
    }
    if (cont.found("defineProduction")) {
        productionName_ = cont.findTextStr("defineProduction");
    }

    if (cont.found("chemistryOption")) {
        auto optionNameW = cont.findWord("chemistryOption");
        auto optionNameStr = cont.findWord("chemistryOption");
        stringToUpperCase(optionNameStr);

        auto iter = validOptions.find(optionNameStr);

        if (iter != validOptions.end()) {
            chemOptions_ = iter->second;
        } else {
            checkWarningStr(("Unknown chemistry option: " + optionNameW + " in " + cont.name() +
                             ". And the coupled option will be used."));
        }
    }

    if (isCoupled()) {
        Pout << "    Info: setting chemistry option as \"coupled\"" << std::endl;
    } else if (isIntegrated()) {
        Pout << "    Info: setting chemistry option as \"integrated\"" << std::endl;
    } else if (isStrangSplitted()) {
        Pout << "    Info: setting chemistry option as \"StrangSplitted\"" << std::endl;
    }
}

OpenHurricane::uniquePtr<OpenHurricane::combustionModel>
OpenHurricane::combustionModel::creator(flowModel &flows, const controller &cont,
                                        turbulenceModel &turb) {
    string model = cont.subController("combustionModel").findWord("type");

    Pout << "    Setting species combustion model: " << model << std::endl;
    defineInObjCreator(combustionModel, static_cast<std::string>(model), controller,
                       (flows, cont, turb));
}

OpenHurricane::realArray OpenHurricane::combustionModel::calcDamkohler() {
    // Species table
    const auto &spt = species();
    // The number of species
    const integer nsp = spt.size();
    // The number of reactions
    const integer nrc = reactions().size();
    // The temperature array
    const auto &T = flows_.T();
    // The pressure array
    const auto &p = flows_.p();
    // The density array
    const auto &rho = flows_.rho();

    const auto &vol = flows_.mesh().cellVol();
    auto &yi = mixtures().Yi();
    const auto &mu = flows_.mul();
    const auto k = turbulence_.k();
    const auto epsilon = turbulence_.epsilon();
    realArray Da(T.mesh().nTotalCells(), Zero);
    realArray c(nsp, Zero);
    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        // The computation of reaction must be utilized with dimensions.
        real Ti = T[cellI];
        real pi = p[cellI];
        real rhoi = rho[cellI];
        const real mui = mu[cellI];
        const real nui = mui / rhoi;
        const real ki = k[cellI];
        const real ei = epsilon[cellI];

        reactions().updateG0(Ti);

        // Transfer the mass fraction: yi to molar concentration: c
        for (integer i = 0; i < nsp; ++i) {
            c[i] = rhoi * yi[i][cellI] / spt.W(i);
        }

        // integral scale
        real tmix = ki / max(ei, tiny);

        real tch = chemistryPtr_->tc(pi, Ti, c);
        //real kappa = 1.0;
        if (!isinf(tch) && !isnan(tch)) {
            /*if (tmix > tiny)
            {
                    kappa = tch / (tch + tmix);
            }*/
            Da[cellI] = tmix / tch;
        } else {
            Da[cellI] = Zero;
        }
    }
    return Da;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::combustionModel::flameIndex() const {
    realArray f(flows_.mesh().nTotalCells(), Zero);
    const auto ome = omegai();
    real maxW = Zero;
    for (integer n = 0; n < flows_.mesh().nCells(); ++n) {
        f[n] = Zero;

        for (integer isp = 0; isp < species().size(); ++isp) {
            f[n] = max(f[n], mag(ome[isp][n]) / species().W(isp));
        }

        maxW = max(maxW, f[n]);
    }
    HurMPI::allReduce(maxW, MPI_MAX);
    f /= max(maxW, veryTiny);
    return f;
}

hur_nodiscard void OpenHurricane::combustionModel::tcFRR(realArray &tc) {
    if (tc.size() < flows_.mesh().nCells()) {
        tc.resize(flows_.mesh().nTotalCells(), Zero);
    }
    // Species table
    const auto &spt = species();
    // The number of species
    const integer nsp = spt.size();

    // The temperature array
    const auto &T = flows_.T();
    // The pressure array
    const auto &p = flows_.p();
    // The density array
    const auto &rho = flows_.rho();
    auto &yi = mixtures().Yi();

    realArray c(nsp, Zero);
    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        // The computation of reaction must be utilized with dimensions.
        real Ti = T[cellI];
        real pi = p[cellI];
        real rhoi = rho[cellI];
        reactions().updateG0(Ti);

        // Transfer the mass fraction: yi to molar concentration: c
        for (integer i = 0; i < nsp; ++i) {
            c[i] = rhoi * yi[i][cellI] / spt.W(i);
        }

        tc[cellI] = chemistryPtr_->tc(pi, Ti, c);
    }
}

hur_nodiscard void OpenHurricane::combustionModel::tcSFR(realArray &tc) const {
    if (tc.size() < flows_.mesh().nCells()) {
        tc.resize(flows_.mesh().nTotalCells(), Zero);
    }

    // Species table
    const auto &spt = species();
    // The number of species
    const integer nsp = spt.size();

    // The temperature array
    const auto &T = flows_.T();
    // The pressure array
    const auto &p = flows_.p();
    // The density array
    const auto &rho = flows_.rho();
    auto &yi = mixtures().Yi();

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        real Ti = T[cellI];
        real pi = p[cellI];
        real rhoi = rho[cellI];

        tc[cellI] = 0;
        bool hasGot = false;
        for (integer isp = 0; isp < species().size(); ++isp) {
            if (yi[isp][cellI] > tiny) {
                real tcc = rho[cellI] * yi[isp][cellI] /
                           (max(mag(chemistryPtr_->Ri()[isp][cellI]), veryTiny));
                if (tcc > 1e5) {
                    continue;
                }
                if (!isnan(tcc) && !isinf(tcc)) {
                    hasGot = true;
                    tcc = min(tcc, real(200));
                    tc[cellI] = max(tc[cellI], tcc);
                }
            }
        }
        if (!hasGot) {
            tc[cellI] = 200;
        }
    }
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::combustionModel::tcGSPR(const integer isp) {
    realArray tc(flows_.mesh().nTotalCells(), Zero);

    const auto &rho = flows_.rho();
    const auto &yi = mixtures().Yi();

    for (integer n = 0; n < flows_.mesh().nCells(); ++n) {
        tc[n] = rho[n] * yi[isp][n] / (max(mag(chemistryPtr_->Ri()[isp][n]), veryTiny));
        tc[n] = min(tc[n], real(200));
    }
    return tc;
}

void OpenHurricane::combustionModel::tcGSPR(realArray &tc) const {
    if (tc.size() < flows_.mesh().nCells()) {
        tc.resize(flows_.mesh().nTotalCells(), Zero);
    }

    const auto &spId = fuelOxygenProduId();

    const auto &rho = flows_.rho();
    const auto &yi = mixtures().Yi();
    if (spId.size() != 0) {
        for (integer n = 0; n < flows_.mesh().nCells(); ++n) {
            tc[n] = 0;
            bool hasGot = false;
            for (integer i = 0; i < spId.size(); ++i) {
                const auto isp = spId[i];
                if (yi[isp][n] > tiny) {
                    hasGot = true;
                    real tcc =
                        rho[n] * yi[isp][n] / (max(mag(chemistryPtr_->Ri()[isp][n]), veryTiny));

                    tcc = min(tcc, real(200));
                    if (!isnan(tcc) && !isinf(tcc)) {
                        tc[n] = max(tc[n], tcc);
                    } else {
                        tc[n] = 200;
                    }
                }
            }
            if (!hasGot) {
                tc[n] = 200;
            }
        }
    }
}

void OpenHurricane::combustionModel::tcGSPR(realArray &tc, const integerList &spId,
                                            const bool maxOrMin) const {
    if (tc.size() < flows_.mesh().nCells()) {
        tc.resize(flows_.mesh().nTotalCells(), Zero);
    }

    const auto &rho = flows_.rho();
    const auto &yi = mixtures().Yi();
    if (spId.size() != 0) {
        if (maxOrMin) {
            for (integer n = 0; n < flows_.mesh().nCells(); ++n) {
                tc[n] = 0;
                bool hasGot = false;
                for (integer i = 0; i < spId.size(); ++i) {
                    const auto isp = spId[i];
                    if (yi[isp][n] > tiny) {
                        hasGot = true;
                        real tcc =
                            rho[n] * yi[isp][n] / (max(mag(chemistryPtr_->Ri()[isp][n]), veryTiny));

                        tcc = min(tcc, real(200));
                        if (!isnan(tcc) && !isinf(tcc)) {
                            tc[n] = max(tc[n], tcc);
                        } else {
                            tc[n] = 200;
                        }
                    }
                }
                if (!hasGot) {
                    tc[n] = 200;
                }
            }
        } else {
            for (integer n = 0; n < flows_.mesh().nCells(); ++n) {
                tc[n] = 200;
                bool hasGot = false;
                for (integer i = 0; i < spId.size(); ++i) {
                    const auto isp = spId[i];
                    if (yi[isp][n] > tiny) {
                        hasGot = true;
                        real tcc =
                            rho[n] * yi[isp][n] / (max(mag(chemistryPtr_->Ri()[isp][n]), veryTiny));

                        tcc = min(tcc, real(200));
                        if (!isnan(tcc) && !isinf(tcc)) {
                            tc[n] = min(tc[n], tcc);
                        }
                    }
                }
                if (!hasGot) {
                    tc[n] = 200;
                }
            }
        }
    }
}

void OpenHurricane::combustionModel::tcJacDT(realArray &tc) {
    if (tc.size() < flows_.mesh().nCells()) {
        tc.resize(flows_.mesh().nTotalCells(), Zero);
    }
    // Species table
    const auto &spt = species();
    // The number of species
    const integer nsp = spt.size();
    // The number of reactions
    const integer nrc = reactions().size();
    // The temperature array
    const auto &T = flows_.T();
    // The pressure array
    const auto &p = flows_.p();
    // The density array
    const auto &rho = flows_.rho();
    auto &yi = mixtures().Yi();

    realArray c(nsp, Zero);
    realArray yyi(nsp, Zero);
    realArray tci(nsp, Zero);
    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        real Ti = T[cellI];
        real pi = p[cellI];
        real rhoi = rho[cellI];
        reactions().updateG0(Ti);

        // Transfer the mass fraction: yi to molar concentration: c
        for (integer i = 0; i < nsp; ++i) {
            yyi[i] = yi[i][cellI];
            c[i] = rhoi * yi[i][cellI] / spt.W(i);
        }

        chemistryPtr_->timescale2(rhoi, pi, Ti, c, yyi, tci);

        tc[cellI] = 0;
        integer count = 0;
        for (integer i = 0; i < nsp; ++i) {
            if (tci[i] <= 200) {
                tc[cellI] = max(tc[cellI], tci[i]);
                tc[cellI] = min(tc[cellI], real(200));
                count++;
            }
        }
        if (count == 0) {
            tc[cellI] = 200;
        }
    }
}

hur_nodiscard OpenHurricane::integerList OpenHurricane::combustionModel::fuelOxygenProduId() const {
    if (fuelOxygenProduIdPtr_ == nullptr) {
        integerList spId;
        for (integer i = 0; i < fuelName_.size(); ++i) {
            integer id;
            if (!species().contains(fuelName_[i], id)) {
                checkWarningStr(("Cannot found fuel: " + fuelName_[i] + " in species list"));
            } else {
                spId.append(id);
            }
        }

        for (integer i = 0; i < oxygenName_.size(); ++i) {
            integer id;
            if (!species().contains(oxygenName_[i], id)) {
                checkWarningStr(("Cannot found oxygen: " + oxygenName_[i] + " in species list"));
            } else {
                spId.append(id);
            }
        }

        for (integer i = 0; i < productionName_.size(); ++i) {
            integer id;
            if (!species().contains(productionName_[i], id)) {
                checkWarningStr(
                    ("Cannot found production: " + productionName_[i] + " in species list"));
            } else {
                spId.append(id);
            }
        }
        fuelOxygenProduIdPtr_ = new integerList(spId.size());
        *fuelOxygenProduIdPtr_ = spId;
    }
    return *fuelOxygenProduIdPtr_;
}

void OpenHurricane::combustionModel::evaluateSource(timeMarching &times, const integer rhoId,
                                                    const integer rhouId, const integer rhoEId,
                                                    const integer rhoYi0Id) {
    if (chemOptions_ == chemistryOption::coupled) {
        chemistrySourceCoupled(times, rhoId, rhouId, rhoEId, rhoYi0Id);
    } else if (chemOptions_ == chemistryOption::integrated) {
        getIntegratedChemSource();
    }
}

void OpenHurricane::combustionModel::integratedChemistrySource() {
    if (chemOptions_ == chemistryOption::integrated) {
        chemistrySourceIntegrated();
    }
}

void OpenHurricane::combustionModel::strangSplittedChemistrySource(timeMarching &times) {
    if (isStrangSplitted()) {
        const auto &T = flows_.T();
        if (flows_.mesh().Iteration().hasPhysicalTimeStep()) {
            auto &tau = tauInt();
            for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
                tau[cellI] = flows_.mesh().Iteration().pTStep().pTimeStep();
            }
            chemistrySourceStrangSplitted(tau, real(0.5));
        } else {
            chemistrySourceStrangSplitted(times.dt(), real(0.5));
        }
    }
}

void OpenHurricane::combustionModel::chemistrySourceCoupled(timeMarching &times,
                                                            const integer rhoId,
                                                            const integer rhouId,
                                                            const integer rhoEId,
                                                            const integer rhoYi0Id) {
    if (chemOptions_ != chemistryOption::coupled) {
        return;
    }

    if (times.explicitSource()) {
        expChemistrySourceCoupled();
    } else if (times.diagonalImpSource()) {
        diagImpChemistrySourceCoupled();
    } else if (times.fullJacobianSource() || times.fullJacobianSourceTable()) {
        fullImpChemistrySourceCoupled(times.Jacobian(), rhoId, rhouId, rhoEId, rhoYi0Id);
    }
}

void OpenHurricane::combustionModel::getIntegratedChemSource() {
    const auto &T = flows_.T();
    auto &yi = mixtures().Yi();
    const auto &vol = flows_.mesh().cellVol();

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        for (integer i = 0; i < species().size() - 1; ++i) {
            real omega = chemistryPtr_->Ri(i)[cellI];
            if (isnan(omega) || isinf(omega)) {
                continue;
            }
            yi[i].rhs()[cellI] += omega * vol[cellI];
        }
    }
}

bool OpenHurricane::combustionModel::calcMixtureFraction() {
    if (fuelMassFraction_.size() == 0) {
        checkWarning("The mass fraction of fuel is not given for computing mixture fraction");
        return false;
    }
    if (oxygenMassFraction_.size() == 0) {
        checkWarning("The mass fraction of oxygen is not given for computing mixture fraction");
        return false;
    }

    auto &Z = mixZ();

    real YF1 = 0;
    for (integer j = 0; j < fuelId_.size(); ++j) {
        if (fuelId_[j] != -1) {
            YF1 += fuelMassFraction_[j];
        }
    }

    real YO22 = 0;
    for (integer j = 0; j < oxygenId_.size(); ++j) {
        if (oxygenId_[j] != -1) {
            YO22 += oxygenMassFraction_[j];
        }
    }
    const real s = stoichiometricRatio_;

    for (integer n = 0; n < Z.size(); ++n) {
        real YF = 0;
        for (integer j = 0; j < fuelId_.size(); ++j) {
            if (fuelId_[j] != -1) {
                integer i = fuelId_[j];
                YF += mixtures().Yi()[i][n];
            }
        }

        real YO2 = 0;
        for (integer j = 0; j < oxygenId_.size(); ++j) {
            if (oxygenId_[j] != -1) {
                integer i = oxygenId_[j];
                YO2 += mixtures().Yi()[i][n];
            }
        }

        Z[n] = max(real(0), (s * YF - YO2 + YO22) / (s * YF1 + YO22));
    }

    Pout.unsetReal();
    return true;
}

OpenHurricane::realArray
OpenHurricane::combustionModel::mixtureDisspationRate(const spatialScheme &sps) {
    realArray xi(flows_.mesh().nTotalCells(), Zero);

    if (calcMixtureFraction()) {
        auto &Z = mixZ();
        sps.grad(Z);

        for (integer n = 0; n < flows_.mesh().nCells(); ++n) {
            /*real Dm = 0;
            for (integer i = 0; i < species().size(); ++i)
            {
                    Dm += mixtures().Yi()[i][n] * mixtures().Dim(i)[n];
            }*/
            real Dm = mixtures().transTable().DiffMix(flows_.p()[n], flows_.T()[n], mixtures().Yi(),
                                                      mixtures().Dim(), n);

            xi[n] = (2.0 * Dm) * Z.grad()[n].magSqr();
        }
    }

    return xi;
}
