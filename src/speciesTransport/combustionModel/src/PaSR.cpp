/*!
 * \file PaSR.cpp
 * \brief Main subroutines for PaSR combustion model.
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

#include "PaSR.hpp"

namespace OpenHurricane {
    createClassNameStr(PaSR, "PaSR");
}

namespace OpenHurricane {
    registerObjFty(combustionModel, PaSR, controller);
}

OpenHurricane::real OpenHurricane::PaSR::TauMixing(const real nu, const real k,
                                                   const real eps) const {
    switch (mixingTimeScale_) {
    case (mixingTimeScale::Kolmogorov):
        return sqrt(nu / max(eps, tiny));
        break;
    case (mixingTimeScale::Integral):
        return cMix_ * k / max(eps, tiny);
        break;
    case (mixingTimeScale::Mean):
        return sqrt(k * sqrt(nu) / max(eps * sqrt(eps), tiny));
        break;
    default:
        break;
    }
    return 0.0;
}
OpenHurricane::PaSR::PaSR(flowModel &flows, const controller &cont, turbulenceModel &turb)
    : finiteRate(flows, cont, turb), mixingTimeScale_(mixingTimeScale::Kolmogorov),
      chemicalTimeScale_(chemicalTimeScale::SpeciesFormationRates), cMix_(1.0),
      dtInit_(flows.mesh().nCells(), large),
      kappaField_(object("PaSRKappa", flows.mesh()), flows.mesh()),
      tchField_(object("PaSRTauch", flows.mesh()), flows.mesh()),
      tmixField_(object("PaSRTaumix", flows.mesh()), flows.mesh()), maxTaustar_(0.10), DRet_(4.0),
      givSpId_(), limitTauStar_(false), tauLimitePtr_(nullptr) {
    if (cont.subController("combustionModel").found("PaSR")) {
        const auto &PaSRCont = cont.subController("combustionModel").subController("PaSR");
        cMix_ = PaSRCont.findOrDefault<real>("cMix", 1.0);
        maxTaustar_ = PaSRCont.findOrDefault<real>("maxTau", 0.10);
        minDtFactor_ = PaSRCont.findOrDefault<real>("minDtFactor", 0.1);
        if (PaSRCont.found("mixingTimeScale")) {
            const auto mw = PaSRCont.findWord("mixingTimeScale");
            if (mw == "Kolmogorov") {
                mixingTimeScale_ = mixingTimeScale::Kolmogorov;
            } else if (mw == "Integral") {
                mixingTimeScale_ = mixingTimeScale::Integral;
            } else if (mw == "Mean") {
                mixingTimeScale_ = mixingTimeScale::Mean;
            } else if (mw == "basedLocalRet") {
                mixingTimeScale_ = mixingTimeScale::basedLocalRet;
                DRet_ = PaSRCont.findOrDefault<real>("DRet", DRet_);
            } else {
                errorAbortStr(("Unknown mixing time scale type: " + mw));
            }
        }
        controllerSwitch myContS(PaSRCont);
        limitTauStar_ = myContS("limitTauStar", limitTauStar_);
        if (limitTauStar_) {
            if (PaSRCont.found("tauLimiteFactor")) {
                real tlf = PaSRCont.findOrDefault<real>("tauLimiteFactor", 0);
                tauLimitePtr_ = new real(tlf);
            }
        }
        if (PaSRCont.found("chemicalTimeScale")) {
            const auto mw = PaSRCont.findWord("chemicalTimeScale");
            if (mw == "ForwardReactionRates") {
                chemicalTimeScale_ = chemicalTimeScale::ForwardReactionRates;
            } else if (mw == "SpeciesFormationRates") {
                chemicalTimeScale_ = chemicalTimeScale::SpeciesFormationRates;
            } else if (mw == "GivenSpeciesProductionRates" ||
                       mw == "GivenSpeciesProductionRatesMax" ||
                       mw == "GivenSpeciesProductionRatesMin") {
                chemicalTimeScale_ = chemicalTimeScale::GivenSpeciesProductionRatesMax;
                if (mw == "GivenSpeciesProductionRatesMin") {
                    chemicalTimeScale_ = chemicalTimeScale::GivenSpeciesProductionRatesMin;
                }
                if (PaSRCont.found("GivenSpeciesProductionRates")) {
                    auto spN = PaSRCont.findTextStr("GivenSpeciesProductionRates");
                    for (integer i = 0; i < spN.size(); ++i) {
                        integer id;
                        if (!species().contains(spN[i], id)) {
                            checkWarningStr(("Cannot found fuel: " + spN[i] + " in species list"));
                        } else {
                            givSpId_.append(id);
                        }
                    }
                } else {
                    integer id;
                    if (species().contains("CH4", id)) {
                        givSpId_.append(id);
                    }
                    if (species().contains("H2", id)) {
                        givSpId_.append(id);
                    }
                    if (species().contains("O2", id)) {
                        givSpId_.append(id);
                    }
                    if (species().contains("H2O", id)) {
                        givSpId_.append(id);
                    }
                    if (species().contains("CO2", id)) {
                        givSpId_.append(id);
                    }
                }
            } else if (mw == "JacobianDiagonal") {
                chemicalTimeScale_ = chemicalTimeScale::JacobianDiagonal;
            } else {
                errorAbortStr(("Unknown chemical time scale type: " + mw));
            }
        }
    }
    kappaField_.setWriteResult();
    tchField_.setWriteResult();
    tmixField_.setWriteResult();
}

void OpenHurricane::PaSR::calcSourceTerms(const realArray &dt) {
    finiteRate::calcSourceTerms(dt);
    //computeFineStructure();
}

void OpenHurricane::PaSR::expChemistrySource(realArray &dt, const bool isModifiedDt) {
    calcSourceTerms(dt);
    getChemistrySource();
    /*const auto& T = flows_.T();
    auto& yi = mixtures().Yi();
    const auto& vol = flows_.mesh().cellVol();

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI)
    {
            for (integer i = 0; i < species().size() - 1; ++i)
            {
                    real omega = chemistryPtr_->Ri(i)[cellI];
                    if (isnan(omega) || isinf(omega))
                    {
                            continue;
                    }
                    yi[i].rhs()[cellI] += kappaField_[cellI] * omega * vol[cellI];
            }
    }*/
}

void OpenHurricane::PaSR::impCalcSourceTerms(const realArray &dt) {
    finiteRate::impCalcSourceTerms(dt);
    //computeFineStructure();
}

void OpenHurricane::PaSR::impChemistrySource(realArray &dt, const bool isModifiedDt) {
    impCalcSourceTerms(dt);
    getImpChemistrySource(dt, isModifiedDt);
    //const auto& T = flows_.T();
    //auto& yi = mixtures().Yi();
    //const auto& vol = flows_.mesh().cellVol();

    //for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI)
    //{
    //	real minDt = veryLarge;
    //	for (integer i = 0; i < species().size() - 1; ++i)
    //	{
    //		real omega = chemistryPtr_->Ri(i)[cellI];
    //		if (isnan(omega) || isinf(omega))
    //		{
    //			yi[i].diagSource()[cellI] = 0;
    //			continue;
    //		}

    //		yi[i].rhs()[cellI] += kappaField_[cellI] * omega * vol[cellI];
    //		if (this->instantaneousReactionRate_)
    //		{
    //			if (isnan(yi[i].diagSource()[cellI]) || isinf(yi[i].diagSource()[cellI]))
    //			{
    //				yi[i].diagSource()[cellI] = 0;
    //				continue;
    //			}
    //			else
    //			{
    //				if (isModifiedDt)
    //				{
    //					minDt = min(minDt, mag(inv(kappaField_[cellI] * yi[i].diagSource()[cellI])));
    //				}
    //			}
    //			yi[i].diagSource()[cellI] *= kappaField_[cellI];
    //			yi[i].diagSource()[cellI] *= vol[cellI];
    //		}
    //		else
    //		{
    //			yi[i].diagSource()[cellI] = 0;
    //		}
    //	}
    //	if (isModifiedDt)
    //	{
    //		if (minDt < dt[cellI])
    //		{
    //			// Reduce the time step and limit to 0.5dt. Do not reduce too much.
    //			dt[cellI] = max(minDt, minDtFactor_ * dt[cellI]);
    //		}
    //	}
    //}
}

void OpenHurricane::PaSR::getChemistrySource() {
    computeFineStructure();
    const auto &T = flows_.T();
    auto &yi = mixtures().Yi();
    const auto &vol = flows_.mesh().cellVol();

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        for (integer i = 0; i < species().size() - 1; ++i) {
            real omega = chemistryPtr_->Ri(i)[cellI];
            if (isnan(omega) || isinf(omega)) {
                continue;
            }
            yi[i].rhs()[cellI] += kappaField_[cellI] * omega * vol[cellI];
        }
    }
}

void OpenHurricane::PaSR::getImpChemistrySource(realArray &dt, const bool isModifiedDt) {
    computeFineStructure();
    const auto &T = flows_.T();
    auto &yi = mixtures().Yi();
    const auto &vol = flows_.mesh().cellVol();

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        real minDt = veryLarge;
        for (integer i = 0; i < species().size() - 1; ++i) {
            real omega = chemistryPtr_->Ri(i)[cellI];
            if (isnan(omega) || isinf(omega)) {
                yi[i].diagSource()[cellI] = 0;
                continue;
            }

            yi[i].rhs()[cellI] += kappaField_[cellI] * omega * vol[cellI];

            if (isnan(yi[i].diagSource()[cellI]) || isinf(yi[i].diagSource()[cellI])) {
                yi[i].diagSource()[cellI] = 0;
                continue;
            } else {
                if (isModifiedDt) {
                    minDt = min(minDt, mag(inv(kappaField_[cellI] * yi[i].diagSource()[cellI])));
                }
            }
            yi[i].diagSource()[cellI] *= kappaField_[cellI];
            yi[i].diagSource()[cellI] *= vol[cellI];
        }
        if (isModifiedDt) {
            if (minDt < dt[cellI]) {
                // Reduce the time step and limit to 0.5dt. Do not reduce too much.
                dt[cellI] = max(minDt, minDtFactor_ * dt[cellI]);
            }
        }
    }
}

void OpenHurricane::PaSR::fullPointImpChemistrySource(realArray &dt, cellRealSquareMatrixArray &Jac,
                                                      const integer rhoId, const integer rhouId,
                                                      const integer rhoEId,
                                                      const integer rhoYi0Id) {
    if (chemicalTimeScale_ == chemicalTimeScale::GivenSpeciesProductionRatesMax ||
        chemicalTimeScale_ == chemicalTimeScale::GivenSpeciesProductionRatesMin ||
        chemicalTimeScale_ == chemicalTimeScale::SpeciesFormationRates) {
        finiteRate::calcSourceTerms(dt);
    }
    computeFineStructure();
    const auto &rho = flows_.rho();
    const auto &T = flows_.T();
    const auto &p = flows_.p();
    const auto &v = flows_.v();
    auto &yi = mixtures().Yi();
    auto &hi = mixtures().hi();
    realArray ei(species().size(), Zero);
    realArray dTdrhoYi(species().size() - 1, Zero);
    const integer nsp = species().size();
    const integer Tid = nsp;
    const auto &vol = flows_.mesh().cellVol();

    realSquareMatrix dwdciT(nsp + 1, Zero);
    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        for (integer i = 0; i < species().size(); ++i) {
            ei[i] = hi[i][cellI] - T[cellI] * species().Ri(i);
        }
        const auto cv = mixtures().thermalTable().cv_p(p[cellI], T[cellI], yi, cellI);
        const real dTdrhoE = real(1) / max(rho[cellI] * cv, veryTiny);
        const real dTdrho = dTdrhoE * (0.5 * v[cellI].magSqr() - ei[nsp - 1]);
        const vector dTdrhov = -dTdrhoE * v[cellI];
        for (integer i = 0; i < nsp - 1; ++i) {
            dTdrhoYi[i] = -dTdrhoE * (ei[i] - ei[nsp - 1]);
        }

        chemistryPtr_->sourceTermsForCellI(cellI, dwdciT);
        for (integer i = 0; i < species().size() - 1; ++i) {
            real omega = chemistryPtr_->Ri(i)[cellI];

            // If omega for i_th species is not a number or infinite,
            // set its Jacobian row to Zero and continue to next species.
            if (isnan(omega) || isinf(omega)) {
                Jac[cellI].row(rhoYi0Id + i).setZero();
                continue;
            }

            yi[i].rhs()[cellI] += kappaField_[cellI] * omega * vol[cellI];

            const auto Wi = species().W(i);
            const auto WNs = species().W(nsp - 1);
            const real dRiddT = Wi * dwdciT(i, Tid);
            const real dRiddrhoY_Ns = Wi / WNs * dwdciT(i, nsp - 1);

            const real factor = kappaField_[cellI] * vol[cellI];

            Jac[cellI](rhoYi0Id + i, rhoId) = (dRiddT * dTdrho + dRiddrhoY_Ns) * factor;
            Jac[cellI](rhoYi0Id + i, rhouId) = dRiddT * dTdrhov[0] * factor;
            Jac[cellI](rhoYi0Id + i, rhouId + 1) = dRiddT * dTdrhov[1] * factor;
            Jac[cellI](rhoYi0Id + i, rhouId + 2) = dRiddT * dTdrhov[2] * factor;
            Jac[cellI](rhoYi0Id + i, rhoEId) = dRiddT * dTdrhoE * factor;
            for (integer j = 0; j < nsp - 1; ++j) {
                Jac[cellI](rhoYi0Id + i, rhoYi0Id + j) =
                    (dRiddT * dTdrhoYi[j] - dRiddrhoY_Ns + Wi / species().W(j) * dwdciT(i, j)) *
                    factor;
            }
        }
    }
}

OpenHurricane::realArray OpenHurricane::PaSR::heatReleaseRate() {
    auto hr = finiteRate::heatReleaseRate();
    hr *= kappaField_;
    return hr;
}

void OpenHurricane::PaSR::computeFineStructure() {
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
    const auto &mut = flows_.mut();

    const auto k = turbulence_.k();
    const auto epsilon = turbulence_.epsilon();
    tchField_ = Zero;
    if (chemicalTimeScale_ == chemicalTimeScale::GivenSpeciesProductionRatesMax) {
        combustionModel::tcGSPR(tchField_, givSpId_, true);
    } else if (chemicalTimeScale_ == chemicalTimeScale::GivenSpeciesProductionRatesMin) {
        combustionModel::tcGSPR(tchField_, givSpId_, false);
    } else if (chemicalTimeScale_ == chemicalTimeScale::SpeciesFormationRates) {
        combustionModel::tcSFR(tchField_);
    } else if (chemicalTimeScale_ == chemicalTimeScale::ForwardReactionRates) {
        combustionModel::tcFRR(tchField_);
    } else if (chemicalTimeScale_ == chemicalTimeScale::JacobianDiagonal) {
        combustionModel::tcJacDT(tchField_);
    }

    //realArray c(nsp, Zero);
    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        // The computation of reaction must be utilized with dimensions.
        //real Ti = T[cellI];
        //real pi = p[cellI];
        const real rhoi = rho[cellI];
        const real mui = mu[cellI];
        const real nui = mui / rhoi;
        const real ki = k[cellI];
        const real ei = epsilon[cellI];

        real tmix = 0;
        switch (mixingTimeScale_) {
        case (mixingTimeScale::Kolmogorov):
            tmix = sqrt(nui / max(ei, tiny));
            break;
        case (mixingTimeScale::Integral):
            tmix = cMix_ * ki / max(ei, tiny);
            break;
        case (mixingTimeScale::basedLocalRet):
            tmix = pow(real(0.09) * mui / mut[cellI], alphaRetExp()) * ki / max(ei, tiny);
            break;
        case (mixingTimeScale::Mean):
        default:
            tmix = sqrt(ki * sqrt(nui) / max(ei * sqrt(ei), tiny));
            break;
        }

        real kappa = 1.0;
        if (!isinf(tchField_[cellI]) && !isnan(tchField_[cellI])) {
            if (tmix > tiny) {
                kappa = tchField_[cellI] / (tchField_[cellI] + tmix);
            }
        }
        tmixField_[cellI] = tmix;
        kappaField_[cellI] = kappa;
    }
}

void OpenHurricane::PaSR::constVolReactor(const realArray &dt, const real dtFactor,
                                          const bool firstCall) {
    computeFineStructure();
    auto &spt = species();
    const integer nsp = spt.size();         // The size of species
    const integer nrc = reactions().size(); // The size of reactions
    auto &T = flows_.T();
    auto &p = flows_.p();
    auto &E = flows_.E();
    auto &v = flows_.v();
    const auto &rho = flows_.rho();

    auto &yi = mixtures().Yi();

    //chemistryPtr_->setConstantVolume();
    chemistryPtr_->setConstantPressure();
    realArray y0(nsp);

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        // The computation of reaction must be utilized with dimensions.
        real Ti = T[cellI];
        real pi = p[cellI];
        real rhoi = rho[cellI];
        integer n = maxSubStep_;
        //chemistryPtr_->setODEFactor(kappaField_[cellI]);

        //const real e0 = (E[cellI] - 0.5 * v[cellI].magSqr())*E.refValue();
        for (integer i = 0; i < nsp; ++i) {
            y0[i] = yi[i][cellI];
        }

        for (integer ii = 0; ii < n; ++ii) {
            chemistryPtr_->solve(dt[cellI] * dtFactor / real(n), dt[cellI] * dtFactor / real(n), pi,
                                 Ti, rhoi, y0);
            pi = chemistryPtr_->p();
            Ti = chemistryPtr_->T();
        }
        for (integer i = 0; i < nsp; ++i) {
            yi[i][cellI] = y0[i];
        }
        p[cellI] = chemistryPtr_->p();
        T[cellI] = chemistryPtr_->T();
    }
    //chemistryPtr_->unsetODEFactor();
}

OpenHurricane::realArrayArray OpenHurricane::PaSR::omegai() const {
    realArrayArray ome;
    ome.resize(species().size());
    for (integer i = 0; i < species().size(); ++i) {
        ome[i].resize(flows_.mesh().nCells(), Zero);
    }
    for (integer cellI = 0; cellI < flows_.mesh().nCells(); ++cellI) {
        for (integer i = 0; i < species().size(); ++i) {
            real omega = chemistryPtr_->Ri(i)[cellI];

            ome[i][cellI] = kappaField_[cellI] * omega;
        }
    }
    return ome;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::PaSR::tcSFR() {
    //chemistryPtr_->calculateSourceTerms();
    return combustionModel::tcSFR();
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::PaSR::tcGSPR() {
    //chemistryPtr_->calculateSourceTerms();
    return combustionModel::tcGSPR();
}

void OpenHurricane::PaSR::expChemistrySourceCoupled() {
    chemistryPtr_->calculateSourceTerms();
    computeFineStructure();
    const auto &T = flows_.T();
    auto &yi = mixtures().Yi();
    const auto &vol = flows_.mesh().cellVol();

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        for (integer i = 0; i < species().size() - 1; ++i) {
            real omega = chemistryPtr_->Ri(i)[cellI];
            if (isnan(omega) || isinf(omega)) {
                continue;
            }
            yi[i].rhs()[cellI] += kappaField_[cellI] * omega * vol[cellI];
        }
    }
}

void OpenHurricane::PaSR::diagImpChemistrySourceCoupled() {
    chemistryPtr_->calculateSourceTermsImp();

    computeFineStructure();
    const auto &T = flows_.T();
    auto &yi = mixtures().Yi();
    const auto &vol = flows_.mesh().cellVol();

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        for (integer i = 0; i < species().size() - 1; ++i) {
            real omega = chemistryPtr_->Ri(i)[cellI];
            if (isnan(omega) || isinf(omega)) {
                yi[i].diagSource()[cellI] = 0;
                continue;
            }

            yi[i].rhs()[cellI] += kappaField_[cellI] * omega * vol[cellI];
            if (isnan(yi[i].diagSource()[cellI]) || isinf(yi[i].diagSource()[cellI])) {
                yi[i].diagSource()[cellI] = 0;
                continue;
            }
            yi[i].diagSource()[cellI] *= kappaField_[cellI];
            yi[i].diagSource()[cellI] *= vol[cellI];
        }
    }
}

void OpenHurricane::PaSR::fullImpChemistrySourceCoupled(cellRealSquareMatrixArray &Jac,
                                                        const integer rhoId, const integer rhouId,
                                                        const integer rhoEId,
                                                        const integer rhoYi0Id) {
    if (chemicalTimeScale_ == chemicalTimeScale::GivenSpeciesProductionRatesMax ||
        chemicalTimeScale_ == chemicalTimeScale::GivenSpeciesProductionRatesMin ||
        chemicalTimeScale_ == chemicalTimeScale::SpeciesFormationRates) {
        chemistryPtr_->calculateSourceTerms();
    }
    computeFineStructure();

    const auto &rho = flows_.rho();
    const auto &T = flows_.T();
    const auto &p = flows_.p();
    const auto &v = flows_.v();
    auto &yi = mixtures().Yi();
    auto &hi = mixtures().hi();
    realArray ei(species().size(), Zero);
    realArray dTdrhoYi(species().size() - 1, Zero);
    const integer nsp = species().size();
    const integer Tid = nsp;
    const auto &vol = flows_.mesh().cellVol();

    realSquareMatrix dwdciT(nsp + 1, Zero);

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        for (integer i = 0; i < species().size(); ++i) {
            ei[i] = hi[i][cellI] - T[cellI] * species().Ri(i);
        }
        const auto cv = mixtures().thermalTable().cv_p(p[cellI], T[cellI], yi, cellI);
        const real dTdrhoE = real(1) / max(rho[cellI] * cv, veryTiny);
        const real dTdrho = dTdrhoE * (0.5 * v[cellI].magSqr() - ei[nsp - 1]);
        const vector dTdrhov = -dTdrhoE * v[cellI];
        for (integer i = 0; i < nsp - 1; ++i) {
            dTdrhoYi[i] = -dTdrhoE * (ei[i] - ei[nsp - 1]);
        }

        chemistryPtr_->sourceTermsForCellI(cellI, dwdciT);
        for (integer i = 0; i < species().size() - 1; ++i) {
            real omega = chemistryPtr_->Ri(i)[cellI];

            // If omega for i_th species is not a number or infinite,
            // set its Jacobian row to Zero and continue to next species.
            if (isnan(omega) || isinf(omega)) {
                Jac[cellI].row(rhoYi0Id + i).setZero();
                continue;
            }

            yi[i].rhs()[cellI] += kappaField_[cellI] * omega * vol[cellI];

            const auto Wi = species().W(i);
            const auto WNs = species().W(nsp - 1);
            const real dRiddT = Wi * dwdciT(i, Tid);
            const real dRiddrhoY_Ns = Wi / WNs * dwdciT(i, nsp - 1);

            const real factor = kappaField_[cellI] * vol[cellI];

            Jac[cellI](rhoYi0Id + i, rhoId) = (dRiddT * dTdrho + dRiddrhoY_Ns) * factor;
            Jac[cellI](rhoYi0Id + i, rhouId) = dRiddT * dTdrhov[0] * factor;
            Jac[cellI](rhoYi0Id + i, rhouId + 1) = dRiddT * dTdrhov[1] * factor;
            Jac[cellI](rhoYi0Id + i, rhouId + 2) = dRiddT * dTdrhov[2] * factor;
            Jac[cellI](rhoYi0Id + i, rhoEId) = dRiddT * dTdrhoE * factor;
            for (integer j = 0; j < nsp - 1; ++j) {
                Jac[cellI](rhoYi0Id + i, rhoYi0Id + j) =
                    (dRiddT * dTdrhoYi[j] - dRiddrhoY_Ns + Wi / species().W(j) * dwdciT(i, j)) *
                    factor;
            }
        }
    }
}

void OpenHurricane::PaSR::chemistrySourceIntegrated() {
    if (chemicalTimeScale_ == chemicalTimeScale::GivenSpeciesProductionRatesMax ||
        chemicalTimeScale_ == chemicalTimeScale::GivenSpeciesProductionRatesMin ||
        chemicalTimeScale_ == chemicalTimeScale::SpeciesFormationRates) {
        chemistryPtr_->calculateSourceTerms();
    }
    computeFineStructure();

    const auto &T = flows_.T();
    const auto &iter = T.mesh().Iteration();
    auto &tau = tauInt();
    if (limitTauStar_ && tauLimitePtr_ != nullptr) {
        for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
            tau[cellI] = min(min(tchField_[cellI], tmixField_[cellI]), *tauLimitePtr_);
        }
    } else if (limitTauStar_ && iter.hasPhysicalTimeStep()) {
        for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
            tau[cellI] = min(min(tchField_[cellI], tmixField_[cellI]), iter.pTStep().pTimeStep());
        }
    } else if (limitTauStar_ && T.mesh().foundOnlyObject("timeStep")) {
        const auto &dt = T.mesh().findObjectRef<cellRealArray>("timeStep");
        for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
            tau[cellI] = min(min(tchField_[cellI], tmixField_[cellI]), dt[cellI]);
        }
    } else {
        for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
            tau[cellI] = min(tchField_[cellI], tmixField_[cellI]);
        }
    }

    chemistryPtr_->solve(tau, 1.0);
}

void OpenHurricane::PaSR::getIntegratedChemSource() {
    const auto &T = flows_.T();
    auto &yi = mixtures().Yi();
    const auto &vol = flows_.mesh().cellVol();

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        for (integer i = 0; i < species().size() - 1; ++i) {
            real omega = chemistryPtr_->Ri(i)[cellI];
            if (isnan(omega) || isinf(omega)) {
                yi[i].diagSource()[cellI] = 0;
                continue;
            }

            yi[i].rhs()[cellI] += kappaField_[cellI] * omega * vol[cellI];
        }
    }
}

void OpenHurricane::PaSR::chemistrySourceStrangSplitted(const realArray &dt, const real dtFactor) {
    LFatal("Cannot used PaSR model in StrangSplitted scheme");
}
