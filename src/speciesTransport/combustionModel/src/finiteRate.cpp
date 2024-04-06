/*!
 * \file finiteRate.cpp
 * \brief Main subroutines for finite-rate combustion model.
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

#include "finiteRate.hpp"

namespace OpenHurricane {
    createClassNameStr(finiteRate, "finiteRate");
}

namespace OpenHurricane {
    registerObjFty(combustionModel, finiteRate, controller);
}
OpenHurricane::finiteRate::finiteRate(flowModel &flows, const controller &cont,
                                      turbulenceModel &turb)
    : combustionModel(flows, cont, turb), tchField_(object("tche", flows.mesh()), flows.mesh()),
      tmpNField_(object("tmpN", flows.mesh()), flows.mesh()),
      tmpYoField_(object("tmpYo", flows.mesh()), flows.mesh()),
      tmpYfField_(object("tmpYf", flows.mesh()), flows.mesh()), oIdx_(1), fIdx_(0),
      tauIntFactor_(0.1) {
    if (cont.subController("combustionModel").found("finiteRate")) {
        minDtFactor_ = cont.subController("combustionModel")
                           .subController("finiteRate")
                           .findOrDefault<real>("minDtFactor", 0.1);
        const auto &ftrcont = cont.subController("combustionModel").subController("finiteRate");
        /*if (cont.subController("combustionModel").subController("finiteRate").found("instantaneousReactionRate"))
        {
                auto w = cont.subController("combustionModel").subController("finiteRate").findWord("instantaneousReactionRate");
                trim(w);
                stringToUpperCase(w);
                if (w == "OFF")
                {
                        instantaneousReactionRate_ = false;
                }
                else if (w != "ON")
                {
                        errorAbortStr("Unkown string : " + w + " for instantaneousReactionRate entry in " + cont.subController("combustionModel").subController("finiteRate").name(), HUR_FUNCTION);
                }
        }*/

        tauIntFactor_ = ftrcont.findOrDefault<real>("tauIntFactor", tauIntFactor_);
    }
}

void OpenHurricane::finiteRate::calcSourceTerms(const realArray &dt) {
    chemistryPtr_->calculateSourceTerms();
}

void OpenHurricane::finiteRate::expChemistrySource(realArray &dt, const bool isModifiedDt) {
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
                    yi[i].rhs()[cellI] += omega * vol[cellI];
            }
    }*/
}

void OpenHurricane::finiteRate::impCalcSourceTerms(const realArray &dt) {
    chemistryPtr_->calculateSourceTermsImp();
}

void OpenHurricane::finiteRate::impChemistrySource(realArray &dt, const bool isModifiedDt) {
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
    //		yi[i].rhs()[cellI] += omega * vol[cellI];
    //		if (instantaneousReactionRate_)
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
    //					minDt = min(minDt, mag(inv(yi[i].diagSource()[cellI])));
    //				}
    //			}
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
    //		{	// Reduce the time step and limit to 0.5dt. Do not reduce too much.
    //			dt[cellI] = max(minDt, minDtFactor_ * dt[cellI]);
    //		}
    //	}
    //}
}

void OpenHurricane::finiteRate::getChemistrySource() {
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

void OpenHurricane::finiteRate::getImpChemistrySource(realArray &dt, const bool isModifiedDt) {
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
            yi[i].rhs()[cellI] += omega * vol[cellI];

            if (isnan(yi[i].diagSource()[cellI]) || isinf(yi[i].diagSource()[cellI])) {
                yi[i].diagSource()[cellI] = 0;
                continue;
            } else {
                if (isModifiedDt) {
                    minDt = min(minDt, mag(inv(yi[i].diagSource()[cellI])));
                }
            }
            yi[i].diagSource()[cellI] *= vol[cellI];
        }
        if (isModifiedDt) {
            if (minDt <
                dt[cellI]) { // Reduce the time step and limit to 0.5dt. Do not reduce too much.
                dt[cellI] = max(minDt, minDtFactor_ * dt[cellI]);
            }
        }
    }
}

void OpenHurricane::finiteRate::fullPointImpChemistrySource(
    realArray &dt, cellRealSquareMatrixArray &Jac, const integer rhoId, const integer rhouId,
    const integer rhoEId, const integer rhoYi0Id) {
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
            };
            yi[i].rhs()[cellI] += omega * vol[cellI];

            const auto Wi = species().W(i);
            const auto WNs = species().W(nsp - 1);
            const real dRiddT = Wi * dwdciT(i, Tid);
            const real dRiddrhoY_Ns = Wi / WNs * dwdciT(i, nsp - 1);

            const real factor = vol[cellI];

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

OpenHurricane::realArray OpenHurricane::finiteRate::heatReleaseRate() {
    if (chemOptions_ == chemistryOption::strangSplitted) {
        chemistryPtr_->calculateSourceTerms();
    }
    return chemistryPtr_->heatReleaseRate();
}

void OpenHurricane::finiteRate::constVolReactor(const realArray &dt, const real dtFactor,
                                                const bool firstCall) {
    auto &spt = species();
    const integer nsp = spt.size();         // The size of species
    const integer nrc = reactions().size(); // The size of reactions
    auto &T = flows_.T();
    auto &p = flows_.p();
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
}

OpenHurricane::realArrayArray OpenHurricane::finiteRate::omegai() const {
    if (chemOptions_ == chemistryOption::strangSplitted) {
        chemistryPtr_->calculateSourceTerms();
    }
    return chemistryPtr_->Ri();
}

void OpenHurricane::finiteRate::expChemistrySourceCoupled() {
    chemistryPtr_->calculateSourceTerms();

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

void OpenHurricane::finiteRate::diagImpChemistrySourceCoupled() {
    chemistryPtr_->calculateSourceTermsImp();

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
            yi[i].rhs()[cellI] += omega * vol[cellI];

            if (isnan(yi[i].diagSource()[cellI]) || isinf(yi[i].diagSource()[cellI])) {
                yi[i].diagSource()[cellI] = 0;
                continue;
            }
            yi[i].diagSource()[cellI] *= vol[cellI];
        }
    }
}

void OpenHurricane::finiteRate::fullImpChemistrySourceCoupled(cellRealSquareMatrixArray &Jac,
                                                              const integer rhoId,
                                                              const integer rhouId,
                                                              const integer rhoEId,
                                                              const integer rhoYi0Id) {
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
            };
            yi[i].rhs()[cellI] += omega * vol[cellI];

            const auto Wi = species().W(i);
            const auto WNs = species().W(nsp - 1);
            const real dRiddT = Wi * dwdciT(i, Tid);
            const real dRiddrhoY_Ns = Wi / WNs * dwdciT(i, nsp - 1);

            const real factor = vol[cellI];

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

void OpenHurricane::finiteRate::chemistrySourceIntegrated() {
    auto &tau = tauInt();

    const auto &T = flows_.T();
    const auto &rho = flows_.rho();
    auto &yi = mixtures().Yi();
    const auto &vol = flows_.mesh().cellVol();
    const auto &v = flows_.v();
    const auto &mul = flows_.mul();

    for (integer cellI = 0; cellI < T.mesh().nCells(); ++cellI) {
        if (flows_.mesh().Iteration().hasPhysicalTimeStep()) {
            tau[cellI] = flows_.mesh().Iteration().pTStep().pTimeStep();
        } else {
            const auto dv = v[cellI].magnitude();

            const auto dx = pow(vol[cellI], real(1.0 / 3.0));

            tau[cellI] = dx / max(dv, tiny);

            if (mul.size() != 0) {
                tau[cellI] = min(tau[cellI], rho[cellI] * sqr(dx) / max(veryTiny, mul[cellI]));
            }
            tau[cellI] *= tauIntFactor_;
        }
    }

    chemistryPtr_->solve(tau, real(1.0));
}

void OpenHurricane::finiteRate::chemistrySourceStrangSplitted(const realArray &dt,
                                                              const real dtFactor) {
    for (integer n = 0; n < maxSubStep_; ++n) {
        chemistryPtr_->solveUpdate(dt, dtFactor / real(maxSubStep_));
    }
}
