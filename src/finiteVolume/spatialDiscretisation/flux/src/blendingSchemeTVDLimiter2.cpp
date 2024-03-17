/*!
 * \file blendingSchemeTVDLimiter2.cpp
 * \brief Main subroutines for the mix of centered and upwind-biased Riemann flux scheme with a TVD limiter 2.
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

#include "blendingSchemeTVDLimiter2.hpp"
#include "HLLC.hpp"

namespace OpenHurricane {
    createClassNameStr(blendingSchemeTVDLimiter2, "blendingSchemeTVDLimiter2");
    registerObjFty(spatialScheme, blendingSchemeTVDLimiter2, controller);
} // namespace OpenHurricane

void OpenHurricane::blendingSchemeTVDLimiter2::makeAlphaf() const {
    HurDelete(alphafPtr_);

    alphafPtr_ = new realArray(mesh_.nFaces(), minAlphaf_);
    auto &apf = *alphafPtr_;
    const auto fsk = meshCheckers::faceSkewness(mesh_, mesh_.points(), mesh_.faceCentre(),
                                                mesh_.faceArea(), mesh_.cellCentre());

    for (integer fi = 0; fi < mesh_.nFaces(); ++fi) {
        if (fsk[fi] >= 1.0) {
            apf[fi] = 1.0;
        } else if (fsk[fi] <= 1e-2) {
            apf[fi] = minAlphaf_;
        } else {
            apf[fi] = (1 - minAlphaf_) * sqrt(fsk[fi]) + minAlphaf_;
        }
    }

    const auto &aspr = mesh_.aspectRatio();

    for (integer n = 0; n < aspr.size(); ++n) {
        if (aspr[n] > 100.0) {
            for (integer fi = 0; fi < mesh_.cells()[n].faceSize(); ++fi) {
                const auto facei = mesh_.cells()[n].facei(fi);
                apf[facei] = min(maxAlphaf_, max(minAlphaf_, max(real(0.3), apf[facei])));
            }
        } else if (aspr[n] >= 10.0) {
            for (integer fi = 0; fi < mesh_.cells()[n].faceSize(); ++fi) {
                const auto facei = mesh_.cells()[n].facei(fi);
                apf[facei] = min(maxAlphaf_, max(minAlphaf_, max(real(0.1), apf[facei])));
            }
        }
    }

    const auto fot = meshCheckers::faceOrthogonality(mesh_, mesh_.faceArea(), mesh_.faceCentre(),
                                                     mesh_.cellCentre());

    for (integer fi = 0; fi < mesh_.nFaces(); ++fi) {
        if (fot[fi] <= 0.5) {
            apf[fi] = 1.0;
        } else if (fot[fi] >= 1) {
            apf[fi] = max(apf[fi], minAlphaf_);
        } else {
            apf[fi] = max(apf[fi],
                          real(4) * (real(1) - minAlphaf_) * sqr(real(1) - fot[fi]) + minAlphaf_);
        }
    }

#ifdef HUR_DEBUG
    real maxAlphaf = minAlphaf_;
#endif // HUR_DEBUG
    for (integer fi = 0; fi < mesh_.nFaces(); ++fi) {
        apf[fi] = min(maxAlphaf_, max(apf[fi], minAlphaf_));

#ifdef HUR_DEBUG
        maxAlphaf = max(maxAlphaf, apf[fi]);
#endif // HUR_DEBUG
    }

#ifdef HUR_DEBUG
    HurMPI::reduce(maxAlphaf, MPI_MAX);
    Pout << "    Info: The maximum alphaf is " << toString(maxAlphaf) << std::endl;
#endif // HUR_DEBUG
}

void OpenHurricane::blendingSchemeTVDLimiter2::cntUpwdBlendFct() {
    if (mesh().Iteration().hasSubIteration()) {
        if (mesh().Iteration().subIter().cSubStep() >= 2) {
            return;
        }
    }
    shocks_.calcShockSensor();
    upwindFactor_ = Zero;
    extendShockSensor();
    const auto &fZ = mesh().faceZones();
    const auto &fL = mesh().faces();
    const auto &v = v_;
    const auto &cv = flows_.mesh().cellVolume();
    shocks_.checkTemperatureDiffer(upwindFactor_);

    const auto &TFlag = flows_.temperatureFlag();
    const auto &T = flows_.T();
    const auto &pFlag = flows_.pressureFlag();

    for (integer n = 0; n < mesh().nCells(); ++n) {
        if (TFlag[n] != 1 || pFlag[n] != 1) {
            upwindFactor_[n] = 1.0;
            continue;
        }
        if (limitTemperature_) {
            if (T[n] >= THighLimit_ || T[n] <= TLowLimit_) {
                upwindFactor_[n] = maxAlphaf_;
                //continue;
            }
        }
        const real S = skewMagnitude(v.grad()[n]);

        const real dil = v.grad()[n].xx() + v.grad()[n].yy() + v.grad()[n].zz();

        const real dil2 = sqr(dil);
        const real CSS = dil2 / (dil2 + sqr(S) + tiny);

        upwindFactor_[n] *= CSS;
    }

    realTransfer myTransfer(mesh(), upwindFactor_, false, true);
    myTransfer.transferInit();
    for (integer fZI = 0; fZI < fZ.size(); ++fZI) {
        if (fZ[fZI].isInterior()) {
            continue;
        } else if (fZ[fZI].isCutFace()) {
            continue;
        } else if (fZ[fZI].isPeriodic() || fZ[fZI].isPeriodicShadow()) {
            continue;
        } else {
            for (integer fI = fZ[fZI].firstIndex(); fI <= fZ[fZI].lastIndex(); ++fI) {
                const auto &cl = fL[fI].leftCell();
                const auto &cr = fL[fI].rightCell();
                upwindFactor_[cr] = upwindFactor_[cl];
            }
        }
    }
    myTransfer.transferring();
}

void OpenHurricane::blendingSchemeTVDLimiter2::calcFlux(
    const real rhol, const real rhor, const vector &VL, const vector &VR, const real pl,
    const real pr, const real gl, const real gr, const real el, const real er, const real cl,
    const real cr, const vector &faceArea, const real blend, realArray &flux) const {
    LFatal("Attempt to access an empty function");
}

void OpenHurricane::blendingSchemeTVDLimiter2::upwindFlux(
    const real rhol, const real rhor, const vector &VL, const vector &VR, const real pl,
    const real pr, const real gl, const real gr, const real el, const real er, const real cl,
    const real cr, const vector &faceArea, const real blend, realArray &flux) const {
    upwindMethodPtr_->calcFlux(rhol, rhor, VL, VR, pl, pr, gl, gr, el, er, cl, cr, faceArea, blend,
                               flux);
}

void OpenHurricane::blendingSchemeTVDLimiter2::pureCentralFlux(
    const real rhol, const real rhor, const vector &VL, const vector &VR, const real pl,
    const real pr, const real el, const real er, const vector &faceArea, const integer fi,
    const real fw, realArray &flux) const {
    const real rhof = fw * rhol + (1.0 - fw) * rhor;
    const vector rhovf = fw * rhol * VL + (real(1.0) - fw) * rhor * VR;
    const vector vf = fw * VL + (real(1.0) - fw) * VR;
    const real pf = fw * pl + (real(1.0) - fw) * pr;

    const real rhoef = fw * rhol * el + (real(1.0) - fw) * rhor * er;

    const auto n = faceArea.normalized();
    const auto Sf = faceArea.magnitude();
    const real Vn = vf * n;
    vnSf_[fi] = Vn * Sf;
    flux[0] = rhof * Vn * Sf;

    flux[1] = (rhovf.x() * Vn + pf * n.x()) * Sf;
    flux[2] = (rhovf.y() * Vn + pf * n.y()) * Sf;
    flux[3] = (rhovf.z() * Vn + pf * n.z()) * Sf;
    flux[4] = (rhoef + pf) * Vn * Sf;
}

OpenHurricane::blendingSchemeTVDLimiter2::blendingSchemeTVDLimiter2(const controller &cont,
                                                                    const runtimeMesh &mesh,
                                                                    flowModel &flow)
    : spatialScheme(cont, mesh, flow), upwindMethodPtr_(nullptr), fluxC_(5, Zero), fluxU_(5, Zero),
      vnSf_(object("VnSF", mesh, object::NOT_WRITE), mesh),
      shocks_(flow, cont.subController("blendingSchemeTVDLimiter2")),
      upwindFactor_(object("upwindFactor", mesh, object::NOT_WRITE), mesh),
      minAlphaf_(
          cont.subController("blendingSchemeTVDLimiter2").findOrDefault<real>("minAlphaf", 0.1)),
      maxAlphaf_(
          cont.subController("blendingSchemeTVDLimiter2").findOrDefault<real>("maxAlphaf", 1.0)),
      alphafPtr_(nullptr), limitTemperature_(true), THighLimit_(3000.0), TLowLimit_(100.0) {
    const auto &blendCont = cont.subController("blendingSchemeTVDLimiter2");
    controllerSwitch myCont(blendCont);
    limitTemperature_ = myCont("checkTemperature", limitTemperature_);
    THighLimit_ = blendCont.findOrDefault<real>("THighLimit", THighLimit_);
    TLowLimit_ = blendCont.findOrDefault<real>("TLowLimit", TLowLimit_);

    if (blendCont.found("upwindScheme")) {
        upwindMethodPtr_ = upwindScheme::creator(blendCont, mesh, flow);
    } else {
        LFatal("Upwind scheme method is not set");
    }
    if (minAlphaf_ < 0 || minAlphaf_ > 1) {
        errorAbortStr(("The value of minAlpha must be within [0,1], while it is set to " +
                       toString(minAlphaf_)));
    }
    if (maxAlphaf_ < 0 || maxAlphaf_ > 1) {
        errorAbortStr(("The value of maxAlphaf must be within [0,1], while it is set to " +
                       toString(maxAlphaf_)));
    }

    makeAlphaf();
}

OpenHurricane::blendingSchemeTVDLimiter2::~blendingSchemeTVDLimiter2() noexcept {
    HurDelete(alphafPtr_);
}

void OpenHurricane::blendingSchemeTVDLimiter2::basicFlux() {
    const auto &mesh = v_.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();

    const auto &fA = mesh.faceArea();
    const auto &fC = mesh.faceCntr();
    const auto &fW = mesh.faceWeight();
    const auto &cC = mesh.cellCentre();

    auto &v = v_;
    auto &rho = const_cast<rhoThermo &>(thermo_).rho();
    auto &p = const_cast<rhoThermo &>(thermo_).p();
    auto &T = const_cast<rhoThermo &>(thermo_).T();
    auto &E = const_cast<rhoThermo &>(thermo_).E();
    auto &gama = const_cast<rhoThermo &>(thermo_).gamma();

    rho.rhs() = Zero;
    v.rhs() = Zero;
    T.rhs() = Zero;
    E.rhs() = Zero;

    reconstrPtr_->calcGrad(rho);
    reconstrPtr_->calcGrad(v);
    reconstrPtr_->calcGrad(p);
    reconstrPtr_->calcGrad(T);
    reconstrPtr_->calcGrad(E);
    reconstrPtr_->calcGrad(gama);

    reconstrPtr_->calcLimiter(rho);
    reconstrPtr_->calcLimiter(v);
    reconstrPtr_->calcLimiter(p);
    //reconstrPtr_->calcLimiter(T);
    reconstrPtr_->calcLimiter(E);
    reconstrPtr_->calcLimiter(gama);

    mixture &mixtures = const_cast<rhoThermo &>(thermo_).mixtures();
    realArray psia(5, Zero);
    cntUpwdBlendFct();

    //spectralRadius();

    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        // Symmetric plane
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const integer &cl = fl[fi].leftCell();

                const vector vFlux = p[cl] * fA[fi];
                v.rhs()[cl] += vFlux;
                rhoFlux_[fi] = Zero;
                vnSf_[fi] = Zero;
            }
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                real rhol, rhor;
                real pl, pr;
                vector vl, vr;
                real gaml, gamr;
                real el, er;
                real al;
                real ar;
                if (!switchToUpwind(fi)) {
                    reconstrPtr_->calcReconstruction(rho, fi, rhol, rhor);
                    reconstrPtr_->calcReconstruction(p, fi, pl, pr);
                    if (rhol <= 0.0 || rhor <= 0.0 || pl <= 0.0 || pr <= 0.0) {
                        pl = p[cl];
                        pr = p[cr];
                        rhol = rho[cl];
                        rhor = rho[cr];
                        vl = v[cl];
                        vr = v[cr];
                        gaml = gama[cl];
                        gamr = gama[cr];
                        el = E[cl];
                        er = E[cr];
                    } else {
                        reconstrPtr_->calcReconstruction(v, fi, vl, vr);
                        reconstrPtr_->calcReconstruction(gama, fi, gaml, gamr);
                        reconstrPtr_->calcReconstruction(E, fi, el, er);
                        if (isLowMachCorr_) {
                            lowMachCorrection(rhol, rhor, gaml, gamr, vl, vr, pl, pr, fA[fi]);
                        }
                    }
                    const real upwdFct = actualUpwdFct(fi);
                    al = sqrt(gaml * pl / rhol);
                    ar = sqrt(gamr * pr / rhor);

                    real blend = min(shockFactor_[cl], shockFactor_[cr]);

                    upwindFlux(rhor, rhol, vr, vl, pr, pl, gamr, gaml, er, el, ar, al, fA[fi],
                               blend, fluxU_);

                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];

                    const real rholl = rho[cl] + rho.grad()[cl] * rL;
                    const real rhorr = rho[cr] + rho.grad()[cr] * rR;

                    const vector vll = v[cl] + v.grad()[cl] * rL;
                    const vector vrr = v[cr] + v.grad()[cr] * rR;

                    const real pll = p[cl] + p.grad()[cl] * rL;
                    const real prr = p[cr] + p.grad()[cr] * rR;

                    const real ell = E[cl] + E.grad()[cl] * rL;
                    const real err = E[cr] + E.grad()[cr] * rR;
                    pureCentralFlux(rhorr, rholl, vrr, vll, prr, pll, err, ell, fA[fi], fi, 0.5,
                                    fluxC_);
                    flux_ = (real(1) - upwdFct) * fluxC_ + upwdFct * fluxU_;

                    /*flux_[0] += (real(1) - upwdFct) * artificialDissipation(fi, rhorr - rholl);
                    flux_[1] += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * vrr.x() - rholl * vll.x());
                    flux_[2] += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * vrr.y() - rholl * vll.y());
                    flux_[3] += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * vrr.z() - rholl * vll.z());
                    flux_[4] += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * err - rholl * ell);*/

                    rhoFlux_[fi] = fluxU_[0];
                } else // Upwind scheme
                {
                    reconstrPtr_->calcReconstruction(rho, fi, rhol, rhor);
                    reconstrPtr_->calcReconstruction(p, fi, pl, pr);
                    if (rhol <= 0.0 || rhor <= 0.0 || pl <= 0.0 || pr <= 0.0) {
                        pl = p[cl];
                        pr = p[cr];
                        rhol = rho[cl];
                        rhor = rho[cr];
                        vl = v[cl];
                        vr = v[cr];
                        gaml = gama[cl];
                        gamr = gama[cr];
                        el = E[cl];
                        er = E[cr];
                    } else {
                        reconstrPtr_->calcReconstruction(v, fi, vl, vr);
                        reconstrPtr_->calcReconstruction(gama, fi, gaml, gamr);
                        reconstrPtr_->calcReconstruction(E, fi, el, er);
                        if (isLowMachCorr_) {
                            lowMachCorrection(rhol, rhor, gaml, gamr, vl, vr, pl, pr, fA[fi]);
                        }
                    }
                    vnSf_[fi] = Zero;
                    al = sqrt(gaml * pl / rhol);
                    ar = sqrt(gamr * pr / rhor);
                    flux_ = Zero;
                    real blend = min(shockFactor_[cl], shockFactor_[cr]);
                    upwindFactor_[cl] = 1;
                    upwindFactor_[cr] = 1;
                    upwindFlux(rhor, rhol, vr, vl, pr, pl, gamr, gaml, er, el, ar, al, fA[fi],
                               blend, flux_);
                    rhoFlux_[fi] = flux_[0];
                }
                if (std::isnan(flux_[0]) || isinf(flux_[0])) {
                    std::string errMsg;
                    errMsg = "density inv_flux is not a number\n";
                    errMsg += " rhol = " + toString(rhol);
                    errMsg += " rhor = " + toString(rhor) + "\n";
                    errMsg += " pl = " + toString(pl);
                    errMsg += " pr = " + toString(pr) + "\n";
                    errMsg += " In face center: " + toString(fC[fi]);
                    errorAbortStr(errMsg);
                }
                if (std::isnan(flux_[4]) || isinf(flux_[4])) {
                    std::string errMsg;
                    errMsg = "energy inv_flux is not a number\n";
                    errMsg += " el = " + toString(el);
                    errMsg += " er = " + toString(er) + "\n";
                    errMsg += " In face center: " + toString(fC[fi]);
                    errorAbortStr(errMsg);
                }

                rho.rhs()[cl] += flux_[0];
                v.rhs()[cl][0] += flux_[1];
                v.rhs()[cl][1] += flux_[2];
                v.rhs()[cl][2] += flux_[3];
                E.rhs()[cl] += flux_[4];

                if (cr < mesh.nCells()) {
                    rho.rhs()[cr] -= flux_[0];
                    v.rhs()[cr][0] -= flux_[1];
                    v.rhs()[cr][1] -= flux_[2];
                    v.rhs()[cr][2] -= flux_[3];
                    E.rhs()[cr] -= flux_[4];
                }
            }
        }
    }
}

void OpenHurricane::blendingSchemeTVDLimiter2::invFlux(cellRealArray &cellQ) const {
    cellQ.rhs() = Zero;

    const auto &rho = thermo_.rho();
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);

    const auto &mesh = cellQ.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fW = mesh.faceWeight();

    const auto &fA = mesh.faceArea();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            continue;
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();

                real flux = Zero;
                if (!switchToUpwind(fi)) {
                    real ql, qr;
                    reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);

                    const real upwdFct = actualUpwdFct(fi);

                    // upwind flux
                    real fluxU = Zero;
                    if (rhoFlux_[fi] >= 0.0) {
                        fluxU = rhoFlux_[fi] * qr;
                    } else {
                        fluxU = rhoFlux_[fi] * ql;
                    }

                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];
                    const real rholl = rho[cl] + rho.grad()[cl] * rL;
                    const real rhorr = rho[cr] + rho.grad()[cr] * rR;
                    const auto qll = cellQ[cl] + cellQ.grad()[cl] * rL;
                    const auto qrr = cellQ[cr] + cellQ.grad()[cr] * rR;

                    const real rhoQf = 0.5 * (rholl * qll + rhorr * qrr);
                    //const real rhoQf = fW[fi] * rho[cl] * cellQ[cl] + (1.0 - fW[fi]) * rho[cr] * cellQ[cr];

                    const real fluxC = rhoQf * vnSf_[fi];

                    flux = (1.0 - upwdFct) * fluxC + upwdFct * fluxU;

                    //flux += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * qrr - rholl * qll);
                } else // upwind scheme
                {
                    real ql, qr;
                    reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                    if (rhoFlux_[fi] >= 0.0) {
                        flux = rhoFlux_[fi] * qr;
                    } else {
                        flux = rhoFlux_[fi] * ql;
                    }
                }
                //Add the numerical flux to the right-hand-term of both cells sharing the face
                cellQ.rhs()[cl] += flux;
                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= flux;
                }
            }
        }
    }
}

void OpenHurricane::blendingSchemeTVDLimiter2::invFlux(cellVectorArray &cellQ) const {
    const auto &rho = thermo_.rho();
    cellQ.rhs() = Zero;
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);
    const auto &mesh = cellQ.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fW = mesh.faceWeight();

    const auto &fA = mesh.faceArea();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            continue;
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();

                vector flux = Zero;
                if (!switchToUpwind(fi)) {
                    vector ql, qr;
                    reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);

                    const real upwdFct = actualUpwdFct(fi);

                    // upwind flux
                    vector fluxU = Zero;
                    if (rhoFlux_[fi] >= 0.0) {
                        fluxU = rhoFlux_[fi] * qr;
                    } else {
                        fluxU = rhoFlux_[fi] * ql;
                    }
                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];
                    const real rholl = rho[cl] + rho.grad()[cl] * rL;
                    const real rhorr = rho[cr] + rho.grad()[cr] * rR;
                    const auto qll = cellQ[cl] + cellQ.grad()[cl] * rL;
                    const auto qrr = cellQ[cr] + cellQ.grad()[cr] * rR;

                    const vector rhoQf = real(0.5) * (rholl * qll + rhorr * qrr);
                    //const vector rhoQf = fW[fi] * rho[cl] * cellQ[cl] + (1.0 - fW[fi]) * rho[cr] * cellQ[cr];

                    const vector fluxC = rhoQf * vnSf_[fi];

                    flux = (real(1.0) - upwdFct) * fluxC + upwdFct * fluxU;
                    //flux += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * qrr - rholl * qll);
                } else // upwind scheme
                {
                    vector ql, qr;
                    reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                    if (rhoFlux_[fi] >= 0.0) {
                        flux = rhoFlux_[fi] * qr;
                    } else {
                        flux = rhoFlux_[fi] * ql;
                    }
                }
                //Add the numerical flux to the right-hand-term of both cells sharing the face
                cellQ.rhs()[cl] += flux;
                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= flux;
                }
            }
        }
    }
}

void OpenHurricane::blendingSchemeTVDLimiter2::invFlux(cellRealArray &cellQ,
                                                       const realBounded &bound) const {
    const auto &rho = thermo_.rho();
    cellQ.rhs() = Zero;
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);
    const auto &mesh = cellQ.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fW = mesh.faceWeight();

    const auto &fA = mesh.faceArea();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            continue;
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                real flux = Zero;
                if (!switchToUpwind(fi)) {
                    real ql, qr;
                    reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                    ql = bound.bounding(ql, cellQ[cl]);
                    qr = bound.bounding(qr, cellQ[cr]);

                    const real upwdFct = actualUpwdFct(fi);

                    // upwind flux
                    real fluxU = Zero;
                    if (rhoFlux_[fi] >= 0.0) {
                        fluxU = rhoFlux_[fi] * qr;
                    } else {
                        fluxU = rhoFlux_[fi] * ql;
                    }
                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];
                    const real rholl = rho[cl] + rho.grad()[cl] * rL;
                    const real rhorr = rho[cr] + rho.grad()[cr] * rR;
                    auto qll = cellQ[cl] + cellQ.grad()[cl] * rL;
                    auto qrr = cellQ[cr] + cellQ.grad()[cr] * rR;
                    qll = bound.bounding(qll, cellQ[cl]);
                    qrr = bound.bounding(qrr, cellQ[cr]);
                    const real rhoQf = 0.5 * (rholl * qll + rhorr * qrr);
                    //const real rhoQf = fW[fi] * rho[cl] * cellQ[cl] + (1.0 - fW[fi]) * rho[cr] * cellQ[cr];

                    const real fluxC = rhoQf * vnSf_[fi];

                    flux = (1.0 - upwdFct) * fluxC + upwdFct * fluxU;
                    // flux += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * qrr - rholl * qll);
                } else // upwind scheme
                {
                    real ql, qr;
                    reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                    ql = bound.bounding(ql, cellQ[cl]);
                    qr = bound.bounding(qr, cellQ[cr]);
                    if (rhoFlux_[fi] >= 0.0) {
                        flux = rhoFlux_[fi] * qr;
                    } else {
                        flux = rhoFlux_[fi] * ql;
                    }
                }
                //Add the numerical flux to the right-hand-term of both cells sharing the face
                cellQ.rhs()[cl] += flux;
                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= flux;
                }
            }
        }
    }
}

void OpenHurricane::blendingSchemeTVDLimiter2::invFlux(cellVectorArray &cellQ,
                                                       const vectorBounded &bound) const {
    const auto &rho = thermo_.rho();
    cellQ.rhs() = Zero;
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);
    const auto &mesh = cellQ.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fW = mesh.faceWeight();

    const auto &fA = mesh.faceArea();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            continue;
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                vector flux = Zero;
                if (!switchToUpwind(fi)) {
                    vector ql, qr;
                    reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                    ql = bound.bounding(ql, cellQ[cl]);
                    qr = bound.bounding(qr, cellQ[cr]);

                    const real upwdFct = actualUpwdFct(fi);

                    // upwind flux
                    vector fluxU = Zero;
                    if (rhoFlux_[fi] >= 0.0) {
                        fluxU = rhoFlux_[fi] * qr;
                    } else {
                        fluxU = rhoFlux_[fi] * ql;
                    }
                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];
                    const real rholl = rho[cl] + rho.grad()[cl] * rL;
                    const real rhorr = rho[cr] + rho.grad()[cr] * rR;
                    auto qll = cellQ[cl] + cellQ.grad()[cl] * rL;
                    auto qrr = cellQ[cr] + cellQ.grad()[cr] * rR;
                    qll = bound.bounding(qll, cellQ[cl]);
                    qrr = bound.bounding(qrr, cellQ[cr]);
                    const vector rhoQf = 0.5 * (rholl * qll + rhorr * qrr);
                    //const vector rhoQf = fW[fi] * rho[cl] * cellQ[cl] + (1.0 - fW[fi]) * rho[cr] * cellQ[cr];

                    const vector fluxC = rhoQf * vnSf_[fi];

                    flux = (real(1.0) - upwdFct) * fluxC + upwdFct * fluxU;
                    //flux += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * qrr - rholl * qll);
                } else // upwind scheme
                {
                    vector ql, qr;
                    reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                    ql = bound.bounding(ql, cellQ[cl]);
                    qr = bound.bounding(qr, cellQ[cr]);
                    if (rhoFlux_[fi] >= 0.0) {
                        flux = rhoFlux_[fi] * qr;
                    } else {
                        flux = rhoFlux_[fi] * ql;
                    }
                }
                //Add the numerical flux to the right-hand-term of both cells sharing the face
                cellQ.rhs()[cl] += flux;
                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= flux;
                }
            }
        }
    }
}

void OpenHurricane::blendingSchemeTVDLimiter2::invFluxSpecies(PtrList<cellRealArray> &yi,
                                                              const bool withLastSpc) const {
    const auto &rho = thermo_.rho();
    if (withLastSpc) {
        for (integer i = 0; i < yi.size(); ++i) {
            yi[i].rhs() = Zero;
            reconstrPtr_->calcGrad(yi[i]);
            reconstrPtr_->calcLimiter(yi[i]);
        }
    } else {
        for (integer i = 0; i < yi.size() - 1; ++i) {
            yi[i].rhs() = Zero;
            reconstrPtr_->calcGrad(yi[i]);
            reconstrPtr_->calcLimiter(yi[i]);
        }
    }
    const auto &mesh = yi[0].mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fW = mesh.faceWeight();

    const auto &fA = mesh.faceArea();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    realArray ql;
    realArray qr;
    if (withLastSpc) {
        ql.resize(yi.size());
        qr.resize(yi.size());
    } else {
        ql.resize(yi.size() - 1);
        qr.resize(yi.size() - 1);
    }

    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            continue;
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                if (withLastSpc) {
                    real yilt = Zero;
                    real yirt = Zero;
                    real flux = Zero;
                    if (!switchToUpwind(fi)) {
                        const real upwdFct = actualUpwdFct(fi);
                        for (integer i = 0; i < yi.size(); ++i) {
                            reconstrPtr_->calcReconstruction(yi[i], fi, ql[i], qr[i]);
                            ql[i] = min(real(1), max(real(0), ql[i]));
                            qr[i] = min(real(1), max(real(0), qr[i]));
                            yilt += ql[i];
                            yirt += qr[i];
                        }
                        for (integer i = 0; i < yi.size(); ++i) {
                            ql[i] /= max(tiny, yilt);
                            qr[i] /= max(tiny, yirt);

                            // upwind flux
                            real fluxU = Zero;
                            if (rhoFlux_[fi] >= 0.0) {
                                fluxU = rhoFlux_[fi] * qr[i];
                            } else {
                                fluxU = rhoFlux_[fi] * ql[i];
                            }
                            const vector &rL = faceLeftCellCtr[fi];
                            const vector &rR = faceRightCellCtr[fi];
                            const real rholl = rho[cl] + rho.grad()[cl] * rL;
                            const real rhorr = rho[cr] + rho.grad()[cr] * rR;
                            auto qll = yi[i][cl] + yi[i].grad()[cl] * rL;
                            auto qrr = yi[i][cr] + yi[i].grad()[cr] * rR;
                            qll = min(real(1), max(real(0), qll));
                            qrr = min(real(1), max(real(0), qrr));
                            const real rhoQf = 0.5 * (rholl * qll + rhorr * qrr);
                            //const real rhoQf = fW[fi] * rho[cl] * yi[i][cl] + (1.0 - fW[fi]) * rho[cr] * yi[i][cr];

                            const real fluxC = rhoQf * vnSf_[fi];

                            flux = (1.0 - upwdFct) * fluxC + upwdFct * fluxU;
                            //flux += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * qrr - rholl * qll);

                            //Add the numerical flux to the right-hand-term of both cells sharing the face
                            yi[i].rhs()[cl] += flux;
                            if (cr < mesh.nCells()) {
                                yi[i].rhs()[cr] -= flux;
                            }
                        }
                    } else {
                        for (integer i = 0; i < yi.size(); ++i) {
                            reconstrPtr_->calcReconstruction(yi[i], fi, ql[i], qr[i]);
                            ql[i] = min(real(1), max(real(0), ql[i]));
                            qr[i] = min(real(1), max(real(0), qr[i]));
                            yilt += ql[i];
                            yirt += qr[i];
                        }
                        for (integer i = 0; i < yi.size(); ++i) {
                            ql[i] /= max(tiny, yilt);
                            qr[i] /= max(tiny, yirt);
                            if (rhoFlux_[fi] >= 0.0) {
                                flux = rhoFlux_[fi] * qr[i];
                            } else {
                                flux = rhoFlux_[fi] * ql[i];
                            }

                            //Add the numerical flux to the right-hand-term of both cells sharing the face
                            yi[i].rhs()[cl] += flux;
                            if (cr < mesh.nCells()) {
                                yi[i].rhs()[cr] -= flux;
                            }
                        }
                    }
                } else {
                    real yilt = Zero;
                    real yirt = Zero;
                    real flux = Zero;
                    if (!switchToUpwind(fi)) {
                        const real upwdFct = actualUpwdFct(fi);
                        for (integer i = 0; i < yi.size() - 1; ++i) {
                            reconstrPtr_->calcReconstruction(yi[i], fi, ql[i], qr[i]);
                            ql[i] = min(real(1), max(real(0), ql[i]));
                            qr[i] = min(real(1), max(real(0), qr[i]));
                            yilt += ql[i];
                            yirt += qr[i];
                        }
                        real yiLastL = 1.0 - yilt;
                        real yiLastR = 1.0 - yirt;

                        for (integer i = 0; i < yi.size() - 1; ++i) {
                            if (yiLastL < 0) {
                                ql[i] /= yilt;
                            }
                            if (yiLastR < 0) {
                                qr[i] /= yirt;
                            }

                            // upwind flux
                            real fluxU = Zero;
                            if (rhoFlux_[fi] >= 0.0) {
                                fluxU = rhoFlux_[fi] * qr[i];
                            } else {
                                fluxU = rhoFlux_[fi] * ql[i];
                            }
                            const vector &rL = faceLeftCellCtr[fi];
                            const vector &rR = faceRightCellCtr[fi];
                            const real rholl = rho[cl] + rho.grad()[cl] * rL;
                            const real rhorr = rho[cr] + rho.grad()[cr] * rR;
                            auto qll = yi[i][cl] + yi[i].grad()[cl] * rL;
                            auto qrr = yi[i][cr] + yi[i].grad()[cr] * rR;
                            qll = min(real(1), max(real(0), qll));
                            qrr = min(real(1), max(real(0), qrr));
                            const real rhoQf = 0.5 * (rholl * qll + rhorr * qrr);
                            //const real rhoQf = fW[fi] * rho[cl] * yi[i][cl] + (1.0 - fW[fi]) * rho[cr] * yi[i][cr];

                            const real fluxC = rhoQf * vnSf_[fi];

                            flux = (1.0 - upwdFct) * fluxC + upwdFct * fluxU;
                            //flux += (real(1) - upwdFct) * artificialDissipation(fi, rhorr * qrr - rholl * qll);

                            //Add the numerical flux to the right-hand-term of both cells sharing the face
                            yi[i].rhs()[cl] += flux;
                            if (cr < mesh.nCells()) {
                                yi[i].rhs()[cr] -= flux;
                            }
                        }
                    } else {
                        for (integer i = 0; i < yi.size() - 1; ++i) {
                            reconstrPtr_->calcReconstruction(yi[i], fi, ql[i], qr[i]);
                            ql[i] = min(real(1), max(real(0), ql[i]));
                            qr[i] = min(real(1), max(real(0), qr[i]));
                            yilt += ql[i];
                            yirt += qr[i];
                        }
                        real yiLastL = 1.0 - yilt;
                        real yiLastR = 1.0 - yirt;

                        for (integer i = 0; i < yi.size() - 1; ++i) {
                            if (yiLastL < 0) {
                                ql[i] /= yilt;
                            }
                            if (yiLastR < 0) {
                                qr[i] /= yirt;
                            }

                            if (rhoFlux_[fi] >= 0.0) {
                                flux = rhoFlux_[fi] * qr[i];
                            } else {
                                flux = rhoFlux_[fi] * ql[i];
                            }
                            //Add the numerical flux to the right-hand-term of both cells sharing the face
                            yi[i].rhs()[cl] += flux;
                            if (cr < mesh.nCells()) {
                                yi[i].rhs()[cr] -= flux;
                            }
                        }
                    }
                }
            }
        }
    }
}

void OpenHurricane::blendingSchemeTVDLimiter2::extendShockSensor() {
    const auto &fl = mesh().faces();
    for (integer fi = 0; fi < mesh().nFaces(); ++fi) {
        const auto &cl = fl[fi].leftCell();
        const auto &cr = fl[fi].rightCell();

        if (shocks_.sensor()[cl] == 1.0) {
            if (shocks_.sensor()[cr] != 1.0) {
                upwindFactor_[cr] = 0.85;
            }
        } else if (shocks_.sensor()[cr] == 1.0) {
            if (shocks_.sensor()[cl] != 1.0) {
                upwindFactor_[cl] = 0.85;
            }
        } else if (upwindFactor_[cl] == 0.85) {
            if (shocks_.sensor()[cr] != 1.0 && upwindFactor_[cr] != 0.85) {
                upwindFactor_[cr] = 0.65;
            }
        } else if (upwindFactor_[cr] == 0.85) {
            if (shocks_.sensor()[cl] != 1.0 && upwindFactor_[cl] != 0.85) {
                upwindFactor_[cl] = 0.65;
            }
        } else if (upwindFactor_[cl] == 0.65) {
            if (shocks_.sensor()[cr] != 1.0 && upwindFactor_[cr] != 0.85 &&
                upwindFactor_[cr] != 0.65) {
                upwindFactor_[cr] = 0.45;
            }
        } else if (upwindFactor_[cr] == 0.65) {
            if (shocks_.sensor()[cl] != 1.0 && upwindFactor_[cl] != 0.85 &&
                upwindFactor_[cl] != 0.65) {
                upwindFactor_[cl] = 0.45;
            }
        }
    }
}
