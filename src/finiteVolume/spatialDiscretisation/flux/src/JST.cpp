/*!
 * \file JST.cpp
 * \brief Main subroutines for JST central scheme..
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

#include "JST.hpp"
#include <cmath>

namespace OpenHurricane {
    createClassNameStr(JST, "JST");
    registerObjFty(spatialScheme, JST, controller);
} // namespace OpenHurricane

void OpenHurricane::JST::calcFlux(const real rhol, const real rhor, const vector &VL,
                                  const vector &VR, const real pl, const real pr, const real gl,
                                  const real gr, const real el, const real er, const real cl,
                                  const real cr, const vector &faceArea, const real blend,
                                  realArray &flux) const {
    LFatal("Not implemented");
}

void OpenHurricane::JST::spectralRadius() {
    const auto &mesh = v_.mesh();
    const auto &fl = mesh.faces();
    const auto &fA = mesh.faceArea();
    const auto &fW = mesh.faceWeight();
    const auto nCells = mesh.nCells();
    const auto &v = v_;
    const auto &rho = thermo_.rho();
    const auto &p = thermo_.p();
    const auto &gama = thermo_.gamma();
    cellRealArray Ai(object("Ai_JST", mesh, object::NOT_WRITE, object::TEMPORARY), mesh, Zero);

    for (integer fi = 0; fi < mesh.nFaces(); ++fi) {
        const auto &cl = fl[fi].leftCell();
        const auto &cr = fl[fi].rightCell();
        const real wl = fW[fi];
        const real wr = 1 - wl;

        const real rhof = wl * rho[cl] + wr * rho[cr];
        const real pf = wl * p[cl] + wr * p[cr];
        const real gamf = wl * gama[cl] + wr * gama[cr];
        const vector vf = wl * v[cl] + wr * v[cr];

        const auto n = fA[fi].normalized();
        const auto Sf = fA[fi].magnitude();
        const real Vn = mag(vf * n);

        const real cf = sqrt(gamf * pf / rhof);
        const real Af = (Vn + cf) * Sf;
        Ai[cl] += Af;
        if (cr < nCells) {
            Ai[cr] += Af;
        }
    }

    // Update boundary
    realTransfer myTransfer(mesh, Ai, false, true);
    myTransfer.transferInit();
    const auto &fzl = mesh.faceZones();
    for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
        if (fzl[fzi].isInterior() || fzl[fzi].isCutFace() || fzl[fzi].isPeriodic() ||
            fzl[fzi].isPeriodicShadow()) {
            continue;
        } else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
                const auto cl = fl(fi).leftCell();
                const auto cr = fl(fi).rightCell();
                Ai[cr] = Ai[cl];
            }
        }
    }
    //fv::transfer(Ai, true);
    myTransfer.transferring();

    for (integer fi = 0; fi < mesh.nFaces(); ++fi) {
        const auto &cl = fl[fi].leftCell();
        const auto &cr = fl[fi].rightCell();
        const real wl = fW[fi];
        const real wr = 1 - wl;
        AIJ_[fi] = wl * Ai[cl] + wr * Ai[cr];
    }
}

void OpenHurricane::JST::calcModifyFactor() {
    if (!useModifyFactor_) {
        return;
    }
    auto &mdf = *modifyFactorPtr_;
    auto &v = v_;

    for (integer n = 0; n < mesh().nCells(); ++n) {
        const real Du = tr(v.grad()[n]);

        const real omega2 = skewMagSqr(v.grad()[n]);

        mdf[n] = sqr(max(real(0), -Du)) / (sqr(Du) + omega2 + tiny);
    }

    const auto &fZ = mesh().faceZones();
    const auto &fL = mesh().faces();

    realTransfer myTransfer(mesh(), mdf, false, true);
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
                mdf[cr] = mdf[cl];
            }
        }
    }
    myTransfer.transferring();
}

void OpenHurricane::JST::getGeometricalWeights() {
    const auto &mesh = v_.mesh();
    const auto &fl = mesh.faces();
    const auto &cells = mesh.cells();
    const auto &fA = mesh.faceArea();
    //const auto& fC = mesh.faceCntr();
    //const auto& fW = mesh.faceWeight();
    const auto &cC = mesh.cellCentre();

    const auto fsk = meshCheckers::faceSkewness(mesh_, mesh_.points(), mesh_.faceCentre(),
                                                mesh_.faceArea(), mesh_.cellCentre());

    for (integer i = 0; i < mesh.nCells(); ++i) {
        const auto &ci = cC[i];
        tensor I = Zero;
        vector Rs = Zero;

        for (integer j = 0; j < cells[i].faceSize(); ++j) {
            const auto fi = cells[i].facei(j);
            const auto cl = fl[fi].leftCell();
            const auto cr = fl[fi].rightCell();
            const integer jj = cl + cr - i;

            vector tempV = (cC[jj] - ci);
            I += inv(tempV.magSqr()) * (tempV & tempV);

            Rs += tempV;
        }
        vector lamb = Rs / I;

        for (integer j = 0; j < cells[i].faceSize(); ++j) {
            const auto fi = cells[i].facei(j);
            const auto cl = fl[fi].leftCell();
            const auto cr = fl[fi].rightCell();
            const integer jj = cl + cr - i;
            vector tempV = (cC[jj] - ci);
            real ski = max(real(1.0 - fsk[fi]), real(0.01));
            real thet = real(1.0) + inv(tempV.magSqr()) * lamb * tempV;
            thet = min(real(2), max(thet, real(0)));
            if (ski == 0.01) {
                thet *= ski;
            }
            if (i == cl) {
                thetaIJ_[fi][0] = thet;
            } else {
                thetaIJ_[fi][1] = thet;
            }
        }
    }
}

void OpenHurricane::JST::pressureSensor() {
    const auto &p = thermo_.p();
    const integer nCells = mesh().nCells();
    const auto &faces = mesh().faces();
    const auto &cells = mesh().cells();
    for (integer celli = 0; celli < nCells; ++celli) {
        // numerator
        real num = Zero;

        // denominator
        real den = Zero;

        for (integer i = 0; i < cells[celli].faceSize(); ++i) {
            const integer fi = cells[celli].facei(i);
            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const integer I = celli;
            const integer J = (cl == celli) ? cr : cl;
            if (cl == celli) {
                num += (thetaIJ_[fi][0] * (p[J] - p[I]));
            } else {
                num += (thetaIJ_[fi][1] * (p[J] - p[I]));
            }
            den += (p[J] + p[I]);
        }
        pgi_[celli] = mag(num) / max(den, tiny);
    }

    const auto &fZ = mesh().faceZones();
    const auto &fL = mesh().faces();

    realTransfer myTransfer(mesh(), pgi_, false, true);
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
                pgi_[cr] = pgi_[cl];
            }
        }
    }
    myTransfer.transferring();
}

void OpenHurricane::JST::calcEsp() {
    const auto &faces = mesh().faces();
    const auto nFaces = mesh().nFaces();
    pressureSensor();
    if (useModifyFactor_) {
        calcModifyFactor();
    }
    for (integer fi = 0; fi < nFaces; ++fi) {
        const auto I = faces[fi].leftCell();
        const auto J = faces[fi].rightCell();
        if (useModifyFactor_) {
            espIJ2_[fi] =
                k2_ * max((*modifyFactorPtr_)[I] * pgi_[I], (*modifyFactorPtr_)[J] * pgi_[J]);
        } else {
            espIJ2_[fi] = k2_ * max(pgi_[I], pgi_[J]);
        }

        espIJ4_[fi] = max(real(0), (k4_ - espIJ2_[fi]));
        //espIJ4_[fi] = 0;
    }
}

OpenHurricane::real OpenHurricane::JST::artificialDissipation(const integer fi, const real QJI,
                                                              const real WRL) const {
    const real Ae2 = AIJ_[fi] * espIJ2_[fi];
    const real Ae4 = 4.0 * AIJ_[fi] * espIJ4_[fi];
    const real DIJ2 = Ae2 * QJI;
    const real DIJ4 = Ae4 * WRL;
    return (DIJ2 + DIJ4);
}

OpenHurricane::vector OpenHurricane::JST::artificialDissipation(const integer fi, const vector &QJI,
                                                                const vector &WRL) const {
    const real Ae2 = AIJ_[fi] * espIJ2_[fi];
    const real Ae4 = 4.0 * AIJ_[fi] * espIJ4_[fi];
    const vector DIJ2 = Ae2 * QJI;
    const vector DIJ4 = Ae4 * WRL;
    return (DIJ2 + DIJ4);
}

void OpenHurricane::JST::addArtificialDissipation(const integer fi, const real QJI, const real WRL,
                                                  const integer cl, const integer cr,
                                                  cellRealArray &rhs) const {
    const real Ae2 = AIJ_[fi] * espIJ2_[fi];
    const real Ae4 = 4.0 * AIJ_[fi] * espIJ4_[fi];
    const real DIJ2 = Ae2 * QJI;
    const real DIJ4 = Ae4 * WRL;
    rhs[cl] += (DIJ2 * thetaIJ_[fi][0] + DIJ4);
    if (cr < mesh_.nCells()) {
        rhs[cr] -= (DIJ2 * thetaIJ_[fi][1] + DIJ4);
    }
}

void OpenHurricane::JST::addArtificialDissipation(const integer fi, const vector QJI,
                                                  const vector WRL, const integer cl,
                                                  const integer cr, cellVectorArray &rhs) const {
    const real Ae2 = AIJ_[fi] * espIJ2_[fi];
    const real Ae4 = 4.0 * AIJ_[fi] * espIJ4_[fi];
    const vector DIJ2 = Ae2 * QJI;
    const vector DIJ4 = Ae4 * WRL;
    rhs[cl] += (DIJ2 * thetaIJ_[fi][0] + DIJ4);
    if (cr < mesh_.nCells()) {
        rhs[cr] -= (DIJ2 * thetaIJ_[fi][1] + DIJ4);
    }
}

OpenHurricane::JST::JST(const controller &cont, const runtimeMesh &mesh, flowModel &flow)
    : spatialScheme(cont, mesh, flow), AIJ_(mesh.nFaces(), Zero), thetaIJ_(mesh.nFaces(), Zero),
      espIJ2_(mesh.nFaces(), Zero), espIJ4_(mesh.nFaces(), Zero), rhoLR_(mesh.nFaces(), Zero),
      pgi_(object("pressureSensorI", mesh, object::NOT_WRITE), mesh), useModifyFactor_(false),
      modifyFactorPtr_(nullptr), shocks_(flow, cont.subController("JST")),
      vnf_(mesh.nFaces(), Zero), k2_(cont.subController("JST").findOrDefault<real>("k2", 0.5)),
      k4_(cont.subController("JST").findOrDefault<real>("k4", 1.0 / 128.0)) {
    controllerSwitch myContS(cont.subController("JST"));
    useModifyFactor_ = myContS("useModifyFactor", useModifyFactor_);

    if (useModifyFactor_) {
        modifyFactorPtr_ = new cellRealArray(object("pSModFct", mesh, object::NOT_WRITE), mesh);
    }
    getGeometricalWeights();
}

void OpenHurricane::JST::basicFlux() {
    const auto &mesh = v_.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();

    const auto &fA = mesh.faceArea();
    //const auto& fC = mesh.faceCntr();
    const auto &fW = mesh.faceWeight();
    //const auto& cC = mesh.cellCentre();

    auto &v = v_;
    auto &rho = const_cast<rhoThermo &>(thermo_).rho();
    auto &p = const_cast<rhoThermo &>(thermo_).p();
    auto &T = const_cast<rhoThermo &>(thermo_).T();
    auto &E = const_cast<rhoThermo &>(thermo_).E();
    auto &gama = const_cast<rhoThermo &>(thermo_).gamma();

    rho.rhs() = Zero;
    v.rhs() = Zero;
    //T.rhs() = Zero;
    E.rhs() = Zero;

    shocks_.calcShockSensor();
    spectralRadius();
    calcEsp();
    reconstrPtr_->calcGrad(rho);
    reconstrPtr_->calcGrad(v);
    reconstrPtr_->calcGrad(p);
    reconstrPtr_->calcGrad(E);
    reconstrPtr_->calcLimiter(rho);
    reconstrPtr_->calcLimiter(v);
    reconstrPtr_->calcLimiter(p);
    reconstrPtr_->calcLimiter(E);

    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const integer &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();

                const vector vFlux = p[cl] * fA[fi];
                v.rhs()[cl] += vFlux;
                rhoFlux_[fi] = Zero;
                real rhol, rhor;
                reconstrPtr_->calcReconstruction(rho, fi, rhol, rhor);
                real pl, pr;
                vector vl, vr;
                real el, er;
                reconstrPtr_->calcReconstruction(p, fi, pl, pr);
                if (rhol <= 0.0 || rhor <= 0.0 || pl <= 0.0 || pr <= 0.0) {
                    pl = p[cl];
                    pr = p[cr];
                    rhol = rho[cl];
                    rhor = rho[cr];
                    vl = v[cl];
                    vr = v[cr];
                    el = E[cl];
                    er = E[cr];
                } else {
                    reconstrPtr_->calcReconstruction(v, fi, vl, vr);
                    reconstrPtr_->calcReconstruction(E, fi, el, er);
                }
                rhoLR_[fi][0] = rhol;
                rhoLR_[fi][1] = rhor;

                addArtificialDissipation(fi, rho[cr] - rho[cl], rhor - rhol, cl, cr, rho.rhs());
                addArtificialDissipation(fi, rho[cr] * v[cr] - rho[cl] * v[cl],
                                         rhor * vr - rhol * vl, cl, cr, v.rhs());

                addArtificialDissipation(fi, rho[cr] * E[cr] - rho[cl] * E[cl],
                                         rhor * er - rhol * el, cl, cr, E.rhs());
            }
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();

                real rhol, rhor;
                reconstrPtr_->calcReconstruction(rho, fi, rhol, rhor);
                real pl, pr;
                vector vl, vr;
                real el, er;
                reconstrPtr_->calcReconstruction(p, fi, pl, pr);
                if (rhol <= 0.0 || rhor <= 0.0 || pl <= 0.0 || pr <= 0.0) {
                    pl = p[cl];
                    pr = p[cr];
                    rhol = rho[cl];
                    rhor = rho[cr];
                    vl = v[cl];
                    vr = v[cr];
                    el = E[cl];
                    er = E[cr];
                } else {
                    reconstrPtr_->calcReconstruction(v, fi, vl, vr);
                    reconstrPtr_->calcReconstruction(E, fi, el, er);
                }
                rhoLR_[fi][0] = rhol;
                rhoLR_[fi][1] = rhor;

                real rholl = rhol;
                real rhorr = rhor;
                real pll = pl;
                real prr = pr;
                real ell = el;
                real err = er;
                vector vll = vl;
                vector vrr = vr;

                if (shocks_.sensor()[cl] == 0 && shocks_.sensor()[cr] == 0) {
                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];

                    rholl = rho[cl] + rho.grad()[cl] * rL;
                    rhorr = rho[cr] + rho.grad()[cr] * rR;

                    vll = v[cl] + v.grad()[cl] * rL;
                    vrr = v[cr] + v.grad()[cr] * rR;

                    pll = p[cl] + p.grad()[cl] * rL;
                    prr = p[cr] + p.grad()[cr] * rR;

                    ell = E[cl] + E.grad()[cl] * rL;
                    err = E[cr] + E.grad()[cr] * rR;
                }

                const real rhof = 0.5 * (rholl + rhorr);
                const vector rhovf = real(0.5) * (rholl * vll + rhorr * vrr);
                const vector vf = 0.5 * (vll + vrr);
                const real pf = 0.5 * (pll + prr);
                const real rhoef = 0.5 * (rholl * ell + rhorr * err);

                const auto n = fA[fi].normalized();
                const auto Sf = fA[fi].magnitude();
                const real Vn = vf * n;
                vnf_[fi] = Vn;
                flux_[0] = rhof * Vn * Sf;
                flux_[1] = (rhovf.x() * Vn + pf * n.x()) * Sf;
                flux_[2] = (rhovf.y() * Vn + pf * n.y()) * Sf;
                flux_[3] = (rhovf.z() * Vn + pf * n.z()) * Sf;
                flux_[4] = (rhoef + pf) * Vn * Sf;

                rhoFlux_[fi] = flux_[0];

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
                addArtificialDissipation(fi, rho[cr] - rho[cl], rhor - rhol, cl, cr, rho.rhs());
                addArtificialDissipation(fi, rho[cr] * v[cr] - rho[cl] * v[cl],
                                         rhor * vr - rhol * vl, cl, cr, v.rhs());

                addArtificialDissipation(fi, rho[cr] * E[cr] - rho[cl] * E[cl],
                                         rhor * er - rhol * el, cl, cr, E.rhs());
            }
        }
    }

    reconstrPtr_->calcGrad(T);
    reconstrPtr_->calcGrad(gama);
}

void OpenHurricane::JST::invFlux(cellRealArray &cellQ) const {
    const auto &mesh = v_.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();

    const auto &fA = mesh.faceArea();
    //const auto& fC = mesh.faceCntr();
    const auto &fW = mesh.faceWeight();
    //const auto& cC = mesh.cellCentre();

    const auto &v = v_;
    const auto &rho = thermo_.rho();

    cellQ.rhs() = Zero;

    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);

    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                real ql, qr;
                reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                addArtificialDissipation(fi, rho[cr] * cellQ[cr] - rho[cl] * cellQ[cl],
                                         rhoLR_[fi][1] * qr - rhoLR_[fi][0] * ql, cl, cr,
                                         cellQ.rhs());
            }
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();

                real ql, qr;
                reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);

                real qll = ql;
                real qrr = qr;
                real rholl = rhoLR_[fi][0];
                real rhorr = rhoLR_[fi][1];
                if (shocks_.sensor()[cl] == 0 && shocks_.sensor()[cr] == 0) {
                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];

                    rholl = rho[cl] + rho.grad()[cl] * rL;
                    rhorr = rho[cr] + rho.grad()[cr] * rR;

                    qll = cellQ[cl] + cellQ.grad()[cl] * rL;
                    qrr = cellQ[cr] + cellQ.grad()[cr] * rR;
                }
                const real rhoQf = 0.5 * (rholl * qll + rhorr * qrr);

                const auto Sf = fA[fi].magnitude();
                const real fluxQ = rhoQf * vnf_[fi] * Sf;

                cellQ.rhs()[cl] += fluxQ;

                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= fluxQ;
                }

                addArtificialDissipation(fi, rho[cr] * cellQ[cr] - rho[cl] * cellQ[cl],
                                         rhoLR_[fi][1] * qr - rhoLR_[fi][0] * ql, cl, cr,
                                         cellQ.rhs());
            }
        }
    }
}

void OpenHurricane::JST::invFlux(cellVectorArray &cellQ) const {
    const auto &mesh = v_.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();

    const auto &fA = mesh.faceArea();
    //const auto& fC = mesh.faceCntr();
    const auto &fW = mesh.faceWeight();
    //const auto& cC = mesh.cellCentre();

    const auto &v = v_;
    const auto &rho = thermo_.rho();

    cellQ.rhs() = Zero;
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                vector ql, qr;
                reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                addArtificialDissipation(fi, rho[cr] * cellQ[cr] - rho[cl] * cellQ[cl],
                                         rhoLR_[fi][1] * qr - rhoLR_[fi][0] * ql, cl, cr,
                                         cellQ.rhs());
            }
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                vector ql, qr;
                reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);

                vector qll = ql;
                vector qrr = qr;
                real rholl = rhoLR_[fi][0];
                real rhorr = rhoLR_[fi][1];
                if (shocks_.sensor()[cl] == 0 && shocks_.sensor()[cr] == 0) {
                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];

                    rholl = rho[cl] + rho.grad()[cl] * rL;
                    rhorr = rho[cr] + rho.grad()[cr] * rR;

                    qll = cellQ[cl] + cellQ.grad()[cl] * rL;
                    qrr = cellQ[cr] + cellQ.grad()[cr] * rR;
                }
                const vector rhoQf = 0.5 * (rholl * qll + rhorr * qrr);

                const auto Sf = fA[fi].magnitude();
                const vector fluxQ = rhoQf * vnf_[fi] * Sf;

                cellQ.rhs()[cl] += fluxQ;

                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= fluxQ;
                }

                addArtificialDissipation(fi, rho[cr] * cellQ[cr] - rho[cl] * cellQ[cl],
                                         rhoLR_[fi][1] * qr - rhoLR_[fi][0] * ql, cl, cr,
                                         cellQ.rhs());
            }
        }
    }
}

void OpenHurricane::JST::invFlux(cellRealArray &cellQ, const realBounded &bound) const {
    const auto &mesh = v_.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();

    const auto &fA = mesh.faceArea();
    //const auto& fC = mesh.faceCntr();
    const auto &fW = mesh.faceWeight();
    //const auto& cC = mesh.cellCentre();

    const auto &v = v_;
    const auto &rho = thermo_.rho();

    cellQ.rhs() = Zero;
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);

    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                real ql, qr;
                reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                ql = bound.bounding(ql, cellQ[cl]);
                qr = bound.bounding(qr, cellQ[cr]);
                addArtificialDissipation(fi, rho[cr] * cellQ[cr] - rho[cl] * cellQ[cl],
                                         rhoLR_[fi][1] * qr - rhoLR_[fi][0] * ql, cl, cr,
                                         cellQ.rhs());
            }
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                real ql, qr;
                reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                ql = bound.bounding(ql, cellQ[cl]);
                qr = bound.bounding(qr, cellQ[cr]);

                real qll = ql;
                real qrr = qr;
                real rholl = rhoLR_[fi][0];
                real rhorr = rhoLR_[fi][1];
                if (shocks_.sensor()[cl] == 0 && shocks_.sensor()[cr] == 0) {
                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];

                    rholl = rho[cl] + rho.grad()[cl] * rL;
                    rhorr = rho[cr] + rho.grad()[cr] * rR;

                    qll = cellQ[cl] + cellQ.grad()[cl] * rL;
                    qrr = cellQ[cr] + cellQ.grad()[cr] * rR;
                    qll = bound.bounding(qll, cellQ[cl]);
                    qrr = bound.bounding(qrr, cellQ[cr]);
                }
                const real rhoQf = 0.5 * (rholl * qll + rhorr * qrr);

                const auto Sf = fA[fi].magnitude();
                const real fluxQ = rhoQf * vnf_[fi] * Sf;

                cellQ.rhs()[cl] += fluxQ;

                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= fluxQ;
                }

                addArtificialDissipation(fi, rho[cr] * cellQ[cr] - rho[cl] * cellQ[cl],
                                         rhoLR_[fi][1] * qr - rhoLR_[fi][0] * ql, cl, cr,
                                         cellQ.rhs());
            }
        }
    }
}

void OpenHurricane::JST::invFlux(cellVectorArray &cellQ, const vectorBounded &bound) const {
    const auto &mesh = v_.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();

    const auto &fA = mesh.faceArea();
    //const auto& fC = mesh.faceCntr();
    const auto &fW = mesh.faceWeight();
    //const auto& cC = mesh.cellCentre();

    const auto &v = v_;
    const auto &rho = thermo_.rho();

    cellQ.rhs() = Zero;
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                vector ql, qr;
                reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                ql = bound.bounding(ql, cellQ[cl]);
                qr = bound.bounding(qr, cellQ[cr]);
                addArtificialDissipation(fi, rho[cr] * cellQ[cr] - rho[cl] * cellQ[cl],
                                         rhoLR_[fi][1] * qr - rhoLR_[fi][0] * ql, cl, cr,
                                         cellQ.rhs());
            }
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                vector ql, qr;
                reconstrPtr_->calcReconstruction(cellQ, fi, ql, qr);
                ql = bound.bounding(ql, cellQ[cl]);
                qr = bound.bounding(qr, cellQ[cr]);

                vector qll = ql;
                vector qrr = qr;
                real rholl = rhoLR_[fi][0];
                real rhorr = rhoLR_[fi][1];
                if (shocks_.sensor()[cl] == 0 && shocks_.sensor()[cr] == 0) {
                    const vector &rL = faceLeftCellCtr[fi];
                    const vector &rR = faceRightCellCtr[fi];

                    rholl = rho[cl] + rho.grad()[cl] * rL;
                    rhorr = rho[cr] + rho.grad()[cr] * rR;

                    qll = cellQ[cl] + cellQ.grad()[cl] * rL;
                    qrr = cellQ[cr] + cellQ.grad()[cr] * rR;
                    qll = bound.bounding(qll, cellQ[cl]);
                    qrr = bound.bounding(qrr, cellQ[cr]);
                }
                const vector rhoQf = 0.5 * (rholl * qll + rhorr * qrr);

                const auto Sf = fA[fi].magnitude();
                const vector fluxQ = rhoQf * vnf_[fi] * Sf;

                cellQ.rhs()[cl] += fluxQ;

                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= fluxQ;
                }

                addArtificialDissipation(fi, rho[cr] * cellQ[cr] - rho[cl] * cellQ[cl],
                                         rhoLR_[fi][1] * qr - rhoLR_[fi][0] * ql, cl, cr,
                                         cellQ.rhs());
            }
        }
    }
}

void OpenHurricane::JST::invFluxSpecies(PtrList<cellRealArray> &yi, const bool withLastSpc) const {
    if (withLastSpc) {
        for (integer i = 0; i < yi.size(); ++i) {
            this->invFlux(yi[i], realBounded::lowerUpperBound(real(0), real(1),
                                                              realBounded::BOUNDED_WITH_VALUE));
        }
    } else {
        for (integer i = 0; i < yi.size() - 1; ++i) {
            this->invFlux(yi[i], realBounded::lowerUpperBound(real(0), real(1),
                                                              realBounded::BOUNDED_WITH_VALUE));
        }
    }
}
