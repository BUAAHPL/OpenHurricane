/*!
 * \file SpalartAllmaras.cpp
 * \brief Main subroutines for the SpalartAllmaras turbulence model.
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

#include "SpalartAllmaras.hpp"

namespace OpenHurricane {
    createClassNameStr(SpalartAllmaras, "SpalartAllmaras");
    registerObjFty(turbulenceModel, SpalartAllmaras, controller);
} // namespace OpenHurricane

OpenHurricane::SpalartAllmaras::SpalartAllmaras(const controller &cont, flowModel &ev)
    : RANSModel(cont, ev),
      nut_(object("nut", mesh(), object::WRITE_RELAY_OUTPUT, object::PRIMITIVE), mesh()),
      nutLastTime_(mesh().size(), Zero), cb1_(cont.findOrDefault<real>("cb1", 0.1355)),
      cb2_(cont.findOrDefault<real>("cb2", 0.622)),
      sigma_(cont.findOrDefault<real>("sigma", 2.0 / 3.0)),
      cv1_(cont.findOrDefault<real>("cv1", 7.1)), kappa_(cont.findOrDefault<real>("kappa", 0.41)),
      cw2_(cont.findOrDefault<real>("cw2", 0.3)), cw3_(cont.findOrDefault<real>("cw3", 2.0)),
      ct3_(cont.findOrDefault<real>("ct3", 1.2)), ct4_(cont.findOrDefault<real>("ct4", 0.5)),
      cv2_(cont.findOrDefault<real>("cv2", 0.7)), cv3_(cont.findOrDefault<real>("cv3", 0.9)),
      SAVer_(SA), nutRelax_(1.0) {
    if (cont.found("SpalartAllmaras")) {
        if (cont.isController("SpalartAllmaras")) {
            const auto &SAcont = cont.subController("SpalartAllmaras");
            cb1_ = SAcont.findOrDefault<real>("cb1", 0.1355);
            cb2_ = SAcont.findOrDefault<real>("cb2", 0.622);
            sigma_ = SAcont.findOrDefault<real>("sigma", 2.0 / 3.0);
            cv1_ = SAcont.findOrDefault<real>("cv1", 7.1);
            kappa_ = SAcont.findOrDefault<real>("kappa", 0.41);
            cw2_ = SAcont.findOrDefault<real>("cw2", 0.3);
            cw3_ = SAcont.findOrDefault<real>("cw3", 2.0);
            ct3_ = SAcont.findOrDefault<real>("ct3", 1.2);
            ct4_ = SAcont.findOrDefault<real>("ct4", 0.5);
            cv2_ = SAcont.findOrDefault<real>("cv2", 0.7);
            cv3_ = SAcont.findOrDefault<real>("cv3", 0.9);
            nutRelax_ = SAcont.findOrDefault<real>("nutRelax", 1.0);
            if (SAcont.found("SAVersion")) {
                string sstw = SAcont.findWord("SAVersion");
                trim(sstw);
                stringToUpperCase(sstw);
                if (sstw == "SA") {
                    SAVer_ = SA;
                } else if (sstw == "SANEG") {
                    SAVer_ = SANEG;
                } else {
                    errorAbortStr(("Unknown SA version: " + sstw));
                }
            }
        }
    }
    cw1_ = cb1_ / sqr(kappa_) + (1.0 + cb2_) / sigma_;
    RANSModel::nEq_ = 1;

    bndValueSetting(cont.topCont());
    const auto &tcont = cont.topCont();
    if (tcont.found("flow")) {
        const auto &fcont = tcont.subController("flow");
        if (fcont.found("limits")) {
            const auto &lmtcont = fcont.subController("limits");
            if (lmtcont.found("limitMethod")) {
                const auto lmtw = lmtcont.findWord("limitMethod");
                if (lmtw == "directly") {
                    averageLimit_ = false;
                } else if (lmtw == "average") {
                    averageLimit_ = true;
                } else {
                    LFatal("Unknown limit method: %s in %s", lmtw.c_str(), lmtcont.name().c_str());
                }
            }
        }
    }
}

void OpenHurricane::SpalartAllmaras::turbParaInitialize() {
    nut_.initialize();
    nut_.updateBoundary();
    real nut0 = nut_.initialize();
    const real nu = mul()[0] / rho()[0];
    const real chi = nut0 / nu;
    const real fv1 = pow3(chi) / (pow3(chi) + pow3(cv1_));
    mut() = rho()[0] * nut0 * fv1;

    mutBoundaryUpdate();
}

void OpenHurricane::SpalartAllmaras::expSource() {
    if (isSplitting()) {
        return;
    }
    const auto &meshVol = mesh().cellVol();
    const auto nCell = mesh().nCells();
    const auto &grdV = v().grad();

    for (integer n = 0; n < nCell; ++n) {
        const real nul = mul()[n] / rho()[n];
        const real nut = nut_[n];
        nutLastTime_[n] = nut_[n];

        const real dist = wallDist()[n];
        const real vor = skewMagnitude(grdV[n]);
        const real chi = nut / nul;
        const real chi3 = pow3(chi);
        const real fv1 = chi3 / (chi3 + pow3(cv1_));
        const real fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const real mid1 = 1.0 / (sqr(dist));
        const real mid2 = mid1 / sqr(kappa_);
        const real ssline = nut * fv2 * mid2;
        real ss = vor + ssline;
        switch (SAVer_) {
        case (SA):
            ss = max(real(0.3) * vor, ss);
            break;
        case (SANEG):
            if (ssline < -cv2_ * vor) {
                real temp1 = sqr(cv2_) * vor + cv3_ * fv2 * nut * mid2;
                real temp2 = (cv3_ - 2 * cv2_) * vor - nut * mid2;
                ss = vor + vor * temp1 / temp2;
            }
            break;
        default:
            break;
        }
        /*real ss = vor + nut * fv2 * mid2;
        ss = max(ss, real(0.3) * vor); */
        const real r = min(nut / ss * mid2, real(10.0));
        const real gg = g(r);
        const real fww = fw(gg);
        const real ftt2 = ft2(chi);
        const real p = cb1_ * (real(1.0) - ftt2) * ss;
        const real d = (cw1_ * fww - cb1_ * ftt2 / sqr(kappa_)) * nut * mid1;
        real sp = p * nut * rho()[n];
        real sd = d * nut * rho()[n];
        real NC = cb2_ / sigma_ * nut_.grad()[n].magSqr() * rho()[n];
        NC -= real(1.0) / sigma_ * (nul + nut) * nut_.grad()[n] * rho().grad()[n];
        nut_.rhs()[n] += (sp - sd + NC) * meshVol[n];
    }
}

void OpenHurricane::SpalartAllmaras::impSource() {
    if (isSplitting()) {
        return;
    }
    const auto &meshVol = mesh().cellVol();
    const auto nCell = mesh().nCells();
    const auto &grdV = v().grad();
    for (integer n = 0; n < nCell; ++n) {
        const real nul = mul()[n] / rho()[n];
        const real nut = nut_[n];
        nutLastTime_[n] = nut_[n];

        const real dist = wallDist()[n];
        const real vor = skewMagnitude(grdV[n]);
        const real chi = nut / nul;
        const real chi3 = pow3(chi);
        const real fv1 = chi3 / (chi3 + pow3(cv1_));
        const real fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const real mid1 = 1.0 / (sqr(dist));
        const real mid2 = mid1 / sqr(kappa_);
        const real ssline = nut * fv2 * mid2;
        real ss = vor + ssline;
        switch (SAVer_) {
        case (SA):
            ss = max(real(0.3) * vor, ss);
            break;
        case (SANEG):
            if (ssline < -cv2_ * vor) {
                real temp1 = sqr(cv2_) * vor + cv3_ * fv2 * nut * mid2;
                real temp2 = (cv3_ - 2 * cv2_) * vor - nut * mid2;
                ss = vor + vor * temp1 / temp2;
            }
            break;
        default:
            break;
        }
        /*real ss = vor + nut * fv2 * mid2;
        ss = max(ss, real(0.3) * vor); */
        real r = min(nut / ss * mid2, real(10.0));
        real gg = g(r);
        real fww = fw(gg);
        const real ftt2 = ft2(chi);
        const real p = cb1_ * (real(1.0) - ftt2) * ss;
        const real d = (cw1_ * fww - cb1_ * ftt2 / sqr(kappa_)) * nut * mid1;
        real sp = p * nut * rho()[n];
        real sd = d * nut * rho()[n];
        real NC = cb2_ / sigma_ * nut_.grad()[n].magSqr() * rho()[n];
        NC -= real(1.0) / sigma_ * (nul + nut) * nut_.grad()[n] * rho().grad()[n];
        nut_.rhs()[n] += (sp - sd + NC) * meshVol[n];

        real dfv1 = 0.0;
        real dfv2 = 0.0;
        real dss = 0.0;
        real dp = cb1_ * (1.0 - ftt2) * ss;
        real drr = 0.0;
        real dgg = 0.0;
        real dfw = 0.0;
        real dd = (cw1_ * fww - cb1_ * ftt2 / sqr(kappa_)) * nut / sqr(dist);
        nut_.diagSource()[n] = min((dp - dd), real(0.0)) * meshVol[n];
    }
}

void OpenHurricane::SpalartAllmaras::fullImpSource(cellRealSquareMatrixArray &Jac, const integer rhoId,
                                               const integer rhoTurb0) {
    if (isSplitting()) {
        return;
    }
    const auto &meshVol = mesh().cellVol();
    const auto nCell = mesh().nCells();
    const auto &grdV = v().grad();

    for (integer n = 0; n < nCell; ++n) {
        const real nul = mul()[n] / rho()[n];
        const real nut = nut_[n];
        nutLastTime_[n] = nut_[n];

        const real dist = wallDist()[n];
        const real vor = skewMagnitude(grdV[n]);
        const real chi = nut / nul;
        const real chi3 = pow3(chi);
        const real fv1 = chi3 / (chi3 + pow3(cv1_));
        const real fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const real mid1 = 1.0 / (sqr(dist));
        const real mid2 = mid1 / sqr(kappa_);
        const real ssline = nut * fv2 * mid2;
        real ss = vor + ssline;
        switch (SAVer_) {
        case (SA):
            ss = max(real(0.3) * vor, ss);
            break;
        case (SANEG):
            if (ssline < -cv2_ * vor) {
                real temp1 = sqr(cv2_) * vor + cv3_ * fv2 * nut * mid2;
                real temp2 = (cv3_ - 2 * cv2_) * vor - nut * mid2;
                ss = vor + vor * temp1 / temp2;
            }
            break;
        default:
            break;
        }
        /*real ss = vor + nut * fv2 * mid2;
        ss = max(ss, real(0.3) * vor); */
        real r = min(nut / ss * mid2, real(10.0));
        real gg = g(r);
        real fww = fw(gg);
        const real ftt2 = ft2(chi);
        const real p = cb1_ * (real(1.0) - ftt2) * ss;
        const real d = (cw1_ * fww - cb1_ * ftt2 / sqr(kappa_)) * nut * mid1;
        real sp = p * nut * rho()[n];
        real sd = d * nut * rho()[n];
        real NC = cb2_ / sigma_ * nut_.grad()[n].magSqr() * rho()[n];
        NC -= real(1.0) / sigma_ * (nul + nut) * nut_.grad()[n] * rho().grad()[n];
        nut_.rhs()[n] += (sp - sd + NC) * meshVol[n];

        real dfv1 = 0.0;
        real dfv2 = 0.0;
        real dss = 0.0;
        real dp = cb1_ * (1.0 - ftt2) * ss;
        real drr = 0.0;
        real dgg = 0.0;
        real dfw = 0.0;
        real dd = (cw1_ * fww - cb1_ * ftt2 / sqr(kappa_)) * nut / sqr(dist);
        //nut_.diagSource()[n] = min((dp - dd), real(0.0)) * meshVol[n];
        Jac[n](rhoTurb0, rhoId) = 0;
        Jac[n](rhoTurb0, rhoTurb0) = min((dp - dd), real(0.0)) * meshVol[n];
    }
}

void OpenHurricane::SpalartAllmaras::visFlux(const faceRealArray &rhof, const faceRealArray &mulf,
                                         const faceRealArray &mutf, const cellRealArray &mul,
                                         const cellRealArray &mut, const cellRealArray &rho) {
    if (isSplitting()) {
        return;
    }
    const auto &fzl = mesh().faceZones();
    const auto &fcl = mesh().faces();
    const auto &fA = mesh().fA();

    auto gnutf = fv::gradf(nut_);
    auto nutf = fv::interpolate(nut_);

    for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
        if (fzl[fzi].isWall()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
                const integer cl = fcl[fi].leftCell();
                const integer cr = fcl[fi].rightCell();
                const real coef_turb = mulf[fi] / (sigma_);
                const real flux = coef_turb * (gnutf[fi] * fA[fi]);
                nut_.rhs()[cl] -= flux;
            }
        } else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
                const integer cl = fcl[fi].leftCell();
                const integer cr = fcl[fi].rightCell();
                const real coef_turb = (mulf[fi] + rhof[fi] * nutf[fi]) / (sigma_);
                const real flux = coef_turb * (gnutf[fi] * fA[fi]);
                nut_.rhs()[cl] -= flux;
                nut_.rhs()[cr] += flux;
            }
        }
    }
}

void OpenHurricane::SpalartAllmaras::update() {
    const integer nCells = mesh().nCells();

    for (integer n = 0; n < nCells; ++n) {
        // calculate viscosity  based on Sparlat-Allmaras Model
        const real nul = mul()[n] / rho()[n];
        const real nut = nut_[n];
        const real chi0 = nut / nul;
        const real chi03 = pow3(chi0);
        const real fv1 = chi03 / (chi03 + pow3(cv1_));
        if (nut < 0.0) {
            mut()[n] = 0.0;
        } else {
            mut()[n] = rho()[n] * nut * fv1;
        }
    }
    limit();
    mut().updateBoundary();
}

void OpenHurricane::SpalartAllmaras::mutBoundaryUpdate() {
    const faceZoneList &fZ = mesh().faceZones();
    const faceList &fL = mesh().faces();
    realTransfer myTransfer(mesh(), mut(), false, true);
    myTransfer.transferInit();
    for (integer fZI = 0; fZI < fZ.size(); ++fZI) {
        if (fZ[fZI].isInterior() || fZ[fZI].isCutFace() || fZ[fZI].isPeriodic() ||
            fZ[fZI].isPeriodicShadow()) {
            continue;
        } else {
            if (fZ[fZI].isWall()) {
                for (integer fi = fZ[fZI].firstIndex(); fi <= fZ[fZI].lastIndex(); ++fi) {
                    const integer cl = fL[fi].leftCell();
                    const integer cr = fL[fi].rightCell();

                    mut()[cr] = -mut()[cl];
                }
            } else if (fZ[fZI].bcType() == faceBCType::bcTypes::PRESSUREFARFIELD ||
                       fZ[fZI].bcType() == faceBCType::bcTypes::PRESSUREINLET ||
                       fZ[fZI].bcType() == faceBCType::bcTypes::VELOCITYINLET ||
                       fZ[fZI].bcType() == faceBCType::bcTypes::MASSFLOWINLET) {
                for (integer fi = fZ[fZI].firstIndex(); fi <= fZ[fZI].lastIndex(); ++fi) {
                    const integer cl = fL[fi].leftCell();
                    const integer cr = fL[fi].rightCell();
                    const real nul = mul()[cl] / rho()[cl];
                    const real nut = nut_[cl];
                    const real chi0 = nut / nul;
                    const real chi03 = pow3(chi0);
                    const real fv1 = chi03 / (chi03 + pow3(cv1_));
                    mut()[cr] = rho()[cr] * nut * fv1;
                }
            } else {
                for (integer fi = fZ[fZI].firstIndex(); fi <= fZ[fZI].lastIndex(); ++fi) {
                    const integer cl = fL[fi].leftCell();
                    const integer cr = fL[fi].rightCell();

                    mut()[cr] = mut()[cl];
                }
            }
        }
    }

    myTransfer.transferring();
}

void OpenHurricane::SpalartAllmaras::limit() {
    const integer nCells = mesh().nCells();

    for (integer n = 0; n < nCells; ++n) {
        mut()[n] = max(mut()[n], mutLow() * mul()[n]);
        mut()[n] = min(mut()[n], mutHigh() * mul()[n]);
    }
}

OpenHurricane::realArray OpenHurricane::SpalartAllmaras::k() const {
    return realArray();
}

OpenHurricane::realArray OpenHurricane::SpalartAllmaras::epsilon() const {
    return realArray();
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::SpalartAllmaras::Ret() const {
    return realArray();
}

OpenHurricane::cellRealArray &OpenHurricane::SpalartAllmaras::var(const integer i) {
    return nut_;
}

const OpenHurricane::cellRealArray &OpenHurricane::SpalartAllmaras::var(const integer i) const {
    return nut_;
}

void OpenHurricane::SpalartAllmaras::solving(const realArray &dt) {
    if (isCoupled()) {
        return;
    }
    const auto &meshVol = mesh().cellVol();
    const auto nCell = mesh().nCells();
    const auto nFace = mesh().nFaces();
    const auto &grdV = v().grad();

    // Source terms
    for (integer n = 0; n < nCell; ++n) {
        const real nul = mul()[n] / rho()[n];
        const real nut = nut_[n];
        const real dist = wallDist()[n];
        const real vol = meshVol[n];
        nutLastTime_[n] = nut_[n];

        const real ux = v().grad()[n].xx();
        const real uy = v().grad()[n].xy();
        const real uz = v().grad()[n].xz();

        const real vx = v().grad()[n].yx();
        const real vy = v().grad()[n].yy();
        const real vz = v().grad()[n].yz();

        const real wx = v().grad()[n].zx();
        const real wy = v().grad()[n].zy();
        const real wz = v().grad()[n].zz();
        const real w1 = wy - vz;
        const real w2 = uz - wx;
        const real w3 = vx - uy;
        const real vor = sqrt(sqr(w1) + sqr(w2) + sqr(w3));
        const real chi = nut / nul;
        const real chi3 = pow3(chi);
        const real fv1 = chi3 / (chi3 + pow3(cv1_));
        const real fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const real mid1 = 1.0 / sqr(dist);
        const real mid2 = mid1 / sqr(kappa_);
        const real ssline = nut * fv2 * mid2;
        real ss = vor + ssline;
        switch (SAVer_) {
        case (SA):
            ss = max(real(0.3) * vor, ss);
            break;
        case (SANEG):
            if (ssline < -cv2_ * vor) {
                real temp1 = sqr(cv2_) * vor + cv3_ * fv2 * nut * mid2;
                real temp2 = (cv3_ - 2 * cv2_) * vor - nut * mid2;
                ss = vor + vor * temp1 / temp2;
            }
            break;
        default:
            break;
        }
        /*real ss = vor + nut * fv2 * mid2;
        ss = max(ss, real(0.3) * vor);*/
        const real r = min(nut / ss * mid2, real(10.0));
        const real gg = g(r);
        const real fww = fw(gg);

        const real p = cb1_ * ss;
        const real d = cw1_ * fww * nut * mid1;

        const real sp = p * nut;
        const real sd = d * nut;
        nut_.rhs()[n] = (sp - sd) * vol;

        const real dfv1 = (fv1 - sqr(fv1)) * 3.0 / nut;
        const real dfv2 = (fv2 - 1.0) / nut + (sqr(1.0 - fv2)) * (fv1 / nut + dfv1);
        const real dss = (fv2 + nut * dfv2) * mid2;
        const real dp = cb1_ * dss;

        const real drr = r / nut - r * r * (fv2 / nut + dfv2);
        const real dgg = (1.0 - cw2_ + 6.0 * cw2_ * pow(r, real(5.0))) * drr;
        const real dfw =
            pow(gg, real(-7)) / (1.0 + (pow(cw3_, real(-6)))) * pow(fww, real(7.0)) * dgg;
        const real dd = cw1_ * mid1 * (fww + nut * dfw);

        nut_.diagSource()[n] = -(min((p - d), real(0.0)) + min((dp - dd), real(0.0)) * nut) * vol;
    }

    vector2DArray ra(nFace);
    const auto &fcl = mesh().faces();
    const auto &fA = mesh().fA();
    const auto &cC = mesh().cellCentre();

    // Flux
    for (integer fi = 0; fi < nFace; ++fi) {
        const integer cl = fcl[fi].leftCell();
        const integer cr = fcl[fi].rightCell();

        const real uul = v()[cl] * fA[fi];
        const real uum = 0.5 * (uul - fabs(uul));

        const real uur = v()[cr] * fA[fi];
        const real uup = 0.5 * (uur + fabs(uur));

        const real flux = nut_[cl] * uum + nut_[cr] * uup;

        // Inviscous flux
        nut_.rhs()[cl] += flux;
        nut_.rhs()[cr] -= flux;

        nut_.diagSource()[cl] -= uum;
        nut_.diagSource()[cr] += uup;

        ra[fi][0] = -uup;
        ra[fi][1] = uum;

        // Viscous flux
        vector rlr = cC[cl] - cC[cr];
        const real nnij = rlr.magnitude();

        const real null = mul()[cl] / rho()[cl];
        const real nulr = mul()[cr] / rho()[cr];

        const real nulf = 0.5 * (null + nulr);
        const real nutf = 0.5 * (nut_[cl] + nut_[cr]);
        const real nn = fA[fi].magnitude();
        const real coe1 = (nulf + (1.0 + cb2_) * nutf - cb2_ * nut_[cl]) / nnij / sigma_ * nn;
        const real coe2 = (nulf + (1.0 + cb2_) * nutf - cb2_ * nut_[cr]) / nnij / sigma_ * nn;

        const real fluxl = coe1 * (nut_[cl] - nut_[cr]);
        const real fluxr = coe2 * (nut_[cl] - nut_[cr]);

        nut_.rhs()[cl] -= fluxl;
        nut_.rhs()[cr] += fluxr;

        nut_.diagSource()[cl] += coe1;
        nut_.diagSource()[cr] += coe2;

        ra[fi][0] -= coe1;
        ra[fi][1] -= coe2;
    }
    realArray dq(mesh().nTotalCells(), Zero);
    for (integer i = 0; i < dq.size(); ++i) {
        dq[i] = 0.0;
    }
    const auto &cells = mesh().cells();
    const auto &vol = mesh().cellVolume();
    // LU-SGS lower loop
    LUSGS::lowerLoop(nut_, ra, dt, dq);

    // LU-SGS upper loop
    LUSGS::upperLoop(nut_, ra, dt, dq);

    // LU-SGS update
    for (integer n = 0; n < nCell; ++n) {
        nut_[n] = nut_[n] + dq[n] * nutRelax_;
    }
}

void OpenHurricane::SpalartAllmaras::updateBoundary() {
    nut_.updateBoundary();
}

void OpenHurricane::SpalartAllmaras::limitAndUpdateBoundary() {
    if (mesh().Iteration().isBeginLimit()) {
        limitTurbVarIncrease(nut_, nutLastTime_);
    }
    limitNegative(nut_);
    nut_.updateBoundary();
}

OpenHurricane::symmTensorArray OpenHurricane::SpalartAllmaras::ReynoldsStressTensor() const {
    symmTensorArray tt(
        -mut() / rho() *
        (twoSymm(v().grad()) - (real(2.0 / 3.0) * (div(diagToVector(v().grad()))) * I)));

    return symmTensorArray();
}

OpenHurricane::faceSymmTensorArray
OpenHurricane::SpalartAllmaras::tauEff(const faceRealArray &rhof, const faceRealArray &muf,
                                   const faceRealArray &mutf,
                                   const faceTensorArray &deltafV) const {

    faceSymmTensorArray tau(
        object("tau", mesh(), object::NOT_WRITE, object::TEMPORARY), mesh(),
        (muf + mutf) * (twoSymm(deltafV) - (real(2.0 / 3.0) * (div(diagToVector(deltafV))) * I)));
    return tau;
}

void OpenHurricane::SpalartAllmaras::calcGrad(const spatialScheme &sps) {
    if (isCoupled()) {
        return;
    }
    sps.grad(nut_);
}
