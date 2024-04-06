/*!
 * \file SST.cpp
 * \brief Main subroutines for the SST turbulence model.
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

#include "SST.hpp"



 namespace OpenHurricane{
	createClassNameStr(SST,"SST");
}
namespace OpenHurricane {
     registerObjFty(turbulenceModel,SST,controller);
}

OpenHurricane::real OpenHurricane::SST::F1(const integer celli, const real cdsst) const {
    const auto dist = wallDist()[celli];

    const real CDkw = CDkOmega(cdsst);

    const real arg11 = sqrt(k_[celli]) / (betas_ * w_[celli] * dist);
    const real dist2 = sqr(dist);
    const real arg12 = 500.0 * mul()[celli] / (rho()[celli] * dist2 * w_[celli]);
    const real arg13 = 4.0 * rho()[celli] * sigmaw2_ * k_[celli] / (CDkw * dist2);
    const real arg1 = min(max(arg11, arg12), arg13);
    return tanh(pow4(min(arg1, real(100))));
}

OpenHurricane::real OpenHurricane::SST::F2(const integer celli) const {
    const auto dist = wallDist()[celli];

    const real arg21 = 2.0 * sqrt(k_[celli]) / (betas_ * w_[celli] * dist);
    const real arg22 = 500.0 * mul()[celli] / (rho()[celli] * sqr(dist) * w_[celli]);
    const real arg2 = max(arg21, arg22);
    return tanh(sqr(min(arg2, real(1000))));
}

OpenHurricane::real OpenHurricane::SST::yPlusLam(const real kappa, const real E) {
    real yl = 11.0;
    for (uinteger i = 0; i < 10; ++i) {
        yl = log(max(E * yl, real(1))) / kappa;
    }
    return yl;
}

OpenHurricane::real OpenHurricane::SST::CDkOmega(const real cdsst) const {
    switch (SSTVer_) {
    case (standardSST):
    case (SSTSUST):
    case (SSTV):
        return max(cdsst, real(1.0e-20));
        break;
    case (SST2003):
        return max(cdsst, real(1.0e-10));
        break;
    default:
        return max(cdsst, real(1.0e-20));
        break;
    }
}

OpenHurricane::real OpenHurricane::SST::minP(const real P, const real brwk) const {
    switch (SSTVer_) {
    case (standardSST):
    case (SSTSUST):
    case (SSTV):
        return min(P, real(20.0) * brwk);
        break;
    case (SST2003):
        return min(P, real(10.0) * brwk);
        break;
    default:
        return min(P, real(20.0) * brwk);
        break;
    }
}

OpenHurricane::real OpenHurricane::SST::dpkdrk(const real P, const real brwk, const real rhok) const {
    switch (SSTVer_) {
    case (standardSST):
    case (SSTSUST):
    case (SSTV):
        if (P > 20.0 * brwk) {
            return 20.0 * brwk / rhok;
        } else {
            return P / rhok;
        }
        break;
    case (SST2003):
        if (P > 10.0 * brwk) {
            return 10.0 * brwk / rhok;
        } else {
            return P / rhok;
        }
        break;
    default:
        return P / rhok;
        break;
    }
}

OpenHurricane::real OpenHurricane::SST::dpkdrk(const real pk, const real brwk, const real rhok,
                                       const real dpkdrkt) const {
    switch (SSTVer_) {
    case (standardSST):
    case (SSTSUST):
    case (SSTV):
        if (pk > 20.0 * brwk) {
            return 20.0 * brwk / rhok;
        } else {
            return dpkdrkt;
        }
        break;
    case (SST2003):
        if (pk > 10.0 * brwk) {
            return 10.0 * brwk / rhok;
        } else {
            return dpkdrkt;
        }
        break;
    default:
        return dpkdrkt;
        break;
    }
}

OpenHurricane::real OpenHurricane::SST::dpkdrw(const real P, const real brwk, const real rhow) const {
    switch (SSTVer_) {
    case (standardSST):
    case (SSTSUST):
    case (SSTV):
        /*if (P > 20.0 * brwk)
        {
            return 20.0 * brwk / rhow;
        }
        else
        {
            return  0.0;
        }*/
        return 0.0;
        break;
    case (SST2003):
        /* if (P > 10.0 * brwk)
         {
             return 10.0 * brwk / rhow;
         }
         else
         {
             return 0.0;
         }*/
        return 0.0;
        break;
    default:
        return 0.0;
        break;
    }
}

OpenHurricane::SST::SST(const controller &cont, flowModel &ev)
    : RANSModel(cont, ev),
      k_(object("kt", mesh(), object::WRITE_RELAY_OUTPUT, object::PRIMITIVE), mesh()),
      w_(object("wt", mesh(), object::WRITE_RELAY_OUTPUT, object::PRIMITIVE), mesh()),
      kLastTime_(mesh().size(), Zero), wLastTime_(mesh().size(), Zero),
      F1_(object("F1", mesh(), object::WRITE_OUTPUT), mesh(), Zero),
      kSource(object("kSource", mesh(), object::WRITE_OUTPUT), mesh(), Zero),
      wSource(object("wSource", mesh(), object::WRITE_OUTPUT), mesh(), Zero),
      kInvFlux(object("kInvFlux", mesh(), object::WRITE_OUTPUT), mesh(), Zero),
      wInvFlux(object("wInvFlux", mesh(), object::WRITE_OUTPUT), mesh(), Zero),
      kVisFluxTmp(object("kVisFluxTmp", mesh(), object::WRITE_OUTPUT), mesh(), Zero),
      wVisFluxTmp(object("wVisFluxTmp", mesh(), object::WRITE_OUTPUT), mesh(), Zero),
      F2_(object("F2", mesh(), object::WRITE_OUTPUT), mesh(), Zero),
      SS_(object("vorticityInSST", mesh(), object::WRITE_OUTPUT), mesh(), Zero),
      a1_(cont.findOrDefault<real>("a1", 0.31)), a2_(cont.findOrDefault<real>("a2", 0.15)),
      a3_(cont.findOrDefault<real>("a3", 0.2)), betas_(cont.findOrDefault<real>("betas", 0.09)),
      kappa_(cont.findOrDefault<real>("kappa", 0.41)),
      sigmak1_(cont.findOrDefault<real>("sigmak1", 0.85)),
      sigmaw1_(cont.findOrDefault<real>("sigmaw1", 0.5)),
      beta1_(cont.findOrDefault<real>("beta1", 0.075)),
      sigmak2_(cont.findOrDefault<real>("sigmak2", 1.0)),
      sigmaw2_(cont.findOrDefault<real>("sigmaw2", 0.856)),
      beta2_(cont.findOrDefault<real>("beta2", 0.0828)), gam1_(0.0), gam2_(0.0), Cmu_(0.09),
      E_(9.8), kamb_(0.0), wamb_(0.0), SSTVer_(standardSST), kRelax_(1.0), wRelax_(1.0), rak_(),
      raw_(), kVisFlux_(mesh().nFaces(), Zero), Pk_(mesh().nCells(), Zero), dqk_(), dqw_(),
      yplusLam_(11.0), omegaWallFunction_(false), omegaWallFunctionFaceZoneList_(), minK_(1e-14),
      minW_(1e-20) {
    if (cont.found("SST")) {
        if (cont.isController("SST")) {
            const auto &SSTcont = cont.subController("SST");
            a1_ = SSTcont.findOrDefault<real>("a1", 0.31);
            a2_ = SSTcont.findOrDefault<real>("a2", 0.15);
            a3_ = SSTcont.findOrDefault<real>("a3", 0.2);
            betas_ = SSTcont.findOrDefault<real>("betas", 0.09);
            kappa_ = SSTcont.findOrDefault<real>("kappa", 0.41);
            sigmak1_ = SSTcont.findOrDefault<real>("sigmak1", 0.85);
            sigmaw1_ = SSTcont.findOrDefault<real>("sigmaw1", 0.5);
            beta1_ = SSTcont.findOrDefault<real>("beta1", 0.075);
            sigmak2_ = SSTcont.findOrDefault<real>("sigmak2", 1.0);
            sigmaw2_ = SSTcont.findOrDefault<real>("sigmaw2", 0.856);
            beta2_ = SSTcont.findOrDefault<real>("beta2", 0.0828);
            kRelax_ = SSTcont.findOrDefault<real>("kRelax", 1.0);
            wRelax_ = SSTcont.findOrDefault<real>("wRelax", 1.0);
            if (SSTcont.found("SSTVersion")) {
                string sstw = SSTcont.findWord("SSTVersion");
                trim(sstw);
                stringToUpperCase(sstw);
                if (sstw == "STANDARDSST") {
                    SSTVer_ = standardSST;
                } else if (sstw == "SST-2003") {
                    SSTVer_ = SST2003;
                } else if (sstw == "SST-SUST") {
                    SSTVer_ = SSTSUST;
                } else if (sstw == "SST-V") {
                    SSTVer_ = SSTV;
                } else {
                    errorAbortStr(("Unknown SST version: " + sstw));
                }
            }
        }
    }
    switch (SSTVer_) {
    case (standardSST):
    case (SSTSUST):
    case (SSTV):
        gam1_ = beta1_ / betas_ - sigmaw1_ * sqr(kappa_) / sqrt(betas_);
        gam2_ = beta2_ / betas_ - sigmaw2_ * sqr(kappa_) / sqrt(betas_);
        break;
    case (SST2003):
        gam1_ = 5.0 / 9.0;
        gam2_ = 0.44;
        break;
    default:
        break;
    }

    RANSModel::nEq_ = 2;

    k_.rhs().setWriteResult();
    w_.rhs().setWriteResult();
    k_.diagSource().setWriteResult();
    w_.diagSource().setWriteResult();
    bndValueSetting(cont.topCont());

    if (isSplitting()) {
        rak_.resize(mesh().nFaces());
        raw_.resize(mesh().nFaces());
        dqk_.resize(mesh().nTotalCells());
        dqw_.resize(mesh().nTotalCells());
    }
    yplusLam_ = yPlusLam(kappa_, E_);

    const auto &tcont = cont.topCont();
    if (tcont.found("flow")) {
        const auto &fcont = tcont.subController("flow");
        if (fcont.found("limits")) {
            const auto &lmtcont = fcont.subController("limits");
            minK_ = lmtcont.findOrDefault<real>("minKt", minK_);
            minW_ = lmtcont.findOrDefault<real>("minWt", minW_);
            if (lmtcont.found("limitMethod")) {
                const auto lmtw = lmtcont.findWord("limitMethod");
                if (lmtw == "directly") {
                    averageLimit_ = false;
                } else if (lmtw == "average") {
                    averageLimit_ = true;
                } else {
                    errorAbortStr(("Unknown limit method: " + lmtw + " in " + lmtcont.name()));
                }
            }
        }
    }
}

void OpenHurricane::SST::turbParaInitialize() {
    //k_.initialize();
    //w_.initialize();
    real k0 = k_.initialize();
    real w0 = w_.initialize();
    /* k0 = 1.0e-14 / ref().e();
     k_ = k0;
     w0 = 1.0e-20 * ref().time();
     w_ = w0;*/
    k_.updateBoundary();
    w_.updateBoundary();

    mut() = rho()[0] * k0 / w0;

    //20210318 ����mut�߽����
    mut().updateBoundary();
}

void OpenHurricane::SST::expSource() {
    if (isSplitting()) {
        return;
    }
    const auto &meshVol = mesh().cellVol();
    const auto nCell = mesh().nCells();
    const auto &grdV = v().grad();

    kInvFlux = k_.rhs();
    wInvFlux = w_.rhs();
    // The production of k
    getPk(Pk_);
    for (integer n = 0; n < nCell; ++n) {
        const real dist = wallDist()[n];
        const real nut = mut()[n] / rho()[n];
        const real rhok = rho()[n] * k_[n];
        const real rhow = rho()[n] * w_[n];

        // 20210603 ��˼��
        kLastTime_[n] = k_[n];
        wLastTime_[n] = w_[n];

        //Blending function
        const real divkw = k_.grad()[n] * w_.grad()[n];
        const real cdsst = 2.0 * sigmaw2_ * rho()[n] * divkw / w_[n];
        /*const real CDkw = CDkOmega(cdsst);
        const real arg11 = sqrt(k_[n]) / (betas_ * w_[n] * dist);
        const real dist2 = sqr(dist);
        const real arg12 = 500.0 * mul()[n] / (rho()[n] * dist2 * w_[n] * Reref);
        const real arg13 = 4.0 * rho()[n] * sigmaw2_ * k_[n] / (CDkw * dist2);
        const real arg1 = min(max(arg11, arg12), arg13);
        const real blend = tanh(pow4(arg1));*/
        const auto blend = F1(n, cdsst);
        F1_[n] = blend;

        real pk = Pk_[n];
        const real dk = betas_ * rhow * k_[n];
        pk = minP(pk, dk);
        real sk = pk - dk;
        if (SSTVer_ == SSTSUST) {
            sk += betas_ * rho()[n] * wamb_ * kamb_;
        }

        const real w1 = blend;
        const real w2 = 1.0 - w1;
        const real betat = w1 * beta1_ + w2 * beta2_;
        const real gamt = w1 * gam1_ + w2 * gam2_;
        real pw;
        if (SSTVer_ == SST2003) {
            pw = gamt * pk / nut;
        } else {
            pw = gamt * Pk_[n] / nut;
        }
        const real dw = betat * rhow * w_[n];
        const real pd = w2 * cdsst;
        real sw = pw - dw + pd;
        if (SSTVer_ == SSTSUST) {
            sw += betat * rho()[n] * wamb_ * wamb_;
        }
        k_.rhs()[n] += sk * meshVol[n];
        w_.rhs()[n] += sw * meshVol[n];
        kSource[n] = sk * meshVol[n];
        wSource[n] = sw * meshVol[n];
    }
    correctOmegaSource();
    correctOmegaRHS();
    updateF1();
}

void OpenHurricane::SST::impSource() {
    if (isSplitting()) {
        return;
    }
    const auto &meshVol = mesh().cellVol();
    const auto nCell = mesh().nCells();
    const auto &grdV = v().grad();

    realArray dpdrkt(nCell, Zero);

    kInvFlux = k_.rhs();
    wInvFlux = w_.rhs();

    // The production of k
    getPk(Pk_, &dpdrkt);

    for (integer n = 0; n < nCell; ++n) {
        k_.diagSource()[n] = 0.0;
        w_.diagSource()[n] = 0.0;
        const real dist = wallDist()[n];
        const real nut = mut()[n] / rho()[n];
        const real rhok = rho()[n] * k_[n];
        const real rhow = rho()[n] * w_[n];

        // 20210603 ��˼��
        kLastTime_[n] = k_[n];
        wLastTime_[n] = w_[n];

        //Blending function
        const real divkw = k_.grad()[n] * w_.grad()[n];
        const real cdsst = 2.0 * sigmaw2_ * rho()[n] * divkw / w_[n];
        /*const real CDkw = CDkOmega(cdsst);
        const real arg11 = sqrt(k_[n]) / (betas_ * w_[n] * dist);
        const real dist2 = sqr(dist);
        const real arg12 = 500.0 * mul()[n] / (rho()[n] * dist2 * w_[n] * Reref);
        const real arg13 = 4.0 * rho()[n] * sigmaw2_ * k_[n] / (CDkw * dist2);
        const real arg1 = min(max(arg11, arg12), arg13);
        const real blend = tanh(pow4(arg1));*/
        const auto blend = F1(n, cdsst);
        F1_[n] = blend;

        real pk = Pk_[n];
        const real pk0 = pk;
        const real dk = betas_ * rhow * k_[n];
        pk = minP(pk, dk);
        real sk = pk - dk;
        if (SSTVer_ == SSTSUST) {
            sk += betas_ * rho()[n] * wamb_ * kamb_;
        }
        const real w1 = blend;
        const real w2 = 1.0 - w1;
        const real betat = w1 * beta1_ + w2 * beta2_;
        const real gamt = w1 * gam1_ + w2 * gam2_;
        real pw;
        if (SSTVer_ == SST2003) {
            pw = gamt * pk / nut;
        } else {
            pw = gamt * pk0 / nut;
        }
        const real dw = betat * rhow * w_[n];
        const real pd = w2 * cdsst;
        real sw = pw - dw + pd;

        if (SSTVer_ == SSTSUST) {
            sw += betat * rho()[n] * wamb_ * wamb_;
        }
        k_.rhs()[n] += sk * meshVol[n];
        w_.rhs()[n] += sw * meshVol[n];
        kSource[n] = sk * meshVol[n];
        wSource[n] = sw * meshVol[n];

        //k_.diagSource()[n] = (dpkdrk(pk0, dk, rhok, dpdrkt[n]) - betas_ * w_[n]) * meshVol[n];
        //w_.diagSource()[n] = (gamt / nut * dpkdrw(pk0, dk, rhow) - betat * w_[n]) * meshVol[n];
        k_.diagSource()[n] = (-betas_ * w_[n]) * meshVol[n];
        w_.diagSource()[n] = (-2.0 * betat * w_[n] - pd / rhow) * meshVol[n];
    }

    correctOmegaSourceImp();
    correctOmegaRHSDiagImp();
    updateF1();
}

void OpenHurricane::SST::fullImpSource(cellRealSquareMatrixArray &Jac, const integer rhoId,
                                   const integer rhoTurb0) {
    if (isSplitting()) {
        return;
    }
    const auto &meshVol = mesh().cellVol();
    const auto nCell = mesh().nCells();
    const auto &grdV = v().grad();

    realArray dpdrkt(nCell, Zero);

    kInvFlux = k_.rhs();
    wInvFlux = w_.rhs();

    // The production of k
    getPk(Pk_, &dpdrkt);

    for (integer n = 0; n < nCell; ++n) {
        const real dist = wallDist()[n];
        const real nut = mut()[n] / rho()[n];
        const real rhok = rho()[n] * k_[n];
        const real rhow = rho()[n] * w_[n];

        // 20210603 ��˼��
        kLastTime_[n] = k_[n];
        wLastTime_[n] = w_[n];

        //Blending function
        const real divkw = k_.grad()[n] * w_.grad()[n];
        const real cdsst = 2.0 * sigmaw2_ * rho()[n] * divkw / w_[n];
        const auto blend = F1(n, cdsst);
        F1_[n] = blend;

        real pk = Pk_[n];
        const real pk0 = pk;
        const real dk = betas_ * rhow * k_[n];
        pk = minP(pk, dk);
        real sk = pk - dk;
        if (SSTVer_ == SSTSUST) {
            sk += betas_ * rho()[n] * wamb_ * kamb_;
        }
        const real w1 = blend;
        const real w2 = 1.0 - w1;
        const real betat = w1 * beta1_ + w2 * beta2_;
        const real gamt = w1 * gam1_ + w2 * gam2_;
        real pw;
        if (SSTVer_ == SST2003) {
            pw = gamt * pk / nut;
        } else {
            pw = gamt * pk0 / nut;
        }
        const real dw = betat * rhow * w_[n];
        const real pd = w2 * cdsst;
        real sw = pw - dw + pd;

        if (SSTVer_ == SSTSUST) {
            sw += betat * rho()[n] * wamb_ * wamb_;
        }
        k_.rhs()[n] += sk * meshVol[n];
        w_.rhs()[n] += sw * meshVol[n];
        kSource[n] = sk * meshVol[n];
        wSource[n] = sw * meshVol[n];

        Jac[n](rhoTurb0, rhoId) = betas_ * w_[n] * k_[n] * meshVol[n];
        Jac[n](rhoTurb0, rhoTurb0) =
            (dpkdrk(pk0, dk, rhok, dpdrkt[n]) - betas_ * w_[n]) * meshVol[n];
        Jac[n](rhoTurb0, rhoTurb0 + 1) = (-betas_ * k_[n]) * meshVol[n];
        Jac[n](rhoTurb0 + 1, rhoId) = (betat * w_[n] * w_[n] + 2 * pd / rho()[n]) * meshVol[n];
        Jac[n](rhoTurb0 + 1, rhoTurb0) =
            (gamt * dpkdrk(pk0, dk, rhok, dpdrkt[n]) / nut) * meshVol[n];
        Jac[n](rhoTurb0 + 1, rhoTurb0 + 1) = (-2.0 * betat * w_[n] - pd / rhow) * meshVol[n];
    }

    correctOmegaSourceImp();
    correctOmegaRHSImp(Jac, rhoId, rhoTurb0);
    updateF1();
}

void OpenHurricane::SST::visFlux(const faceRealArray &rhof, const faceRealArray &mulf,
                             const faceRealArray &mutf, const cellRealArray &mul,
                             const cellRealArray &mut, const cellRealArray &rho) {
    if (isSplitting()) {
        return;
    }
    const auto &fzl = mesh().faceZones();
    const auto &fcl = mesh().faces();
    const auto &fA = mesh().faceArea();

    faceVectorArray gkf(fv::gradf(k_));
    faceVectorArray gwf(fv::gradf(w_));
    faceRealArray F1f(fv::interpolate(F1_));

    kVisFluxTmp = Zero;
    wVisFluxTmp = Zero;
    for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
        if (fzl[fzi].isWall()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
                const auto &cl = fcl[fi].leftCell();
                const auto &cr = fcl[fi].rightCell();
                const real coef_turb = mulf[fi];
                const real fluxk = coef_turb * (gkf[fi] * fA[fi]);
                const real fluxw = coef_turb * (gwf[fi] * fA[fi]);
                k_.rhs()[cl] -= fluxk;
                w_.rhs()[cl] -= fluxw;
                kVisFluxTmp[cl] -= fluxk;
                wVisFluxTmp[cl] -= fluxw;
                // 20210322 ��˼�� ���� �洢k����ճ��ͨ��
                kVisFlux_[fi] = fluxk;
            }
        } else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
                const auto &cl = fcl[fi].leftCell();
                const auto &cr = fcl[fi].rightCell();
                const real sigmak = F1f[fi] * sigmak1_ + (1.0 - F1f[fi]) * sigmak2_;
                real coef_turb = (mulf[fi] + sigmak * mutf[fi]);
                const real fluxk = coef_turb * (gkf[fi] * fA[fi]);
                k_.rhs()[cl] -= fluxk;
                k_.rhs()[cr] += fluxk;
                kVisFluxTmp[cl] -= fluxk;
                kVisFluxTmp[cr] += fluxk;

                const real sigmaw = F1f[fi] * sigmaw1_ + (1.0 - F1f[fi]) * sigmaw2_;
                coef_turb = (mulf[fi] + sigmaw * mutf[fi]);
                const real fluxw = coef_turb * (gwf[fi] * fA[fi]);
                w_.rhs()[cl] -= fluxw;
                w_.rhs()[cr] += fluxw;
                wVisFluxTmp[cl] -= fluxw;
                wVisFluxTmp[cr] += fluxw;

                // 20210322 ��˼�� ���� �洢k����ճ��ͨ��
                kVisFlux_[fi] = fluxk;
            }
        }
    }
}

void OpenHurricane::SST::update() {
    const integer nCells = mesh().nCells();

    correctOmegaRHS();
    for (integer n = 0; n < nCells; ++n) {
        // calculate viscosity based on SST Model
        const real dist = wallDist()[n];
        const real rhok = rho()[n] * k_[n];
        /*const real arg21 = 2.0 * sqrt(k_[n]) / (betas_ * w_[n] * dist);
        const real arg22 = 500.0 * mul()[n] / (rho()[n] * sqr(dist) * w_[n] * Reref);
        const real arg2 = max(arg21, arg22);
        const real F2 = tanh(arg2 * arg2);*/
        const auto F22 = F2(n);
        F2_[n] = F22;
        real S = 0.0;
        if (SSTVer_ == standardSST || SSTVer_ == SSTSUST || SSTVer_ == SSTV) {
            // Use the vorticity magnitude
            S = skewMagnitude(v().grad()[n]);
        } else if (SSTVer_ == SST2003) {
            // Use the strain invariant
            const symmTensor Sij = symm(v().grad()[n]);
            S = 2.0 * aijbij(Sij, Sij);
            S = sqrt(S);
        }
        SS_[n] = S;
        mut()[n] = rhok * a1_ / max(a1_ * w_[n], S * F22);

        if (std::isnan(mut()[n]) || std::isinf(mut()[n])) {
            errorAbortStr(("Mut is nan: k = " + toString(k_[n]) + ", w = " + toString(w_[n])));
        }
    }
    limit();
    mut().updateBoundary();
}

void OpenHurricane::SST::mutBoundaryUpdate() {
    const faceZoneList &fZ = mesh().faceZones();
    const faceList &fL = mesh().faces();

    const real Cmu025 = pow(Cmu_, real(0.25));
    realTransfer myTransfer(mesh(), mut(), false, true);
    myTransfer.transferInit();
    for (integer fZI = 0; fZI < fZ.size(); ++fZI) {
        if (fZ[fZI].isInterior() || fZ[fZI].isCutFace() || fZ[fZI].isPeriodic() ||
            fZ[fZI].isPeriodicShadow()) {
            continue;
        } else {
            if (fZ[fZI].isWall()) {
                for (integer fi = fZ[fZI].firstIndex(); fi <= fZ[fZI].lastIndex(); ++fi) {
                    const auto &cl = fL[fi].leftCell();
                    const auto &cr = fL[fi].rightCell();
                    if (omegaWallFunction_) {
                        real rhow = (rho()[cl] + rho()[cr]) * real(0.50);
                        real muw = (mul()[cl] + mul()[cr]) * real(0.50);
                        const real yPlus = Cmu025 * y().wallDist()[cl] * sqrt(k()[cl]) / (muw / rhow);
                        if (yPlus > yplusLam_) {
                            const real mutw = muw * (yPlus * kappa_ / log(E_ * yPlus) - 1.0);
                            mut()[cr] = 2.0 * mutw - mut()[cl];
                        } else {
                            mut()[cr] = -mut()[cl];
                        }
                    } else {
                        mut()[cr] = -mut()[cl];
                    }
                }
            } else if (fZ[fZI].bcType() == faceBCType::bcTypes::PRESSUREFARFIELD ||
                       fZ[fZI].bcType() == faceBCType::bcTypes::PRESSUREINLET ||
                       fZ[fZI].bcType() == faceBCType::bcTypes::VELOCITYINLET ||
                       fZ[fZI].bcType() == faceBCType::bcTypes::MASSFLOWINLET) {
                for (integer fi = fZ[fZI].firstIndex(); fi <= fZ[fZI].lastIndex(); ++fi) {
                    const auto &cl = fL[fi].leftCell();
                    const auto &cr = fL[fi].rightCell();

                    mut()[cr] = rho()[cr] * k_[cr] / w_[cr];
                }
            } else {
                for (integer fi = fZ[fZI].firstIndex(); fi <= fZ[fZI].lastIndex(); ++fi) {
                    const auto &cl = fL[fi].leftCell();
                    const auto &cr = fL[fi].rightCell();

                    mut()[cr] = mut()[cl];
                }
            }
        }
    }
    //fv::transfer(mut(), true);
    myTransfer.transferring();
}

void OpenHurricane::SST::updateF1() {
    const faceZoneList &fZ = mesh().faceZones();

    const faceList &fL = mesh().faces();
    realTransfer myTransfer(mesh(), F1_, false, true);
    myTransfer.transferInit();
    for (integer fZI = 0; fZI < fZ.size(); ++fZI) {
        if (fZ[fZI].isInterior() || fZ[fZI].isCutFace() || fZ[fZI].isPeriodic() ||
            fZ[fZI].isPeriodicShadow()) {
            continue;
        } else {
            for (integer fI = fZ[fZI].firstIndex(); fI <= fZ[fZI].lastIndex(); ++fI) {
                const auto &cl = fL[fI].leftCell();
                const auto &cr = fL[fI].rightCell();
                F1_(cr) = F1_(cl);
            }
        }
    }
    //fv::transfer(F1_, true);
    myTransfer.transferring();
}

void OpenHurricane::SST::limit() {
    const integer nCells = mesh().nCells();
    for (integer n = 0; n < nCells; ++n) {
        mut()[n] = max(mut()[n], mutLow() * mul()[n]);
        mut()[n] = min(mut()[n], mutHigh() * mul()[n]);
    }
}

OpenHurricane::realArray OpenHurricane::SST::k() const {
    return k_.array_ref();
}

OpenHurricane::realArray OpenHurricane::SST::epsilon() const {
    return betas_ * w_.array_ref() * k_.array_ref();
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::SST::Ret() const {
    return k_ * k_ / (flowM_.nu() * epsilon());
}

OpenHurricane::cellRealArray &OpenHurricane::SST::var(const integer i) {
    if (i == 0) {
        return k_;
    } else if (i == 1) {
        return w_;
    } else {
        errorAbortStr(("Unknow equation index:\"" + toString(i) + "\" of SST moddel"));
    }
    return k_;
}

const OpenHurricane::cellRealArray &OpenHurricane::SST::var(const integer i) const {
    if (i == 0) {
        return k_;
    } else if (i == 1) {
        return w_;
    } else {
        errorAbortStr(("Unknow equation index:\"" + toString(i) + "\" of SST moddel"));
    }
    return k_;
}

void OpenHurricane::SST::solving(const realArray &dt) {
    //return;
    if (isCoupled()) {
        return;
    }

    // First step: Calculate turbulence source
    const auto &meshVol = mesh().cellVol();
    const auto nCell = mesh().nCells();
    const auto &grdV = v().grad();

    const auto nFace = mesh().nFaces();
    const auto &faceWeight = mesh().faceWgt();
    realArray dpdrkt(nCell, Zero);
    getPk(Pk_, &dpdrkt);

    for (integer n = 0; n < nCell; ++n) {
        k_.diagSource()[n] = 0.0;
        w_.diagSource()[n] = 0.0;
        const real dist = wallDist()[n];
        const real nut = mut()[n] / rho()[n];
        const real rhok = rho()[n] * k_[n];
        const real rhow = rho()[n] * w_[n];

        // 20210603 ��˼��
        kLastTime_[n] = k_[n];
        wLastTime_[n] = w_[n];

        //Blending function
        const real divkw = k_.grad()[n] * w_.grad()[n];
        const real cdsst = 2.0 * sigmaw2_ * rho()[n] * divkw / w_[n];
        /*const real CDkw = CDkOmega(cdsst);
        const real arg11 = sqrt(k_[n]) / (betas_ * w_[n] * dist);
        const real dist2 = sqr(dist);
        const real arg12 = 500.0 * mul()[n] / (rho()[n] * dist2 * w_[n] * Reref);
        const real arg13 = 4.0 * rho()[n] * sigmaw2_ * k_[n] / (CDkw * dist2);
        const real arg1 = min(max(arg11, arg12), arg13);
        const real blend = tanh(pow4(arg1));*/
        const auto blend = F1(n, cdsst);
        F1_[n] = blend;

        real pk = Pk_[n];
        const real pk0 = pk;

        const real dk = betas_ * rhow * k_[n];
        pk = minP(pk, dk);
        real sk = pk - dk;
        if (SSTVer_ == SSTSUST) {
            sk += betas_ * rho()[n] * wamb_ * kamb_;
        }
        const real w1 = blend;
        const real w2 = 1.0 - w1;
        const real betat = w1 * beta1_ + w2 * beta2_;
        const real gamt = w1 * gam1_ + w2 * gam2_;
        real pw;
        if (SSTVer_ == SST2003) {
            pw = gamt * pk / nut;
        } else {
            pw = gamt * pk0 / nut;
        }
        const real dw = betat * rhow * w_[n];
        const real pd = w2 * cdsst;
        real sw = pw - dw + pd;
        if (SSTVer_ == SSTSUST) {
            sw += betat * rho()[n] * wamb_ * wamb_;
        }

        k_.rhs()[n] = sk * meshVol[n];

        w_.rhs()[n] = sw * meshVol[n];
        kSource[n] = sk * meshVol[n];
        wSource[n] = sw * meshVol[n];

        //k_.diagSource()[n] = -(dpkdrk(pk0, dk, rhok, dpdrkt[n]) - betas_ * w_[n]) * meshVol[n];
        //w_.diagSource()[n] = -(gamt / nut * dpkdrw(pk0, dk, rhow) - betat * w_[n]) * meshVol[n];
        k_.diagSource()[n] = -(-betas_ * w_[n]) * meshVol[n];
        w_.diagSource()[n] = -(-2.0 * betat * w_[n] - pd / rhow) * meshVol[n];
    }
    correctOmegaSourceImp();
    updateF1();

    // Second step: Calculate the flux of k and w
    const auto &fcl = mesh().faces();
    const auto &fA = mesh().faceArea();
    const auto &cC = mesh().cellCentre();

    kInvFlux = Zero;
    wInvFlux = Zero;
    kVisFluxTmp = Zero;
    wVisFluxTmp = Zero;
    const faceRealArray F1f(fv::interpolate(F1_));

    for (integer fi = 0; fi < nFace; ++fi) {
        const auto &cl = fcl[fi].leftCell();
        const auto &cr = fcl[fi].rightCell();

        const real uul = v()[cl] * fA[fi];
        const real uum = 0.5 * (uul - fabs(uul));

        const real uur = v()[cr] * fA[fi];
        const real uup = 0.5 * (uur + fabs(uur));

        // Inviscous flux
        real fluxk = rho()[cl] * k_[cl] * uum + rho()[cr] * k_[cr] * uup;
        real fluxw = rho()[cl] * w_[cl] * uum + rho()[cr] * w_[cr] * uup;

        k_.rhs()[cl] += fluxk;
        k_.rhs()[cr] -= fluxk;

        //----------------------------
        kInvFlux[cl] += fluxk;
        kInvFlux[cr] -= fluxk;
        //------------------------
        w_.rhs()[cl] += fluxw;
        w_.rhs()[cr] -= fluxw;

        //===============================
        wInvFlux[cl] += fluxw;
        wInvFlux[cr] -= fluxw;
        //---------------------------

        k_.diagSource()[cl] -= uum;
        k_.diagSource()[cr] += uup;

        w_.diagSource()[cl] -= uum;
        w_.diagSource()[cr] += uup;
        rak_[fi][0] = -uup;
        rak_[fi][1] = uum;
        raw_[fi][0] = -uup;
        raw_[fi][1] = uum;

        // viscous flux
        vector rlr = cC[cl] - cC[cr];
        const real nnij = rlr.magnitude();

        const real mull = mul()[cl];
        const real mulr = mul()[cr];
        const real mulf = faceWeight[fi] * mull + (1.0 - faceWeight[fi]) * mulr;

        const real mutl = mut()[cl];
        const real mutr = mut()[cr];
        const real mutf = faceWeight[fi] * mutl + (1.0 - faceWeight[fi]) * mutr;
        const real sigmak = F1f[fi] * sigmak1_ + (1.0 - F1f[fi]) * sigmak2_;
        const real sigmaw = F1f[fi] * sigmaw1_ + (1.0 - F1f[fi]) * sigmaw2_;
        const real nn = fA[fi].magnitude();

        const real coek1 = (mulf + mutf * sigmak) / nnij / rho()[cl] * nn;
        const real coek2 = (mulf + mutf * sigmak) / nnij / rho()[cr] * nn;
        fluxk = coek1 * rho()[cl] * k_[cl] - coek2 * rho()[cr] * k_[cr];

        const real coew1 = (mulf + mutf * sigmaw) / nnij / rho()[cl] * nn;
        const real coew2 = (mulf + mutf * sigmaw) / nnij / rho()[cr] * nn;
        fluxw = coew1 * rho()[cl] * w_[cl] - coew2 * rho()[cr] * w_[cr];

        // 20210322 ��˼�� ���� �洢k����ճ��ͨ��
        kVisFlux_[fi] = fluxk;

        k_.rhs()[cl] -= fluxk;
        k_.rhs()[cr] += fluxk;

        w_.rhs()[cl] -= fluxw;
        w_.rhs()[cr] += fluxw;

        //----------------
        kVisFluxTmp[cl] -= fluxk;
        kVisFluxTmp[cr] += fluxk;
        wVisFluxTmp[cl] -= fluxw;
        wVisFluxTmp[cr] += fluxw;
        //-----------

        k_.diagSource()[cl] += coek1;
        k_.diagSource()[cr] += coek2;

        w_.diagSource()[cl] += coew1;
        w_.diagSource()[cr] += coew2;

        rak_[fi][0] -= coek2;
        rak_[fi][1] -= coek1;
        raw_[fi][0] -= coew2;
        raw_[fi][1] -= coew1;
    }

    for (integer i = 0; i < dqk_.size(); ++i) {
        dqk_[i] = 0.0;
        dqw_[i] = 0.0;
    }
    //correctOmegaRHS();
    // LU-SGS lower loop
    LUSGS::lowerLoop(k_, rak_, dt, dqk_);
    LUSGS::lowerLoop(w_, raw_, dt, dqw_);

    //20210409 ��˼�� ���Ӵ�������߽��k��w��dqֵ
    /*fv::cutTransfer(mesh(), dqk_, true);
    fv::cutTransfer(mesh(), dqw_, true);*/

    // LU-SGS upper loop
    LUSGS::upperLoop(k_, rak_, dt, dqk_);
    LUSGS::upperLoop(w_, raw_, dt, dqw_);

    const auto &ccl = mesh().cells();
    // LU-SGS update
    for (integer n = 0; n < nCell; ++n) {
        const real rhok = rho()[n] * k_[n];
        const real rhow = rho()[n] * w_[n];

        const real newRhok = rhok + dqk_[n] * kRelax_;
        const real newRhow = rhow + dqw_[n] * wRelax_;

        k_[n] = newRhok / rho()[n];
        w_[n] = newRhow / rho()[n];
    }
}

void OpenHurricane::SST::updateBoundary() {
    k_.updateBoundary();
    w_.updateBoundary();
}

void OpenHurricane::SST::limitAndUpdateBoundary() {
    // 20210408 ��˼�� ���Ӷ�k��w�����ƣ�����ĳ����Ϊ��ֵʱ��ȡ�����������ƽ��ֵ��
    limitKAndW();
    k_.updateBoundary();
    w_.updateBoundary();
}

void OpenHurricane::SST::limitKAndW() {
    /* if (mesh().Iteration().isBeginLimit())
     {*/
    limitTurbVarIncrease(k_, kLastTime_);
    limitTurbVarIncrease(w_, wLastTime_);
    //}
    limitNegative(k_, minK_, averageLimit_, true);
    limitNegative(w_, minW_, averageLimit_, true);
}

OpenHurricane::symmTensorArray OpenHurricane::SST::ReynoldsStressTensor() const {
    return -mut() / rho() *
               (twoSymm(v().grad()) - (real(2.0 / 3.0) * (div(diagToVector(v().grad()))) * I)) +
           real(2.0 / 3.0) * k_ * I;
}

OpenHurricane::faceSymmTensorArray OpenHurricane::SST::tauEff(const faceRealArray &rhof,
                                                      const faceRealArray &mulf,
                                                      const faceRealArray &mutf,
                                                      const faceTensorArray &deltafV) const {
    faceSymmTensorArray tau(
        object("tau", mesh(), object::NOT_WRITE, object::TEMPORARY), mesh(),
        (mulf + mutf) * (twoSymm(deltafV) - (real(2.0 / 3.0) * (div(diagToVector(deltafV))) * I)) -
            real(2.0 / 3.0) * rhof * fv::interpolate(k_) * I);
    return tau;
}

void OpenHurricane::SST::calcGrad(const spatialScheme &sps) {
    if (isCoupled()) {
        return;
    }
    sps.grad(k_);
    sps.grad(w_);
}

void OpenHurricane::SST::correctEnergyEqVisFlux(cellRealArray &E) const {
    const auto &fzl = mesh().faceZones();
    const auto &fcl = mesh().faces();
    const auto &fA = mesh().faceArea();

    for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            const auto &cl = fcl[fi].leftCell();
            const auto &cr = fcl[fi].rightCell();
            E.rhs()[cl] -= kVisFlux_[fi];
            E.rhs()[cr] += kVisFlux_[fi];
        }
    }
}

void OpenHurricane::SST::getPk(realArray &Pk, realArray *dPkdrhokPtr) {
    const auto nCell = mesh().nCells();

    for (integer n = 0; n < nCell; ++n) {
        const real dist = wallDist()[n];
        const real rhok = rho()[n] * k_[n];

        //const tensor tau = mut()[n] * (twoSymm(v().grad()[n]) - real(2.0 / 3.0) * (diag(v().grad()[n]))) / Reref - real(2.0 / 3.0) * rhok * I;
        if (SSTVer_ == SSTV) {
            Pk[n] = (mut()[n] * skewMagSqr(v().grad()[n]) -
                     real(2.0 / 3.0) * rhok * div(diagToVector(v().grad()[n])));
        } else {
            const tensor tau = mut()[n] * (twoSymm(v().grad()[n]) -
                                           real(2.0 / 3.0) * div(diagToVector(v().grad()[n])) * I) -
                               real(2.0 / 3.0) * rhok * I;
            Pk[n] = aijbij(tau, v().grad()[n]);
        }
        if (dPkdrhokPtr != nullptr) {
            (*dPkdrhokPtr)[n] = Pk[n] / rhok;
        }
    }
    if (omegaWallFunction_) {
        const real Cmu025 = pow(Cmu_, real(0.25));
        const auto &fZL = mesh().faceZones();
        const auto &faces = mesh().faces();
        const auto &fA = mesh().faceArea();
        const auto &y = mesh().findObjectRef<cellRealArray>("wallDist");
        for (integer i = 0; i < omegaWallFunctionFaceZoneList_.size(); ++i) {
            const integer fzi = omegaWallFunctionFaceZoneList_[i];
            for (integer facei = fZL[fzi].firstIndex(); facei <= fZL[fzi].lastIndex(); ++facei) {
                const integer cl = faces[facei].leftCell();
                const integer cr = faces[facei].rightCell();
                real uu = v()[cl] * fA[facei].normalized();
                real ut = sqrt(max(v()[cl].magSqr() - sqr(uu), real(0)));
                real rhow = (rho()[cl] + rho()[cr]) * real(0.50);
                real muw = (mul()[cl] + mul()[cr]) * real(0.50);
                real mutw = (mut()[cl] + mut()[cr]) * real(0.50);
                //real tauw = muw * ut / (y[cl] * Reref);
                const real yPlus = Cmu025 * y[cl] * sqrt(k()[cl]) / (muw / rhow);
                if (yPlus > yplusLam_) {
                    Pk[cl] =
                        (mutw + muw) * (ut / y[cl]) * Cmu025 * sqrt(k()[cl]) / (kappa_ * y[cl]);
                    //Pk[cl] = muw * sqr(ut / y[cl]) / Reref;
                    //Pk[cl] = sqr(tauw) / (kappa_ * rhow * Cmu025 * sqrt(k()[cl]) * y[cl]);

                } else {
                    Pk[cl] = 0;
                }
                if (dPkdrhokPtr != nullptr) {
                    (*dPkdrhokPtr)[cl] = 0.0;
                }
                if (isnan(Pk[cl])) {
                    std::stringstream sstr;
                    sstr << muw << ", ut = " << ut << ", v = " << v()[cl].magnitude()
                         << ", uu = " << uu << ", y = " << y[cl] << ", k = " << k()[cl]
                         << std::endl;
                    sstr << toString(v()[cl]) << ", n = " << toString(fA[facei].normalized())
                         << std::endl;
                    errorAbortStr(sstr.str());
                }
            }
        }
    }
}

void OpenHurricane::SST::correctOmegaRHS() {
    if (omegaWallFunction_) {
        const auto &fZL = mesh().faceZones();
        const auto &faces = mesh().faces();
        for (integer i = 0; i < omegaWallFunctionFaceZoneList_.size(); ++i) {
            const integer fzi = omegaWallFunctionFaceZoneList_[i];
            for (integer facei = fZL[fzi].firstIndex(); facei <= fZL[fzi].lastIndex(); ++facei) {
                const integer cl = faces[facei].leftCell();
                w().rhs()[cl] = 0;
            }
        }
    }
}

void OpenHurricane::SST::correctOmegaRHSDiagImp() {
    if (omegaWallFunction_) {
        const auto &fZL = mesh().faceZones();
        const auto &faces = mesh().faces();
        for (integer i = 0; i < omegaWallFunctionFaceZoneList_.size(); ++i) {
            const integer fzi = omegaWallFunctionFaceZoneList_[i];
            for (integer facei = fZL[fzi].firstIndex(); facei <= fZL[fzi].lastIndex(); ++facei) {
                const integer cl = faces[facei].leftCell();
                w().rhs()[cl] = 0;
                w().diagSource()[cl] = 0;
            }
        }
    }
}

void OpenHurricane::SST::correctOmegaRHSImp(cellRealSquareMatrixArray &Jac, const integer rhoId,
                                        const integer rhoTurb0) {
    if (omegaWallFunction_) {
        const auto &fZL = mesh().faceZones();
        const auto &faces = mesh().faces();
        for (integer i = 0; i < omegaWallFunctionFaceZoneList_.size(); ++i) {
            const integer fzi = omegaWallFunctionFaceZoneList_[i];
            for (integer facei = fZL[fzi].firstIndex(); facei <= fZL[fzi].lastIndex(); ++facei) {
                const integer cl = faces[facei].leftCell();
                w().rhs()[cl] = 0;
                Jac[cl](rhoTurb0 + 1, rhoId) = 0;
                Jac[cl](rhoTurb0 + 1, rhoTurb0) = 0;
                Jac[cl](rhoTurb0 + 1, rhoTurb0 + 1) = 0;
            }
        }
    }
}
