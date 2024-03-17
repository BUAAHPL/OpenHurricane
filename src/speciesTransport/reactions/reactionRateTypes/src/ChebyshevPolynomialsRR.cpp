/*!
 * \file ChebyshevPolynomialsRR.hpp
 * \brief Main subroutines for using Chebyshev polynomials reaction rate.
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

#include "ChebyshevPolynomialsRR.hpp"

namespace OpenHurricane {
    createClassName(ChebyshevPolynomialsRR);
         registerObjFty(reactionRateTypes, ChebyshevPolynomialsRR, controller);
} // namespace OpenHurricane

OpenHurricane::ChebyshevPolynomialsRR::ChebyshevPolynomialsRR(const realArrayArray &ANM, const real min,
                                                          const real max, const bool defaultTorP)
    : reactionRateTypes(), Anm_(ANM), thirdBodyEfficiency_(nullptr) {
    if (defaultTorP) {
        Pmin_ = min;
        Pmax_ = max;
        Tmin_ = ChebyshevRRLimit::temperatureLimit::Tmin;
        Tmax_ = ChebyshevRRLimit::temperatureLimit::Tmax;
    } else {
        Tmin_ = min;
        Tmax_ = max;
        Pmin_ = ChebyshevRRLimit::pressureLimit::Pmin;
        Pmax_ = ChebyshevRRLimit::pressureLimit::Pmax;
    }
    powTminusOneMin_ = 1.0 / Tmin_;
    powTminusOneMax_ = 1.0 / Tmax_;
    log10Pmin_ = std::log10(Pmin_);
    log10Pmax_ = std::log10(Pmax_);
}

OpenHurricane::ChebyshevPolynomialsRR::ChebyshevPolynomialsRR(const realArrayArray &ANM, const real min,
                                                          const real max, const bool defaultTorP,
                                                          const thirdBodyEfficiency &tbe)
    : reactionRateTypes(), Anm_(ANM), thirdBodyEfficiency_(nullptr) {
    if (defaultTorP) {
        Pmin_ = min;
        Pmax_ = max;
        Tmin_ = ChebyshevRRLimit::temperatureLimit::Tmin;
        Tmax_ = ChebyshevRRLimit::temperatureLimit::Tmax;
    } else {
        Tmin_ = min;
        Tmax_ = max;
        Pmin_ = ChebyshevRRLimit::pressureLimit::Pmin;
        Pmax_ = ChebyshevRRLimit::pressureLimit::Pmax;
    }
    powTminusOneMin_ = 1.0 / Tmin_;
    powTminusOneMax_ = 1.0 / Tmax_;
    log10Pmin_ = std::log10(Pmin_);
    log10Pmax_ = std::log10(Pmax_);

    thirdBodyEfficiency_.reset(new thirdBodyEfficiency(tbe));
}

OpenHurricane::ChebyshevPolynomialsRR::ChebyshevPolynomialsRR(const realArrayArray &ANM,
                                                          const real Tmin, const real Tmax,
                                                          const real Pmin, const real Pmax)
    : reactionRateTypes(), Anm_(ANM), Tmin_(Tmin), Tmax_(Tmax), Pmin_(Pmin), Pmax_(Pmax),
      thirdBodyEfficiency_(nullptr) {
    powTminusOneMin_ = 1.0 / Tmin_;
    powTminusOneMax_ = 1.0 / Tmax_;
    log10Pmin_ = std::log10(Pmin_);
    log10Pmax_ = std::log10(Pmax_);
}

OpenHurricane::ChebyshevPolynomialsRR::ChebyshevPolynomialsRR(const realArrayArray &ANM,
                                                          const real Tmin, const real Tmax,
                                                          const real Pmin, const real Pmax,
                                                          const thirdBodyEfficiency &tbe)
    : reactionRateTypes(), Anm_(ANM), Tmin_(Tmin), Tmax_(Tmax), Pmin_(Pmin), Pmax_(Pmax),
      thirdBodyEfficiency_(nullptr) {
    powTminusOneMin_ = 1.0 / Tmin_;
    powTminusOneMax_ = 1.0 / Tmax_;
    log10Pmin_ = std::log10(Pmin_);
    log10Pmax_ = std::log10(Pmax_);
    thirdBodyEfficiency_.reset(new thirdBodyEfficiency(tbe));
}

OpenHurricane::ChebyshevPolynomialsRR::ChebyshevPolynomialsRR(const realArray &ANM,
                                                          const realArray &TMinMax,
                                                          const realArray &PMinMax,
                                                          const speciesList &spt,
                                                          const realArray &efficiences)
    : reactionRateTypes(), Anm_(), Tmin_(ChebyshevRRLimit::temperatureLimit::Tmin),
      Tmax_(ChebyshevRRLimit::temperatureLimit::Tmax), Pmin_(ChebyshevRRLimit::pressureLimit::Pmin),
      Pmax_(ChebyshevRRLimit::pressureLimit::Pmax), thirdBodyEfficiency_(nullptr) {
    integer NN = 0;
    integer MM = 0;
    if (ANM.size() != 0) {
        NN = integer(ANM[0]);
        MM = integer(ANM[1]);
        if (ANM.size() != (NN * MM + 2)) {
            errorAbortStr(("The size of input list is wrong. For n = " + toString(NN) +
                           ", m = " + toString(MM) + ", and the size should be " +
                           toString(NN * MM + 2) + ", but it is " + toString(ANM.size())));
        }
        Anm_.resize(NN);
        for (integer i = 0; i < NN; i++) {
            Anm_[i].resize(MM);
            for (integer j = 0; j < MM; ++j) {
                Anm_[i][j] = ANM[i * MM + j + 2];
            }
        }
    }

    if (TMinMax.size() == 2) {
        Tmin_ = TMinMax[0];
        Tmax_ = TMinMax[1];
    }
    if (PMinMax.size() == 2) {
        Pmin_ = PMinMax[0];
        Pmax_ = PMinMax[1];
    }
    if (efficiences.size() != 0) {
        thirdBodyEfficiency_.reset(new thirdBodyEfficiency(spt, efficiences));
    }
    powTminusOneMin_ = 1.0 / Tmin_;
    powTminusOneMax_ = 1.0 / Tmax_;
    log10Pmin_ = std::log10(Pmin_);
    log10Pmax_ = std::log10(Pmax_);
}

OpenHurricane::ChebyshevPolynomialsRR::ChebyshevPolynomialsRR(const speciesList &sp,
                                                          const controller &cont)
    : reactionRateTypes(sp, cont), Anm_(),
      Tmin_(cont.findOrDefault<real>("Tmin", ChebyshevRRLimit::temperatureLimit::Tmin)),
      Tmax_(cont.findOrDefault<real>("Tmax", ChebyshevRRLimit::temperatureLimit::Tmax)),
      Pmin_(cont.findOrDefault<real>("Pmin", ChebyshevRRLimit::pressureLimit::Pmin)),
      Pmax_(cont.findOrDefault<real>("Pmax", ChebyshevRRLimit::pressureLimit::Pmax)),
      thirdBodyEfficiency_(nullptr) {
    integer ASize = cont.findType<integer>("ASize", integer());
    Anm_.resize(ASize);

    for (integer i = 0; i < ASize; i++) {
        Anm_[i] = cont.findType<realList>(toString(i), realList());
    }

    if (cont.found("thirdBody")) {
        thirdBodyEfficiency_.reset(new thirdBodyEfficiency(sp, cont));
    }

    powTminusOneMin_ = 1.0 / Tmin_;
    powTminusOneMax_ = 1.0 / Tmax_;
    log10Pmin_ = std::log10(Pmin_);
    log10Pmax_ = std::log10(Pmax_);
}

hur_nodiscard OpenHurricane::real OpenHurricane::ChebyshevPolynomialsRR::k(const real p, const real T,
                                                                   const realArray &c) const {
    const real powTminusOne = 1.0 / max(T, tiny);
    const real log10P = std::log10(p);
    const real Ttransform = transformation(powTminusOne, powTminusOneMin_, powTminusOneMax_);

    const real Ptransform = transformation(log10P, log10Pmin_, log10Pmax_);

#ifdef HUR_DEBUG
    bool checkOk = true;
    if (Anm_.size() == 0) {
        checkOk = false;
    } else {
        for (integer n = 0; n < Anm_.size(); n++) {
            if (Anm_[n].size() == 0) {
                checkOk = false;
                break;
            }
        }
    }
    if (!checkOk) {
        LFatal("The size of Anm_ is zero. Please check!");
    }
#endif // HUR_DEBUG

    real logk = Zero;
    for (integer n = 0; n < Anm_.size(); n++) {
        for (integer m = 0; m < Anm_[n].size(); m++) {
            logk += Anm_[n][m] * phi(n, Ttransform) * phi(m, Ptransform);
        }
    }

    // Change unit change from mole/cm^3-s to kmole/m^3-s
    real kj = real(0.001) * std::log10(logk);

    if (thirdBodyEfficiency_) {
        kj *= thirdBodyEfficiency_->M(c);
    }

    return kj;
}

hur_nodiscard OpenHurricane::real OpenHurricane::ChebyshevPolynomialsRR::DkDT(const real kj, const real p,
                                                                      const real T,
                                                                      const realArray &c) const {
    const real ln10 = log(real(10));
    const real dTPdT =
        -real(2) / (sqr(max(T, tiny)) * max(powTminusOneMax_ - powTminusOneMin_, tiny));
    const real powTminusOne = 1.0 / max(T, tiny);
    const real log10P = std::log10(p);
    const real Ttransform = transformation(powTminusOne, powTminusOneMin_, powTminusOneMax_);

    const real Ptransform = transformation(log10P, log10Pmin_, log10Pmax_);
    real dlogkdTP = real(0);
    for (integer n = 0; n < Anm_.size(); n++) {
        for (integer m = 0; m < Anm_[n].size(); m++) {
            dlogkdTP += Anm_[n][m] * DPhiDX(n, Ttransform) * phi(m, Ptransform);
        }
    }
    // The third-body term has been included in k.
    return kj * ln10 * dlogkdTP * dTPdT;
}

hur_nodiscard OpenHurricane::real OpenHurricane::ChebyshevPolynomialsRR::DkDP(const real kj, const real p,
                                                                      const real T,
                                                                      const realArray &c) const {
    const real ln10 = log(real(10));
    const real dPPdP = real(2) / (p * ln10 * max(powTminusOneMax_ - powTminusOneMin_, tiny));
    const real powTminusOne = 1.0 / max(T, tiny);
    const real log10P = std::log10(p);
    const real Ttransform = transformation(powTminusOne, powTminusOneMin_, powTminusOneMax_);

    const real Ptransform = transformation(log10P, log10Pmin_, log10Pmax_);
    real dlogkdPP = real(0);
    for (integer n = 0; n < Anm_.size(); n++) {
        for (integer m = 0; m < Anm_[n].size(); m++) {
            dlogkdPP += Anm_[n][m] * phi(n, Ttransform) * DPhiDX(m, Ptransform);
        }
    }
    // The third-body term has been included in k.
    return kj * ln10 * dlogkdPP * dPPdP;
}

void OpenHurricane::ChebyshevPolynomialsRR::gamDGamDCi(const real P, const real T, const realArray &c,
                                                   realArray &gdgdci) const {
    if (isModefiedByThirdBody()) {
        if (c.size() != thirdBodyEfficiency_->size()) {
            LFatal("Size not equal");
        }
        if (gdgdci.size() != c.size()) {
            gdgdci.resize(c.size(), Zero);
        }
        for (integer i = 0; i < thirdBodyEfficiency_->size(); ++i) {
            gdgdci[i] = (*thirdBodyEfficiency_)[i] / max(thirdBodyEfficiency_->M(c), tiny);
        }
    }
}

void OpenHurricane::ChebyshevPolynomialsRR::restructFromTmpAnm() {
    integer NN = 0;
    integer MM = 0;
    if (tmpAnm_.size() != 0) {
        NN = integer(tmpAnm_[0]);
        MM = integer(tmpAnm_[1]);
        if (tmpAnm_.size() != (NN * MM + 2)) {
            errorAbortStr(("The size of input list is wrong. For n = " + toString(NN) +
                           ", m = " + toString(MM) + ", and the size should be " +
                           toString(NN * MM + 2) + ", but it is " + toString(tmpAnm_.size())));
        }
        Anm_.resize(NN);
        for (integer i = 0; i < NN; i++) {
            Anm_[i].resize(MM);
            for (integer j = 0; j < MM; ++j) {
                Anm_[i][j] = tmpAnm_[i * MM + j + 2];
            }
        }
    }
    tmpAnm_.clear();
}

void OpenHurricane::ChebyshevPolynomialsRR::resetThirdBodyEff(const thirdBodyEfficiency &tbe) {
    thirdBodyEfficiency_.clear();
    thirdBodyEfficiency_.reset(new thirdBodyEfficiency(tbe));
}
