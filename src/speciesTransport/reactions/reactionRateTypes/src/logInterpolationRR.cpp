/*!
 * \file logInterpolationRR.cpp
 * \brief Main subroutines for general pressure dependence reaction rate using logarithmic interpolation.
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

#include "logInterpolationRR.hpp"

namespace OpenHurricane {
    createClassName(logInterpolationRR);
         registerObjFty(reactionRateTypes, logInterpolationRR, controller);
} // namespace OpenHurricane

void OpenHurricane::logInterpolationRR::checkSizeAndPresureRank() {
    if (p_.size() != k_.size()) {
        std::string errMsg;
        errMsg = "The size of the pressure: ";
        errMsg += toString(p_.size());
        errMsg += " is not equal to that of the Arrhenius parameters: ";
        errMsg += toString(k_.size());
        errMsg += " in reaction: " + type();
        errorAbortStr(errMsg);
    }

    // bubble sorting
    for (integer pi = 0; pi < p_.size() - 1; pi++) {
        for (integer pj = 0; pj < p_.size() - 1 - pi; pj++) {
            if (p_[pj] > p_[pj + 1]) {
                Swap(p_[pj], p_[pj + 1]);
                Swap(k_[pj], k_[pj + 1]);
            }
        }
    }
}
OpenHurricane::logInterpolationRR::logInterpolationRR(const speciesList &sp, const controller &cont)
    : reactionRateTypes(sp, cont), p_(cont.findType<realList>("p", realList())), k_(),
      thirdBodyEfficiencyPtr_(nullptr) {
    integer kSize = cont.findOrDefault<integer>("kSize", integer());
    k_.resize(kSize);

    for (integer i = 0; i < kSize; ++i) {
        if (!cont.found(toString(i))) {
            errorAbortStr(("Arrhenius Reaction Rate list in No." + toString(i) + " is missing."));
        }
        k_[i] = ArrheniusRR(sp, cont.subController(toString(i)));
    }

    if (cont.found("thirdBody")) {
        thirdBodyEfficiencyPtr_.reset(new thirdBodyEfficiency(sp, cont));
    }

    checkSizeAndPresureRank();
}

hur_nodiscard OpenHurricane::real OpenHurricane::logInterpolationRR::k(const real p, const real T,
                                                               const realArray &c) const {
    real k = Zero;
    if (p <= p_[0]) {
        k = k_[0].k(p_[0], T, c);
        for (integer pi = 0; pi < p_.size() - 1; pi++) {
            if (p_[pi + 1] > p_[pi]) {
                break;
            } else // In case of multiple minimum pressure coefficients.
            {
                k += k_[pi + 1].k(p_[pi + 1], T, c);
            }
        }
    }

    else if (p >= p_.last()) {
        k = k_.last().k(p_.last(), T, c);

        for (integer pi = p_.size() - 1; pi > 1; pi--) {
            if (p_[pi - 1] < p_[pi]) {
                break;
            } else // In case of multiple maximum pressure coefficients.
            {
                k += k_[pi - 1].k(p_[pi - 1], T, c);
            }
        }

    }

    else {
        for (integer pi = 0; pi < p_.size() - 1; pi++) {
            if (p < p_[pi + 1]) {
                real k0 = k_[pi].k(p_[pi], T, c);
                for (integer pj = pi - 1; pj >= 0; pj--) {
                    if (p_[pj + 1] > p_[pj]) {
                        break;
                    } else {
                        k0 += k_[pj].k(p_[pj], T, c);
                    }
                }

                real k1 = k_[pi + 1].k(p_[pi + 1], T, c);

                for (integer pj = pi + 1; pj < p_.size() - 1; pj++) {
                    if (p_[pj + 1] > p_[pj]) {
                        break;
                    } else {
                        k1 += k_[pj + 1].k(p_[pj + 1], T, c);
                    }
                }

                real lnk0 = std::log(k0);
                real lnk1 = std::log(k1);

                real lnp0 = std::log(p_[pi]);
                real lnp1 = std::log(p_[pi + 1]);
                real lnp = std::log(p);

                real lnk = lnk0 + (lnk1 - lnk0) * (lnp - lnp0) / (lnp1 - lnp0);
                k = std::exp(lnk);
            }
        }
    }

    if (thirdBodyEfficiencyPtr_) {
        k *= thirdBodyEfficiencyPtr_->M(c);
    }

    return k;
}

hur_nodiscard OpenHurricane::real OpenHurricane::logInterpolationRR::DkDT(const real kj, const real p,
                                                                  const real T,
                                                                  const realArray &c) const {
    real dkdT = Zero;
    if (p <= p_[0]) {
        real k0 = k_[0].k(p_[0], T, c);
        dkdT = k_[0].DkDT(k0, p, T, c);
        for (integer pi = 0; pi < p_.size() - 1; pi++) {
            if (p_[pi + 1] > p_[pi]) {
                break;
            } else // In case of multiple minimum pressure coefficients.
            {
                real kpi = k_[pi + 1].k(p_[pi + 1], T, c);
                dkdT += k_[pi + 1].DkDT(kpi, p, T, c);
            }
        }
    }

    else if (p >= p_.last()) {
        real kl = k_.last().k(p_.last(), T, c);
        dkdT = k_.last().DkDT(kl, p, T, c);

        for (integer pi = p_.size() - 1; pi > 1; pi--) {
            if (p_[pi - 1] < p_[pi]) {
                break;
            } else // In case of multiple maximum pressure coefficients.
            {
                real kpi = k_[pi - 1].k(p_[pi - 1], T, c);
                dkdT += k_[pi + 1].DkDT(kpi, p, T, c);
            }
        }

    }

    else {
        for (integer pi = 0; pi < p_.size() - 1; pi++) {
            if (p < p_[pi + 1]) {
                real k0 = k_[pi].k(p_[pi], T, c);
                real dk0dT = k_[pi].DkDT(k0, p_[pi], T, c);
                for (integer pj = pi - 1; pj >= 0; pj--) {
                    if (p_[pj + 1] > p_[pj]) {
                        break;
                    } else {
                        real kpj = k_[pj].k(p_[pj], T, c);
                        k0 += kpj;
                        dk0dT += k_[pj].DkDT(kpj, p_[pj], T, c);
                    }
                }

                real k1 = k_[pi + 1].k(p_[pi + 1], T, c);
                real dk1dT = k_[pi + 1].DkDT(k1, p_[pi + 1], T, c);
                for (integer pj = pi + 1; pj < p_.size() - 1; pj++) {
                    if (p_[pj + 1] > p_[pj]) {
                        break;
                    } else {
                        real kpj = k_[pj + 1].k(p_[pj + 1], T, c);
                        k1 += kpj;
                        dk1dT += k_[pj + 1].DkDT(kpj, p_[pj + 1], T, c);
                    }
                }

                real lnp0 = std::log(p_[pi]);
                real lnp1 = std::log(p_[pi + 1]);
                real lnp = std::log(p);

                real dlnkdT = 1.0 / k0 * dk0dT +
                              (1.0 / k1 * dk1dT - 1.0 / k0 * dk0dT) * (lnp - lnp0) / (lnp1 - lnp0);
                dkdT = kj * dlnkdT;
            }
        }
    }

    if (thirdBodyEfficiencyPtr_) {
        dkdT *= thirdBodyEfficiencyPtr_->M(c);
    }

    return dkdT;
}

hur_nodiscard OpenHurricane::real OpenHurricane::logInterpolationRR::DkDP(const real kj, const real p,
                                                                  const real T,
                                                                  const realArray &c) const {
    real dlnkdp = Zero;
    if (p <= p_[0]) {
        dlnkdp = 0;
    } else if (p >= p_.last()) {
        dlnkdp = 0;
    } else {
        for (integer pi = 0; pi < p_.size() - 1; pi++) {
            if (p < p_[pi + 1]) {
                real k0 = k_[pi].k(p_[pi], T, c);
                for (integer pj = pi - 1; pj >= 0; pj--) {
                    if (p_[pj + 1] > p_[pj]) {
                        break;
                    } else {
                        k0 += k_[pj].k(p_[pj], T, c);
                    }
                }

                real k1 = k_[pi + 1].k(p_[pi + 1], T, c);

                for (integer pj = pi + 1; pj < p_.size() - 1; pj++) {
                    if (p_[pj + 1] > p_[pj]) {
                        break;
                    } else {
                        k1 += k_[pj + 1].k(p_[pj + 1], T, c);
                    }
                }

                real lnk0 = std::log(k0);
                real lnk1 = std::log(k1);

                real lnp0 = std::log(p_[pi]);
                real lnp1 = std::log(p_[pi + 1]);
                real lnp = std::log(p);

                dlnkdp = (lnk1 - lnk0) / (p * (lnp1 - lnp0));
            }
        }
    }
    real dkdp = kj * dlnkdp;

    if (thirdBodyEfficiencyPtr_) {
        dkdp *= thirdBodyEfficiencyPtr_->M(c);
    }

    return dkdp;
}

void OpenHurricane::logInterpolationRR::gamDGamDCi(const real P, const real T, const realArray &c,
                                               realArray &gdgdci) const {
    if (isModefiedByThirdBody()) {
        if (c.size() != thirdBodyEfficiencyPtr_->size()) {
            LFatal("Size not equal");
        }
        if (gdgdci.size() != c.size()) {
            gdgdci.resize(c.size(), Zero);
        }
        for (integer i = 0; i < thirdBodyEfficiencyPtr_->size(); ++i) {
            gdgdci[i] = (*thirdBodyEfficiencyPtr_)[i] / max(thirdBodyEfficiencyPtr_->M(c), tiny);
        }
    }
}

void OpenHurricane::logInterpolationRR::resetThirdBodyEff(const thirdBodyEfficiency &tbe) {
    thirdBodyEfficiencyPtr_.reset(new thirdBodyEfficiency(tbe));
}
