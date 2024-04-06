/*!
 * \file fallOffFunctions.cpp
 * \brief Main subroutines for fall-off functions.
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

#include "fallOffFunctions.hpp"

namespace OpenHurricane {
    namespace fallOffFunctions {
        createClassName(fallOffFunction);
                 createObjFty(fallOffFunction, controller);

        createClassName(Lindemann);
                 registerObjFty(fallOffFunction, Lindemann, controller);

        createClassName(SRI);
                 registerObjFty(fallOffFunction, SRI, controller);

        createClassName(Troe);
                 registerObjFty(fallOffFunction, Troe, controller);
    } // namespace fallOffFunctions
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::fallOffFunctions::fallOffFunction>
OpenHurricane::fallOffFunctions::fallOffFunction::creator(const controller &cont) {
    string fofType = cont.findWord(fallOffFunction::className_);
    defineInObjCreator(fallOffFunction, fofType, controller, (cont));
}

hur_nodiscard OpenHurricane::real OpenHurricane::fallOffFunctions::SRI::F(const real T,
                                                                  const real Pr) const {
    real logPr = log10(max(Pr, tiny));

    real X = 1.0 / (1.0 + sqr(logPr));

    return d_ * std::pow(a_ * exp(-b_ / T) + exp(-T / c_), X) * std::pow(T, e_);
}

hur_nodiscard OpenHurricane::real OpenHurricane::fallOffFunctions::SRI::DFDci(const real T, const real Pr,
                                                                      const real F,
                                                                      const real DPrDci) const {
    const real X = real(1.0) / (real(1) + sqr(log10(max(Pr, tiny))));
    const real dXdc = -X * X * real(2) * log10(Pr) / Pr / log(real(10.0)) * DPrDci;
    return (F * dXdc * log(a_ * exp(-b_ / T) + exp(-T / c_)));
}

hur_nodiscard OpenHurricane::real OpenHurricane::fallOffFunctions::SRI::DFDT(const real T, const real Pr,
                                                                     const real F,
                                                                     const real DPrDT) const {
    const real X = real(1.0) / (real(1) + sqr(log10(max(Pr, tiny))));
    const real dXdPr = -X * X * real(2) * log10(Pr) / Pr / log(real(10.0));
    return (F * (e_ / T + dXdPr * DPrDT * log(a_ * exp(-b_ / T) + exp(-T / c_)) +
                 X * (a_ * b_ * exp(-b_ / T) / sqr(T) - exp(-T / c_) / c_) /
                     (a_ * exp(-b_ / T) + exp(-T / c_))));
}

hur_nodiscard OpenHurricane::real OpenHurricane::fallOffFunctions::Troe::F(const real T,
                                                                   const real Pr) const {
    real logFcent;
    if (isTssOmitted_) {
        logFcent =
            log10(max((real(1.0) - alpha_) * exp(-T / Tsss_) + alpha_ * exp(-T / Ts_), tiny));
    } else {
        logFcent = log10(
            max((real(1.0) - alpha_) * exp(-T / Tsss_) + alpha_ * exp(-T / Ts_) + exp(-Tss_ / T),
                tiny));
    }

    real c = -0.4 - 0.67 * logFcent;
    real n = 0.75 - 1.27 * logFcent;
    const real d = 0.14;

    real logPr = log10(max(Pr, tiny));

    return std::pow(real(10.0), logFcent / (real(1.0) + sqr((logPr + c) / (n - d * (logPr + c)))));
}

hur_nodiscard OpenHurricane::real OpenHurricane::fallOffFunctions::Troe::DFDci(const real T, const real Pr,
                                                                       const real F,
                                                                       const real DPrDci) const {
    const real logPr = log10(max(Pr, tiny));

    real logFcent;
    if (isTssOmitted_) {
        logFcent =
            log10(max((real(1.0) - alpha_) * exp(-T / Tsss_) + alpha_ * exp(-T / Ts_), tiny));
    } else {
        logFcent = log10(
            max((real(1.0) - alpha_) * exp(-T / Tsss_) + alpha_ * exp(-T / Ts_) + exp(-Tss_ / T),
                tiny));
    }

    real c = -0.4 - 0.67 * logFcent;
    real n = 0.75 - 1.27 * logFcent;
    const real d = 0.14;

    const real logPrpc = logPr + c;
    const real dlogPrpc = d * (logPrpc);
    const real nmdlogPrpc = n - dlogPrpc;
    const real invnmdlogPrpc = logPrpc / nmdlogPrpc;

    return (F * ((-logFcent) / sqr(real(1) + sqr(invnmdlogPrpc))) *
            (real(2) * logPrpc / sqr(nmdlogPrpc)) * (real(1) + d * invnmdlogPrpc) *
            (real(1) / (Pr)*DPrDci));
}

hur_nodiscard OpenHurricane::real OpenHurricane::fallOffFunctions::Troe::DFDT(const real T, const real Pr,
                                                                      const real F,
                                                                      const real DPrDT) const {
    const real ln10 = log(real(10));
    const real logPr = log10(max(Pr, tiny));

    real Fcent;
    real DFcentDT;
    if (isTssOmitted_) {
        Fcent = max((real(1.0) - alpha_) * exp(-T / Tsss_) + alpha_ * exp(-T / Ts_), tiny);
        DFcentDT =
            (real(1) - alpha_) * exp(-T / Tsss_) / (-Tsss_) + alpha_ * exp(-T / Ts_) / (-Ts_);
    } else {
        Fcent = max(
            (real(1.0) - alpha_) * exp(-T / Tsss_) + alpha_ * exp(-T / Ts_) + exp(-Tss_ / T), tiny);
        DFcentDT = (real(1) - alpha_) * exp(-T / Tsss_) / (-Tsss_) +
                   alpha_ * exp(-T / Ts_) / (-Ts_) + exp(-Tss_ / T) * (Tss_ / sqr(T));
    }
    const real logFcent = log10(Fcent);
    const real DlogFcentDT = real(1) / (Fcent * ln10) * DFcentDT;
    const real c = -0.4 - 0.67 * logFcent;
    const real DcDT = -0.67 * DlogFcentDT;
    const real n = 0.75 - 1.27 * logFcent;
    const real DnDT = -1.27 * DlogFcentDT;
    const real d = 0.14;
    const real dlogPrDT = DPrDT / Pr / ln10;
    const real logPrpc = logPr + c;
    const real dlogPrpc = d * (logPrpc);
    const real nmdlogPrpc = n - dlogPrpc;
    const real invnmdlogPrpc = logPrpc / nmdlogPrpc;

    const real dPPDT = real(2) * logPrpc / sqr(nmdlogPrpc) *
                       ((dlogPrDT + DcDT) - logPrpc / nmdlogPrpc * (DnDT - d * (dlogPrDT + DcDT)));

    return (F * ln10 *
            (DlogFcentDT / (real(1) + sqr(invnmdlogPrpc)) -
             logFcent / sqr(real(1) + sqr(invnmdlogPrpc)) * dPPDT));
}