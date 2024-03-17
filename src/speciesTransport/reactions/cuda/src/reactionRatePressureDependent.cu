/*!
 * \file reactionRatePressureDependent.cu
 * \brief The subroutines and functions of pressure dependent reaction in CUDA platform.
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

#include "reactionRatePressureDependent.hpp"
#include <cmath>
#ifdef CUDA_PARALLEL

cu_dual cu_real
OpenHurricane::cuChem::pressureDependentCoefficients::fallOffFunctions::
    SRI(const cu_real a, const cu_real b, const cu_real c,
        const cu_real d, const cu_real e, const cu_real T,
        const cu_real Pr) {
    cu_real logPr = std::log10(fmax(Pr, cu_tiny));

    cu_real X = 1.0 / (1.0 + logPr * logPr);

    return d * std::pow(a * std::exp(-b / T) + std::exp(-T / c), X) *
           std::pow(T, e);
}

cu_dual cu_real OpenHurricane::cuChem::
    pressureDependentCoefficients::fallOffFunctions::Troe(const cu_real alpha,
                                                          const cu_real Tsss,
                                                          const cu_real Ts,
                                                          const cu_real T,
                                                          const cu_real Pr) {
    cu_real logFcent;
    logFcent = std::log10(fmax((cu_real(1) - alpha) * std::exp(-T / Tsss) +
                                   alpha * std::exp(-T / Ts),
                               cu_tiny));
    cu_real c = cu_real(-0.4) - cu_real(0.67) * logFcent;
    cu_real n = cu_real(0.75) - cu_real(1.27) * logFcent;
    const cu_real d = cu_real(0.14);

    cu_real logPr = std::log10(fmax(Pr, cu_tiny));
    cu_real tmp = ((logPr + c) / (n - d * (logPr + c)));
    return std::pow(cu_real(10), logFcent / (cu_real(1) + tmp * tmp));
}

cu_dual cu_real
OpenHurricane::cuChem::pressureDependentCoefficients::fallOffFunctions::
    Troe(const cu_real alpha, const cu_real Tsss, const cu_real Ts,
         const cu_real Tss, const cu_real T, const cu_real Pr) {
    cu_real logFcent;
    logFcent =
        std::log10(fmax((cu_real(1) - alpha) * std::exp(-T / Tsss) +
                            alpha * std::exp(-T / Ts) + std::exp(-Tss / T),
                        cu_tiny));
    cu_real c = cu_real(-0.4) - cu_real(0.67) * logFcent;
    cu_real n = cu_real(0.75) - cu_real(1.27) * logFcent;
    const cu_real d = cu_real(0.14);

    cu_real logPr = std::log10(fmax(Pr, cu_tiny));
    cu_real tmp = ((logPr + c) / (n - d * (logPr + c)));
    return std::pow(cu_real(10), logFcent / (cu_real(1) + tmp * tmp));
}

cu_device cu_real
OpenHurricane::cuChem::pressureDependentCoefficients::F(
    const cu_ushort ri, const cu_real T, const cu_real Pr) const {
    const auto rii = index_(ri);
    if (fallOffType_(rii) == Lindmann) {
        return fallOffFunctions::Lindemann(T, Pr);
    } else if (fallOffType_(rii) == TroeWithoutLastTerm) {
        const auto alpha = fallOffCoeff_(0, rii);
        const auto Tsss = fallOffCoeff_(1, rii);
        const auto Ts = fallOffCoeff_(2, rii);
        return fallOffFunctions::Troe(alpha, Tsss, Ts, T, Pr);
    } else if (fallOffType_(rii) == TroeWithLastTerm) {
        const auto alpha = fallOffCoeff_(0, rii);
        const auto Tsss = fallOffCoeff_(1, rii);
        const auto Ts = fallOffCoeff_(2, rii);
        const auto Tss = fallOffCoeff_(3, rii);
        return fallOffFunctions::Troe(alpha, Tsss, Ts, Tss, T, Pr);
    } else if (fallOffType_(rii) == SRI) {
        const auto a = fallOffCoeff_(0, rii);
        const auto b = fallOffCoeff_(1, rii);
        const auto c = fallOffCoeff_(2, rii);
        const auto d = fallOffCoeff_(3, rii);
        const auto e = fallOffCoeff_(4, rii);

        return fallOffFunctions::SRI(a, b, c, d, e, T, Pr);
    }
    return cu_real(1);
}

#endif // CUDA_PARALLEL