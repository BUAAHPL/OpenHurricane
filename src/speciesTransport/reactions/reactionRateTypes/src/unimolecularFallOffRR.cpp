/*!
 * \file unimolecularFallOffRR.hpp
 * \brief Main subroutines for unimolecular fall-off reaction rate.
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

#include "unimolecularFallOffRR.hpp"

namespace OpenHurricane {
    createClassName(unimolecularFallOffRR);
         registerObjFty(reactionRateTypes, unimolecularFallOffRR, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::unimolecularFallOffRR::k(const real p, const real T,
                                                                  const realArray &c) const {
    real k0 = k0_.k(p, T, c);
    real kinf = kinf_.k(p, T, c);
    real Pr = k0 / notZeroD(kinf) * thirdBodyEfficiency_.M(c);

    return kinf * (Pr / (1.0 + Pr)) * F_->F(T, Pr);
}

hur_nodiscard OpenHurricane::real OpenHurricane::unimolecularFallOffRR::DkDT(const real kj, const real p,
                                                                     const real T,
                                                                     const realArray &c) const {
    const real k0 = k0_.k(p, T, c);
    const real kinf = kinf_.k(p, T, c);
    const real Pr = k0 / notZeroD(kinf) * thirdBodyEfficiency_.M(c);

    return (Pr / (1.0 + Pr)) * F_->F(T, Pr) * kinf_.DkDT(kinf, p, T, c);
}

void OpenHurricane::unimolecularFallOffRR::gamDGamDCi(const real P, const real T, const realArray &c,
                                                  realArray &gdgdci) const {
    const real M = thirdBodyEfficiency_.M(c);

    if (M > tiny) {
        const real k0 = k0_.k(P, T, c);
        const real kinf = kinf_.k(P, T, c);
        const real Pr = k0 * M / notZeroD(kinf);
        const real F = F_->F(T, Pr);

        for (int i = 0; i < thirdBodyEfficiency_.size(); ++i) {
            const real dPrdci = thirdBodyEfficiency_[i] * k0 / notZeroD(kinf);
            const real dFdci = F_->DFDci(T, Pr, F, dPrdci);
            gdgdci[i] = (dPrdci / (Pr * (1 + Pr)) + dFdci / F);
        }
    } else {
        gdgdci = Zero;
    }
}

hur_nodiscard OpenHurricane::real
OpenHurricane::unimolecularFallOffRR::gamDGamDT(const real p, const real T, const realArray &c) const {
    const real M = thirdBodyEfficiency_.M(c);

    if (M > tiny) {
        const real k0 = k0_.k(p, T, c);
        const real kinf = kinf_.k(p, T, c);

        const real Pr = k0 * M / notZeroD(kinf);
        const real F = F_->F(T, Pr);
        const real dPrdT = Pr * (k0_.DkDT(k0, p, T, c) / notZeroD(k0) -
                                 kinf_.DkDT(kinf, p, T, c) / notZeroD(kinf));
        const real dFdT = F_->DFDT(T, Pr, F, dPrdT);

        return (dPrdT / (Pr * (1 + Pr)) + dFdT / F);
    } else {
        return 0;
    }
}

void OpenHurricane::unimolecularFallOffRR::resetThirdBodyEff(const thirdBodyEfficiency &tbe) {
    thirdBodyEfficiency_.resize(tbe.size());
    for (integer i = 0; i < tbe.size(); ++i) {
        thirdBodyEfficiency_[i] = tbe[i];
    }
}
