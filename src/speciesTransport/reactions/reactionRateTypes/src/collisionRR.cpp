/*!
 * \file collisionRR.cpp
 * \brief Main subroutines for Collision Frequency Efficiency Expression Rate.
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

#include "collisionRR.hpp"
namespace OpenHurricane {
    createClassName(collisionRR);
         registerObjFty(reactionRateTypes, collisionRR, controller);
} // namespace OpenHurricane

void OpenHurricane::collisionFrequencies::calcAB() {
    real WA = species_[index1_].W();
    real WB = species_[index2_].W();
    WAB_ = (WA * WB) / (WA + WB);

    RAB_ = constant::physicalConstant::Ru / WAB_;
}

hur_nodiscard OpenHurricane::real OpenHurricane::collisionRR::k(const real p, const real T,
                                                        const realArray &c) const {
    return gammai(p, T, c) * collFre_.Zb(T);
}