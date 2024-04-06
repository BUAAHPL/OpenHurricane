/*!
 * \file variableHardSphereModel.cpp
 * \brief Main subroutines for transport properties by the variable hard sphere model.
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

#include "variableHardSphereModel.hpp"

namespace OpenHurricane {
    createClassNameStr(variableHardSphereModel, "variableHardSphereModel");
    registerObjFty(transport, variableHardSphereModel, controller);
} // namespace OpenHurricane

void OpenHurricane::variableHardSphereModel::calcMuRef() {
    const auto pi = constant::mathConstants::pi;
    const auto kB = constant::physicalConstant::k;
    muref_ = 15.0 * sqrt(pi * m_ * kB * Tref_) /
             (2.0 * pi * sqr(dref_) * (5.0 - 2.0 * omega_) * (7.0 - 2.0 * omega_));
}

OpenHurricane::variableHardSphereModel &
OpenHurricane::variableHardSphereModel::operator=(const variableHardSphereModel &tra) {
    if (this != std::addressof(tra)) {
        transport::operator=(tra);

        omega_ = tra.omega_;
        Tref_ = tra.Tref_;
        dref_ = tra.dref_;
        m_ = tra.m_;
        muref_ = tra.muref_;
    }
    return *this;
}