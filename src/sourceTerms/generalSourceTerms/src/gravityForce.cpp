/*!
 * \file gravityForce.cpp
 * \brief Main subroutines for gravity force source terms.
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

#include "gravityForce.hpp"
 namespace OpenHurricane{
	createClassNameStr(gravityForce,"gravityForce");
}

namespace OpenHurricane {
     registerObjFty(sourceTerm,gravityForce,controller);
}

OpenHurricane::gravityForce::gravityForce(const flowModel &flows, const iteration &iter,
                                      const controller &cont)
    : generalSourceTerms(flows, iter, cont),
      g_(cont.findOrDefault<vector>("gravity",
                                    vector(Zero, -constant::physicalConstant::gn, Zero))) {}

void OpenHurricane::gravityForce::addSourceTerms(const cellRealArray &rho, cellVectorArray &U) const {
    const auto &mesh = rho.mesh();

    const auto &cV = mesh.cellVol();
    for (integer n = 0; n < mesh.nCells(); ++n) {
        U.rhs()[n] += rho[n] * cV[n] * g_;
    }
}
