/*!
 * \file surfaceIntegrals.cpp
 * \brief Main subroutines for base class of computing surface integrals.
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

#include "surfaceIntegrals.hpp"

namespace OpenHurricane {
    createClassName(surfaceIntegrals);
    createObjFty(surfaceIntegrals, controller);
} // namespace OpenHurricane

OpenHurricane::surfaceIntegrals::surfaceIntegrals(const iteration &iter, const runtimeMesh &mesh,
                                                  const controller &cont,
                                                  const integerList &zoneIdList)
    : iter_(iter), mesh_(mesh), zoneIdList_(zoneIdList) {}

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::surfaceIntegrals>
OpenHurricane::surfaceIntegrals::creator(const iteration &iter, const runtimeMesh &mesh,
                                         const controller &cont, const integerList &zd) {
    string reportType = cont.findWord("reportType");
    defineInObjCreator(surfaceIntegrals, reportType, controller, (iter, mesh, cont, zd));
}