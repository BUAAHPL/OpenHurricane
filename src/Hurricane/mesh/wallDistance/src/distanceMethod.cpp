/*!
 * \file distanceMethod.cpp
 * \brief Main subroutines for distance method.
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

#include "distanceMethod.hpp"

namespace OpenHurricane {
    createClassNameStr(distanceMethod, "distanceMethod");
}

namespace OpenHurricane {
    createObjFty(distanceMethod, controller);
}

OpenHurricane::distanceMethod::distanceMethod(const controller &cont, const runtimeMesh &mesh,
                                          const integerList zoneIdL)
    : mesh_(mesh), zoneIdList_(zoneIdL) {}

OpenHurricane::distanceMethod::distanceMethod(const runtimeMesh &mesh, const integerList zoneIdL)
    : mesh_(mesh), zoneIdList_(zoneIdL) {}

OpenHurricane::uniquePtr<OpenHurricane::distanceMethod>
OpenHurricane::distanceMethod::creator(const controller &cont, const runtimeMesh &mesh,
                                   const integerList zoneIdL) {
    string distanceMethodType = cont.findWord(distanceMethod::className_);
#ifdef HUR_DEBUG
    PLInfo("    Info: creating %s distance method...\n", distanceMethodType.c_str());
#endif // HUR_DEBUG
    defineInObjCreator(distanceMethod, distanceMethodType, controller, (cont, mesh, zoneIdL));
}