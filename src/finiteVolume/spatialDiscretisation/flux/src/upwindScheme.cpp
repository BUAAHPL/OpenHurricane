/*!
 * \file upwindScheme.cpp
 * \brief Main subroutines for upwind scheme.
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

#include "upwindScheme.hpp"

namespace OpenHurricane {
    createClassNameStr(upwindScheme, "upwindScheme");
    createObjFty(upwindScheme, controller);
} // namespace OpenHurricane

OpenHurricane::upwindScheme::upwindScheme(const controller &cont, const runtimeMesh &mesh,
                                          flowModel &flow)
    : mesh_(mesh), flows_(flow), flux_(5, Zero) {}

OpenHurricane::uniquePtr<OpenHurricane::upwindScheme>
OpenHurricane::upwindScheme::creator(const controller &cont, const runtimeMesh &mesh,
                                     flowModel &flow) {
    string spschemeType = cont.findWord(upwindScheme::className_);

    Pout << "    Info: setting upwind scheme: " << spschemeType << std::endl;
    defineInObjCreator(upwindScheme, static_cast<std::string>(spschemeType), controller,
                       (cont, mesh, flow));
}

OpenHurricane::upwindScheme::~upwindScheme() noexcept {}
