/*!
 * \file noSlipWall.cpp
 * \brief Main subroutines for noSlipWall.
 * \author Yang Hongzhen
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
#include "noSlipWall.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::noSlipWall::className_ = "noSlipWall";

namespace OpenHurricane {
    registerObjFty(vectorBoundary, noSlipWall, controller);
} // namespace OpenHurricane

OpenHurricane::noSlipWall::noSlipWall(const faceZone &fZ, geometryArray<vector, cellMesh> &gf,
                                            const controller &cont)
    : Base(fZ, gf, cont) {
    Base::setSpecified();
}

void OpenHurricane::noSlipWall::updateBoundary() {
    const faceList &fL = this->varField().mesh().faces();
    auto &v = this->varArray_;

    for (integer fi = boundaryZone_.firstIndex(); fi < boundaryZone_.lastIndex() + 1; fi++) {
        const integer cl = fL[fi].leftCell();
        const integer cr = fL[fi].rightCell();
        v[cr] = -v[cl];
    }
}
