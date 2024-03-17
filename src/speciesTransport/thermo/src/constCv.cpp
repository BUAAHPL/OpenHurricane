﻿/*!
 * \file constCv.cpp
 * \brief Main subroutines for constCv.
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
#include "constCv.hpp"

namespace OpenHurricane {
    createClassNameStr(constCv, "constCv");
    registerObjFty(thermo, constCv, controller);
} // namespace OpenHurricane

namespace OpenHurricane {}

OpenHurricane::constCv::constCv(const controller &cont, const equationOfState &st, const integer id)
    : thermo(cont, st, id), cv_(cont.findType<real>("cv", real(0.0))),
      hc_(cont.findOrDefault<real>("hc", real(0.0))) {}

OpenHurricane::constCv &OpenHurricane::constCv::operator=(const constCv &cC) {
    if (this != std::addressof(cC)) {
        thermo::operator=(cC);
        cv_ = cC.cv_;
        hc_ = cC.hc_;
    }
    return *this;
}
