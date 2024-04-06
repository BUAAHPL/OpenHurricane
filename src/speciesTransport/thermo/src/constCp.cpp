/*!
 * \file constCp.cpp
 * \brief Main subroutines for constCp.
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
#include "constCp.hpp"

namespace OpenHurricane {
    createClassNameStr(constCp, "constCp");
    registerObjFty(thermo, constCp, controller);
} // namespace OpenHurricane


OpenHurricane::constCp::constCp(const controller &cont, const equationOfState &st, const integer id)
    : thermo(cont, st, id), cp_(cont.findOrDefault<real>("cp", real(1004.0))),
      hc_(cont.findOrDefault<real>("hc", real(0.0))) {}

OpenHurricane::constCp& OpenHurricane::constCp::operator=(const constCp &cC) {
    if (this != std::addressof(cC)) {
        thermo::operator=(cC);
        cp_ = cC.cp_;
        hc_ = cC.hc_;
    }
    return *this;
}
