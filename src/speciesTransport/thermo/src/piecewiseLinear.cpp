/*!
 * \file piecewiseLinear.cpp
 * \brief Main subroutines for piecewise-linear.
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
#include "piecewiseLinear.hpp"

namespace OpenHurricane {
    createClassNameStr(piecewiseLinear, "piecewiseLinear");
}

namespace OpenHurricane {
    registerObjFty(thermo, piecewiseLinear, controller);
}

OpenHurricane::piecewiseLinear::piecewiseLinear(const controller &cont, const equationOfState &st,
                                            const integer id)
    : thermo(cont, st, id), cp_(), T_(), hc_(cont.findOrDefault<real>("hc", real(0.0))) {
    const auto cpl = cont.findTextStr("cp");
    cp_.resize(cpl.size());
    for (integer i = 0; i < cpl.size(); ++i) {
        std::stringstream sstr(cpl[i]);
        sstr >> cp_[i];
    }
    const auto Tl = cont.findTextStr("T");
    if (Tl.size() != cpl.size()) {
        errorAbortStr(("The size of T does not equal to cp in " + cont.name()));
    }
    T_.resize(cpl.size());
    for (integer i = 0; i < Tl.size(); ++i) {
        std::stringstream sstr(Tl[i]);
        sstr >> T_[i];
    }

    for (integer i = 0; i < T_.size() - 1; ++i) {
        if (T_[i] >= T_[i + 1]) {
            errorAbortStr(("The sort of T must be in ascending order in " + cont.name()));
        }
    }
}

OpenHurricane::piecewiseLinear &OpenHurricane::piecewiseLinear::operator=(const piecewiseLinear &cC) {
    if (this != std::addressof(cC)) {
        thermo::operator=(cC);
        cp_ = cC.cp_;
        T_ = cC.T_;
        hc_ = cC.hc_;
    }
    return *this;
}