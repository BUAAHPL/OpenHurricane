/*!
 * \file linearCFL.cpp
 * \brief The subroutines and functions of linear adapt CFL number.
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

#include "linearCFL.hpp"

namespace OpenHurricane {
    createClassName(linearCFL);
    registerObjFty(CFL, linearCFL, controller);
} // namespace OpenHurricane

OpenHurricane::real OpenHurricane::linearCFL::getCFLFromStep(const integer cstep,
                                                             const bool isPrintCFL) const {
    if (cstep <= CFLConst_) {
        cCFL_ = CFL::CFLMin_;
    } else if (cstep <= CFLLinear_) {
        cCFL_ = CFL::CFLMin_ + 0.9 * CFL::CFLMax_ * fabs(real(cstep) - real(CFLConst_)) /
                                   max(real(1.0), fabs(real(CFLLinear_) - real(CFLConst_)));
        cCFL_ = min(cCFL_, CFLMax_);
    } else {
        cCFL_ = CFL::CFLMax_;
    }
    if (isPrintCFL) {
        printCFL();
    }
    return cCFL_;
}

OpenHurricane::linearCFL::linearCFL(const iteration &iter, const runtimeMesh &mesh,
                                    const controller &cont)
    : CFL(iter, mesh, cont), CFLConst_(), CFLLinear_() {
    if (iter.hasSubIteration()) {
        CFLConst_ = cont.subController("linearCFL").findOrDefault<integer>("cflConst", 5);
        CFLLinear_ = cont.subController("linearCFL").findOrDefault<integer>("cflLinear", 10);
        CFLLinear_ = min(CFLLinear_, iter.subIter().minSubStep());
        CFLConst_ = min(CFLConst_, CFLLinear_);
    } else {
        CFLConst_ = cont.subController("linearCFL").findOrDefault<integer>("cflConst", 100);
        CFLLinear_ = cont.subController("linearCFL").findOrDefault<integer>("cflLinear", 200);
    }
    if (cont.found("cfl")) {
        real cfl = cont.findOrDefault<real>("cfl", 10.0);
        CFL::CFLMax_ = cfl;
        CFL::CFLMin_ = 0.1 * cfl;
    }

    if (!cont.found("stepForPrintCFL")) {
        stepForPrintCFL_ = 1000;
    }
}

OpenHurricane::integer OpenHurricane::linearCFL::CFLConst() const noexcept {
    return CFLConst_;
}

OpenHurricane::integer OpenHurricane::linearCFL::CFLLinear() const noexcept {
    return CFLLinear_;
}

OpenHurricane::real OpenHurricane::linearCFL::getCFL() const {
    if (iter().hasSubIteration()) {
        return getCFLFromStep(iter().subIter().cSubStep());
    } else {
        return getCFLFromStep(iter().cStep(), true);
    }
}
