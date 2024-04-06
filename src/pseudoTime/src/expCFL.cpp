/*!
 * \file expCFL.cpp
 * \brief The subroutines and functions of exponential law adapt CFL number.
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

#include "expCFL.hpp"

namespace OpenHurricane {
    createClassName(expCFL);
    registerObjFty(CFL, expCFL, controller);
} // namespace OpenHurricane

OpenHurricane::real OpenHurricane::expCFL::getCFLFromStep(const integer cstep,
                                                          const bool isPrintCFL) const {
    if (cstep <= CFLConst_) {
        cCFL_ = CFLMin_;
    } else {
        cCFL_ = min(CFLMax_, pow(CFLMin_, cstep - CFLConst_));
    }
    if (isPrintCFL) {
        printCFL();
    }
    return cCFL_;
}

OpenHurricane::expCFL::expCFL(const iteration &iter, const runtimeMesh &mesh,
                              const controller &cont)
    : CFL(iter, mesh, cont), CFLConst_() {
    if (iter.hasSubIteration()) {
        CFLConst_ = cont.subController("expCFL").findOrDefault<integer>("cflConst", 5);
        CFLConst_ = min(CFLConst_, iter.subIter().minSubStep());
    } else {
        CFLConst_ = cont.subController("expCFL").findOrDefault<integer>("cflConst", 100);
    }
    CFL::CFLMin_ = cont.findOrDefault<real>("cflMin", 1.1);
    CFL::CFLMax_ = cont.findOrDefault<real>("cflMax", 10.0);
    if (CFL::CFLMax_ < CFL::CFLMin_) {
        Pout << "    Info: The maximum CFL: " << CFL::CFLMax_
             << " is smaller than the minmum CFL: " << CFL::CFLMin_
             << ". And it is reset to: " << CFL::CFLMin_ << std::endl;
        CFL::CFLMax_ = CFL::CFLMin_;
    }
    if (!cont.found("stepForPrintCFL")) {
        stepForPrintCFL_ = 1000;
    }
}

OpenHurricane::integer OpenHurricane::expCFL::CFLConst() const noexcept {
    return CFLConst_;
}

OpenHurricane::real OpenHurricane::expCFL::getCFL() const {
    if (iter().hasSubIteration()) {
        return getCFLFromStep(iter().subIter().cSubStep());
    } else {
        return getCFLFromStep(iter().cStep(), true);
    }
}
