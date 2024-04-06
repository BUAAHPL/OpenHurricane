/*!
 * \file CFL.cpp
 * \brief The subroutines and functions of CFL number.
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
#include "CFL.hpp"
#include "cellArrays.hpp"
#include "iteration.hpp"
#include "linearCFL.hpp"

namespace OpenHurricane {
    createClassName(CFL);
    createObjFty(CFL, controller);
} // namespace OpenHurricane

OpenHurricane::CFL::CFL(const iteration &iter, const runtimeMesh &mesh, const controller &cont)
    : iter_(iter), mesh_(mesh), CFLMax_(cont.findOrDefault<real>("cflMax", 10.0)),
      CFLMin_(cont.findOrDefault<real>("cflMin", 0.1)), cCFL_(1.0), notChange_(false),
      stepForPrintCFL_(cont.findOrDefault<integer>("stepForPrintCFL", 100)),
      initialStage_(cont.findOrDefault<integer>("initialStage", 100)) {
    if (CFLMax_ < CFLMin_) {
        Pout << "    Info: The maximum CFL: " << CFLMax_
             << " is smaller than the minmum CFL: " << CFLMin_
             << ". And it is reset to: " << CFLMin_ << std::endl;
        CFLMax_ = CFLMin_;
    }
    cCFL_ = CFLMax_;
}

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::CFL>
OpenHurricane::CFL::creator(const iteration &iter, const runtimeMesh &mesh,
                            const controller &cont) {
    string CFLType = cont.findWord("CFLType");
    Pout << "    Info: setting CFL type: " << CFLType << std::endl;
    defineInObjCreator(CFL, CFLType, controller, (iter, mesh, cont));
}

hur_nodiscard const OpenHurricane::iteration &OpenHurricane::CFL::iter() const noexcept {
    return iter_;
}

hur_nodiscard const OpenHurricane::runtimeMesh &OpenHurricane::CFL::mesh() const noexcept {
    return mesh_;
}

OpenHurricane::integer OpenHurricane::CFL::CFLConst() const noexcept {
    return 0;
}

OpenHurricane::integer OpenHurricane::CFL::CFLLinear() const noexcept {
    return 0;
}

void OpenHurricane::CFL::setRelativeCorrection(const real relCorre) const {}

void OpenHurricane::CFL::setBreakdowns() const {}

void OpenHurricane::CFL::setNotChange() const {
    notChange_ = true;
}

void OpenHurricane::CFL::printCFL() const {
    if (HurMPI::master()) {
        if (iter_.totalStep() % stepForPrintCFL_ == 0) {
            Pout << std::endl;
            Pout << "    Info: The global CFL number is " << cCFL_ << std::endl;
            Pout << std::endl;
        }
    }
}

hur_nodiscard bool OpenHurricane::CFL::isInitialStage() const noexcept {
    if (iter().isSteadyFlow()) {
        return iter().totalStep() <= initialStage_;
    } else {
        if (iter().hasSubIteration()) {
            return iter().subIter().cSubStep() <= iter().subIter().maxSubStep();
        }
    }
    return false;
}
