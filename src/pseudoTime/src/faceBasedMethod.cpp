/*!
 * \file faceBasedMethod.cpp
 * \brief Main subroutines for face-based method for computing time step.
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

#include "faceBasedMethod.hpp"
#include "cellArrays.hpp"

namespace OpenHurricane {
    createClassNameStr(faceBasedMethod, "faceBasedMethod");
    registerObjFty(pseudoTime, faceBasedMethod, controller);

} // namespace OpenHurricane
OpenHurricane::faceBasedMethod::faceBasedMethod(
    const controller &cont, const runtimeMesh &mesh, const iteration &iter, const flowModel &flows,
    realArray &dt, const cellRealArray &shockFactor, const faceVector2DArray &rai,
    const faceVector2DArray &rav, const integerArray &temperatureFlag,
    const integerArray &pressureFlag, const cellIntegerArray &CFLFlag)
    : pseudoTime(cont, mesh, iter, flows, dt, shockFactor, rai, rav, temperatureFlag, pressureFlag,
                 CFLFlag) {}

void OpenHurricane::faceBasedMethod::computingTimeStep() {
    computingTimeStep(dt_, cflPtr_->getCFL());
}

void OpenHurricane::faceBasedMethod::computingTimeStep(realArray &dt, const real cfl0) {
    const real C = CForTimeStep_;

    const faceList &faces = mesh_.faces();
    const auto &cells = mesh_.cells();
    const auto &cV = mesh_.cellVolume();

    for (integer n = 0; n < mesh_.nCells(); ++n) {
        real ra_sum = 0.0;
        for (integer i = 0; i < cells[n].faceSize(); ++i) {
            const integer j = cells[n].facei(i);
            const auto &cl = faces[j].leftCell();
            const auto &cr = faces[j].rightCell();

            real sl = fabs(real(cl) - real(n));
            real sr = fabs(real(cr) - real(n));

            ra_sum += (rai_[j][0] + C * rav_[j][0]) * (sr / max(real(1.0), sr)) +
                      (rai_[j][1] + C * rav_[j][1]) * (sl / max(real(1.0), sl));
        }
        dt[n] = 2.0 * cfl0 * cV[n] / ra_sum;
    }
}
