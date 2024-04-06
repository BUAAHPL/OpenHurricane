/*!
 * \file unsteadyTVDRK.cpp
 * \brief Main subroutines for unsteady TVDRK schemes.
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
#include "unsteadyTVDRK.hpp"
#include "solver.hpp"
namespace OpenHurricane {
    createClassNameStr(unsteadyTVDRK, "unsteadyTVDRK");
    registerObjFty(timeMarching, unsteadyTVDRK, controller);
} // namespace OpenHurricane

OpenHurricane::unsteadyTVDRK::unsteadyTVDRK(const controller &cont, const runtimeMesh &mesh,
                                            const flowModel &flowM, solver &_solver,
                                            cellVectorArray &v)
    : TVDRK(cont, mesh, flowM, _solver, v), unsteady(mesh, flowM, _solver) {
    // Set the time type to global time-step manually
    const_cast<iteration &>(iter()).setGlobalTimeStep();
}

void OpenHurricane::unsteadyTVDRK::timeStep() {
    if (iter().pTStep().isDynamicTimeStep()) {
        calcFluxSpectRadius();
        if (iter().pTStep().isDynamicCFLSet()) {
            const_cast<iteration &>(iter()).pTStep().setTimeStep(
                pseudoTimes().getGlobalTimeStep(iter().pTStep().dyCFL(), true));
        } else {
            const_cast<iteration &>(iter()).pTStep().setTimeStep(pseudoTimes().timeStep());
        }
    }

    getTimeStep(dt_);

    getShockFactor(shockFactor_);
}
