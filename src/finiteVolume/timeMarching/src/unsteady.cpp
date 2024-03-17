/*!
 * \file unsteady.cpp
 * \brief Main subroutines for unsteady time schemes.
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
#include "unsteady.hpp"
#include "solver.hpp"

namespace OpenHurricane {
    createClassNameStr(unsteady, "unsteady");
}

OpenHurricane::unsteady::unsteady(const runtimeMesh &mesh, const flowModel &flowM, solver &_solver)
    : mesh_(mesh), flowM_(flowM), solverM_(_solver) {}

void OpenHurricane::unsteady::getTimeStep(realArray &dt) const {
    dt = solverM_.iter().pTStep().pTimeStep();
}

void OpenHurricane::unsteady::getShockFactor(cellRealArray &shockFactor_) const {
    const auto &cells = mesh_.cells();
    const auto &faces = mesh_.faces();
    auto &p = const_cast<flowModel &>(flowM_).p();
    for (integer n = 0; n < mesh_.nCells(); ++n) {
        real ff = 1.0;
        for (integer i = 0; i < cells[n].faceSize(); ++i) {
            const integer j = cells[n].facei(i);
            const auto &cl = faces[j].leftCell();
            const auto &cr = faces[j].rightCell();
            ff = min(ff, min(p[cl], p[cr]) / (fabs(p[cr] - p[cl]) + tiny));
        }
        shockFactor_[n] = ff;
    }

    realTransfer myTransfer(mesh_, shockFactor_, false, true);
    myTransfer.transferInit();
    myTransfer.transferring();
}
