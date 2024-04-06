/*!
 * \file drag.cpp
 * \brief Main subroutines for drag.
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

#include "drag.hpp"
#include "fVArraysInclude.hpp"

namespace OpenHurricane {
    createClassName(drag);
    registerObjFty(forces, drag, controller);
} // namespace OpenHurricane

void OpenHurricane::drag::setReportName() {
    integer count = Zero;
    reportName_[count++] = "Drag";
    if (reportCoordinateComponent_) {
        reportName_[count++] = "Dx";
        reportName_[count++] = "Dy";
        reportName_[count++] = "Dz";
    }
    if (reportPressureComponent_) {
        reportName_[count++] = "Dp";
        if (reportCoordinateComponent_) {
            reportName_[count++] = "Dpx";
            reportName_[count++] = "Dpy";
            reportName_[count++] = "Dpz";
        }
    }
    if (reportViscousComponent_) {
        reportName_[count++] = "Dv";
        if (reportCoordinateComponent_) {
            reportName_[count++] = "Dvx";
            reportName_[count++] = "Dvy";
            reportName_[count++] = "Dvz";
        }
    }
}

OpenHurricane::drag::drag(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                          const integerList &zd)
    : force(iter, mesh, cont, zd) {
    setReportName();
}
