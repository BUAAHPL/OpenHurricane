/*!
 * \file solverKinds.hpp
 * \brief Header of kind of solver
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
#pragma once

#include "Map.hpp"
#include "string.hpp"

namespace OpenHurricane {
namespace Options {
enum solverKinds {
    NO_SOLVER = 0,     /*!< \brief Definition of no solver. */
    EULER = 1,         /*!< \brief Definition of Euler's solver. */
    NAVIER_STOKES = 2, /*!< \brief Definition of Navier-Stokes' solver. */
    RANS =
        3 /*!< \brief Definition of Reynolds-averaged Navier-Stokes' solver. */
};

static const std::map<std::string, solverKinds> solverKindsMap =
    createMap<std::string, solverKinds>("NO_SOLVER", NO_SOLVER)("EULER", EULER)(
        "NAVIER_STOKES", NAVIER_STOKES)("RANS", RANS);

} // End namespace Options
} // namespace OpenHurricane