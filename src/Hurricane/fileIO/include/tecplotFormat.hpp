/*!
 * \file tecplotFormat.hpp
 * \brief Headers of tecplotFormat.
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

#pragma once
#include "stdMaps.hpp"

namespace OpenHurricane {
    namespace tecplotFormat {
        enum FACENEIGHBORMODE : short {
            LOCALONETOONE = 0,
            LOCALONETOMANY = 1,
            GLOBALONETOONE = 2,
            GLOBALONETOMANY = 3,
            NONENEIGHBOUR
        };

        static std::map<short, std::string> FACENEIGHBORMODEMap = createMap<short, std::string>(
            LOCALONETOONE, "LOCALONETOONE")(LOCALONETOMANY, "LOCALONETOMANY")(
            GLOBALONETOONE, "GLOBALONETOONE")(GLOBALONETOMANY, "GLOBALONETOMANY");

        enum ZONETYPE : short {
            FELINESEG = 0,
            FETRIANGLE = 1,
            FEQUADRILATERAL = 2,
            FEPOLYGON = 3,
            FETETRAHEDRON,
            FEBRICK,
            FEPOLYHEDRAL
        };

        static std::map<short, std::string> ZONETYPEMap =
            createMap<short, std::string>(FELINESEG, "FELINESEG")(FETRIANGLE, "FETRIANGLE")(
                FEQUADRILATERAL, "FEQUADRILATERAL")(FEPOLYGON, "FEPOLYGON")(
                FETETRAHEDRON, "FETETRAHEDRON")(FEBRICK, "FEBRICK")(FEPOLYHEDRAL, "FEPOLYHEDRAL");

        enum FILETYPE : short { FULL = 0, GRID = 1, SOLUTION = 2 };

    } // namespace tecplotFormat
} // namespace OpenHurricane