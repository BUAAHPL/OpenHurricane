/*!
 * \file geometryModel.hpp
 * \brief Header of geometry model
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
#include <string>

namespace OpenHurricane {
    namespace Options {
        
        namespace ghostCellLayers {
            enum ghostCellLayer : short {
                NO_GHOST_CELL = 0,
                ONE_LAYER_GHOST_CELL = 1, // One layer of ghost cells.
                TWO_LAYER_GHOST_CELL = 2, // Two layers of ghost cells.
                THREE_LAYER_GHOST_CELL = 3
            };

            static std::map<short, ghostCellLayer> ghostCellLayerMap =
                createMap<short, ghostCellLayer>(0, NO_GHOST_CELL)(1, ONE_LAYER_GHOST_CELL)(
                    2, TWO_LAYER_GHOST_CELL)(3, THREE_LAYER_GHOST_CELL);

        } // namespace ghostCellLayers

        

    } // End namespace Options
} // namespace OpenHurricane