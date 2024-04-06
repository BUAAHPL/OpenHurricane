/*!
 * \file cellGradientOptions.hpp
 * \brief Header of cell gradient options
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
namespace Gradients {
enum cellGradient : short {
    greenGauss = 0,
    weightedLeastSquares = 1,
    hybridGreenGaussLeastSquares = 2
};

static const std::map<std::string, cellGradient> cellGradientMap =
    createMap<std::string, cellGradient>("greenGauss", greenGauss)(
        "weightedLeastSquares", weightedLeastSquares)(
        "hybridGreenGaussLeastSquares", hybridGreenGaussLeastSquares);

enum WLSQ : short { WLSQG = 0, WLSQ0 = 1, WLSQ1 = 2, WLSQ2 = 3, WLSQ3 = 4 };

static const std::map<std::string, WLSQ> WLSQMap =
    createMap<std::string, WLSQ>("WLSQG", WLSQG)("WLSQ0", WLSQ0)(
        "WLSQ1", WLSQ1)("WLSQ2", WLSQ2)("WLSQ3", WLSQ3);

enum neighbour : short {
    faceNeighbour = 0,
    nodeNeighbour = 1,
    twoLayersOfFaceNeighbour = 2
};

static const std::map<std::string, neighbour> neighbourMap =
    createMap<std::string, neighbour>("faceNeighbour", faceNeighbour)(
        "nodeNeighbour", nodeNeighbour)("twoLayersOfFaceNeighbour",
                                        twoLayersOfFaceNeighbour);
} // End namespace Gradients
} // End namespace Options
} // namespace OpenHurricane