/*!
 * \file spaceSchemeOptions.hpp
 * \brief Header of space scheme options
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
namespace SpaceSchemes {
/*!
 * \brief types of spatial discretizations
 */
enum spaceType : short {
    SPACE_CENTERED =
        0,            /*!< \brief Space centered convective numerical method. */
    SPACE_UPWIND = 1, /*!< \brief Upwind convective numerical method. */
    FINITE_ELEMENT =
        2 /*!< \brief Finite element convective numerical method. */
};
static const std::map<std::string, spaceType> spaceMap =
    createMap<std::string, spaceType>("SPACE_CENTERED", SPACE_CENTERED)(
        "SPACE_UPWIND", SPACE_UPWIND)("FINITE_ELEMENT", FINITE_ELEMENT);

enum spatialScheme : short {
    HLLC = 0,
    Roe = 1,
    HLL = 2,
    HLLC_HLL = 3,
    Vanleer = 4,
    Steger_warming = 5,
    AUSM = 6,
    AUSMPWPLUS = 7,
    AUSMPLUSUP = 8,
    AUSMPUP2 = 9,
    JST = 10, /*!< \brief Jameson-Smith-Turkel centered numerical method. */
    LAX = 11, /*!< \brief Lax-Friedrich centered numerical method. */
    JST_KE =
        12 /*!< \brief Kinetic Energy preserving Jameson-Smith-Turkel centered numerical method. */
};

static const std::map<std::string, spatialScheme> spatialSchemeMap =
    createMap<std::string, spatialScheme>("HLLC", HLLC)("Roe", Roe)("HLL", HLL)(
        "HLLC_HLL", HLLC_HLL)("Vanleer", Vanleer)("Steger_warming",
                                                  Steger_warming)("AUSM", AUSM)(
        "AUSMPWPLUS", AUSMPWPLUS)("AUSMPLUSUP", AUSMPLUSUP)(
        "AUSMPUP2", AUSMPUP2)("JST", JST)("JST_KE", JST_KE)("LAX-FRIEDRICH",
                                                            LAX);

//enum centeredScheme :short
//{
//	JST = 0,   /*!< \brief Jameson-Smith-Turkel centered numerical method. */
//	LAX = 1,   /*!< \brief Lax-Friedrich centered numerical method. */
//	JST_KE = 2 /*!< \brief Kinetic Energy preserving Jameson-Smith-Turkel centered numerical method. */
//};

//static const std::map<std::string, centeredScheme> centeredSchemeMap =
//	createMap<std::string, centeredScheme>
//	("JST", JST)
//	("JST_KE", JST_KE)
//	("LAX-FRIEDRICH", LAX);

/*!
 * \brief types of FEM spatial discretizations
 */
enum FEM : short {
    NO_FEM = 0, /*!< \brief No finite element scheme is used. */
    DG = 1      /*!< \brief Discontinuous Galerkin numerical method. */
};
static const std::map<std::string, FEM> FEMMap =
    createMap<std::string, FEM>("NONE", NO_FEM)("DG", DG);
} // End namespace SpaceSchemes
} // End namespace Options
} // namespace OpenHurricane