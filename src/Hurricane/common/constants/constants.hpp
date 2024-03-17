﻿/*!
 * \file constants.hpp
 * \brief Headers of all constants file.
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

#include "physicalConstant.hpp"

#include "mathConstants.hpp"

namespace OpenHurricane {
    /**\brief Conversion from atm to Pa.*/
   hur_nodiscard inline constexpr real atmToPa(const real atm) {
        return (atm * constant::physicalConstant::Patm);
    }

    /**\brief Conversion from atm to Pa.*/
    hur_nodiscard inline constexpr real paToAtm(const real pa) {
        return (pa / constant::physicalConstant::Patm);
    }

    /**\brief Conversion from degrees to radians.*/
    hur_nodiscard inline constexpr real degToRad(const real deg) {
        return (deg * constant::mathConstants::pi / real(180.0));
    }

    /**\brief Conversion from radians to degrees.*/
    hur_nodiscard inline constexpr real radToDeg(const real rad) {
        return (rad / constant::mathConstants::pi * real(180.0));
    }
} // namespace OpenHurricane