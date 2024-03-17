/*!
 * \file sutherlandThreeCoeffs.cpp
 * \brief Main subroutines for transport properties by Sutherland three coeffs.
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

#include "sutherlandThreeCoeffs.hpp"

namespace OpenHurricane {
    createClassNameStr(sutherlandThreeCoeffs, "sutherlandThree");
    registerObjFty(transport, sutherlandThreeCoeffs, controller);
} // namespace OpenHurricane

namespace OpenHurricane {}

OpenHurricane::sutherlandThreeCoeffs &
OpenHurricane::sutherlandThreeCoeffs::operator=(const sutherlandThreeCoeffs &tra) {
    if (this != std::addressof(tra)) {
        transport::operator=(tra);
        mu0_ = tra.mu0_;
        T0_ = tra.T0_;
        S_ = tra.S_;
    }
    return *this;
}