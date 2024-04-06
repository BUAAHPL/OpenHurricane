#include "upwindSchemes.hpp"
/*!
 * \file upwindSchemes.inl
 * \brief In-Line subroutines of the <i>upwindSchemes.hpp</i> file.
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

inline void OpenHurricane::upwindSchemes::calcFlux(
    const real rhol, const real rhor, const vector &VL, const vector &VR, const real pl,
    const real pr, const real gl, const real gr, const real el, const real er, const real cl,
    const real cr, const vector &faceArea, const real blend, realArray &flux) const {
    upwindMethodPtr_->calcFlux(rhol, rhor, VL, VR, pl, pr, gl, gr, el, er, cl, cr, faceArea, blend,
                               flux);
}