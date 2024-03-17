/*!
 * \file vectors.cpp
 * \brief The subroutines and functions of vectors.
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

#include "vectors.hpp"

template <> const std::streamsize OpenHurricane::vector::vsType::precision(feature<real>::precision);

template <> const std::streamsize OpenHurricane::vector2D::vsType::precision(feature<real>::precision);

template <>
const std::streamsize OpenHurricane::integerVector::vsType::precision(feature<integer>::precision);

template <>
const std::streamsize OpenHurricane::integerVector2D::vsType::precision(feature<integer>::precision);

#ifdef HURRICANE_DP
#include "floatReal.hpp"

template <>
const std::streamsize OpenHurricane::floatVector::vsType::precision(feature<floatReal>::precision);

template <>
const std::streamsize OpenHurricane::floatVector2D::vsType::precision(feature<floatReal>::precision);

#endif // HURRICANE_SP

template <>
const std::streamsize OpenHurricane::complexVector::vsType::precision(feature<complex>::precision);

template <>
const std::streamsize OpenHurricane::complexVector2D::vsType::precision(feature<complex>::precision);
