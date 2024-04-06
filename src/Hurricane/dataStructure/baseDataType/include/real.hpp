/*!
 * \file real.hpp
 * \brief Header of real
 *       The subroutines and functions are in the <i>real.cpp</i> file.
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

#include "doubleReal.hpp"
#include "floatReal.hpp"
#include <string>

#if defined(HURRICANE_SP)

// Define real as a float

namespace OpenHurricane {
    using real = floatReal;

    constexpr hur_inline_var real large = floatRealLarge;
    constexpr hur_inline_var real veryLarge = floatRealVeryLarge;
    constexpr hur_inline_var real rootVeryLarge = floatRealRootVeryLarge;
    constexpr hur_inline_var real tiny = floatRealTiny;
    constexpr hur_inline_var real rootTiny = floatRealRootTiny;
    constexpr hur_inline_var real veryTiny = floatRealVeryTiny;
    constexpr hur_inline_var real rootVeryTiny = floatRealRootVeryTiny;

} // namespace OpenHurricane

#elif defined(HURRICANE_DP)

// Define real as a double

namespace OpenHurricane {
    using real = doubleReal;

    constexpr hur_inline_var real large = doubleRealLarge;
    constexpr hur_inline_var real veryLarge = doubleRealVeryLarge;
    constexpr hur_inline_var real rootVeryLarge = doubleRealRootVeryLarge;
    constexpr hur_inline_var real tiny = doubleRealTiny;
    constexpr hur_inline_var real rootTiny = doubleRealRootTiny;
    constexpr hur_inline_var real veryTiny = doubleRealVeryTiny;
    constexpr hur_inline_var real rootVeryTiny = doubleRealRootVeryTiny;

} // namespace OpenHurricane

#endif

namespace OpenHurricane {
    bool readReal(const char *_str, real &val);

    bool readReal(const std::string &str, real &val);
}