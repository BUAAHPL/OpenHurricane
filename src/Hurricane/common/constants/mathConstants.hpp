/*!
 * \file mathematicals.hpp
 * \brief All the information about the definition of the mathematicals constants.
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

#include "real.hpp"

namespace OpenHurricane {
    namespace constant {
        /**
         * \brief The namespace of mathematical conatsnt.
         */
        namespace mathConstants {
            /**\brief Natural base: e.*/
            constexpr hur_inline_var real e = 2.7182818284590452353602874;

            /**\brief Constant Pi.*/
            constexpr hur_inline_var real pi = 3.1415926535897932384626;
            constexpr hur_inline_var real twoPi(2.0 * pi);
            constexpr hur_inline_var real piBytwo(0.5 * pi);

            /**\brief Euler's constant.*/
            constexpr hur_inline_var real Eu(0.57721566490153286060651209);

            constexpr hur_inline_var real OneThird(1.0 / 3.0);
            constexpr hur_inline_var real TwoThird(2.0 / 3.0);
            constexpr hur_inline_var real FourThird(4.0 / 3.0);
        } // namespace mathConstants
    }     // namespace constant

} // namespace OpenHurricane