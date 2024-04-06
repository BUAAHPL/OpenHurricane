/*!
 * \file Lapacians.hpp
 * \brief Headers of the Lapacians.
 *        The subroutines and functions are in the <i>Lapacians.inl</i> file.
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

#include "Lapacian.hpp"

namespace OpenHurricane {

    namespace fv {
        /*!\brief Compute the lapacian of array: cellQ.
         * \param[in] cellQ - The array.
         * \param[out] lap - The lapacian of cellQ.
         */
        template <class Type>
        static void lapacian(const geometryArray<Type, cellMesh> &cellQ,
                             geometryArray<Type, cellMesh> &lap);

        /*!\brief Compute the lapacian of  array: cellQ.
         * \param[in] cellQ - The array.
         * \return The lapacian of cellQ.
         */
        template <class Type>
        static geometryArray<Type, cellMesh> lapacian(const geometryArray<Type, cellMesh> &cellQ);

    } // namespace fv

} // namespace OpenHurricane

#include "Lapacians.inl"