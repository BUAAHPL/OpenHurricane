/*!
 * \file Lapacian.hpp
 * \brief Headers of the Lapacian.
 *        The subroutines and functions are in the <i>Lapacian.cpp</i> file.
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
#include "GaussDiv.hpp"
#include "cellArrays.hpp"
#include "controller.hpp"
#include "faceArrays.hpp"
#include "objectFactory.hpp"
#include "runtimeMesh.hpp"
#include "transfer.hpp"

namespace OpenHurricane {

    /*!\brief The class of Lapacian.*/
    class Lapacian {
    public:
        declareClassName(Lapacian);

        Lapacian() {}

        virtual ~Lapacian() noexcept {}

        /*!\brief Compute the lapacian of real array: cellQ.
         * \param[in] cellQ - The real array.
         * \param[out] lap - The lapacian of cellQ.
         */
        static void calcLapacian(const geometryArray<real, cellMesh> &cellQ,
                                 geometryArray<real, cellMesh> &lap);

        /*!\brief Compute the lapacian of vector array: cellQ.
         * \param[in] cellQ - The vector array.
         * \param[out] lap - The lapacian of cellQ.
         */
        static void calcLapacian(const geometryArray<vector, cellMesh> &cellQ,
                                 geometryArray<vector, cellMesh> &lap);

        /*!\brief Compute the lapacian of real array: cellQ.
         * \param[in] cellQ - The real array.
         * \return The lapacian of cellQ.
         */
        static geometryArray<real, cellMesh> lapacian(const geometryArray<real, cellMesh> &cellQ);

        /*!\brief Compute the lapacian of vector array: cellQ.
         * \param[in] cellQ - The vector array.
         * \return The lapacian of cellQ.
         */
        static geometryArray<vector, cellMesh>
        lapacian(const geometryArray<vector, cellMesh> &cellQ);

    protected:
    };

} // namespace OpenHurricane
