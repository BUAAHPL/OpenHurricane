/*!
 * \file divs.hpp
 * \brief Headers of the divergence.
 *        The subroutines and functions are in the <i>divs.inl</i> file.
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

namespace OpenHurricane {
    namespace fv {
        /*!\brief Compute divergence of vector array: cellQ, with Gauss theorem.
         * \param[in] cellQ - The vector array.
         * \param[out] div - The divergence of cellQ.
         */
        static void div(const geometryArray<vector, cellMesh> &cellQ,
                        geometryArray<real, cellMesh> &divQ);

        /*!\brief Compute divergence of tensor array: cellQ, with Gauss theorem.
         * \param[in] cellQ - The tensor array.
         * \param[out] div - The divergence of cellQ.
         */
        static void div(const geometryArray<tensor, cellMesh> &cellQ,
                        geometryArray<vector, cellMesh> &divs);

        /*!\brief Compute divergence of vector array: cellQ, with Gauss theorem.
         * \param[in] cellQ - The vector array.
         * \param[in] faceQ - The vector array interpolated into the face center.
         * \param[out] div - The divergence of cellQ.
         */
        static void div(const geometryArray<vector, cellMesh> &cellQ,
                        const geometryArray<vector, faceMesh> &faceQ,
                        geometryArray<real, cellMesh> &divs);

        /*!\brief Compute divergence of tensor array: cellQ, with Gauss theorem.
         * \param[in] cellQ - The tensor array.
         * \param[in] faceQ - The tensor array interpolated into the face center.
         * \param[out] div - The divergence of cellQ.
         */
        static void div(const geometryArray<tensor, cellMesh> &cellQ,
                        const geometryArray<tensor, faceMesh> &faceQ,
                        geometryArray<vector, cellMesh> &divs);

        /*!\brief Compute divergence of vector array: cellQ, with Gauss theorem.
         * \param[in] cellQ - The vector array.
         * \return The divergence of cellQ.
         */
        static geometryArray<real, cellMesh> div(const geometryArray<vector, cellMesh> &cellQ);

        /*!\brief Compute divergence of tensor array: cellQ, with Gauss theorem.
         * \param[in] cellQ - The tensor array.
         * \return The divergence of cellQ.
         */
        static geometryArray<vector, cellMesh> div(const geometryArray<tensor, cellMesh> &cellQ);

        /*!\brief Compute divergence of vector array: cellQ, with Gauss theorem.
         * \param[in] cellQ - The vector array.
         * \param[in] faceQ - The vector array interpolated into the face center.
         * \return The divergence of cellQ.
         */
        static geometryArray<real, cellMesh> div(const geometryArray<vector, cellMesh> &cellQ,
                                                 const geometryArray<vector, faceMesh> &faceQ);

        /*!\brief Compute divergence of tensor array: cellQ, with Gauss theorem.
         * \param[in] cellQ - The tensor array.
         * \param[in] faceQ - The tensor array interpolated into the face center.
         * \return The divergence of cellQ.
         */
        static geometryArray<vector, cellMesh> div(const geometryArray<tensor, cellMesh> &cellQ,
                                                   const geometryArray<tensor, faceMesh> &faceQ);
    } // namespace fv

} // namespace OpenHurricane

#include "divs.inl"