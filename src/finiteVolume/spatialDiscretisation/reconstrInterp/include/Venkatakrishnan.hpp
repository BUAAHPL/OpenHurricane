/*!
 * \file Venkatakrishnan.hpp
 * \brief Header of Venkatakrishnan's limiter for piecewise linear reconstruction.
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

#include "limitersForLinear.hpp"

namespace OpenHurricane {
    /*!\brief The class of Venkatakrishnan's limiters for piecewise linear reconstruction.*/
    class Venkatakrishnan : public limitersForLinear {
    private:

        /*!\brief The constant coefficient of O(1) for Venkatakrishnan's limiters*/
        real K_;

        /*!\brief K_^3*/
        real K3_;

        /*!\brief Disallow copy constructor.*/
        Venkatakrishnan(const Venkatakrishnan &) = delete;

        /*!\brief Disallow bitwise assignment.*/
        void operator=(const Venkatakrishnan &) = delete;

    public:
        declareClassNames;

        /*!\brief Construct from controller.*/
        Venkatakrishnan(const controller &cont);

        /*!\brief Destructor.*/
        virtual ~Venkatakrishnan() noexcept {}

        /*!
         * \brief Calculating the limiters for real cellQ.
         * \param[in] cellQ - The variable stored in cell-centroid.
         * \param[in] gradQ - The gradient of cellQ stored in cell-centroid.
         * \param[out] limiters - The limiters for cellQ.
         */
        virtual void calcLimiters(
            const geometryArray<real, cellMesh> &cellQ,
            const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &gradQ,
            geometryArray<real, cellMesh> &limiters) const;

        /*!
         * \brief Calculating the limiters for vector cellQ.
         * \param[in] cellQ - The variable stored in cell-centroid.
         * \param[in] gradQ - The gradient of cellQ stored in cell-centroid.
         * \param[out] limiters - The limiters for cellQ.
         */
        virtual void calcLimiters(
            const geometryArray<vector, cellMesh> &cellQ,
            const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &grad,
            geometryArray<vector, cellMesh> &limiters) const;
    };

} // namespace OpenHurricane

#include "Venkatakrishnan.inl"
