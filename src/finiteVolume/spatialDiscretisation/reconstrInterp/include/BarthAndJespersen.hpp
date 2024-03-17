/*!
 * \file BarthAndJespersen.hpp
 * \brief Header of Barth and Jespersen limiter for piecewise linear reconstruction.
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

#include "limitersForLinear.hpp"

namespace OpenHurricane {
    /*!\brief The class of Barth and Jespersen limiters for piecewise linear reconstruction.*/
    class BarthAndJespersen : public limitersForLinear {
    private:
    public:
        declareClassNames;

        /*!\brief Construct from controller.*/
        inline BarthAndJespersen(const controller &cont) {}

        /*!\brief Disallow copy constructor.*/
        BarthAndJespersen(const BarthAndJespersen &) = delete;

        /*!\brief Disallow bitwise assignment.*/
        BarthAndJespersen &operator=(const BarthAndJespersen &) = delete;

        /*!\brief Destructor.*/
        virtual inline ~BarthAndJespersen() noexcept {}

        /*!
         * \brief Calculating the limiters for real cellQ.
         * \param[in] cellQ - The variable stored in cell-centroid.
         * \param[in] gradQ - The gradient of cellQ stored in cell-centroid.
         * \param[out] limiters - The limiters for cellQ.
         */
        virtual void
        calcLimiters(const geometryArray<real, cellMesh> &cellQ,
                     const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &grad,
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
