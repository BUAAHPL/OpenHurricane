/*!
 * \file Riemann.hpp
 * \brief Headers of the Riemann.
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
#include "boundaries.hpp"

namespace OpenHurricane {

    class Riemann : public realBoundary {
    public:
        using Base = realBoundary;

    private:
        /*!\brief Free stream velocity.*/
        vector vFree_;

        /*!\brief Free stream density.*/
        real rhoFree_;

        /*!\brief Free stream pressure.*/
        real pFree_;

        /*!\breif Free stream temperature.*/
        real tFree_;

        /*!\brief Free stream sonic speed.*/
        real cFree_;

        /*!\brief Free stream specific heat capicity.*/
        real gFree_;

        /*!\breif Free stream entropy.*/
        real sFree_;

        realArray yiFree_;

    public:
        declareClassNames;

        Riemann(const faceZone &fZ, geometryArray<real, cellMesh> &gf, const controller &cont);

        /**\brief Construct from components.*/
        Riemann(const Riemann &bB)
            : Base(bB), vFree_(bB.vFree_), rhoFree_(bB.rhoFree_), pFree_(bB.pFree_),
              tFree_(bB.tFree_), cFree_(bB.cFree_), gFree_(bB.gFree_), sFree_(bB.sFree_),
              yiFree_(bB.yiFree_) {}

        virtual ~Riemann() noexcept {}

        virtual void updateBoundary();
    };

} // namespace OpenHurricane