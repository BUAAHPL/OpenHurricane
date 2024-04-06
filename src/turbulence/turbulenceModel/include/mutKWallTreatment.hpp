/*!
 * \file mutKWallTreatment.hpp
 * \brief Headers of turbulent mut near-wall treatment.
 *        The subroutines and functions are in the <i>mutKWallTreatment.cpp</i> file.
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

#include "boundaries.hpp"
#include "cellArrays.hpp"

namespace OpenHurricane {
    /**\brief The class of omega near-wall treatment.*/
    class mutKWallTreatment : public realBoundary {
    public:
        using Base = realBoundary;

    private:
        real beta_;

        real Cmiu_;
        real kappa_;

        real E_;

        real yplusLam_;

    public:
        declareClassNames;

        static real yPlusLam(const real kappa, const real E);

        /**\brief Construct from components.*/
        mutKWallTreatment(const faceZone &fZ, geometryArray<real, cellMesh> &gf,
                          const controller &cont);

        /**\brief Construct from components.*/
        mutKWallTreatment(const mutKWallTreatment &bB);

        /*!\brief Destructor.*/
        virtual ~mutKWallTreatment() noexcept {}

        virtual void updateBoundary();
    };
} // namespace OpenHurricane
