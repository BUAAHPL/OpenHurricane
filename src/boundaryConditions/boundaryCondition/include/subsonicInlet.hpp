/*!
 * \file subsonicInlet.hpp
 * \brief Headers of the subsonicInlet.
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

    class subsonicInlet : public realBoundary {
    public:
        using Base = realBoundary;

    private:
        enum directionType : short { normalToBoundary, givenDirect };

        directionType directionType_;

        vector direct_;

        real p0_;

        real p_;

        real T0_;

    public:
        declareClassNames;

        subsonicInlet(const faceZone &fZ, geometryArray<real, cellMesh> &gf,
                       const controller &cont);

        /**\brief Construct from components.*/
        subsonicInlet(const subsonicInlet &bB)
            : Base(bB), directionType_(bB.directionType_), direct_(bB.direct_), p0_(bB.p0_),
              p_(bB.p_), T0_(bB.T0_) {}

        virtual ~subsonicInlet() noexcept {}

        virtual void updateBoundary();
    };

} // namespace OpenHurricane