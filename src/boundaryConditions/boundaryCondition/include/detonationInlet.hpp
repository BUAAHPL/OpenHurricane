/*!
 * \file pressureInlet.hpp
 * \brief Headers of the detonationInlet.
 * \author Chen Zhenyi
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

namespace OpenHurricane {

    class detonationInlet : public vectorBoundary {
    public:
        using Base = vectorBoundary;

    private:
        enum directionType : short { normalToBoundary, givenDirect };

        directionType directionType_;

        enum class inletStateType : short { blocking, subsonic, supersonic };

        List<inletStateType> inletType_;

        vector direct_;

        /**\brief Total temperature.*/
        real Tt_;

        /**\breif total pressure.*/
        real Pt_;

        realArray yiBnd_;
                
        void getValue(real &value, const std::string& name, const controller &cont)const;
    public:
        declareClassNames;

        detonationInlet(const faceZone &fZ, geometryArray<vector, cellMesh> &gf,
                      const controller &cont);

        /**\brief Construct from components.*/
        detonationInlet(const detonationInlet &bB)
            : Base(bB), directionType_(bB.directionType_), direct_(bB.direct_), Tt_(bB.Tt_),
              Pt_(bB.Pt_), yiBnd_(bB.yiBnd_) {}

        virtual ~detonationInlet() noexcept {}

        virtual void updateBoundary();
    };

} // namespace OpenHurricane