﻿/*!
 * \file limitersForMUSCL.hpp
 * \brief Header of limiter for MUSCL reconstruction.
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
#include "cellArrays.hpp"
#include "controller.hpp"
#include "faceArrays.hpp"
#include "objectFactory.hpp"
#include "smartPointer.hpp"

namespace OpenHurricane {
    /*!\brief The base class of limiters for MUSCL reconstruction.*/
    class limitersForMUSCL {
    public:
        declareClassNames;

        declareObjFty(limitersForMUSCL, controller, (const controller &cont), (cont));

        inline limitersForMUSCL() {}

        /*!\brief Select null constructed.*/
        static uniquePtr<limitersForMUSCL> creator(const controller &cont);

        /*!\brief Destructor.*/
        virtual ~limitersForMUSCL() noexcept {}

        /*!\brief Calculating the limiters for real.*/
        hur_nodiscard virtual real limiter(const real r) const = 0;

        /*!\brief Calculating the limiters for vector.*/
        hur_forceinline vector limiter(const vector r) const;
    };

} // namespace OpenHurricane

#include "limitersForMUSCL.inl"
