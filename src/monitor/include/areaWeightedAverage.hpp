﻿/*!
 * \file areaWeightedAverage.hpp
 * \brief Headers of class of computing the area-weighted average of a quantity through a surface.
 *        The subroutines and functions are in the <i>areaWeightedAverage.cpp</i> file.
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
#include "surfaceIntegrals.hpp"

namespace OpenHurricane {
    /*!\brief The class of area of computing the area-weighted average of a quantity through a surface.*/
    class areaWeightedAverage : public surfaceIntegrals {
    public:
        declareClassName(areaWeightedAverage);

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        areaWeightedAverage(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                            const integerList &zd);

        /**
         * \brief Destructor.
         */
        virtual ~areaWeightedAverage() noexcept {}

        hur_nodiscard virtual real surIntegral(const cellRealArray &phi) const;
        hur_nodiscard virtual vector surIntegral(const cellVectorArray &phi) const;

        hur_nodiscard inline virtual std::string printInfo() const {
            return std::string("area-weighted average");
        }
    };
} // namespace OpenHurricane