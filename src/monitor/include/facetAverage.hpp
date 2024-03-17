﻿/*!
 * \file facetAverage.hpp
 * \brief Headers of class of the facet average of a specified field variable on a surface
 *        The subroutines and functions are in the <i>facetAverage.cpp</i> file.
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
    /*!
     * \brief The class of the facet average of a specified field variable on a surface.
     */
    class facetAverage : public surfaceIntegrals {
    public:
        declareClassName(facetAverage);

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        facetAverage(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                     const integerList &zd);

        /**
         * \brief Destructor.
         */
        virtual ~facetAverage() noexcept {}

        hur_nodiscard virtual real surIntegral(const cellRealArray &phi) const;
        hur_nodiscard virtual vector surIntegral(const cellVectorArray &phi) const;

        hur_nodiscard inline virtual std::string printInfo() const {
            return std::string("facet average");
        }
    };
} // namespace OpenHurricane
