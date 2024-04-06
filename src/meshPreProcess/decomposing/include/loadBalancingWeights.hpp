/*!
 * \file loadBalancingWeights.hpp
 * \brief Headers of load balancing weights.
 *        The subroutines and functions are in the <i>loadBalancingWeights.cpp</i> file.
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

#include "Lists.hpp"
#include "geometryMesh.hpp"
#include "logFile.hpp"
#include "originMeshRead.hpp"

namespace OpenHurricane {
    class loadBalancingWeights {
    public:
        enum class weightTypes : short { unweight, facesPerCell, givenInMeshFile };

    private:
        const originMeshRead &originMesh_;

        weightTypes weightType_;

        integer faceAdditionalWeights_;

    public:

        /*!\brief Construct from components.*/
        loadBalancingWeights(const originMeshRead &originMeshes, const controller &cont);

        /*!\brief Destructor.*/
        inline ~loadBalancingWeights() noexcept;

        hur_nodiscard inline weightTypes weightType() const noexcept;
        hur_nodiscard inline bool isUnWeight() const noexcept;
        hur_nodiscard inline bool isFacesPerCell() const noexcept;
        hur_nodiscard inline bool isGivenByMeshFile() const noexcept;

        hur_nodiscard integerArray weights() const;
    };
} // namespace OpenHurricane

#include "loadBalancingWeights.inl"