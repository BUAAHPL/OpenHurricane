/*!
 * \file wallDistance.hpp
 * \brief Headers of class of wallDistance.
 *        The subroutines and functions are in the <i>wallDistance.cpp</i> file.
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
#include "distanceMethod.hpp"

namespace OpenHurricane {
    /*!\brief The class of wall distance.*/
    class wallDistance {
    private:
        /*!\brief Hold const-reference to the runtime mesh.*/
        const runtimeMesh &mesh_;

        /*!\brief Zone id list for distance calcuation.*/
        integerList zoneIdList_;

        /*!\brief The bc bcType of the zone.*/
        faceBCType::bcTypes zoneType_;

        /*!\brief The distance to wall.*/
        mutable cellRealArray wallDist_;

        /*!\brief Distance method pointer.*/
        mutable uniquePtr<distanceMethod> dmPtr_;

    public:
        wallDistance() = delete;

        wallDistance(const runtimeMesh &mesh, const controller &cont,
                     faceBCType::bcTypes tp = faceBCType::bcTypes::WALL);

        wallDistance(const wallDistance &) = delete;
        wallDistance &operator=(const wallDistance &) = delete;

        /*!\brief Destructor.*/
        inline ~wallDistance() noexcept { dmPtr_.clear(); }

        /*!\brief Zone id list for distance calcuation.*/
        hur_nodiscard inline const integerList &zoneIdList() const noexcept;

        /*!\brief The distance to wall.*/
        hur_nodiscard inline const cellRealArray &wallDist() const noexcept;
    };
} // namespace OpenHurricane

#include "wallDistance.inl"