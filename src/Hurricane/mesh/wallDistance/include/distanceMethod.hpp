/*!
 * \file distanceMethod.hpp
 * \brief Headers of class of distance method.
 *        The subroutines and functions are in the <i>distanceMethod.cpp</i> file.
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
#include "fVArraysInclude.hpp"

namespace OpenHurricane {
    /*!\brief The base class of distance method.*/
    class distanceMethod {
    protected:
        /*!\brief Hold const-reference to the runtime mesh.*/
        const runtimeMesh &mesh_;

        /*!\brief Zone id list for distance calcuation.*/
        const integerList zoneIdList_;

    public:
        declareClassNames;
        declareObjFty(distanceMethod, controller,
                      (const controller &cont, const runtimeMesh &mesh, const integerList zoneIdL),
                      (cont, mesh, zoneIdL));

        distanceMethod() = delete;

        distanceMethod(const controller &cont, const runtimeMesh &mesh, const integerList zoneIdL);

        distanceMethod(const runtimeMesh &mesh, const integerList zoneIdL);

        distanceMethod(const distanceMethod &) = delete;
        distanceMethod &operator=(const distanceMethod &) = delete;

        /*!\brief Selector.*/
        static uniquePtr<distanceMethod> creator(const controller &cont, const runtimeMesh &mesh,
                                                 const integerList zoneIdL);

        /*!\brief Destructor.*/
        inline virtual ~distanceMethod() noexcept {}

        /*!\brief Return const reference to the zone id list.*/
        hur_nodiscard inline const integerList &zoneIdList() const noexcept;

        hur_nodiscard inline virtual bool movePoints() noexcept;

        virtual bool getDistance(cellRealArray &y) = 0;

        virtual bool getDistance(cellRealArray &y, cellVectorArray &ny) = 0;
    };
} // namespace OpenHurricane

#include "distanceMethod.inl"