/*!
 * \file relayScatter.hpp
 * \brief Headers of the relayScatter.
 *        The subroutines and functions are in the <i>relayScatter.cpp</i> file.
 * \author Yang Hongzhen
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
#include "cellMesh.hpp"
#include "geometryArrays.hpp"

namespace OpenHurricane {
    namespace relayScatterFunc {
        template <class Type>
        void relayScatter(const integerArrayArray &cellOfProc, const Array<Type> &allArray,
                          Array<Type> &ipArray);

        template <>
        void relayScatter(const integerArrayArray &cellOfProc, const Array<real> &allArray,
                          Array<real> &ipArray);

        template <class Type>
        void relayReorder(const Array<Type> &pArray, Array<Type> &ipArray, const integerList &perm_,
                          const integerList &iperm_);
    } // namespace relayScatterFunc
} // namespace OpenHurricane
#include "relayScatter.inl"