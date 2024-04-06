/*!
 * \file RCMReordering.hpp
 * \brief Headers of the  reversed Cuthill-McKee ordering.
 *        The subroutines and functions are in the <i>RCMReordering.cpp</i> file.
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
#include "reordering.hpp"

namespace OpenHurricane {

class RCMReordering : public reordering {
private:
    // Private functions

    /*!\brief Reversed Cuthill-Mckee ordering.*/
    void RCM();

    /*!\brief Disallow null constructor.*/
    RCMReordering() = delete;

    /*!\brief Disallow construct as copy.*/
    RCMReordering(const RCMReordering &) = delete;

public:
    // Static data
    declareClassNames;

    /*!\brief Constructors.*/

    /*!\brief Construct from components.*/
    inline RCMReordering(decomposingMeshByMetis &dcpMesh)
        : reordering(dcpMesh) {}

    /*!\brief Destructor.*/
    virtual inline ~RCMReordering() noexcept {};

    /*!\brief Reordering.*/
    virtual inline void reorder();
};

} // namespace OpenHurricane

#include "RCMReordering.inl"