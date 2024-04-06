/*!
 * \file ODEsCUDA.hpp
 * \brief Headers of the base class of ordinary differential equations in CUDA.
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
#ifdef CUDA_PARALLEL
#include "CUDAFunctions.hpp"

namespace OpenHurricane {
/**
 * \brief The base lass of a set of ODEs.
 *            dy/dt = f(t,y)
 */
class ODEsCUDA {

public:
    // Constructors

    /**\brief Construct null.*/
    cu_dual ODEsCUDA() {}

    /**\brief Construct null.*/
    cu_dual ODEsCUDA(const ODEsCUDA &) {}

    /**\brief Destructor.*/
    cu_dual virtual ~ODEsCUDA() noexcept {}

    // Member Functions
};

} // namespace OpenHurricane

#endif // CUDA_PARALLEL