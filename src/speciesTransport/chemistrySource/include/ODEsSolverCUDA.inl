#include "ODEsSolverCUDA.hpp"
/*!
 * \file ODEsSolverCUDA.inl
 * \brief In-Line subroutines of the <i>ODEsSolverCUDA.hpp</i> file.
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

inline cu_device cu_integer
OpenHurricane::ODEsSolverCUDA::maxSteps() const noexcept {
    return maxSteps_;
}

inline cu_device cu_real
OpenHurricane::ODEsSolverCUDA::ATOL() const noexcept {
    return ATOL_;
}

inline cu_device cu_real
OpenHurricane::ODEsSolverCUDA::RTOL() const noexcept {
    return RTOL_;
}

inline cu_host OpenHurricane::ODEsSolverCUDA::ODEsSolverCUDA(
    const cu_real atol, const cu_real rtol, const unsigned int maxStep)
    : maxSteps_(maxStep), ATOL_(atol), RTOL_(rtol) {}

inline cu_dual
OpenHurricane::ODEsSolverCUDA::ODEsSolverCUDA(const ODEsSolverCUDA &os)
    : maxSteps_(os.maxSteps_), ATOL_(os.ATOL_), RTOL_(os.RTOL_) {}
#endif // CUDA_PARALLEL