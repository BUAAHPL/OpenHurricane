#include "ODEsSolverCUDA.hpp"
/*!
 * \file ODEsSolverCUDA.cpp
 * \brief Main subroutines for the ODEs solver in CUDA.
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
#ifdef CUDA_PARALLEL
#include "ODEsSolverCUDA.hpp"
#include "cudaReduction.hpp"

cu_device cu_real OpenHurricane::ODEsSolverCUDA::maxError(
    CUDAChemistrySourceODEs &odes, const cu_real *__restrict__ y0,
    const cu_real *__restrict__ y, cu_real *__restrict__ e,
    unsigned int tid, const unsigned int dimX) const {
    auto j = tid;
    while (j < odes.nEqns()) {
        cu_real tolelence = ATOL_ + RTOL_ * cu_max(fabs(y0[j]), fabs(y[j]));
        e[j] = fabs(e[j]) / tolelence;
        j += dimX;
    }
    __syncthreads();
    cudaReduction::reduceMax(e, tid, odes.nEqns(), dimX);
    return e[0];
}

#endif // CUDA_PARALLEL