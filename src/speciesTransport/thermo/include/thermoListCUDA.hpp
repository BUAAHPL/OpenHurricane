/*!
 * \file thermoListCUDA.hpp
 * \brief Header of thermo in CUDA platform.
 *       The subroutines and functions are in the <i>thermoListCUDA.cpp</i> file.
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

#ifdef CUDA_PARALLEL

#include "cudaStreams.hpp"
#include "kineticTheoryCUDA.hpp"
#include "nThreadsAndBlocks.hpp"

namespace OpenHurricane {
    namespace CUDAThermo {
        void calchai(const cu_real *__restrict__ hostYiTP, const cuChem::speciesTable &spcs,
                     const nThreadsAndBlocks &nTB, const cu_integer nsp,
                     const cu_integer nCells, cu_real *__restrict__ host_hai);

        void calchaiAsync(const cu_real *__restrict__ hostYiTP, cu2DArray<cu_real> &dYiTP,
                          const cuChem::speciesTable &spcs, const nThreadsAndBlocks &nTB,
                          const cu_integer nsp, const cu_integer nCells,
                          cu_real *__restrict__ host_hai, cu2DArray<cu_real> &d_hai,
                          const cudaStreams &streams);

        void calchaicp(const cu_real *__restrict__ hostYiTP, const cuChem::speciesTable &spcs,
                       const nThreadsAndBlocks &nTB, const cu_integer nsp,
                       const cu_integer nCells, cu_real *__restrict__ host_haicp);

        void calchaicpAsync(const cu_real *__restrict__ hostYiTP, cu2DArray<cu_real> &dYiTP,
                            const cuChem::speciesTable &spcs, const nThreadsAndBlocks &nTB,
                            const cu_integer nsp, const cu_integer nCells,
                            cu_real *__restrict__ host_haicp, cu2DArray<cu_real> &d_haicp,
                            const cudaStreams &streams);
    } // namespace CUDAThermo

} // namespace OpenHurricane

#ifdef Rgen
#undef Rgen
#endif // !Rgen

#endif // CUDA_PARALLEL
