/*!
 * \file nThreadsAndBlocks.hpp
 * \brief Header of number of threads and blocks for launching CUDA kernel function.
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
     * \brief The class of number of threads and blocks for launching CUDA kernel function.
     */
    class nThreadsAndBlocks {
    private:
        dim3 blockSize_;
        dim3 gridSize_;

        size_t sharedMemSize_;

    public:
        // Constructors

        /**
         * \brief Null constructor.
         */
        inline nThreadsAndBlocks();

        /**
         * \brief Construct from the given number of species and of cells.
         */
        inline nThreadsAndBlocks(const cu_integer nsp, const cu_integer nCells);

        /**
         * \brief Destructor.
         */
        inline ~nThreadsAndBlocks() noexcept;

        /**
         * \brief Return the size of CUDA block.
         */
        hur_nodiscard inline const dim3 &blockSize() const noexcept;
        hur_nodiscard inline dim3 &blockSize() noexcept { return blockSize_; }

        /**
         * \brief Return the size of CUDA grid.
         */
        hur_nodiscard inline const dim3 &gridSize() const noexcept;
        hur_nodiscard inline dim3 &gridSize() noexcept { return gridSize_; }

        inline size_t sharedMemSize() const noexcept;
        inline size_t &sharedMemSize() noexcept;

        inline void setBlockAndGridSize(const cu_integer nsp, const cu_integer nCells);
        inline void setBlockAndGridSize64(const cu_integer nsp, const cu_integer nCells);

        template <typename Type> inline void setSharedMemSize1(const cu_integer num) noexcept;
    };
} // namespace OpenHurricane
#include "nThreadsAndBlocks.inl"

#endif // CUDA_PARALLEL