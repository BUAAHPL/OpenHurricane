#include "nThreadsAndBlocks.hpp"
/*!
 * \file nThreadsAndBlocks.inl
 * \brief The In-Line functions of the <i>nThreadsAndBlocks.hpp</i> file.
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

inline OpenHurricane::nThreadsAndBlocks::nThreadsAndBlocks()
    : blockSize_(), gridSize_(), sharedMemSize_() {}

inline OpenHurricane::nThreadsAndBlocks::nThreadsAndBlocks(const cu_integer nsp,
                                                       const cu_integer nCells)
    : blockSize_(), gridSize_(), sharedMemSize_() {
    setBlockAndGridSize(nsp, nCells);

    setSharedMemSize1<cu_real>(3 * nsp * blockSize_.y);
}

inline OpenHurricane::nThreadsAndBlocks::~nThreadsAndBlocks() noexcept {}

hur_nodiscard inline const dim3 &OpenHurricane::nThreadsAndBlocks::blockSize() const noexcept {
    return blockSize_;
}

hur_nodiscard inline const dim3 &OpenHurricane::nThreadsAndBlocks::gridSize() const noexcept {
    return gridSize_;
}
inline size_t OpenHurricane::nThreadsAndBlocks::sharedMemSize() const noexcept {
    return sharedMemSize_;
}

inline size_t &OpenHurricane::nThreadsAndBlocks::sharedMemSize() noexcept {
    return sharedMemSize_;
}

template <typename Type>
inline void OpenHurricane::nThreadsAndBlocks::setSharedMemSize1(const cu_integer num) noexcept {
    sharedMemSize_ = num * sizeof(Type);
}

inline void OpenHurricane::nThreadsAndBlocks::setBlockAndGridSize(const cu_integer nsp,
                                                              const cu_integer nCells) {
    unsigned int blocksize_x, blocksize_y;
    unsigned int gridsize_x;
    if (nsp < 8) {
        blocksize_x = 4;
        blocksize_y = 32;
        //blocksize_y = 64;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 16) {
        blocksize_x = 8;
        blocksize_y = 16;
        //blocksize_y = 32;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 32) {
        blocksize_x = 16;
        blocksize_y = 8;
        //blocksize_y = 16;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 64) {
        blocksize_x = 32;
        blocksize_y = 4;
        //blocksize_y = 8;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 128) {
        blocksize_x = 64;
        blocksize_y = 2;
        //blocksize_y = 4;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else {
        blocksize_x = 128;
        blocksize_y = 1;
        //blocksize_y = 2;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    }

    dim3 blockSizeTmp(blocksize_x, blocksize_y, 1);
    dim3 gridSizeTmp(gridsize_x, 1, 1);

    blockSize_ = blockSizeTmp;
    gridSize_ = gridSizeTmp;
}

inline void OpenHurricane::nThreadsAndBlocks::setBlockAndGridSize64(const cu_integer nsp,
                                                                const cu_integer nCells) {
    unsigned int blocksize_x, blocksize_y;
    unsigned int gridsize_x;
    if (nsp < 8) {
        blocksize_x = 4;
        blocksize_y = 16;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 16) {
        blocksize_x = 8;
        blocksize_y = 8;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 32) {
        blocksize_x = 16;
        blocksize_y = 4;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 64) {
        blocksize_x = 32;
        blocksize_y = 2;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else {
        blocksize_x = 64;
        blocksize_y = 1;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    }

    dim3 blockSizeTmp(blocksize_x, blocksize_y, 1);
    dim3 gridSizeTmp(gridsize_x, 1, 1);

    blockSize_ = blockSizeTmp;
    gridSize_ = gridSizeTmp;
}

#endif // CUDA_PARALLEL
