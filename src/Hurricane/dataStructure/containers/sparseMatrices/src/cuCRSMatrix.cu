#include "cuCRSMatrix.hpp"
/*!
 * \file cuCRSMatrix.cu
 * \brief Main subroutines of cuCRSMatrix.
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
#include "cuCRSMatrix.hpp"

cu_host OpenHurricane::cuMCRSMatrixAddressing::cuMCRSMatrixAddressing(
    const size_type nRows, const size_type nColumns, const size_type NNZ,
    const size_type *__restrict__ rowPtr, const size_type *__restrict__ col)
    : nRows_(nRows), nColumns_(nColumns), NNZ_(NNZ), rowPtr_(nRows + 1), col_(NNZ),
      diagIndex_(nRows) {
    for (cu_integer i = 0; i < nRows_ + 1; ++i) {
        rowPtr_[i] = *(rowPtr + i);
    }
    for (cu_integer i = 0; i < NNZ; ++i) {
        col_[i] = *(col + i);
    }
        
    cu_integer uEnd = rowPtr[0];
    cu_integer uStart;

    for (cu_integer i = 0; i < nRows_; ++i) {
        // Start and end for this row of matrix A
        uStart = uEnd;
        uEnd = rowPtr[i + 1];

        for (cu_integer j0 = uStart; j0 < uEnd; ++j0) {
            const auto j = col[j0];
            if (i == j) {
                diagIndex_[i] = j0;
            }
        }
    }
}

cu_host void OpenHurricane::cuMCRSMatrixAddressing::clear() noexcept {
    nRows_ = 0;
    nColumns_ = 0;
    NNZ_ = 0;
    rowPtr_.clear();
    col_.clear();
    diagIndex_.clear();
}

cu_host OpenHurricane::cuCRSMatrixAddressing::cuCRSMatrixAddressing(
    const size_type nRows, const size_type nColumns, const size_type NNZ,
    const size_type *__restrict__ rowPtr, const size_type *__restrict__ col)
    : nRows_(nRows), nColumns_(nColumns), NNZ_(NNZ), rowPtr_(nRows + 1), col_(NNZ),
      diagIndex_(nRows) {
    for (cu_integer i = 0; i < nRows_ + 1; ++i) {
        rowPtr_[i] = *(rowPtr + i);
    }
    for (cu_integer i = 0; i < NNZ; ++i) {
        col_[i] = *(col + i);
    }

    cu_integer uEnd = rowPtr[0];
    cu_integer uStart;

    for (cu_integer i = 0; i < nRows_; ++i) {
        // Start and end for this row of matrix A
        uStart = uEnd;
        uEnd = rowPtr[i + 1];

        for (cu_integer j0 = uStart; j0 < uEnd; ++j0) {
            const auto j = col[j0];
            if (i == j) {
                diagIndex_[i] = j0;
            }
        }
    }
}

cu_host void OpenHurricane::cuCRSMatrixAddressing::clear() noexcept {
    nRows_ = 0;
    nColumns_ = 0;
    NNZ_ = 0;
    rowPtr_.clear();
    col_.clear();
    diagIndex_.clear();
}
#endif // CUDA_PARALLEL
