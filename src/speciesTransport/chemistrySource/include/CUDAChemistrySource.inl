#include "CUDAChemistrySource.hpp"
/*!
 * \file CUDAChemistrySource.inl
 * \brief The In-Line functions of the <i>CUDAChemistrySource.hpp</i> file.
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

inline OpenHurricane::CUDAChemistrySourceODEs &
OpenHurricane::CUDAChemistrySource::CUDAChemODEs() {
    if (CUDAChemODEsPtr_ == nullptr) {
        makeCUDAChemODEs();
    }
    return *CUDAChemODEsPtr_;
}

inline const OpenHurricane::CUDAChemistrySourceODEs &
OpenHurricane::CUDAChemistrySource::CUDAChemODEs() const {
    if (CUDAChemODEsPtr_ == nullptr) {
        makeCUDAChemODEs();
    }
    return *CUDAChemODEsPtr_;
}

inline OpenHurricane::EulerCUDA &
OpenHurricane::CUDAChemistrySource::EulerSolver() {
    if (EulerSolverPtr_ == nullptr) {
        makeEulerSolver();
    }
    return *EulerSolverPtr_;
}

inline const OpenHurricane::EulerCUDA &
OpenHurricane::CUDAChemistrySource::EulerSolver() const {
    if (EulerSolverPtr_ == nullptr) {
        makeEulerSolver();
    }
    return *EulerSolverPtr_;
}

#endif // CUDA_PARALLEL