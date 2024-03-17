#include "cudaEvents.hpp"
/*!
 * \file cudaEvents.inl
 * \brief The In-Line functions of the <i>cudaEvents.hpp</i> file.
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
inline OpenHurricane::cudaEvents::cudaEvents() : evt_() {}

inline OpenHurricane::cudaEvents::~cudaEvents() noexcept {}

inline void OpenHurricane::cudaEvents::create() {
    checkCUDAError(cudaEventCreate(&evt_));
}

inline void OpenHurricane::cudaEvents::destroy() const {
    checkCUDAError(cudaEventDestroy(evt_));
}

inline cudaEvent_t &OpenHurricane::cudaEvents::evt() noexcept {
    return evt_;
}

inline const cudaEvent_t &OpenHurricane::cudaEvents::evt() const noexcept {
    return evt_;
}

inline void OpenHurricane::cudaEvents::record(const cudaStreams &str) const {
    checkCUDAError(cudaEventRecord(evt_, str.stream()));
}

inline auto OpenHurricane::cudaEvents::eventQuery() const {
    return cudaEventQuery(evt_);
}

inline bool OpenHurricane::cudaEvents::isReady() const {
    return eventQuery() == cudaSuccess;
}

inline bool OpenHurricane::cudaEvents::isNotReady() const {
    return !isReady();
}

inline const cudaEvent_t &OpenHurricane::cudaEvents::operator()() const noexcept {
    return evt_;
}

#endif // CUDA_PARALLEL
