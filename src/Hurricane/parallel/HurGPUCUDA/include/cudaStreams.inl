#include "cudaStreams.hpp"
/*!
 * \file cudaStreams.inl
 * \brief The In-Line functions of the <i>cudaStreams.hpp</i> file.
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
inline OpenHurricane::cudaStreams::cudaStreams() : stream_() {}

inline OpenHurricane::cudaStreams::~cudaStreams() noexcept {}

inline void OpenHurricane::cudaStreams::create() {
    checkCUDAError(cudaStreamCreate(&stream_));
}

inline void OpenHurricane::cudaStreams::synchronize() const {
    checkCUDAError(cudaStreamSynchronize(stream_));
}

inline void OpenHurricane::cudaStreams::destroy() const {
    checkCUDAError(cudaStreamDestroy(stream_));
}

inline cudaStream_t &OpenHurricane::cudaStreams::stream() noexcept {
    return stream_;
}

inline const cudaStream_t &OpenHurricane::cudaStreams::stream() const noexcept {
    return stream_;
}

inline auto OpenHurricane::cudaStreams::streamQuery() const {
    return cudaStreamQuery(stream_);
}

inline bool OpenHurricane::cudaStreams::isReady() const {
    return streamQuery() == cudaSuccess;
}

inline bool OpenHurricane::cudaStreams::isNotReady() const {
    return !isReady();
}

inline const cudaStream_t &OpenHurricane::cudaStreams::operator()() const noexcept {
    return stream_;
}

#endif // CUDA_PARALLEL
