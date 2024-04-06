/*!
 * \file cudaStreams.hpp
 * \brief Header of CUDA streams.
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
     * \brief The class of CUDA streams.
     */
    class cudaStreams {
    private:
        cudaStream_t stream_;

    public:
        // Constructors

        /**
         * \brief Null constructor.
         */
        inline cudaStreams();

        /**
         * \brief Disallow bitwise copy constructor.
         */
        cudaStreams(const cudaStreams &) = delete;

        /**
         * \brief Destructor.
         */
        inline ~cudaStreams() noexcept;

        inline void create();

        /**
         * \brief Wait for all operations to finish of this stream.
         */
        inline void synchronize() const;

        inline void destroy() const;

        inline cudaStream_t &stream() noexcept;
        inline const cudaStream_t &stream() const noexcept;

        inline auto streamQuery() const;

        inline bool isReady() const;
        inline bool isNotReady() const;

        /**
         * \brief Disallow bitwise assignment.
         */
        cudaStreams &operator=(const cudaStreams &) = delete;

        inline const cudaStream_t &operator()() const noexcept;
    };
} // namespace OpenHurricane
#include "cudaStreams.inl"

#endif // CUDA_PARALLEL