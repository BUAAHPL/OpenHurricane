/*!
 * \file HurGPU.hpp
 * \brief Headers of the GPU parallel computing.
 *        The subroutines and functions are in the <i>HurGPU.cpp</i> file.
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
#include "HurMPI.hpp"
#include "Lists.hpp"
#include "logFile.hpp"

#ifdef CUDA_PARALLEL

#ifndef checkCUDAGPUError
#define checkCUDAGPUError(val)                                                    \
    if ((val) != cudaSuccess) {                                                   \
        errorAbortStr(("CUDA error: " + std::string(cudaGetErrorString((val))))); \
    }

#endif // !checkCUDAGPUError

#endif // CUDA_PARALLEL

namespace OpenHurricane {
    /**
     * \brief The forward declaration of class of argument parsing.
     */
    class argParse;

    class HurGPU {
    private:
        /*!\brief The number of GPUs used in this process.*/
        static integer nGPU_;

        /** \brief The id of GPU devices used in this process. */
        static integerList GPUDev_;

        /*!\brief The total number of GPUs used in this program.*/
        static integer nTotalGPU_;

        static bool useGPU_;

        /**
         * \brief Is using multi-GPU?
         * Multi-GPU is only supported for serial run.
         */
        static bool isMultiGPU_;

        static bool deviceSupportsMemoryPools_;
        static bool deviceCanMapHostMemory_;

        static int cudaDriverVersion_;

        static bool checkCanMapHostmemory(const int idev);

        /**
         * \brief MPI Communicators for the process with the same GPU.
         */
        static HurMPIBase::subCommtor sameGPUCommtor_;

    public:
        /**
         * \brief Null constructor.
         */
        inline HurGPU() {}

        /**
         * \brief Disallow copy constructor.
         */
        HurGPU(const HurGPU &) = delete;
        HurGPU &operator=(const HurGPU &) = delete;

        /**
         * \brief Destructor.
         */
        inline ~HurGPU() noexcept {}

        /**
         * \brief Initializing with argParse.
         */
        static void init(const argParse &arg);

        /**
         * \brief Set CUDA GPU with given device index.
         * \param[in] i - The CUDA GPU index
         */
        static void setDevice(const int i);
        static void resetDevice();

        static void setSharedMenBankSize();
        static void setSharedMemBankSizeEightByte();
        static void setSharedMemBankSizeFourByte();

        static void setCachePreferShared();

        /*!\brief The number of GPUs used in this process.*/
        hur_nodiscard static inline int nGPU() noexcept { return nGPU_; }

        /** \brief The id of GPU devices used in this process. */
        hur_nodiscard static inline const integerList &GPUDev() noexcept { return GPUDev_; }

        /*!\brief The total number of GPUs used in this program.*/
        hur_nodiscard static inline int nTotalGPU() noexcept { return nTotalGPU_; }

        hur_nodiscard static inline bool useGPU() noexcept { return useGPU_; }

        /**
         * \brief Get GPU memory in MB.
         * \param[out] avail - available GPU memory [MB].
         * \param[out] total - total GPU memory [MB].
         */
        static void getMemInfo(integer &avail, integer &total);

        /**
         * \brief Is using multi-GPU?
         */
        hur_nodiscard static inline bool isMultiGPU() noexcept { return isMultiGPU_; }

        hur_nodiscard static inline bool deviceSupportsMemoryPools() noexcept {
            return deviceSupportsMemoryPools_;
        }
        hur_nodiscard static inline bool deviceCanMapHostMemory() noexcept {
            return deviceCanMapHostMemory_;
        }

        hur_nodiscard static inline int cudaDriverVersion() noexcept { return cudaDriverVersion_; }

        /**
         * \brief MPI Communicators for the process with the same GPU.
         */
        hur_nodiscard static inline HurMPIBase::Comm sameGPUComm() noexcept {
            return sameGPUCommtor_.subComm();
        }
        hur_nodiscard static inline int sameGPUProcRank() noexcept {
            return sameGPUCommtor_.subProcRank();
        }
        hur_nodiscard static inline int sameGPUProcSize() noexcept {
            return sameGPUCommtor_.subProcSize();
        }

        /**
         * \brief MPI Communicators for the process with the same GPU.
         */
        hur_nodiscard static inline const HurMPIBase::subCommtor &sameGPUCommtor() noexcept {
            return sameGPUCommtor_;
        }
    };

} // namespace OpenHurricane
