/*!
 * \file HurGPU.cpp
 * \brief Main subroutines for GPU parallel computing.
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
#include "HurGPU.hpp"
#include "CUDAFunctions.hpp"
#include "GPUInfo.hpp"
#include "HurMPI.hpp"
#include "argParse.hpp"

OpenHurricane::integer OpenHurricane::HurGPU::nGPU_ = 0;
OpenHurricane::integerList OpenHurricane::HurGPU::GPUDev_(1);
OpenHurricane::integer OpenHurricane::HurGPU::nTotalGPU_ = 0;
bool OpenHurricane::HurGPU::useGPU_ = false;
bool OpenHurricane::HurGPU::isMultiGPU_ = false;
bool OpenHurricane::HurGPU::deviceSupportsMemoryPools_ = false;
bool OpenHurricane::HurGPU::deviceCanMapHostMemory_ = false;
int OpenHurricane::HurGPU::cudaDriverVersion_ = 0;
OpenHurricane::HurMPI::subCommtor OpenHurricane::HurGPU::sameGPUCommtor_;

bool OpenHurricane::HurGPU::checkCanMapHostmemory(const int idev) {
#ifdef CUDA_PARALLEL
    checkCUDAGPUError(cudaSetDevice(idev));
    cudaDeviceProp deviceProp;
    /* Verify the selected device supports mapped memory and set the device
     flags for mapping host memory. */
    checkCUDAGPUError(cudaGetDeviceProperties(&deviceProp, idev));

#if CUDART_VERSION >= 2020

    if (!deviceProp.canMapHostMemory) {
        return false;
    } else {
        return true;
    }

#else
    return false;
#endif
#else
    return false;
#endif
}

void OpenHurricane::HurGPU::init(const argParse &arg) {
#ifdef CUDA_PARALLEL
    if (arg.hasGPU()) {
        checkCUDAGPUError(cudaDriverGetVersion(&cudaDriverVersion_));

        int GPUN = (int)GPUInfo::getGPUInfo();
        const auto nGPUPerMache = arg.nGPUPerMachine();
        /*int GPUN;
        checkCUDAGPUError(cudaGetDeviceCount(&GPUN));*/
        if (GPUN == 0) {
            checkWarning("No CUDA supported GPU found, and GPU would not be used");
        } else if (nGPUPerMache == 0) {
            nTotalGPU_ = GPUN;
            if (HurMPIBase::parRun()) {
                nTotalGPU_ = min(nTotalGPU_, HurMPIBase::getProcSize());
                nGPU_ = 1;
                isMultiGPU_ = false;
                GPUDev_.resize(nGPU_);

                GPUDev_[0] = HurMPIBase::procRankNode() % nTotalGPU_;

            } else {
                nGPU_ = nTotalGPU_;
                isMultiGPU_ = (nTotalGPU_ > 1);
                GPUDev_.resize(nGPU_);
                for (int i = 0; i < nGPU_; ++i) {
                    GPUDev_[i] = i;
                }
            }
        } else {
            if (nGPUPerMache > GPUN) {
                PLWarning("The specified GPU number is %d which is larger than the number of GPUs "
                          "of this machine: %d",
                          nGPUPerMache, GPUN);
            }
            nTotalGPU_ = min(nGPUPerMache, GPUN);
            if (HurMPIBase::parRun()) {
                nTotalGPU_ = min(nTotalGPU_, HurMPIBase::getProcSize());
                nGPU_ = 1;
                isMultiGPU_ = false;
                GPUDev_.resize(nGPU_);

                GPUDev_[0] = HurMPIBase::procRankNode() % nTotalGPU_;

            } else {
                nGPU_ = nTotalGPU_;
                isMultiGPU_ = (nTotalGPU_ > 1);
                GPUDev_.resize(nGPU_);
                for (int i = 0; i < nGPU_; ++i) {
                    GPUDev_[i] = i;
                }
            }
        }

        if (nGPU_ != 0) {
            useGPU_ = true;
            if (HurMPIBase::parRun()) {
                if (nTotalGPU_ == 1) {
                    sameGPUCommtor_.setSubComm(HurMPIBase::getComm());
                    sameGPUCommtor_.setSubProcRank(HurMPIBase::getProcRank());
                    sameGPUCommtor_.setSubProcSize(HurMPIBase::getProcSize());
                } else if (HurMPIBase::multiNodes()) {
                    sameGPUCommtor_ = HurMPIBase::subCommtor(HurMPIBase::getNodeComm(), GPUDev_[0],
                                                             HurMPIBase::procRankNode());
                } else {
                    sameGPUCommtor_ = HurMPIBase::subCommtor(HurMPIBase::getComm(), GPUDev_[0],
                                                             HurMPIBase::getProcRank());
                }                
            } else {
                sameGPUCommtor_.setSubComm(HurMPIBase::getComm());
                sameGPUCommtor_.setSubProcRank(HurMPIBase::getProcRank());
                sameGPUCommtor_.setSubProcSize(HurMPIBase::getProcSize());
            }
        }
        HurMPIBase::allReduce(useGPU_, MPI_LAND);
        if (useGPU_) {
            if (cudaDriverVersion_ >= 11020) {
                deviceSupportsMemoryPools_ = true;
                for (integer ig = 0; ig < nGPU_; ++ig) {
                    int val;
                    checkCUDAGPUError(
                        cudaDeviceGetAttribute(&val, cudaDevAttrMemoryPoolsSupported, GPUDev_[ig]));
                    if (val == 0) {
                        deviceSupportsMemoryPools_ = false;
                        break;
                    }
                }
            }
            cuAttributes::init();

            deviceCanMapHostMemory_ = true;
            for (integer ig = 0; ig < nGPU_; ++ig) {
                bool canMapHM = checkCanMapHostmemory(ig);
                if (!canMapHM) {
                    deviceCanMapHostMemory_ = false;
                    break;
                }
            }
        }
    }
#else

    if (arg.hasGPU()) {
        checkWarning("This program does not support GPU, and the GPU would not be used");
    }
#endif
}

void OpenHurricane::HurGPU::setDevice(const int i) {
#ifdef CUDA_PARALLEL
    checkCUDAGPUError(cudaSetDevice(i));
    if (argParse::gpuZeroCopy()) {
        if (deviceCanMapHostMemory()) {
            checkCUDAGPUError(cudaSetDeviceFlags(cudaDeviceMapHost));
        }
    }
#endif
}

void OpenHurricane::HurGPU::resetDevice() {
#ifdef CUDA_PARALLEL
    checkCUDAGPUError(cudaDeviceReset());
#endif
}

void OpenHurricane::HurGPU::setSharedMenBankSize() {
#ifdef CUDA_PARALLEL
#if defined(HURRICANE_SP)
    setSharedMemBankSizeFourByte();
#elif defined(HURRICANE_DP)
    setSharedMemBankSizeEightByte();
#endif
#endif
}

void OpenHurricane::HurGPU::setSharedMemBankSizeEightByte() {
#ifdef CUDA_PARALLEL
    //Devices of compute capability 3.x have configurable bank size. Setting the bank size to eight bytes can help
    //avoid shared memory bank conflicts when accessing double precision data.
    checkCUDAGPUError(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
#endif
}

void OpenHurricane::HurGPU::setSharedMemBankSizeFourByte() {
#ifdef CUDA_PARALLEL
    checkCUDAGPUError(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte));
#endif
}

void OpenHurricane::HurGPU::setCachePreferShared() {
#ifdef CUDA_PARALLEL
    //On devices of compute capability of 3.x, each multiprocessor has 64KB of on-chip memory that can be partitioned
    //between L1 cache and shared memory.For devices of compute capability of 3.x, there are three settings, 48KB
    //shared memory / 16KB L1 cahce, 16KB shared memory / 48KB L1 cache, and 32KB shared memory / 32KB L1 cahce.
    checkCUDAGPUError(cudaDeviceSetCacheConfig(cudaFuncCachePreferShared));
#endif
}

void OpenHurricane::HurGPU::getMemInfo(integer &avail, integer &total) {
#ifdef CUDA_PARALLEL
    size_t av;
    size_t to;
    checkCUDAGPUError(cudaMemGetInfo(&av, &to));
    avail = static_cast<integer>(av / 1024 / 1024);
    total = static_cast<integer>(to / 1024 / 1024);
#else
    avail = 0;
    total = 0;
#endif
}
