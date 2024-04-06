/*!
 * \file GPUInfo.cpp
 * \brief Main subroutines for checking GPU information.
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
#include "GPUInfo.hpp"
#include "CUDAFunctions.hpp"
#include "errorAbort.hpp"

/*!\brief To get the information of GPU devices.
        Return 0 if not have GPU device that support CUDA.
        Return the count of the device if succeed.*/
short int OpenHurricane::GPUInfo::getGPUInfo() {
#ifdef CUDA_PARALLEL

    int deviceCount = 0;
    cudaError_t errorId = cudaGetDeviceCount(&deviceCount);

    if (errorId != cudaSuccess) {
        std::string errStr = "Error: cudaGetDeviceCount() returned ";
        errStr += std::to_string(static_cast<int>(errorId));
        errStr += ": ";
        errStr += cudaGetErrorString(errorId);

        LFatal(errStr.c_str());
        //HurMPIBase::abort(HurMPIBase::getComm(), EXIT_FAILURE);
    }

    if (deviceCount == 0) {
        /*Pout<<"WArning: There are no available GPU device(s) that support CUDA" << std::endl
                 << " GPU acceleration has been turned off!" << std::endl;*/
        PLWarning("WArning: There are no available GPU device(s) that support CUDA\n GPU "
                  "acceleration has been turned off!");

        return 0;
    }
    for (integer iNode = 0; iNode < HurMPIBase::getNode(); ++iNode) {
        if (HurMPIBase::procNode() == iNode) {
            if (HurMPIBase::masterInThisNode()) {
                int namelen;
                char processor_name[MPI_MAX_PROCESSOR_NAME];
                MPI_Get_processor_name(processor_name, &namelen);
                std::cout << std::endl;
                std::cout << "     Detected " << deviceCount
                          << " CUDA Capable GPU device(s) on host: " << processor_name << std::endl;
            }
            int cudaDriverVersion = 0, cudaRuntimeVersion = 0;
            for (int dev = 0; dev < deviceCount; dev++) {
                cudaSetDevice(dev);
                cudaDeviceProp deviceProp;
                cudaGetDeviceProperties(&deviceProp, dev);

                cudaDriverGetVersion(&cudaDriverVersion);
                cudaRuntimeGetVersion(&cudaRuntimeVersion);
                if (HurMPIBase::masterInThisNode()) {
                    std::cout << "     "
                                 "---------------------------------------------"
                                 "---------------------------------"
                              << std::endl;
                    std::cout << "       GPU device " << dev << ": " << deviceProp.name
                              << std::endl;
                    std::cout << "     "
                                 "---------------------------------------------"
                                 "---------------------------------"
                              << std::endl;
                    std::cout << "         CUDA Driver Version / Runtime "
                                 "Version:          "
                              << cudaDriverVersion / 1000 << "." << (cudaDriverVersion % 100) / 10
                              << " / " << cudaRuntimeVersion / 1000 << "."
                              << (cudaRuntimeVersion % 100) / 10 << std::endl;
                    std::cout << "         CUDA Capability of GPU:             "
                                 "            "
                              << deviceProp.major << "." << deviceProp.minor << std::endl;
                    /*char msg[256];
                    SPRINTF(msg, "        Total amount of global memory:                  %.0f MBytes (%.0f Gbytes)\n",
                            (float)deviceProp.totalGlobalMem / 1048576.0f, (float)deviceProp.totalGlobalMem / 1048576.0f / 1024.0f);
                    printf("%s", msg);*/
                    std::cout << "         Total amount of global memory:            "
                                 "      "
                              << int(deviceProp.totalGlobalMem / 1024 / 1024) << " MBytes ("
                              << int(deviceProp.totalGlobalMem / 1024 / 1024 / 1024) << " GBytes)"
                              << std::endl
                              << "         " << deviceProp.multiProcessorCount
                              << " Multiprocessors, "
                              << _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor)
                              << " CUDA cores/MP:          "
                              << _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) *
                                     deviceProp.multiProcessorCount
                              << " CUDA cores" << std::endl
                              << "         Total amount of shared memory per block: "
                              << deviceProp.sharedMemPerBlock / 1024 << " KBytes" << std::endl;
                    std::cout << "     "
                                 "---------------------------------------------"
                                 "---------------------------------"
                              << std::endl
                              << std::endl;
                }
            }
        }
        HurMPIBase::barrier();
    }
    return deviceCount;
#else
    return 0;
#endif // CUDA_PARALLEL
}

int OpenHurricane::GPUInfo::getGPUNum() {
#ifdef CUDA_PARALLEL
    int deviceCount = 0;
    cudaError_t errorId = cudaGetDeviceCount(&deviceCount);

    if (errorId != cudaSuccess) {
        std::string errStr = "Error: cudaGetDeviceCount() returned ";
        errStr += std::to_string(static_cast<int>(errorId));
        errStr += ": ";
        errStr += cudaGetErrorString(errorId);

        LFatal(errStr.c_str());
        //HurMPIBase::abort(HurMPIBase::getComm(), EXIT_FAILURE);
    }
    return deviceCount;
#else
    return 0;
#endif // CUDA_PARALLEL
}
