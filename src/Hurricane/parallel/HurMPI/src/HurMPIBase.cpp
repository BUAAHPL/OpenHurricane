/*!
 * \file HurMPIBase.cpp
 * \brief Main subroutines for the mpi structures.
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
#include "HurMPIBase.hpp"
#include "dataStructure.hpp"
#include "logFile.hpp"
#include <iostream>

#define max(a, b) (a) > (b) ? (a) : (b);

#ifdef _WIN32
#include <winsock2.h>
#include <Windows.h>
#include <process.h>
#pragma comment(lib, "ws2_32.lib")
#elif defined(__linux__) || defined(LINUX)
#include <unistd.h>
#endif

namespace OpenHurricane {
    int HurMPIBase::procRank_ = 0;
    int HurMPIBase::procSize_ = 1;
    int HurMPIBase::minRank_ = 0;
    HurMPIBase::Comm HurMPIBase::HurWorldComm_ = MPI_COMM_WORLD;
    bool HurMPIBase::usingWinMinRank_ = false;
    HurMPIBase::Win HurMPIBase::winMinRank_ = 0;
    int HurMPIBase::node_ = 1;
    int HurMPIBase::procNode_ = 0;
    bool HurMPIBase::multiNode_ = false;
    bool HurMPIBase::parRun_ = false;
    bool HurMPIBase::haveThreads_ = false;
    bool HurMPIBase::isInitialized_ = false;
    HurMPIBase::Info HurMPIBase::info_ = MPI_INFO_NULL;
    HurMPIBase::subCommtor HurMPIBase::nodeCommtor_;
} // namespace OpenHurricane

hur_nodiscard const OpenHurricane::HurMPIBase::subCommtor &OpenHurricane::HurMPIBase::subCommtor::null() {
    return NullRefObj::nullRef<subCommtor>();
}

OpenHurricane::HurMPIBase::subCommtor::subCommtor(const Comm comm, const int color, const int key)
    : subComm_(0), subProcRank_(0), subProcSize_(0) {
#ifdef MPI_PARALLEL
    MPI_Comm_split(comm, color, key, &subComm_);
    MPI_Comm_rank(subComm_, &subProcRank_);
    MPI_Comm_size(subComm_, &subProcSize_);
#endif
 }

 void OpenHurricane::HurMPIBase::subCommtor::freeSubComm() {
#ifdef MPI_PARALLEL
    MPI_Comm_free(&subComm_);
    subProcRank_ = 0;
    subProcSize_ = 0;
#endif
 }

#ifdef MPI_PARALLEL
#else // Not define MPI_PARALLEL
void OpenHurricane::HurMPIBase::copyData(const void *sendbuf, void *recvbuf, int size, Datatype datatype) {
    switch (datatype) {
    case MPI_DOUBLE:
        for (int i = 0; i < size; ++i) {
            static_cast<double *>(recvbuf)[i] = static_cast<const double *>(sendbuf)[i];
        }
        break;

    case MPI_FLOAT:
        for (int i = 0; i < size; ++i) {
            static_cast<float *>(recvbuf)[i] = static_cast<const float *>(sendbuf)[i];
        }
        break;

    case MPI_INT:
        for (int i = 0; i < size; ++i) {
            static_cast<int *>(recvbuf)[i] = static_cast<const int *>(sendbuf)[i];
        }
        break;

    case MPI_LONG:
        for (int i = 0; i < size; ++i) {
            static_cast<long *>(recvbuf)[i] = static_cast<const long *>(sendbuf)[i];
        }
        break;

    case MPI_UNSIGNED_LONG:
        for (int i = 0; i < size; ++i) {
            static_cast<unsigned long *>(recvbuf)[i] =
                static_cast<const unsigned long *>(sendbuf)[i];
        }
        break;

    case MPI_CHAR:
        for (int i = 0; i < size; ++i) {
            static_cast<char *>(recvbuf)[i] = static_cast<const char *>(sendbuf)[i];
        }
        break;

    case MPI_SHORT:
        for (int i = 0; i < size; ++i) {
            static_cast<short *>(recvbuf)[i] = static_cast<const short *>(sendbuf)[i];
        }
        break;

    case MPI_UNSIGNED_SHORT:
        for (int i = 0; i < size; ++i) {
            static_cast<unsigned short *>(recvbuf)[i] =
                static_cast<const unsigned short *>(sendbuf)[i];
        }
        break;

    case MPI_CXX_BOOL:
        for (int i = 0; i < size; ++i) {
            static_cast<bool *>(recvbuf)[i] = static_cast<const bool *>(sendbuf)[i];
        }
        break;

    case MPI_INT32_T:
        for (int i = 0; i < size; ++i) {
            static_cast<int32_t *>(recvbuf)[i] = static_cast<const int32_t *>(sendbuf)[i];
        }
        break;

    case MPI_INT64_T:
        for (int i = 0; i < size; ++i) {
            static_cast<int64_t *>(recvbuf)[i] = static_cast<const int64_t *>(sendbuf)[i];
        }
        break;

    default:
        break;
    }
}
#endif
void OpenHurricane::HurMPIBase::error(const std::string &ErrorMsg, const std::string &FunctionName,
                                  bool isAbort) {
#ifdef MPI_PARALLEL
    if (isInitialized_) {
        bool isCollective = isCollectiveCall(1.0);
        minRank_ = procRank_;
        if (!isCollective) {
            minRank_ = getMinRank(winMinRank_);
        }

        if ((!isCollective&&procRank_ == minRank_) || (isCollective && master())) {
            std::cerr << std::endl << std::endl;
            std::cerr << "========================================================================="
                         "=========="
                      << std::endl;
            std::cerr << "  Fatal error in routine: \"" << FunctionName.c_str()
                      << "\": " << std::endl;
            std::cerr << "-------------------------------------------------------------------------"
                         "----------"
                      << std::endl;
            std::cerr << "  Exit reason: " << ErrorMsg.c_str() << std::endl;
            if (isCollective) {
                std::cerr << "    (The error call is collective)" << std::endl;
            } else {
                std::cerr << "    (The error call is not collective)" << std::endl;
            }
            std::cerr << "  By OpenHurricane CFD program" << std::endl;
            std::cerr << "========================================================================="
                         "=========="
                      << std::endl;
            std::cerr << std::endl << std::endl;
            if (isAbort) {
                abort(HurWorldComm_, EXIT_FAILURE);
            }
        }
    } else {
        if (procRank_ == 0) {
            std::cerr << std::endl << std::endl;
            std::cerr << "========================================================================="
                         "=========="
                      << std::endl;
            std::cerr << "  Fatal error in routine: \"" << FunctionName.c_str()
                      << "\": " << std::endl;
            std::cerr << "-------------------------------------------------------------------------"
                         "----------"
                      << std::endl;
            std::cerr << "  Exit reason: " << ErrorMsg.c_str() << std::endl;
            std::cerr << "  By OpenHurricane CFD program" << std::endl;
            std::cerr << "========================================================================="
                         "=========="
                      << std::endl;
            std::cerr << std::endl << std::endl;
        }
        if (isAbort) {
            abort(HurWorldComm_, 0);
        }
    }
#else // Not define MPI_PARALLEL
    if (procRank_ == 0) {
        std::cerr << std::endl << std::endl;
        std::cerr
            << "==================================================================================="
            << std::endl;
        std::cerr << "  Fatal error in routine: \"" << FunctionName.c_str() << "\": " << std::endl;
        std::cerr
            << "-----------------------------------------------------------------------------------"
            << std::endl;
        std::cerr << "  Exit reason: " << ErrorMsg.c_str() << std::endl;
        std::cerr << "  By OpenHurricane CFD program" << std::endl;
        std::cerr
            << "==================================================================================="
            << std::endl;
        std::cerr << std::endl << std::endl;
    }
    if (isAbort) {
        abort(HurWorldComm_, 0);
    }
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::warning(const std::string& WarningMsg,
    const std::string& FunctionName) {
#ifdef MPI_PARALLEL
    if (isInitialized_) {
        bool isCollective = isCollectiveCall(1.0);
        minRank_ = procRank_;
        if (!isCollective) {
            minRank_ = getMinRank(winMinRank_);
        }
        if ((!isCollective && procRank_ == minRank_) || (isCollective && master())) {
            std::cerr << std::endl << std::endl;
            std::cerr << "Warning in \"" << FunctionName.c_str() << "\": " << std::endl;
            std::cerr << "-------------------------------------------------------------------------"
                      << std::endl;
            std::cerr << WarningMsg.c_str() << std::endl;
            if (isCollective) {
                std::cerr << "(The warning call is collective)" << std::endl;
            } else {
                std::cerr << "(The waring call is not collective)" << std::endl;
            }
            std::cerr << "------------------------------ Warning ----------------------------------"
                      << std::endl;
            std::cerr << std::endl << std::endl;
        }
    } else {
        if (procRank_ == 0) {
            std::cerr << std::endl << std::endl;
            std::cerr << "Waring in \"" << FunctionName.c_str() << "\": " << std::endl;
            std::cerr << "-------------------------------------------------------------------------"
                      << std::endl;
            std::cerr << WarningMsg.c_str() << std::endl;
            std::cerr << "-------------------------------- Warning --------------------------------"
                      << std::endl;
            std::cerr << std::endl << std::endl;
        }
    }
#else // Not define MPI_PARALLEL
    if (procRank_ == 0) {
        std::cerr << std::endl << std::endl;
        std::cerr << "Waring in \"" << FunctionName.c_str() << "\": " << std::endl;
        std::cerr << "-------------------------------------------------------------------------"
                  << std::endl;
        std::cerr << WarningMsg.c_str() << std::endl;
        std::cerr << "-------------------------------- Warning --------------------------------"
                  << std::endl;
        std::cerr << std::endl << std::endl;
    }
#endif // MPI_PARALLEL
}

hur_nodiscard bool OpenHurricane::HurMPIBase::isCollectiveCall(double _times) {
#ifdef MPI_PARALLEL
    int flag = 0;

    /* Find out whether the error call is collective via MPI_Ibarrier. */
    Request barrierRequest;
    MPI_Ibarrier(HurWorldComm_, &barrierRequest);

    /* Try to complete the non-blocking barrier call for a second. */
    double startTime = MPI_Wtime();
    while (true) {

        MPI_Test(&barrierRequest, &flag, MPI_STATUS_IGNORE);
        if (flag) {
            break;
        }
        double currentTime = MPI_Wtime();
        if (currentTime > startTime + _times) {
            break;
        }
    }

    if (flag) {
        return true;
    }
    return false;
#else // Not define MPI_PARALLEL
    return true;
#endif // MPI_PARALLEL
}

hur_nodiscard int OpenHurricane::HurMPIBase::getMinRank(Win winMinRank) {
#ifdef MPI_PARALLEL
    auto minRank = procRank_;
    for (int i = 0; i < procRank_; ++i) {
        int MinRankErrorOther = -1;
        MPI_Win_lock(MPI_LOCK_SHARED, i, 0, winMinRank);
        MPI_Get(&MinRankErrorOther, 1, MPI_INT, i, 0, 1, MPI_INT, winMinRank);
        MPI_Win_unlock(i, winMinRank);
        if (MinRankErrorOther < minRank) {
            minRank = MinRankErrorOther;
            break;
        }
    }
    return minRank;
#else  // Not define MPI_PARALLEL
    return procRank_;
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::setThreads(const bool haveThreads) {
#ifdef HUR_DEBUG

    if (haveThreads) {
        std::cout << "Have support for threads! Rank: " << getProcRank() << std::endl;
    } else {
        std::cout << "Have no support for threads! Rank: " << getProcRank() << std::endl;
    }
#endif // HUR_DEBUG

    haveThreads_ = haveThreads;
}

int OpenHurricane::HurMPIBase::init(int* argc, char*** argv) {
    isInitialized_ = true;
#ifdef MPI_PARALLEL
    int provided_thread_support = 0;

    int MPIERR = MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided_thread_support);

    MPI_Comm_rank(HurWorldComm_, &procRank_);
    MPI_Comm_size(HurWorldComm_, &procSize_);
    nodeCommtor_.setSubProcRank(procRank_);
    nodeCommtor_.setSubProcSize(procSize_);
    setParRun();
    setThreads(provided_thread_support == MPI_THREAD_MULTIPLE);
    minRank_ = procSize_;
    int MPIERR2 = MPI_Win_create(&minRank_, sizeof(int), sizeof(int), MPI_INFO_NULL, HurWorldComm_,
                                 &winMinRank_);
    usingWinMinRank_ = true;
    return MPIERR;
#else  // Not define MPI_PARALLEL
    return MPI_SUCCESS;
#endif // MPI_PARALLEL
}

int OpenHurricane::HurMPIBase::finalize() {
#ifdef MPI_PARALLEL
    if (multiNodes()) {
        nodeCommtor_.freeSubComm();
    }
    if (usingWinMinRank_) {
        MPI_Win_free(&winMinRank_);
    }
    return MPI_Finalize();
#else  // Not define MPI_PARALLEL
    return MPI_SUCCESS;
#endif // MPI_PARALLEL
}

int OpenHurricane::HurMPIBase::bufferAttach(void *buffer, int size) {
#ifdef MPI_PARALLEL
    return MPI_Buffer_attach(buffer, size);
#else  // Not define MPI_PARALLEL
    return MPI_SUCCESS;
#endif // MPI_PARALLEL
}

int OpenHurricane::HurMPIBase::bufferDetach(void *buffer, int *size) {
#ifdef MPI_PARALLEL
    return MPI_Buffer_detach(buffer, size);
#else  // Not define MPI_PARALLEL
    return MPI_SUCCESS;
#endif // MPI_PARALLEL
}

int OpenHurricane::HurMPIBase::commRank(Comm comm, int &rank) {
#ifdef MPI_PARALLEL
    return MPI_Comm_rank(comm, &rank);
#else  // Not define MPI_PARALLEL
    rank = 0;
    return MPI_SUCCESS;
#endif // MPI_PARALLEL
}

int OpenHurricane::HurMPIBase::commSize(Comm comm, int &size) {
#ifdef MPI_PARALLEL
    return MPI_Comm_size(comm, &size);
#else  // Not define MPI_PARALLEL
    size = 1;
    return MPI_SUCCESS;
#endif // MPI_PARALLEL
}

int OpenHurricane::HurMPIBase::abort(Comm comm, int error) {
#ifdef MPI_PARALLEL
    if (isInitialized_) {
        return MPI_Abort(comm, error);
    } else {
        exit(error);
        return error;
    }
#else  // Not define MPI_PARALLEL
    exit(error);
    return error;
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::setMultiNode() {
#ifdef MPI_PARALLEL
    int namelen, pcnmlen_max;
    int *pcnamelen, *ipcnm_this, *ipcnm;
    char c;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &namelen);
    pcnamelen = new int[procSize_];
    barrier(HurWorldComm_);
    MPI_Gather(&namelen, 1, MPI_INTEGER, pcnamelen, 1, MPI_INTEGER, masterNo(), HurWorldComm_);
    if (master()) {
        pcnmlen_max = 1;
        for (int i = 0; i < procSize_; i++) {
            pcnmlen_max = max(pcnmlen_max, pcnamelen[i]);
        }
    }
    barrier(HurWorldComm_);
    MPI_Bcast(&pcnmlen_max, 1, MPI_INTEGER, masterNo(), HurWorldComm_);
    ipcnm_this = new int[int(pcnmlen_max + 1)];
    ipcnm = new int[int(procSize_ * (int(pcnmlen_max + 1)))];
    for (int i = 0; i < pcnmlen_max + 1; i++) {
        c = processor_name[i];
        ipcnm_this[i] = int(c);
    }
    MPI_Gather(ipcnm_this, pcnmlen_max + 1, MPI_INTEGER, ipcnm, pcnmlen_max + 1, MPI_INTEGER,
               masterNo(), HurWorldComm_);
    int ihost_new;

    std::string *pcnameStr;
    int *hostName;
    hostName = new int[procSize_]();
    pcnameStr = new std::string[procSize_];
    if (master()) {
        for (int ip = 0; ip < procSize_; ip++) {
            for (int j = 0; j < pcnmlen_max + 1; j++) {
                c = char(ipcnm[ip * (pcnmlen_max + 1) + j]);
                pcnameStr[ip] += c;
            }
        }

        hostName[0] = 0;

        for (int i = 0; i < procSize_; i++) {
            ihost_new = 1;
            for (int iNode = 0; iNode < node_; iNode++) {
                if (pcnameStr[i] == pcnameStr[hostName[iNode]]) {
                    ihost_new = 0;
                    break;
                }
            }
            if (ihost_new == 1) {
                multiNode_ = true;
                node_++;
                hostName[node_ - 1] = i;
            }
        }
    }

    MPI_Bcast(&multiNode_, 1, feature<bool>::MPIType, masterNo(), HurWorldComm_);
    MPI_Bcast(&node_, 1, MPI_INT, masterNo(), HurWorldComm_);

    if (multiNode_) {
        int *proNode;
        proNode = new int[procSize_];
        if (master()) {
            for (int i = 0; i < procSize_; i++) {
                proNode[i] = -1;
                for (int iNode = 0; iNode < node_; iNode++) {
                    if (pcnameStr[i] == pcnameStr[hostName[iNode]]) {
                        proNode[i] = iNode;
                        break;
                    }
                }
            }
        }
        MPI_Scatter(proNode, 1, MPI_INT, &procNode_, 1, MPI_INT, masterNo(), HurWorldComm_);
        delete[] proNode;
    }

    delete[] pcnamelen;
    delete[] ipcnm_this;
    delete[] ipcnm;
    delete[] hostName;
    delete[] pcnameStr;
    setNodeComm();
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::setNodeComm() {
#ifdef MPI_PARALLEL
    if (multiNode_) {
        nodeCommtor_ = subCommtor(HurWorldComm_, procNode_, procRank_);
    }
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::getPcInfo() {
#ifdef MPI_PARALLEL
    int namelen, pcnmlen_max;
    int *pcnamelen, *ipcnm_this, *ipcnm;
    char c;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &namelen);
    pcnamelen = new int[procSize_];
    barrier(HurWorldComm_);
    MPI_Gather(&namelen, 1, MPI_INTEGER, pcnamelen, 1, MPI_INTEGER, masterNo(), HurWorldComm_);
    if (master()) {
        pcnmlen_max = 1;
        for (int i = 0; i < procSize_; i++) {
            pcnmlen_max = max(pcnmlen_max, pcnamelen[i]);
        }
    }
    barrier(HurWorldComm_);
    MPI_Bcast(&pcnmlen_max, 1, MPI_INTEGER, masterNo(), HurWorldComm_);
    ipcnm_this = new int[int(pcnmlen_max + 1)];
    ipcnm = new int[int(procSize_ * (int(pcnmlen_max + 1)))];
    for (int i = 0; i < pcnmlen_max + 1; i++) {
        c = processor_name[i];
        ipcnm_this[i] = int(c);
    }
    MPI_Gather(ipcnm_this, pcnmlen_max + 1, MPI_INTEGER, ipcnm, pcnmlen_max + 1, MPI_INTEGER,
               masterNo(), HurWorldComm_);
    int ihost_new;

    std::string *pcnameStr;
    int *hostName;
    hostName = new int[procSize_]();
    pcnameStr = new std::string[procSize_];
    if (master()) {
        for (int ip = 0; ip < procSize_; ip++) {
            for (int j = 0; j < pcnmlen_max + 1; j++) {
                c = char(ipcnm[ip * (pcnmlen_max + 1) + j]);
                pcnameStr[ip] += c;
            }
        }

        hostName[0] = 0;

        for (int i = 0; i < procSize_; i++) {
            ihost_new = 1;
            for (int iNode = 0; iNode < node_; iNode++) {
                if (pcnameStr[i] == pcnameStr[hostName[iNode]]) {
                    ihost_new = 0;
                    break;
                }
            }
            if (ihost_new == 1) {
                multiNode_ = true;
                node_++;
                hostName[node_ - 1] = i;
            }
        }
    }

    MPI_Bcast(&multiNode_, 1, feature<bool>::MPIType, masterNo(), HurWorldComm_);
    MPI_Bcast(&node_, 1, MPI_INT, masterNo(), HurWorldComm_);

    if (multiNode_) {
        int *proNode;
        proNode = new int[procSize_];
        if (master()) {
            for (int i = 0; i < procSize_; i++) {
                proNode[i] = -1;
                for (int iNode = 0; iNode < node_; iNode++) {
                    if (pcnameStr[i] == pcnameStr[hostName[iNode]]) {
                        proNode[i] = iNode;
                        break;
                    }
                }
            }
        }
        MPI_Scatter(proNode, 1, MPI_INT, &procNode_, 1, MPI_INT, masterNo(), HurWorldComm_);
        delete[] proNode;
    }

    delete[] pcnamelen;
    delete[] ipcnm_this;
    delete[] ipcnm;
    barrier(HurWorldComm_);
    std::string os_n, os_v;
    getOsInfo(os_n, os_v);
    int numP;
    std::string cV, cB;
    getCpuInfo(numP, cV, cB);
#ifdef _WIN32
    int iPID = (int)_getpid();
#elif defined(__linux__) || defined(LINUX)
    int iPID = (int)getpid();
#endif

    if (master()) {
        std::cout << "     System info:" << std::endl;
        std::cout << "         Operating System:   " << os_n.c_str() << std::endl;
        std::cout << "         CPU vendor:  " << cV.c_str() << "    CPU cores: " << numP
                  << std::endl;
        std::cout << "         Using " << procSize_ << " CPU process(es)" << std::endl;
        if (!multiNode_) {
            std::cout << "         Running on one node" << std::endl;
            //<< std::endl;
        } else {
            std::cout << "         Running on multi nodes: " << node_ << " nodes." << std::endl
                      << std::endl;
            std::cout << "             Hostname:" << std::endl;
            for (int iNode = 0; iNode < node_; iNode++) {
                std::cout << "                  " << pcnameStr[hostName[iNode]] << std::endl;
            }
            //std::cout << std::endl;
        }

        std::cout << "  " << std::left << std::setfill('-') << std::setw(108) << "-" << std::endl;
        std::cout << "  " << std::left << std::setfill(' ') << std::setw(7) << "ID";
        std::cout << std::left << std::setfill(' ') << std::setw(21) << "Hostname";
        std::cout << std::left << std::setfill(' ') << std::setw(6) << "Core";
        std::cout << std::left << std::setfill(' ') << std::setw(23) << "O.S.";
        std::cout << std::left << std::setfill(' ') << std::setw(10) << "PID";
        std::cout << std::left << std::setfill(' ') << std::setw(20) << "CPU Brand";
        std::cout << std::endl;
        std::cout << "  " << std::left << std::setfill('-') << std::setw(108) << "-" << std::endl;
    }
    barrier(HurWorldComm_);
    for (int i = 0; i < procSize_; i++) {
        if (procRank_ == i) {
            std::string idStr = "n" + toString(procRank_);
            std::cout << "  " << std::left << std::setfill(' ') << std::setw(7) << idStr.c_str();
            std::cout << std::left << std::setfill(' ') << std::setw(21) << processor_name;
            std::string corStr = toString(i + 1) + "/" + toString(numP);
            std::cout << std::left << std::setfill(' ') << std::setw(6) << corStr.c_str();
            std::cout << std::left << std::setfill(' ') << std::setw(23) << os_n.c_str();
            std::cout << std::left << std::setfill(' ') << std::setw(10) << iPID;
            std::cout << std::left << std::setfill(' ') << std::setw(20) << cB.c_str();
            std::cout << std::endl;
        }
        barrier(HurWorldComm_);
    }
    for (int i = 0; i < procSize_; i++) {
        barrier(HurWorldComm_);
        if ((procRank_ == procSize_ - 1) && procRank_ == i) {
            std::cout << "  " << std::left << std::setfill('-') << std::setw(108) << "-"
                      << std::endl;
        }
    }
    std::cout << std::right << std::setfill(' ');
    barrier(HurWorldComm_);
    delete[] pcnameStr;
    delete[] hostName;
    setNodeComm();
#else // Not define MPI_PARALLEL
    char hname[256];
    int istart = 1;

#ifdef _WIN32
    WSADATA wsaData;
    if (WSAStartup(MAKEWORD(2, 1), &wsaData)) {
        std::cerr << "Winsock can't be initialized!" << std::endl;
        WSACleanup();
        istart = 0;
    }
#endif // _WIN32
    int igeth = 0;
    if (istart == 1) {
        igeth = gethostname(hname, sizeof(hname));
    }
    std::string os_n, os_v;
    getOsInfo(os_n, os_v);
    int numP;
    std::string cV, cB;
    getCpuInfo(numP, cV, cB);
#ifdef _WIN32
    int iPID = (int)_getpid();
#elif defined(__linux__) || defined(LINUX)
    int iPID = (int)getpid();
#endif

    std::cout << "     System info:" << std::endl;
    std::cout << "         Operating System:   " << os_n.c_str() << std::endl;
    std::cout << "         CPU vendor:  " << cV.c_str() << "    CPU cores: " << numP << std::endl;
    std::cout << "         Using single CPU process" << std::endl;

    std::cout << "         Running on one node" << std::endl;

    /*std::cout << "     ------------------------------------------------------------------------------" << std::endl;
    std::cout << "     ID    Hostname        Core    O.S.       PID    CPU Brand" << std::endl;
    std::cout << "     ------------------------------------------------------------------------------" << std::endl;*/
    std::cout << "  " << std::left << std::setfill('-') << std::setw(108) << "-" << std::endl;
    std::cout << "  " << std::left << std::setfill(' ') << std::setw(7) << "ID";
    std::cout << std::left << std::setfill(' ') << std::setw(21) << "Hostname";
    std::cout << std::left << std::setfill(' ') << std::setw(6) << "Core";
    std::cout << std::left << std::setfill(' ') << std::setw(23) << "O.S.";
    std::cout << std::left << std::setfill(' ') << std::setw(10) << "PID";
    std::cout << std::left << std::setfill(' ') << std::setw(20) << "CPU Brand";
    std::cout << std::endl;
    std::cout << "  " << std::left << std::setfill('-') << std::setw(108) << "-" << std::endl;

    if (istart == 0 || igeth != 0) {
        /*std::cout << "     n" << procRank_ << "   Unknown  " << 1 << "/" << numP << "  "
                << os_n.c_str() << "   " << iPID << "   " << cB.c_str() << std::endl;*/
        std::string idStr = "n" + toString(procRank_);
        std::cout << "  " << std::left << std::setfill(' ') << std::setw(7) << idStr.c_str();
        std::cout << std::left << std::setfill(' ') << std::setw(21) << "Unknown";
        std::string corStr = "1/" + toString(numP);
        std::cout << std::left << std::setfill(' ') << std::setw(6) << corStr.c_str();
        std::cout << std::left << std::setfill(' ') << std::setw(23) << os_n.c_str();
        std::cout << std::left << std::setfill(' ') << std::setw(10) << iPID;
        std::cout << std::left << std::setfill(' ') << std::setw(20) << cB.c_str();
        std::cout << std::endl;
    } else {
        /*std::cout << "     n" << procRank_ << "   " << hname << "  " << 1 << "/" << numP << "  "
                << os_n.c_str() << "   " << iPID << "   " << cB.c_str() << std::endl;*/
        std::string idStr = "n" + toString(procRank_);
        std::cout << "  " << std::left << std::setfill(' ') << std::setw(7) << idStr.c_str();
        std::cout << std::left << std::setfill(' ') << std::setw(21) << hname;
        std::string corStr = "1/" + toString(numP);
        std::cout << std::left << std::setfill(' ') << std::setw(6) << corStr.c_str();
        std::cout << std::left << std::setfill(' ') << std::setw(23) << os_n.c_str();
        std::cout << std::left << std::setfill(' ') << std::setw(10) << iPID;
        std::cout << std::left << std::setfill(' ') << std::setw(20) << cB.c_str();
        std::cout << std::endl;
    }

    std::cout << "  " << std::left << std::setfill('-') << std::setw(108) << "-" << std::endl;
    std::cout << std::right << std::setfill(' ');
#ifdef _WIN32

    WSACleanup();

#endif // _WIN32

#endif // MPI_PARALLEL
}

int OpenHurricane::HurMPIBase::getCount(Status *status, Datatype datatype, int *count) {
#ifdef MPI_PARALLEL
    return MPI_Get_count(status, datatype, count);
#else  // Not define MPI_PARALLEL
    return MPI_SUCCESS;
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::errorString(const int code, std::string &errString) {
#ifdef MPI_PARALLEL
    char *errChar = new char[MPI_MAX_ERROR_STRING]();
    int lenErr = 0;
    MPI_Error_string(code, errChar, &lenErr);
    errChar[lenErr] = '\0';
    errString = errChar;
#else  // Not define MPI_PARALLEL
    if (code == MPI_SUCCESS) {
        errString = "Success";
    } else {
        errString = "Unknown error";
    }
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::printErrorString(const int code) {
#ifdef MPI_PARALLEL
    std::string errMsg;
    errorString(code, errMsg);
    std::string nodeInfo = "Node: ";
    nodeInfo += std::to_string(getProcRank());
    error(errMsg, nodeInfo);
#else  // Not define MPI_PARALLEL
    std::string errMsg;
    errorString(code, errMsg);
    std::string nodeInfo = "Node: ";
    nodeInfo += std::to_string(getProcRank());
    error(errMsg, nodeInfo);
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::bcastString(std::string &bufStr, int root, Comm comm) {
#ifdef MPI_PARALLEL
    if (!parRun()) {
        return;
    }

    int sizeBuf = 0;

    if (!bufStr.empty()) {
        sizeBuf = int(bufStr.size());
    }

    bcast(&sizeBuf, 1, MPI_INT, root, comm);

    char *bufC = new char[sizeBuf + 1]();

    if (!bufStr.empty() && sizeBuf <= bufStr.size()) {
        strcpy(bufC, bufStr.c_str());
    }
    bufC[sizeBuf] = '\0';
    bcast(bufC, sizeBuf + 1, MPI_CHAR, root, comm);
    bufStr = bufC;
    delete[] bufC;
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::bcastString(std::string &bufStr, int root,const subCommtor& comm) {
#ifdef MPI_PARALLEL 
    if (!parRun()) {
        return;
    }

    int sizeBuf = 0;

    if (!bufStr.empty()) {
        sizeBuf = int(bufStr.size());
    }

    bcast(&sizeBuf, 1, MPI_INT, root, comm.subComm());

    char *bufC = new char[sizeBuf + 1]();

    if (!bufStr.empty() && sizeBuf <= bufStr.size()) {
        strcpy(bufC, bufStr.c_str());
    }
    bufC[sizeBuf] = '\0';
    bcast(bufC, sizeBuf + 1, MPI_CHAR, root, comm.subComm());
    bufStr = bufC;
    delete[] bufC;
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::gatherString(std::string &bufStr, int root, Comm comm) {
#ifdef MPI_PARALLEL
    if (!parRun()) {
        return;
    }
    int sizeBuf = static_cast<int>(bufStr.size());
    //sizeBuf++;
    //std::cout << "here" << std::endl;
    int *recvBUf = new int[getProcSize()];
    int *displs = new int[getProcSize()]();
    char *bufC = new char[sizeBuf + 1]();
    strcpy(bufC, bufStr.c_str());
    //bufC[sizeBuf - 1] = '\0';
    gather(&sizeBuf, 1, MPI_INT, recvBUf, 1, MPI_INT, root, comm);

    int totalSize = 0;
    if (isThisProc(root)) {
        for (int i = 0; i < getProcSize(); ++i) {
            totalSize += recvBUf[i];
            //std::cout << recvBUf[i] << std::endl;
        }
    }

    char *recvC = nullptr;
    if (isThisProc(root)) {
        recvC = new char[totalSize + 1]();
    }
    displs[0] = 0;
    for (integer ip = 1; ip < HurMPIBase::getProcSize(); ++ip) {
        displs[ip] = displs[ip - 1] + recvBUf[ip - 1];
    }
    gatherv(bufC, sizeBuf, MPI_CHAR, recvC, recvBUf, displs, MPI_CHAR, root, comm);
    if (isThisProc(root)) {
        recvC[totalSize] = '\0';
        bufStr = recvC;
        //std::cout << recvC << std::endl;
    }
    delete[] recvBUf;
    delete[] displs;
    delete[] bufC;

    if (isThisProc(root)) {
        delete[] recvC;
    }
#endif // MPI_PARALLEL
}

void OpenHurricane::HurMPIBase::gatherString(std::string &bufStr, int root, const subCommtor &comm) {
#ifdef MPI_PARALLEL
    if (!parRun()) {
        return;
    }
    int sizeBuf = static_cast<int>(bufStr.size());
    int *recvBUf = new int[getProcSize()];
    int *displs = new int[getProcSize()]();
    char *bufC = new char[sizeBuf + 1]();
    strcpy(bufC, bufStr.c_str());
    //bufC[sizeBuf - 1] = '\0';
    gather(&sizeBuf, 1, MPI_INT, recvBUf, 1, MPI_INT, root, comm.subComm());

    int totalSize = 0;
    if (comm.isThisProc(root)) {
        for (int i = 0; i < getProcSize(); ++i) {
            totalSize += recvBUf[i];
        }
    }

    char *recvC = nullptr;
    if (comm.isThisProc(root)) {
        recvC = new char[totalSize + 1]();
    }
    displs[0] = 0;
    for (integer ip = 1; ip < HurMPIBase::getProcSize(); ++ip) {
        displs[ip] = displs[ip - 1] + recvBUf[ip - 1];
    }
    gatherv(bufC, sizeBuf, MPI_CHAR, recvC, recvBUf, displs, MPI_CHAR, root, comm.subComm());
    if (comm.isThisProc(root)) {
        recvC[totalSize] = '\0';
        bufStr = recvC;
    }
    delete[] recvBUf;
    delete[] displs;
    delete[] bufC;

    if (comm.isThisProc(root)) {
        delete[] recvC;
    }
#endif // MPI_PARALLEL
}
