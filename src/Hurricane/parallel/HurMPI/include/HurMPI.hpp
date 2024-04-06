/*!
 * \file HurMPI.hpp
 * \brief Headers of the mpi interface for generalized datatypes.
 *        The subroutines and functions are in the <i>HurMPI.cpp</i> file.
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

#include "HurMPIBase.hpp"
#include "List.hpp"
#include "preset.hpp"

namespace OpenHurricane {

    class HurMPI : public HurMPIBase {
    public:
        /**
         * \brief Get the process id List of each machines.
         */
        hur_nodiscard List<List<int>> getProIdOnNodes();

        template <class Type>
        static int sendList(List<Type> &sendL, int dest, int tag, Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return send(sendL.data(), sendL.size(), feature<Type>::MPIType, dest, tag, comm);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class Type>
        static int recvList(List<Type> &bufList, int source, int tag, Comm comm = HurWorldComm_,
                            Status *status = MPI_STATUSES_IGNORE) {
#ifdef MPI_PARALLEL
            return recv(bufList.data(), bufList.size(), feature<Type>::MPIType, source, tag, comm,
                        status);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class Type>
        static int bcastList(List<Type> &bufList, int root = MASTER_PROCCESS,
                             Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            if (!HurMPIBase::parRun()) {
                return MPI_SUCCESS;
            }
            if (feature<Type>::nElements_ == 1) {
                return bcast(bufList.data(), bufList.size(),
                             feature<typename feature<Type>::elementType>::MPIType, root, comm);
            } else {
                LFatal("Please use bcasVectortList()");
                return MPI_SUCCESS;
            }
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class Type>
        static int bcastVectorList(List<Type> &bufList, int root = MASTER_PROCCESS,
                                   Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            if (!HurMPIBase::parRun()) {
                return MPI_SUCCESS;
            }
            typename feature<Type>::elementType *send =
                new typename feature<Type>::elementType[bufList.size() * feature<Type>::nElements_];
            if (isThisProc(root)) {
                for (int i = 0; i < bufList.size(); ++i) {
                    for (int j = 0; j < feature<Type>::nElements_; ++j) {
                        send[i * feature<Type>::nElements_ + j] = bufList[i][j];
                    }
                }
            }
            auto iflag = bcast(send, bufList.size(), feature<Type>::MPIType, root, comm);
            if (!isThisProc(root)) {
                for (int i = 0; i < bufList.size(); ++i) {
                    for (int j = 0; j < feature<Type>::nElements_; ++j) {
                        bufList[i][j] = send[i * feature<Type>::nElements_ + j];
                    }
                }
            }
            delete[] send;
            return iflag;
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class Type>
        static int reduceList(List<Type> &sendbuf, Op op, int root = MASTER_PROCCESS,
                              Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            if (!parRun()) {
                return MPI_SUCCESS;
            }

            if (sendbuf.size() == 0) {
                return MPI_SUCCESS;
            }
            if (feature<Type>::nElements_ == 1) {
                List<Type> recv(sendbuf.size());
                int flag = reduce(sendbuf.data(), recv.data(), sendbuf.size(),
                                  feature<Type>::MPIType, op, root, comm);

                if (isThisProc(root)) {
                    sendbuf.transfer(recv);
                }
                return flag;
            } else {
                LFatal("Please use reduceVectorList()");
                return MPI_SUCCESS;
            }
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class Type>
        static int reduceVectorList(List<Type> &sendbuf, Op op, int root = MASTER_PROCCESS,
                                    Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            if (!parRun()) {
                return MPI_SUCCESS;
            }
            if (sendbuf.size() == 0) {
                return MPI_SUCCESS;
            }
            List<typename feature<Type>::value_type> send(sendbuf.size() *
                                                          feature<Type>::nElements_);
            List<typename feature<Type>::value_type> recv(sendbuf.size() *
                                                          feature<Type>::nElements_);

            for (int i = 0; i < sendbuf.size(); ++i) {
                for (int j = 0; j < feature<Type>::nElements_; ++j) {
                    send[i * feature<Type>::nElements_ + j] = sendbuf[i][j];
                }
            }
            int flag = reduce(send.data(), recv.data(), send.size(),
                              feature<typename feature<Type>::value_type>::MPIType, op, root, comm);
            if (isThisProc(root)) {
                for (int i = 0; i < sendbuf.size(); ++i) {
                    for (int j = 0; j < feature<Type>::nElements_; ++j) {
                        sendbuf[i][j] = send[i * feature<Type>::nElements_ + j];
                    }
                }
            }
            return flag;
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class Type>
        static int allReduceList(List<Type> &sendbuf, Op op, int root = MASTER_PROCCESS,
                                 Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            if (!parRun()) {
                return MPI_SUCCESS;
            }

            if (sendbuf.size() == 0) {
                return MPI_SUCCESS;
            }
            if (feature<Type>::nElements_ == 1) {
                List<Type> recv(sendbuf.size());
                int flag = allReduce(sendbuf.data(), recv.data(), sendbuf.size(),
                                     feature<Type>::MPIType, op, comm);

                sendbuf.transfer(recv);
                return flag;
            } else {
                LFatal("Number of elements is greater than 1");
                return MPI_SUCCESS;
            }
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class Type>
        static int gatherList(List<Type> &l, int root = MASTER_PROCCESS,
                              Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            if (!parRun()) {
                return MPI_SUCCESS;
            }
            Type sendBuf = l.operator[](getProcRank());
            return gather(&sendBuf, 1, feature<Type>::MPIType, l.data(), 1, feature<Type>::MPIType,
                          root, comm);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class Type> static int allGatherList(List<Type> &l, Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            if (!parRun()) {
                return MPI_SUCCESS;
            }
            Type sendBuf = l.operator[](getProcRank());
            return allGather(&sendBuf, 1, feature<Type>::MPIType, l.data(), 1,
                             feature<Type>::MPIType, comm);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }
    };
} // namespace OpenHurricane