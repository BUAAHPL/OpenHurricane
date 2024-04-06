/*!
 * \file HurMPIBase.hpp
 * \brief Headers of the mpi interface for generalized datatypes.
 *        The subroutines and functions are in the <i>HurMPIBase.cpp</i> file.
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

#ifdef MPI_PARALLEL
#include <mpi.h>
typedef MPI_Comm HurComm;
#else // MPI_PARALLEL
#define MPI_COMM_WORLD 0
#define MPI_UNSIGNED_LONG 1
#define MPI_LONG 2
#define MPI_UNSIGNED_SHORT 3
#define MPI_DOUBLE 4
#define MPI_ANY_SOURCE 5
#define MPI_SUM 6
#define MPI_CHAR 7
#define MPI_SHORT 8
#define MPI_MIN 9
#define MPI_MAX 10
#define MPI_INT 11
#define MPI_INT32_T 12
#define MPI_INT64_T 13
#define MPI_FLOAT 14
#define MPI_CXX_BOOL 15
#define MPI_UINT32_T 16
#define MPI_UINT64_T 17
#define MPI_SUCCESS 0
#define MPI_INFO_NULL 0

typedef struct MPI_Status {
    int MPI_TAG;
    int MPI_SOURCE;
    MPI_Status() : MPI_TAG(0), MPI_SOURCE(0) {}
} MPI_Status;
typedef int HurComm;
#define MPI_STATUS_IGNORE (HurMPIBase::Status *)1
#define MPI_STATUSES_IGNORE (HurMPIBase::Status *)1

#endif // MPI_PARALLEL

#include "preset.hpp"
#include "systemInfo.hpp"
#include <string>

namespace OpenHurricane {
    //template <class T> class List;
    template <class PrimitiveType> class feature;
    const int MASTER_PROCCESS = 0; /*!< \brief Master node for MPI parallelization. */
    const int SINGLE_NODE = 1;     /*!< \brief There is only a node in the MPI parallelization. */
    class HurMPIBase {
    public:
#ifdef MPI_PARALLEL
        using Request = MPI_Request;
        using Status = MPI_Status;
        using Datatype = MPI_Datatype;
        using Op = MPI_Op;
        using Comm = MPI_Comm;
        using Win = MPI_Win;
        using Info = MPI_Info;

#else  // Not define MPI_PARALLEL
        using Request = int;
        using Status = MPI_Status;
        using Datatype = int;
        using Op = int;
        using Comm = int;
        using Win = int;
        using Info = int;
#endif // MPI_PARALLEL

        class subCommtor {
        public:
            static constexpr int masterInSub_ = 0;
        private:
            /**
             * \brief MPI Communicators for sub-region.
             */
            Comm subComm_;

            /**
             * \brief MPI process rank in the sur-region.
             */
            int subProcRank_;

            /**
             * \brief Number MPI processes in the sur-region.
             */
            int subProcSize_;

        public:
            hur_nodiscard static const subCommtor &null();

            constexpr inline subCommtor(bool nullSet)
                : subComm_(-1), subProcRank_(-1), subProcSize_(-1) {}

            constexpr inline subCommtor() : subComm_(0), subProcRank_(0), subProcSize_(0) {}

            subCommtor(const Comm comm, const int color, const int key);

            constexpr inline subCommtor(const subCommtor &other)
                : subComm_(other.subComm_), subProcRank_(other.subProcRank_),
                  subProcSize_(other.subProcSize_) {}

            constexpr inline subCommtor &operator=(const subCommtor &other) {
                if (this != std::addressof(other)) {
                    subComm_ = other.subComm_;
                    subProcRank_ = other.subProcRank_;
                    subProcSize_ = other.subProcSize_;
                }
                return *this;
            }

            inline ~subCommtor() noexcept {}

            /**
             * \brief MPI Communicators for sub-region.
             */
            hur_nodiscard constexpr inline Comm subComm() const noexcept { return subComm_; }

            /**
             * \brief MPI process rank in the sur-region.
             */
            hur_nodiscard constexpr inline int subProcRank() const noexcept { return subProcRank_; }

            /**
             * \brief Number MPI processes in the sur-region.
             */
            hur_nodiscard constexpr inline int subProcSize() const noexcept { return subProcSize_; }

            /**
             * \brief MPI Communicators for sub-region.
             */
            constexpr inline void setSubComm(const Comm comm) noexcept { subComm_ = comm; }

            /**
             * \brief MPI process rank in the sur-region.
             */
            constexpr inline void setSubProcRank(const int spr) noexcept { subProcRank_ = spr; }

            /**
             * \brief Number MPI processes in the sur-region.
             */
            constexpr inline void setSubProcSize(const int sps) noexcept { subProcSize_ = sps; }

            constexpr inline bool isMasterInSub() const noexcept {
                return subProcRank_ == masterInSub_;
            }

            constexpr inline int masterNoInSub() const noexcept { return masterInSub_; }

            constexpr inline bool isThisProc(const int proc) const noexcept {
                return subProcRank_ == proc;
            }

            void freeSubComm();
        };

    protected:
#ifdef MPI_PARALLEL
#else  // Not define MPI_PARALLEL
        static void copyData(const void *sendbuf, void *recvbuf, int size, Datatype datatype);
#endif // MPI_PARALLEL
        /*!\brief The world rank (id) of this process.*/
        static int procRank_;

        /*!\brief The total number of the parallel processes in this program (world).*/
        static int procSize_;

        static int minRank_;
        static Comm HurWorldComm_;
        static bool usingWinMinRank_;
        static Win winMinRank_;

        /** \brief The total number of machines (nodes). */
        static int node_;

        /** \brief The node (machine) id of this process. */
        static int procNode_;

        static Info info_;

        /*!\brief By default this is not a multi-node (multi-machines) run.*/
        static bool multiNode_;

        /*!\brief By default this is not a parallel run.*/
        static bool parRun_;

        /*!\brief Have support for threads?.*/
        static bool haveThreads_;

        static bool isInitialized_;

        /*!\brief The communicator in the current node.*/
        static subCommtor nodeCommtor_;

        /*!\brief Set data for parallel running.*/
        static inline void setParRun() noexcept {
            if (procSize_ == 1) {
                parRun_ = false;
            } else {
                parRun_ = true;
            }
        }

    public:
        /*!\brief Process index of the master.*/
        static inline int masterNo() noexcept { return MASTER_PROCCESS; }

        /*!\brief Is this process the master process.*/
        static inline bool master() noexcept { return getProcRank() == masterNo(); }
        /*!\brief Is this process the slaver process.*/
        static inline bool isSlaver() noexcept { return getProcRank() != masterNo(); }

        /*!\brief Return the rank (id) of this process.*/
        static hur_nodiscard inline int getProcRank() noexcept { return procRank_; }
        /*!\brief Return the total number of the parallel processes in this program.*/
        static hur_nodiscard inline int getProcSize() noexcept { return procSize_; }
        static hur_nodiscard inline bool isInitialized() noexcept { return isInitialized_; }

        /**
         * \brief Return World communicator.
         */
        static hur_nodiscard inline Comm getComm() noexcept { return HurWorldComm_; }
        /**
         * \brief Return node communicator.
         */
        static hur_nodiscard inline Comm getNodeComm() noexcept { return nodeCommtor_.subComm(); }
        static hur_nodiscard inline Info getInfo() noexcept { return info_; }
        /*!\brief The rank (id) of this process of the current node.*/
        static hur_nodiscard int procRankNode() noexcept { return nodeCommtor_.subProcRank(); }
        /*!\brief The total number of the parallel processes in the current node.*/
        static hur_nodiscard int procSizeNode() noexcept { return nodeCommtor_.subProcSize(); }

        static hur_nodiscard inline const subCommtor &nodeCommtor() noexcept { return nodeCommtor_; }

        static hur_nodiscard bool masterInThisNode() noexcept {
            return procRankNode() == MASTER_PROCCESS;
        }
        static hur_nodiscard bool slaverInThisNode() noexcept {
            return procRankNode() != MASTER_PROCCESS;
        }
        static inline bool isThisProc(const int ip) noexcept { return getProcRank() == ip; }
        /*!\brief Has multi-nodes?*/
        static inline bool multiNodes() noexcept {
#ifdef MPI_PARALLEL
            return node_ > 1;
#else  // Not define MPI_PARALLEL
            return false;
#endif // MPI_PARALLEL
        }
        /** \brief The total number of machines (nodes). */
        static int getNode() noexcept { return node_; }

        /** \brief The node (machine) id of this process. */
        static int procNode() noexcept { return procNode_; }
        static void error(const std::string &ErrorMsg, const std::string &FunctionName,
                          bool isAbort = true);
        static void warning(const std::string &ErrorMsg, const std::string &FunctionName);
        static hur_nodiscard bool isCollectiveCall(double times = 0.5);

    protected:
        static hur_nodiscard int getMinRank(Win winMinRank = minRank_);

        static void setThreads(const bool haveThreads);

    public:
        /*!\brief Have support for threads.*/
        static inline bool haveThreads() noexcept { return haveThreads_; }

        /*!\brief Is this a main thread?*/
        static inline bool mainThread() {
#ifdef MPI_PARALLEL
            int flag;
            MPI_Is_thread_main(&flag);
            return flag;
#else  // Not define MPI_PARALLEL
            return true;
#endif // MPI_PARALLEL
        }
        /*!\brief Is this a parallel run?*/
        static inline bool &parRun() noexcept { return parRun_; }

        static int init(int *argc, char ***argv);
        static int finalize();
        static int bufferAttach(void *buffer, int size);
        static int bufferDetach(void *buffer, int *size);
        static inline int barrier(Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Barrier(comm);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }
        static int commRank(Comm comm, int &rank);
        static int commSize(Comm comm, int &size);

        /*!\brief Abort program.*/
        static int abort(Comm comm = HurWorldComm_, int error = EXIT_FAILURE);

    protected:
        static void setMultiNode();
        static void setNodeComm();

    public:
        static void getPcInfo();
        static int getCount(Status *status, Datatype datatype, int *count);
        static void errorString(const int code, std::string &errString);

        /*!\brief Print the error string to screen and exit the program.*/
        static void printErrorString(const int code);
        static int wait(Request *request, Status *status = MPI_STATUS_IGNORE) {
#ifdef MPI_PARALLEL
            return MPI_Wait(request, status);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int waitall(int nrequests, Request *request,
                                  Status *status = MPI_STATUS_IGNORE) {
#ifdef MPI_PARALLEL
            return MPI_Waitall(nrequests, request, status);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int waitany(int nrequests, Request *request, int *index,
                                  Status *status = MPI_STATUS_IGNORE) {
#ifdef MPI_PARALLEL
            return MPI_Waitany(nrequests, request, index, status);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int testall(int count, Request *array_of_requests, int *flag,
                                  Status *array_of_statuses = MPI_STATUS_IGNORE) {
#ifdef MPI_PARALLEL
            return MPI_Testall(count, array_of_requests, flag, array_of_statuses);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int test(Request *array_of_requests, int *flag,
                               Status *array_of_statuses = MPI_STATUS_IGNORE) {
#ifdef MPI_PARALLEL
            return MPI_Test(array_of_requests, flag, array_of_statuses);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int probe(int source, int tag, Comm comm,
                                Status *status = MPI_STATUS_IGNORE) {
#ifdef MPI_PARALLEL
            return MPI_Probe(source, tag, comm, status);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int iprobe(int source, int tag, Comm comm, int *flag,
                                 Status *status = MPI_STATUS_IGNORE) {
#ifdef MPI_PARALLEL
            return MPI_Iprobe(source, tag, comm, flag, status);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        /**
         * \brief Performs a standard mode send operation and returns when the send buffer can be safely reused.
         */
        static inline int send(void *buf, int count, Datatype datatype, int dest, int tag,
                               Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Send(buf, count, datatype, dest, tag, comm);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        /**
         * \brief Sends data to a specified process in buffered mode. This function returns when the send buffer can be safely reused.
         */
        static inline int bsend(void *buf, int count, Datatype datatype, int dest, int tag,
                                Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Bsend(buf, count, datatype, dest, tag, comm);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        /**
         * \brief Performs a synchronous mode send operation and returns when the send buffer can be safely reused.
         */
        static inline int ssend(void *buf, int count, Datatype datatype, int dest, int tag,
                                Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Ssend(buf, count, datatype, dest, tag, comm);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        /**
         * \brief Performs a ready mode send operation and returns when the send buffer can be safely reused..
         */
        static inline int rsend(void *buf, int count, Datatype datatype, int dest, int tag,
                                Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Rsend(buf, count, datatype, dest, tag, comm);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        /**
         * \brief Initiates a standard mode send operation and returns a handle to the requested communication operation.
         */
        static inline int isend(void *buf, int count, Datatype datatype, int dest, int tag,
                                Comm comm, Request *request) {
#ifdef MPI_PARALLEL
            return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int ibsend(void *buf, int count, Datatype datatype, int dest, int tag,
                                 Comm comm, Request *request) {
#ifdef MPI_PARALLEL
            return MPI_Ibsend(buf, count, datatype, dest, tag, comm, request);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int issend(void *buf, int count, Datatype datatype, int dest, int tag,
                                 Comm comm, Request *request) {
#ifdef MPI_PARALLEL
            return MPI_Issend(buf, count, datatype, dest, tag, comm, request);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int irsend(void *buf, int count, Datatype datatype, int dest, int tag,
                                 Comm comm, Request *request) {
#ifdef MPI_PARALLEL
            return MPI_Irsend(buf, count, datatype, dest, tag, comm, request);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        /**
         * \brief Performs a receive operation and does not return until a matching message is received.
         */
        static inline int recv(void *buf, int count, Datatype datatype, int source, int tag,
                               Comm comm = HurWorldComm_, Status *status = MPI_STATUS_IGNORE) {
#ifdef MPI_PARALLEL
            return MPI_Recv(buf, count, datatype, source, tag, comm, status);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        /**
         * \brief Initiates a receive operation and returns a handle to the requested communication operation..
         */
        static inline int irecv(void *buf, int count, Datatype datatype, int source, int tag,
                                Comm comm, Request *request) {
#ifdef MPI_PARALLEL
            return MPI_Irecv(buf, count, datatype, source, tag, comm, request);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int sendRecv(void *sendbuf, int sendcnt, Datatype sendtype, int dest,
                                   int sendtag, void *recvbuf, int recvcnt, Datatype recvtype,
                                   int source, int recvtag, Comm comm = HurWorldComm_,
                                   Status *status = MPI_STATUSES_IGNORE) {
#ifdef MPI_PARALLEL
            return MPI_Sendrecv(sendbuf, sendcnt, sendtype, dest, sendtag, recvbuf, recvcnt,
                                recvtype, source, recvtag, comm, status);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcnt, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int scan(void *sendbuf, void *recvbuf, int count, Datatype datatype, Op op,
                               Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, count, datatype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int bcast(void *buf, int count, Datatype datatype, int root = MASTER_PROCCESS,
                                Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Bcast(buf, count, datatype, root, comm);
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }
        static void bcastString(std::string &bufStr, int root = MASTER_PROCCESS,
                                Comm comm = HurWorldComm_);

         static void bcastString(std::string &bufStr, int root, const subCommtor &comm);

        static inline int reduce(void *sendbuf, void *recvbuf, int count, Datatype datatype, Op op,
                                 int root = MASTER_PROCCESS, Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, count, datatype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        /*!
         *\brief Reduce to master node in default comm (Only for global communicators).
         * The reduce result would be stored in sendbuf in master node.
         * The value of sendbuf on slaver nodes is remained.
         */
        template <class Type>
        static inline int reduce(Type &sendbuf, Op op, int root = MASTER_PROCCESS,
                                 Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            if (!parRun()) {
                return MPI_SUCCESS;
            }
            Type recv;
            int flag = reduce(&sendbuf, &recv, feature<Type>::nElements_, feature<Type>::MPIType,
                              op, root, HurWorldComm_);

            if (isThisProc(root)) {
                sendbuf = recv;
            }
            return flag;
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class Type>
        static inline int reduce(Type &sendbuf, Op op, int root, const subCommtor &comm) {
#ifdef MPI_PARALLEL
            if (!parRun()) {
                return MPI_SUCCESS;
            }
            Type recv;
            int flag = reduce(&sendbuf, &recv, feature<Type>::nElements_, feature<Type>::MPIType,
                              op, root, comm.subComm());

            if (comm.isThisProc(root)) {
                sendbuf = recv;
            }
            return flag;
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int allReduce(void *sendbuf, void *recvbuf, int count, Datatype datatype,
                                    Op op, Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, count, datatype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        /*!
         *\brief Reduce to master node in default comm.
         * The reduce result would be stored in sendbuf in all nodes.
         * Please be noticed this function would change the value of sendbuf in all nodes.
         */
        template <class Type>
        static inline int allReduce(Type &sendbuf, Op op, Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            if (!parRun()) {
                return MPI_SUCCESS;
            }
            Type recv;
            int flag = allReduce(&sendbuf, &recv, feature<Type>::nElements_, feature<Type>::MPIType,
                                 op, HurWorldComm_);
            sendbuf = recv;
            return flag;
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int gather(void *sendbuf, int sendcnt, Datatype sendtype, void *recvbuf,
                                 int recvcnt, Datatype recvtype, int root = MASTER_PROCCESS,
                                 Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcnt, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int igather(void *sendbuf, int sendcnt, Datatype sendtype, void *recvbuf,
                                  int recvcnt, Datatype recvtype, int root, Comm comm,
                                  Request *request) {
#ifdef MPI_PARALLEL
            return MPI_Igather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm,
                               request);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcnt, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int gatherv(const void *sendbuf, const int sendcnt, Datatype sendtype,
                                  void *recvbuf, const int *recvcnt, const int *displs,
                                  Datatype recvtype, int root = MASTER_PROCCESS,
                                  Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Gatherv(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, displs, recvtype, root,
                               comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcnt, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int igatherv(const void *sendbuf, const int sendcnt, Datatype sendtype,
                                   void *recvbuf, const int *recvcnt, const int *displs,
                                   Datatype recvtype, int root, Comm comm, Request *request) {
#ifdef MPI_PARALLEL
            return MPI_Igatherv(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, displs, recvtype,
                                root, comm, request);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcnt, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        } 

        /**
         * \brief Gather string from global communicators.
         */
        static void gatherString(std::string &bufStr, int root = MASTER_PROCCESS,
                                 Comm comm = HurWorldComm_);
        
        /**
         * \brief Gather string from sub-communicators.
         */
        static void gatherString(std::string &bufStr, int root, const subCommtor &comm);

        static inline int scatter(void *sendbuf, int sendcnt, Datatype sendtype, void *recvbuf,
                                  int recvcnt, Datatype recvtype, int root = MASTER_PROCCESS,
                                  Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcnt, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int scatterv(const void *sendbuf, const int *sendcnt, const int *displs,
                                   Datatype sendtype, void *recvbuf, int recvcnt, Datatype recvtype,
                                   int root = MASTER_PROCCESS, Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Scatterv(sendbuf, sendcnt, displs, sendtype, recvbuf, recvcnt, recvtype,
                                root, comm);
#else  // Not define MPI_PARALLEL
            copyData(const_cast<void *>(sendbuf), recvbuf, *const_cast<int *>(sendcnt), sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int reduceScatter(void *sendbuf, void *recvbuf, int *recvcounts,
                                        Datatype datatype, Op op, Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, *recvcounts, datatype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int allGather(void *sendbuf, int sendcnt, Datatype sendtype, void *recvbuf,
                                    int recvcnt, Datatype recvtype, Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Allgather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcnt, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int iallGather(void *sendbuf, int sendcnt, Datatype sendtype, void *recvbuf,
                                     int recvcnt, Datatype recvtype, Comm comm, Request *request) {
#ifdef MPI_PARALLEL
            return MPI_Iallgather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, comm,
                                  request);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcnt, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int allGatherv(const void *sendbuf, int sendcount, Datatype sendtype,
                                     void *recvbuf, int *recvcounts, int *displs, Datatype recvtype,
                                     Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                                  recvtype, comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcount, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int iallGatherv(const void *sendbuf, int sendcount, Datatype sendtype,
                                      void *recvbuf, int *recvcounts, int *displs,
                                      Datatype recvtype, Comm comm, Request *request) {
#ifdef MPI_PARALLEL
            return MPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                                   recvtype, comm, request);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcount, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        static inline int allToAll(void *sendbuf, int sendcount, Datatype sendtype, void *recvbuf,
                                   int recvcount, Datatype recvtype, Comm comm = HurWorldComm_) {
#ifdef MPI_PARALLEL
            return MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
#else  // Not define MPI_PARALLEL
            copyData(sendbuf, recvbuf, sendcount, sendtype);
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class T> static inline int allReduceVS(T &sendbuf, Op op) {
            if (!parRun()) {
                return MPI_SUCCESS;
            }
            T recv;
            int flag = allReduce(&sendbuf[0], &recv[0], feature<T>::nElements_, feature<T>::MPIType,
                                 op, HurWorldComm_);

            sendbuf = recv;
            return flag;
        }
        /*!\brief Barrier all processes in the current comm (getComm()).*/
        static inline int barrierDefaultComm() {
#ifdef MPI_PARALLEL
            return MPI_Barrier(getComm());
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }

        template <class T> static inline int reduceVS(T &sendbuf, Op op) {
#ifdef MPI_PARALLEL
            if (!parRun()) {
                return MPI_SUCCESS;
            }
            T recv;
            int flag = reduce(&sendbuf[0], &recv[0], feature<T>::nElements_, feature<T>::MPIType,
                              op, masterNo(), HurWorldComm_);

            if (master()) {
                sendbuf = recv;
            }
            return flag;
#else  // Not define MPI_PARALLEL
            return MPI_SUCCESS;
#endif // MPI_PARALLEL
        }
    };
    //using HurMPI = HurMPIBase;
} // namespace OpenHurricane