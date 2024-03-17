/*!
 * \file CRSMatrixProcessTransfer.hpp
 * \brief Headers of the CRS matrix transfer template.
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
#include "processTopologies.hpp"
#include "smartPointerList.hpp"

namespace OpenHurricane {

    class CRSMatrixAddrCutInterface {
    private:
        cutZoneList cutZones_;

        perZoneList perZones_;

        integerList intfCellId_;

    public:
        inline CRSMatrixAddrCutInterface() {}

        inline CRSMatrixAddrCutInterface(const cutZoneList &cutZones, const perZoneList &perZones,
                                         const integerList &intfCellId)
            : cutZones_(cutZones), perZones_(perZones), intfCellId_(intfCellId) {}

        inline CRSMatrixAddrCutInterface(cutZoneList &cutZones, perZoneList &perZones,
                                         integerList &intfCellId, bool isTransfer) {
            if (isTransfer) {
                cutZones_.swap(cutZones);
                perZones_.swap(perZones);
                intfCellId_.swap(intfCellId);
            } else {
                cutZones_ = cutZones;
                perZones_ = perZones;
                intfCellId_ = intfCellId;
            }
        }

        inline CRSMatrixAddrCutInterface(const CRSMatrixAddrCutInterface &other)
            : cutZones_(other.cutZones_), perZones_(other.perZones_),
              intfCellId_(other.intfCellId_) {}

        inline CRSMatrixAddrCutInterface &operator=(const CRSMatrixAddrCutInterface &other) {
            if (this != std::addressof(other)) {
                cutZones_ = other.cutZones_;
                perZones_ = other.perZones_;
                intfCellId_ = other.intfCellId_;
            }
            return *this;
        }

        inline CRSMatrixAddrCutInterface(CRSMatrixAddrCutInterface &&other) noexcept
            : cutZones_(std::move(other.cutZones_)), perZones_(std::move(other.perZones_)),
              intfCellId_(std::move(other.intfCellId_)) {}

        inline CRSMatrixAddrCutInterface &operator=(CRSMatrixAddrCutInterface &&other) noexcept {
            cutZones_.swap(other.cutZones_);
            perZones_.swap(other.perZones_);
            intfCellId_.swap(other.intfCellId_);
            return *this;
        }

        inline ~CRSMatrixAddrCutInterface() noexcept {}

        inline void clear() noexcept {
            cutZones_.clear();
            perZones_.clear();
            intfCellId_.clear();
        }

        inline void transfer(CRSMatrixAddrCutInterface &other) noexcept {
            cutZones_.transfer(other.cutZones_);
            perZones_.transfer(other.perZones_);
            intfCellId_.transfer(other.intfCellId_);
        }

        hur_nodiscard inline const cutZoneList &cutZones() const noexcept { return cutZones_; }

        hur_nodiscard inline const perZoneList &perZones() const noexcept { return perZones_; }

        hur_nodiscard inline const integerList &intfCellId() const noexcept { return intfCellId_; }
    };

    template <class Type> class CRSMatrixProcessTransfer {
    private:
        const CRSMatrixAddrCutInterface &proceIntf_;

        List<HurMPIBase::Request> requestCut_;
        uniquePtrList<cutProcessReceiv<Type>> cutPRecvs_;
        uniquePtrList<cutProcessSend<Type>> cutPSends_;

        List<HurMPIBase::Request> requestPer_;
        uniquePtrList<perProcessReceiv<Type>> perPRecvs_;
        uniquePtrList<perProcessSend<Type>> perPSends_;
        uniquePtrList<periodicTransferSameProc<Type>> perSameProc_;

    public:
        CRSMatrixProcessTransfer() = delete;

        CRSMatrixProcessTransfer(const CRSMatrixAddrCutInterface &proceIntf)
            : proceIntf_(proceIntf) {
            const auto &cutZ = proceIntf_.cutZones();
            requestCut_.resize(cutZ.size());

            integer cntCutRecv = 0;
            integer cntCutSend = 0;
            for (const auto &e : cutZ) {
                if (e.isSendFromThisProc()) {
                    cntCutSend++;
                } else if (e.isThisProcReceiv()) {
                    cntCutRecv++;
                }
            }
            cutPSends_.resize(cntCutSend);
            cutPRecvs_.resize(cntCutRecv);
            cntCutRecv = 0;
            cntCutSend = 0;
            for (const auto &e : cutZ) {
                if (e.isSendFromThisProc()) {
                    cutPSends_.set(cntCutSend++, new cutProcessSend<Type>(e));
                } else if (e.isThisProcReceiv()) {
                    cutPRecvs_.set(cntCutRecv++, new cutProcessReceiv<Type>(e));
                }
            }

            const auto &perZ = proceIntf_.perZones();
            integer cntPPerZ = 0;
            for (integer pi = 0; pi < perZ.size(); ++pi) {
                if (perZ[pi].isSameProc()) {
                    continue;
                }
                cntPPerZ++;
            }
            requestPer_.resize(cntPPerZ);

            integer cntPerRecv = 0;
            integer cntPerSend = 0;
            integer cntPerSameProc = 0;
            for (const auto &pe : perZ) {
                if (pe.isSameProc()) {
                    cntPerSameProc++;
                } else if (pe.isSendFromThisProc()) {
                    cntPerSend++;
                } else if (pe.isThisProcReceiv()) {
                    cntPerRecv++;
                }
            }
            perSameProc_.resize(cntPerSameProc);
            perPSends_.resize(cntPerSend);
            perPRecvs_.resize(cntPerRecv);
            cntPerRecv = 0;
            cntPerSend = 0;
            cntPerSameProc = 0;
            for (const auto &pe : perZ) {
                if (pe.isSameProc()) {
                    perSameProc_.set(cntPerSameProc++, new periodicTransferSameProc<Type>(pe));
                } else if (pe.isSendFromThisProc()) {
                    perPSends_.set(cntPerSend++, new perProcessSend<Type>(pe));
                } else if (pe.isThisProcReceiv()) {
                    perPRecvs_.set(cntPerRecv++, new perProcessReceiv<Type>(pe));
                }
            }
        }

        CRSMatrixProcessTransfer(const CRSMatrixProcessTransfer &other) = delete;
        CRSMatrixProcessTransfer &operator=(const CRSMatrixProcessTransfer &other) = delete;

        inline ~CRSMatrixProcessTransfer() noexcept {}

        void transferData(const Array<Type> &var, Array<Type> &recvVar) {
            if (!HurMPIBase::parRun()) {
                for (integer i = 0; i < perSameProc_.size(); ++i) {
                    perSameProc_[i].transferData(var, recvVar);
                }
            } else {
                integer cntCutRe = 0;
                for (integer i = 0; i < cutPRecvs_.size(); ++i) {
                    cutPRecvs_[i].recvNonBlock(&requestCut_[cntCutRe++]);
                }

                for (integer i = 0; i < cutPSends_.size(); ++i) {
                    cutPSends_[i].writeBUf(var);
                    cutPSends_[i].sendNonBlock(&requestCut_[cntCutRe++]);
                }

                integer cntPerRe = 0;
                for (integer i = 0; i < perPRecvs_.size(); ++i) {
                    perPRecvs_[i].recvNonBlock(&requestPer_[cntPerRe++]);
                }
                for (integer i = 0; i < perPSends_.size(); ++i) {
                    perPSends_[i].writeBUf(var);
                    perPSends_[i].sendNonBlock(&requestPer_[cntPerRe++]);
                }
                for (integer i = 0; i < perSameProc_.size(); ++i) {
                    perSameProc_[i].transferData(var, recvVar);
                }
                HurMPIBase::waitall(requestCut_.size(), requestCut_.data(), MPI_STATUSES_IGNORE);

                for (integer i = 0; i < cutPRecvs_.size(); ++i) {
                    cutPRecvs_[i].readBUf(recvVar);
                }

                HurMPIBase::waitall(requestPer_.size(), requestPer_.data(), MPI_STATUSES_IGNORE);
                for (integer i = 0; i < perPRecvs_.size(); ++i) {
                    perPRecvs_[i].readBUf(recvVar);
                }
            }
        }
    };

    template <>
    void CRSMatrixProcessTransfer<realArray>::transferData(const Array<realArray> &var,
                                                           Array<realArray> &recvVar);
} // namespace OpenHurricane