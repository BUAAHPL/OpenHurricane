/*!
 * \file zoneMesh.cpp
 * \brief The subroutines and functions of zone mesh
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

#include "zoneMesh.hpp"
#include "matrices.hpp"

OpenHurricane::zone &OpenHurricane::zone::operator=(const zone &z) {
    if (this != std::addressof(z)) {
        name_ = z.name_;
        index_ = z.index_;
        firstIndex_ = z.firstIndex_;
        lastIndex_ = z.lastIndex_;
    }
    return *this;
}

OpenHurricane::pointZone::pointZone(const string &zName, const integer id, const integer firstID,
                                const integer lastID, const Options::geometryModel nD)
    : zone(zName, id, firstID, lastID), type_(ANYTYPE), ND_(nD) {
    if ((ND_ != Options::TWOD) && (ND_ != Options::THREED)) {
        std::string errMsg = "Wrong dimension set: ";
        errMsg = errMsg + std::to_string(ND_);
        errorAbortStr(errMsg);
    }
}

OpenHurricane::pointZone::pointZone(const integer id, const integer firstID, const integer lastID,
                                const short nD, const short t)
    : zone(" ", id, firstID, lastID), type_(ANYTYPE), ND_(Options::geometryModel::TWOD) {
    if (nD == Options::geometryModel::TWOD) {
        ND_ = Options::geometryModel::TWOD;
    } else if (nD == Options::geometryModel::THREED) {
        ND_ = Options::geometryModel::THREED;
    } else {
        std::string errMsg = "Wrong dimension set: ";
        errMsg = errMsg + std::to_string(ND_);
        errorAbortStr(errMsg);
    }

    if (t == pointZone::ANYTYPE) {
        type_ = pointZone::ANYTYPE;
    } else if (t == pointZone::VIRTUALNODES) {
        type_ = pointZone::VIRTUALNODES;
    } else if (t == pointZone::BOUNDARYNODES) {
        type_ = pointZone::BOUNDARYNODES;
    } else {
        std::string errMsg = "Wrong nodes type: ";
        errMsg = errMsg + std::to_string(t);
        errorAbortStr(errMsg);
    }
}

OpenHurricane::pointZone::pointZone(const pointZone &pZ) : zone(pZ), type_(pZ.type_), ND_(pZ.ND_) {}

void OpenHurricane::pointZone::setPointType(const short pT) {
    switch (pT) {
    case 0:
        type_ = VIRTUALNODES;
        break;
    case 1:
        type_ = ANYTYPE;
        break;
    case 2:
        type_ = BOUNDARYNODES;
        break;
    default:
        std::string errMsg;
        errMsg = "No pointType of ";
        errMsg += std::to_string(pT);
        errorAbortStr(errMsg);
        break;
    }
}

void OpenHurricane::pointZone::setND(const short ND) {
    switch (ND) {
    case 2:
        ND_ = Options::geometryModel::TWOD;
        break;
    case 3:
        ND_ = Options::geometryModel::THREED;
        break;
    default:
        std::string errMsg;
        errMsg = "No ND of ";
        errMsg += std::to_string(ND);
        errorAbortStr(errMsg);
        break;
    }
}

OpenHurricane::pointZone &OpenHurricane::pointZone::operator=(const pointZone &pz) {
    if (this != std::addressof(pz)) {
        zone::operator=(pz);
        type_ = pz.type_;
        ND_ = pz.ND_;
    }
    return *this;
}

namespace OpenHurricane {
    const std::string feature<pointZone>::className_ = "pointZone";

    void feature<pointZone>::pack(const pointZone &pZ, short *types) {
        *types = pZ.type();
        *(types + 1) = pZ.ND();
    }

    void feature<pointZone>::objConstruct(List<pointZone> &pointZones, const integer *index,
                                          const integer *cSize, const short *types,
                                          const char *names, const integer &nzone) {
        integer start = 0;
        for (integer i = 0; i < nzone; i++) {
            char *c = new char[*(cSize + i) + 1]();
            strncpy(c, names + start, *(cSize + i));
            c[*(cSize + i)] = '\0';
            string name(c);
            pointZones[i].setIndex(*(index + i));
            pointZones[i].setPointType(*(types + 2 * i));
            pointZones[i].setND(*(types + 2 * i + 1));
            pointZones[i].resetName(name);
            start += *(cSize + i);
            HurDeleteDynArray(c);
        }
    }
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::cellZoneTypes::cellZoneType
OpenHurricane::cellZoneTypes::getCellZoneType(const short czt) noexcept {
    switch (czt) {
    case cellZoneType::FLUID:
        return cellZoneType::FLUID;
        break;
    case cellZoneType::SOLID:
        return cellZoneType::SOLID;
        break;
    default:
        return cellZoneType::OTHER;
        break;
    }
}

OpenHurricane::cellZone::cellZone(const string &zName, const integer id, const integer firstID,
                              const integer lastID, const cellZoneTypes::cellZoneType czT,
                              const cellShapeType::shapeTypes cT)
    : zone(zName, id, firstID, lastID), type_(czT), shapeType_(cT) {}

OpenHurricane::cellZone::cellZone(const integer id, const integer firstID, const integer lastID,
                              const short czT, const short cT)
    : zone(" ", id, firstID, lastID) {
    type_ = cellZoneTypes::getCellZoneType(czT);
    shapeType_ = cellShapeType::getshapeType(cT);
}

OpenHurricane::cellZone::cellZone(const cellZone &cz)
    : zone(cz), type_(cz.type_), shapeType_(cz.shapeType_) {}

OpenHurricane::cellZone &OpenHurricane::cellZone::operator=(const cellZone &cz) {
    if (this != std::addressof(cz)) {
        zone::operator=(cz);
        type_ = cz.type_;
        shapeType_ = cz.shapeType_;
    }
    return *this;
}

namespace OpenHurricane {
    const std::string feature<cellZone>::className_ = "cellZone";

    void feature<cellZone>::objConstruct(List<cellZone> &cellZones, const integer *index,
                                         const integer *cSize, const short *types,
                                         const char *names, const integer &nzone) {
        integer start = 0;
        for (integer i = 0; i < nzone; i++) {
            char *c = new char[*(cSize + i) + 1]();
            strncpy(c, names + start, *(cSize + i));
            c[*(cSize + i)] = '\0';
            string name(c);
            cellZones[i].setIndex(*(index + i));
            cellZones[i].setType(*(types + 2 * i));
            cellZones[i].setCellType(*(types + 2 * i + 1));
            cellZones[i].resetName(name);
            start += *(cSize + i);
            HurDeleteDynArray(c);
        }
    }

    void feature<cellZone>::pack(const cellZone &cZ, short *types) {
        *types = cZ.type();
        *(types + 1) = cZ.shapeType();
    }

} // namespace OpenHurricane

OpenHurricane::faceZone::faceZone(const string &zName, const integer id, const integer firstID,
                              const integer lastID, const faceBCType::bcTypes bcT,
                              const faceShapeType::shapeTypes fT)
    : zone(zName, id, firstID, lastID), bcType_(bcT), faceType_(fT) {}

OpenHurricane::faceZone::faceZone(const integer id, const integer firstID, const integer lastID,
                              const short bcT, const short fT)
    : zone(" ", id, firstID, lastID) {
    bcType_ = faceBCType::getBcTypes(bcT);
    faceType_ = faceShapeType::getShapeTypes(fT);
}

OpenHurricane::faceZone::faceZone(const faceZone &fZ)
    : zone(fZ), bcType_(fZ.bcType_), faceType_(fZ.faceType_) {}

OpenHurricane::faceZone &OpenHurricane::faceZone::operator=(const faceZone &z) {
    if (this != std::addressof(z)) {
        zone::operator=(z);
        bcType_ = z.bcType_;
        faceType_ = z.faceType_;
    }
    return *this;
}

namespace OpenHurricane {
    const std::string feature<faceZone>::className_ = "faceZone";

    void feature<faceZone>::pack(const faceZone &fZ, short *types) {
        *types = fZ.bcType();
        *(types + 1) = fZ.faceType();
    }

    void feature<faceZone>::objConstruct(List<faceZone> &faceZones, const integer *index,
                                         const integer *cSize, const short *types,
                                         const char *names, const integer &nzone) {
        integer start = 0;
        for (integer i = 0; i < nzone; i++) {
            char *c = new char[*(cSize + i) + 1]();
            strncpy(c, names + start, *(cSize + i));
            c[*(cSize + i)] = '\0';
            string name(c);
            faceZones[i].setIndex(*(index + i));
            faceZones[i].setBcType(*(types + 2 * i));
            faceZones[i].setFaceType(*(types + 2 * i + 1));
            faceZones[i].resetName(name);
            start += *(cSize + i);
            HurDeleteDynArray(c);
        }
    }
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::integer OpenHurricane::cutZone::cutSize() const noexcept {
#ifdef HUR_DEBUG
    if (sor_.size() != des_.size()) {
        LFatal("The number of cut bc is not equal in zone: %d sor.size() = %d   des.size() = %d",
               id_, sor_.size(), des_.size());
    }
#endif // HUR_DEBUG
    return sor_.size();
}

template <> void OpenHurricane::cutZone::transfer(Array<real> &f) const {
#ifdef HUR_DEBUG
    if (f.size() == 0) {
        LFatal("The size of array f is zero");
    }
#endif // HUR_DEBUG

    if (HurMPI::parRun()) {
        if (HurMPI::getProcRank() == sendProc_) {
            if (realSendBufPtr_ == nullptr) {
                realSendBufPtr_ = new real[sor_.size()];
            }
            for (integer i = 0; i < sor_.size(); i++) {
                realSendBufPtr_[i] = f[sor_[i]];
            }
            HurMPI::send(realSendBufPtr_, sor_.size(), feature<real>::MPIType, receivProc_,
                         sendProc_, HurMPI::getComm());

        } else if (HurMPI::getProcRank() == receivProc_) {
            if (realRecvBufPtr_ == nullptr) {
                realRecvBufPtr_ = new real[des_.size()];
            }
            HurMPI::Status status;
            HurMPI::recv(realRecvBufPtr_, des_.size(), feature<real>::MPIType, sendProc_, sendProc_,
                         HurMPI::getComm(), &status);
            for (integer i = 0; i < sor_.size(); i++) {
                f[des_[i]] = realRecvBufPtr_[i];
            }
        }
        //HurMPI::barrier(HurMPI::getComm());
    }
}

template <> void OpenHurricane::cutZone::transfer(Array<realArray> &f) const {
#ifdef HUR_DEBUG
    if (f.size() == 0) {
        LFatal("The size of array f is zero");
    }
#endif // HUR_DEBUG

    if (HurMPI::parRun()) {
        if (HurMPI::getProcRank() == sendProc_) {
            real *realSendBufPtr = new real[sor_.size() * f[0].size()];

            for (integer i = 0; i < sor_.size(); i++) {
                for (integer j = 0; j < f[0].size(); ++j) {
                    realSendBufPtr[i * f[0].size() + j] = f[sor_[i]][j];
                }
            }
            HurMPI::send(realSendBufPtr, sor_.size() * f[0].size(), feature<real>::MPIType,
                         receivProc_, sendProc_, HurMPI::getComm());
            delete[] realSendBufPtr;
        } else if (HurMPI::getProcRank() == receivProc_) {
            real *realRecvBufPtr = new real[des_.size() * f[0].size()];

            HurMPI::Status status;
            HurMPI::recv(realRecvBufPtr, des_.size() * f[0].size(), feature<real>::MPIType,
                         sendProc_, sendProc_, HurMPI::getComm(), &status);
            for (integer i = 0; i < sor_.size(); i++) {
                for (integer j = 0; j < f[0].size(); ++j) {
                    f[des_[i]][j] = realRecvBufPtr[i * f[0].size() + j];
                }
            }
            delete[] realRecvBufPtr;
        }
    }
}

OpenHurricane::integer OpenHurricane::cutZone::sendNonBlocking(Array<real> &f,
                                                       HurMPI::Request *request) const {
    if (HurMPI::isThisProc(sendProc_)) {
        if (realSendBufPtr_ == nullptr) {
            realSendBufPtr_ = new real[sor_.size()];
        }
        for (integer i = 0; i < sor_.size(); i++) {
            realSendBufPtr_[i] = f[sor_[i]];
        }

        auto ierr = HurMPI::isend(realSendBufPtr_, sor_.size(), feature<real>::MPIType, receivProc_,
                                  sendProc_, HurMPI::getComm(), request);

        if (ierr != MPI_SUCCESS) {
            std::string estr;
            HurMPI::errorString(ierr, estr);
            errorAbortStr(("Send data error: " + estr));
        }
    }
    return sendProc_;
}

void OpenHurricane::cutZone::recvNonBlocking(Array<real> &f) const {
    if (HurMPI::isThisProc(receivProc_)) {
        if (realRecvBufPtr_ == nullptr) {
            realRecvBufPtr_ = new real[des_.size()];
        }
        HurMPI::Request request;
        auto ierr = HurMPI::irecv(realRecvBufPtr_, des_.size(), feature<real>::MPIType, sendProc_,
                                  sendProc_, HurMPI::getComm(), &request);
        if (ierr != MPI_SUCCESS) {
            std::string estr;
            HurMPI::errorString(ierr, estr);
            errorAbortStr(("Receive data error: " + estr));
        }
        HurMPI::Status status;

        // 这里这样用有问题
        HurMPI::wait(&request, &status);
        for (integer i = 0; i < sor_.size(); i++) {
            f[des_[i]] = realRecvBufPtr_[i];
        }
    }
}

void OpenHurricane::periodicPair::gatherPeriodicPair(const periodicPair &pp, periodicPair &gpp,
                                                 const int root) {
    if (!HurMPI::parRun()) {
        gpp = pp;
        return;
    }
    integerList nSizeL(HurMPI::getProcSize(), Zero);
    nSizeL[HurMPI::getProcRank()] = pp.size();

    integerList displs;
    HurMPI::gatherList(nSizeL, root, HurMPI::getComm());
    if (HurMPI::isThisProc(root)) {
        displs.resize(HurMPI::getProcSize(), Zero);
    }
    integer ssi = 0;

    if (HurMPI::isThisProc(root)) {
        for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
            displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
        }
        for (integer i = 0; i < nSizeL.size(); ++i) {
            ssi += nSizeL[i];
        }
    }

    integerList pp1(pp.size(), Zero);
    integerList rootPp;

    if (HurMPI::isThisProc(root)) {
        rootPp.resize(ssi);
        gpp.resize(ssi);

        gpp.periodicZone_ = pp.periodicZone_;
        gpp.shadowZone_ = pp.shadowZone_;
        gpp.type_ = pp.type_;
        gpp.firstIndex_ = 1;
        gpp.lastIndex_ = ssi;
    }
    HurMPI::barrier(HurMPI::getComm());
    for (int i = 0; i < feature<vector2D>::nElements_; ++i) {
        for (integer j = 0; j < pp.size(); ++j) {
            pp1[j] = pp[j][i];
        }

        HurMPI::Request request;
        HurMPI::igatherv(pp1.data(), pp1.size(), feature<integer>::MPIType, rootPp.data(),
                         nSizeL.data(), displs.data(), feature<integer>::MPIType, root,
                         HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::isThisProc(root)) {
            for (integer j = 0; j < rootPp.size(); ++j) {
                gpp[j][i] = rootPp[j];
            }
        }
    }
}

OpenHurricane::periodicPair &OpenHurricane::periodicPair::operator=(const periodicPair &pp) {
    if (this != std::addressof(pp)) {
        integerVector2DList::resize(pp.size());
        integerVector2DList::operator=(pp);
        periodicZone_ = pp.periodicZone_;
        shadowZone_ = pp.shadowZone_;
        firstIndex_ = pp.firstIndex_;
        lastIndex_ = pp.lastIndex_;
        type_ = pp.type_;
    }
    return *this;
}

OpenHurricane::perZone::perZone(const integer idx, const integer sP, const integer dP,
                            const integerList &sorBc, const integerList &desBc, const short t,
                            const tensor &R)
    : id_(idx), pairId_(0), isShadow_(false), sendProc_(sP), receivProc_(dP), sor_(sorBc),
      des_(desBc), sorFace_(), desFace_(), rotationAxis_(X_ROTATION), RPtr_(nullptr),
      tranDistancePtr_(nullptr), phi_(0.0), isTwoLayer_(false) {
    std::map<short, periodicTypes::periodicType>::const_iterator iter;
    iter = periodicTypes::periodicTypeMap.find(t);
    if (iter == periodicTypes::periodicTypeMap.end()) {
        std::string errMsg = "periodic type not found! type_ = ";
        errMsg = errMsg + std::to_string(t);
        errorAbortStr(errMsg);
    }
    type_ = t;
    if (type_ == periodicTypes::TRANSLATIONAL) {
        RPtr_ = nullptr;
    } else {
        RPtr_ = new tensor(Zero);
        *RPtr_ = R;
    }
    if (sor_.size() != des_.size()) {
        LFatal("The number of cut bc is not equal in periodic zone: %d sor.size() = %d   "
               "des.size() = %d",
               id_, sor_.size(), des_.size());
    }
}

void OpenHurricane::perZone::transferTensorTran(List<tensor> &f) const {
    typename feature<tensor>::value_type *send;
    typename feature<tensor>::value_type *recv;
    if (sendProc_ == receivProc_) //Periodic face pair on the same processor
    {
        if (HurMPI::getProcRank() == sendProc_) {
            for (integer i = 0; i < sor_.size(); i++) {
                for (int j = 0; j < feature<tensor>::nElements_; j++) {
                    f[des_[i]][j] = f[sor_[i]][j];
                }
                if (RPtr_ != nullptr) {
                    f[des_[i]] = (*RPtr_) * f[des_[i]] * (RPtr_)->transpose();
                }
            }
        }
    } else // Periodic face pair on different processor
    {
        if (HurMPI::getProcRank() == sendProc_) {
            send =
                new typename feature<tensor>::value_type[sor_.size() * feature<tensor>::nElements_];
            for (integer i = 0; i < sor_.size(); i++) {
                for (int j = 0; j < feature<tensor>::nElements_; j++) {
                    send[i * feature<tensor>::nElements_ + j] = f[sor_[i]][j];
                }
            }
            HurMPI::send(send, sor_.size() * feature<tensor>::nElements_, feature<tensor>::MPIType,
                         receivProc_, sendProc_, HurMPI::getComm());

            delete[] send;
        } else if (HurMPI::getProcRank() == receivProc_) {
            recv =
                new typename feature<tensor>::value_type[des_.size() * feature<tensor>::nElements_];
            HurMPI::Status status;
            HurMPI::recv(recv, des_.size() * feature<tensor>::nElements_, feature<tensor>::MPIType,
                         sendProc_, sendProc_, HurMPI::getComm(), &status);
            for (integer i = 0; i < sor_.size(); i++) {
                for (int j = 0; j < feature<tensor>::nElements_; j++) {
                    f[des_[i]][j] = recv[i * feature<tensor>::nElements_ + j];
                }
                if (RPtr_ != nullptr) {
                    f[des_[i]] = (*RPtr_) * f[des_[i]] * (RPtr_)->transpose();
                }
            }
            delete[] recv;
        }
    }
    HurMPI::barrier(HurMPI::getComm());
}

void OpenHurricane::perZone::transferPoint(List<vector> &pointField, const vector &origin) const {
    if (isTranslational()) {
        if (tranDistancePtr_ == nullptr) {
            LFatal("The translational distance pointer for periodic "
                   "boundary is not specific!");
        }
        if (sendProc_ == receivProc_) //Periodic face pair on the same processor
        {
            if (HurMPI::getProcRank() == sendProc_) {
                for (integer i = 0; i < sor_.size(); i++) {
                    pointField[des_[i]] = pointField[sor_[i]] + tranDistance();
                }
            }
        } else {
            real *send;
            real *recv;
            if (HurMPI::getProcRank() == sendProc_) {
                send = new real[sor_.size() * feature<vector>::nElements_];
                for (integer i = 0; i < sor_.size(); i++) {
                    for (int j = 0; j < feature<vector>::nElements_; j++) {
                        send[i * feature<vector>::nElements_ + j] = pointField[sor_[i]][j];
                    }
                }
                HurMPI::send(send, sor_.size() * feature<vector>::nElements_,
                             feature<vector>::MPIType, receivProc_, sendProc_, HurMPI::getComm());
                delete[] send;
            } else if (HurMPI::getProcRank() == receivProc_) {
                recv = new real[des_.size() * feature<vector>::nElements_];
                HurMPI::Status status;
                HurMPI::recv(recv, des_.size() * feature<vector>::nElements_,
                             feature<vector>::MPIType, sendProc_, sendProc_, HurMPI::getComm(),
                             &status);
                for (integer i = 0; i < sor_.size(); i++) {
                    for (int j = 0; j < feature<vector>::nElements_; j++) {
                        pointField[des_[i]][j] = recv[i * feature<vector>::nElements_ + j];
                    }
                    pointField[des_[i]] += tranDistance();
                }
                delete[] recv;
            }
        }
    } else {
        if (RPtr_ == nullptr) {
            LFatal("The rotational tensor distance pointer for periodic "
                   "boundary is not specific!");
        }
        realSquareMatrix RR(4, Zero);
        RR(0, 0) = 1.0;
        RR(1, 1) = 1.0;
        RR(2, 2) = 1.0;
        RR(3, 3) = 1.0;
        RR(3, 0) = -origin[0];
        RR(3, 1) = -origin[1];
        RR(3, 2) = -origin[2];

        realSquareMatrix Rot(4, Zero);
        Rot(0, 0) = (*RPtr_).xx();
        Rot(0, 1) = (*RPtr_).xy();
        Rot(0, 2) = (*RPtr_).xz();
        Rot(1, 0) = (*RPtr_).yx();
        Rot(1, 1) = (*RPtr_).yy();
        Rot(1, 2) = (*RPtr_).yz();
        Rot(2, 0) = (*RPtr_).zx();
        Rot(2, 1) = (*RPtr_).zy();
        Rot(2, 2) = (*RPtr_).zz();
        Rot(3, 3) = 1.0;

        realSquareMatrix RRm(4, Zero);
        RRm(0, 0) = 1.0;
        RRm(1, 1) = 1.0;
        RRm(2, 2) = 1.0;
        RRm(3, 3) = 1.0;
        RRm(3, 0) = origin[0];
        RRm(3, 1) = origin[1];
        RRm(3, 2) = origin[2];

        const auto RotL = RR * Rot * RRm;
        realArray p(4, Zero);
        p[3] = 1.0;
        if (sendProc_ == receivProc_) //Periodic face pair on the same processor
        {
            if (HurMPI::getProcRank() == sendProc_) {
                for (integer i = 0; i < sor_.size(); i++) {
                    p[0] = pointField[sor_[i]][0];
                    p[1] = pointField[sor_[i]][1];
                    p[2] = pointField[sor_[i]][2];
                    p = RotL * p;
                    pointField[des_[i]][0] = p[0];
                    pointField[des_[i]][1] = p[1];
                    pointField[des_[i]][2] = p[2];
                }
            }
        } else {
            real *send;
            real *recv;
            if (HurMPI::getProcRank() == sendProc_) {
                send = new real[sor_.size() * feature<vector>::nElements_];
                for (integer i = 0; i < sor_.size(); i++) {
                    for (int j = 0; j < feature<vector>::nElements_; j++) {
                        send[i * feature<vector>::nElements_ + j] = pointField[sor_[i]][j];
                    }
                }
                HurMPI::send(send, sor_.size() * feature<vector>::nElements_,
                             feature<vector>::MPIType, receivProc_, sendProc_, HurMPI::getComm());
                delete[] send;
            } else if (HurMPI::getProcRank() == receivProc_) {
                recv = new real[des_.size() * feature<vector>::nElements_];
                HurMPI::Status status;
                HurMPI::recv(recv, des_.size() * feature<vector>::nElements_,
                             feature<vector>::MPIType, sendProc_, sendProc_, HurMPI::getComm(),
                             &status);
                for (integer i = 0; i < sor_.size(); i++) {
                    for (int j = 0; j < feature<vector>::nElements_; j++) {
                        pointField[des_[i]][j] = recv[i * feature<vector>::nElements_ + j];
                    }
                    if (RPtr_ != nullptr) {
                        //pointField[des_[i]] = (*RPtr_) * pointField[des_[i]];
                        p[0] = pointField[des_[i]][0];
                        p[1] = pointField[des_[i]][1];
                        p[2] = pointField[des_[i]][2];
                        p = RotL * p;
                        pointField[des_[i]][0] = p[0];
                        pointField[des_[i]][1] = p[1];
                        pointField[des_[i]][2] = p[2];
                    }
                }
                delete[] recv;
            }
        }
    }
}

template <> void OpenHurricane::perZone::transfer(Array<realArray> &f) const {
    if (sendProc_ == receivProc_) //Periodic face pair on the same processor
    {
        if (HurMPI::getProcRank() == sendProc_) {
            for (integer i = 0; i < sor_.size(); i++) {
                f[des_[i]] = f[sor_[i]];
            }
        }
    } else // Periodic face pair on different processor
    {
        if (HurMPI::getProcRank() == sendProc_) {
            real *sendBufPtr = new real[sor_.size() * f[0].size()];
            for (integer i = 0; i < sor_.size(); i++) {
                for (integer j = 0; j < f[0].size(); ++j) {
                    sendBufPtr[i * f[0].size() + j] = f[sor_[i]][j];
                }
            }
            HurMPI::send(sendBufPtr, sor_.size() * f[0].size(), feature<real>::MPIType, receivProc_,
                         sendProc_, HurMPI::getComm());
            delete[] sendBufPtr;
        } else if (HurMPI::getProcRank() == receivProc_) {
            real *recvBufPtr = new real[des_.size() * f[0].size()];
            HurMPI::Status status;
            HurMPI::recv(recvBufPtr, des_.size() * f[0].size(), feature<real>::MPIType, sendProc_,
                         sendProc_, HurMPI::getComm(), &status);
            for (integer i = 0; i < sor_.size(); i++) {
                for (integer j = 0; j < f[0].size(); ++j) {
                    f[des_[i]][j] = recvBufPtr[i * f[0].size() + j];
                }
            }
            delete[] recvBufPtr;
        }
    }
}
OpenHurricane::perZone &OpenHurricane::perZone::operator=(const perZone &pz) {
    if (this != std::addressof(pz)) {
        id_ = pz.id_;
        pairId_ = pz.pairId_;
        isShadow_ = pz.isShadow_;
        sendProc_ = pz.sendProc_;
        receivProc_ = pz.receivProc_;

        sorFace_.resize(pz.sorFace_.size());
        sorFace_ = pz.sorFace_;
        desFace_.resize(pz.desFace_.size());
        desFace_ = pz.desFace_;
        rotationAxis_ = pz.rotationAxis_;
        phi_ = pz.phi_;
        isTwoLayer_ = pz.isTwoLayer_;

        sendProc_ = pz.sendProc_;
        receivProc_ = pz.receivProc_;
        sor_.resize(pz.sor_.size());
        sor_ = pz.sor_;
        des_.resize(pz.des_.size());
        des_ = pz.des_;
        type_ = pz.type_;
        if (pz.RPtr_ == nullptr) {
            HurDelete(RPtr_);

        } else if (pz.RPtr_ != nullptr) {
            if (RPtr_ != nullptr) {
                *RPtr_ = *pz.RPtr_;
            } else {
                RPtr_ = new tensor;
                *RPtr_ = *pz.RPtr_;
            }
        }

        if (pz.tranDistancePtr_ == nullptr) {
            HurDelete(tranDistancePtr_);

        } else if (pz.tranDistancePtr_ != nullptr) {
            if (tranDistancePtr_ != nullptr) {
                *tranDistancePtr_ = *pz.tranDistancePtr_;
            } else {
                tranDistancePtr_ = new vector;
                *tranDistancePtr_ = *pz.tranDistancePtr_;
            }
        }
    }
    return *this;
}