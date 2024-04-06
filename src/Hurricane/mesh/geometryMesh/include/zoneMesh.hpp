/*!
 * \file zoneMesh.hpp
 * \brief Header of all kind of zone
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
#include "dataStructure.hpp"
#include "geometryModel.hpp"
#include "meshElements.hpp"
#include "integerArray.hpp"

namespace OpenHurricane {
    class pointZone;
    class cellZone;
    class faceZone;
    class cutZone;
    class periodicPair;
    class perZone;
    /**
     * \brief The point zone list.
     */
    using pointZoneList = List<pointZone>;
    /**
     * \brief The cell zone list.
     */
    using cellZoneList = List<cellZone>;

    /**
     * \brief The face zone list.
     */
    using faceZoneList = List<faceZone>;
    /**
     * \brief The cut zone list.
     */
    using cutZoneList = List<cutZone>;
    
    /**
     * \brief The periodic zone list.
     */
    using perZoneList = List<perZone>;

    /**
     * \brief The periodic pair zone list.
     */
    using periodicPairList = List<periodicPair>;
    namespace Options {
        enum geometryModel : short {
            TWOD = 2,  // Two dimensional
            THREED = 3 // Three dimensional
        };

        static std::map<std::string, geometryModel> geometryModelMap =
            createMap<std::string, geometryModel>("TWOD", TWOD)("THREED", THREED);
    } // namespace Options

    /**
     * \brief The base class of the zone.
     */
    class zone {
    protected:
        /**\brief Name of zone.*/
        string name_;

        /**\brief Index of zone.*/
        integer index_;

        /**\brief First index.*/
        integer firstIndex_;

        /**\brief Last index.*/
        integer lastIndex_;

    public:
        /**\brief Construct from null.*/
        inline zone() : name_(), index_(-1), firstIndex_(-1), lastIndex_(-1) {}

        /**\brief Construct from components.*/
        inline zone(const string &zName, const integer id, const integer firstID,
                    const integer lastID)
            : name_(zName), index_(id), firstIndex_(firstID), lastIndex_(lastID) {}

        /**\brief Construct as copy.*/
        inline zone(const zone &myZone)
            : name_(myZone.name_), index_(myZone.index_), firstIndex_(myZone.firstIndex_),
              lastIndex_(myZone.lastIndex_) {}

        /**\brief Destructor.*/
        inline virtual ~zone() noexcept {}

        /**\brief Return name.*/
        hur_nodiscard inline const string &name() const noexcept { return name_; }

        /**\brief Return the index of this zone.*/
        hur_nodiscard inline integer index() const noexcept { return index_; }

        /**\brief Return the first index.*/
        hur_nodiscard inline integer firstIndex() const noexcept { return firstIndex_; }

        /**\brief Return the last index.*/
        hur_nodiscard inline integer lastIndex() const noexcept { return lastIndex_; }

        /**\brief Return the size of this zone.*/
        hur_nodiscard inline integer size() const noexcept {
            return (lastIndex_ - firstIndex_ + 1);
        }

        /**\brief Reset name.*/
        inline void resetName(const string &n) { name_ = n; }

        /**\brief Reset index.*/
        inline void setIndex(const integer id) noexcept { index_ = id; }

        /**\brief Reset firstIndex.*/
        inline void setFirstIndex(const integer id) noexcept { firstIndex_ = id; }

        /**\brief Reset lastIndex.*/
        inline void setLastIndex(const integer id) noexcept { lastIndex_ = id; }

        zone &operator=(const zone &z);
    };

    /**
     * \brief The class of the point zone.
     */
    class pointZone : public zone {
    public:
        enum pointType : short {
            VIRTUALNODES = 0, //!< Zero for "virtual" nodes.
            ANYTYPE = 1,      //!< One for no (any) bcType.
            BOUNDARYNODES = 2 //!< Two for boundary nodes.
        };

    private:
        /**\brief The node bcType of this point zone.*/
        pointType type_;

        /**\brief The dimensionality of the node data.*/
        Options::geometryModel ND_;

    public:
        /**\brief Construct from null.*/
        inline pointZone() : zone(), type_(ANYTYPE), ND_(Options::THREED) {}

        /**\brief Construct from components.*/
        pointZone(const string &zName, const integer id, const integer firstID,
                  const integer lastID, const Options::geometryModel nD);

        /**\brief Construct from components.*/
        pointZone(const integer id, const integer firstID, const integer lastID, const short nD,
                  const short t);

        /**\brief Construct as copy.*/
        pointZone(const pointZone &);

        /**\brief Destructor.*/
        inline virtual ~pointZone() noexcept {}

        /**\brief Return the dimensionality of the node data.*/
        hur_nodiscard inline Options::geometryModel ND() const noexcept { return ND_; }

        /**\brief Return the node bcType of this point zone.*/
        hur_nodiscard inline pointType type() const noexcept { return type_; }

        void setPointType(const short pT);

        void setND(const short ND);

        pointZone &operator=(const pointZone &pz);
    };

    /*!\brief Template specialization for feature<pointZone>.*/
    template <> class feature<pointZone> {
        pointZone p_;

    public:
        static void pack(const pointZone &, short *);

        static void objConstruct(List<pointZone> &, const integer *, const integer *, const short *,
                                 const char *, const integer &);

        declareClassNames;

        /*!\brief Construct from primitive.*/
        explicit feature(const pointZone &p) : p_(p) {}
    };

    namespace cellZoneTypes {
        enum cellZoneType : short {
            FLUID = 1,  // fluid
            SOLID = 17, // solid
            OTHER = 18  // other
        };

        hur_nodiscard cellZoneType getCellZoneType(const short czt) noexcept;

        static const std::map<std::string, cellZoneType> cellZoneTypeStrMap =
            createMap<std::string, cellZoneType>("FULID", FLUID)("SOLID", SOLID)("OTHER", OTHER);

    } // namespace cellZoneTypes

    /**
     * \brief The class of the cell zone.
     */
    class cellZone : public zone {
    private:
        /**\brief The bcType of this cell zone (fluid or solid).*/
        cellZoneTypes::cellZoneType type_;

        /**\brief The geometry bcType of the cell in the zone.*/
        cellShapeType::shapeTypes shapeType_;

    public:
        /**\brief Construct from null.*/
        inline cellZone()
            : zone(), type_(cellZoneTypes::FLUID),
              shapeType_(cellShapeType::shapeTypes::hexahedral) {}

        /**\brief Construct from components.*/
        cellZone(const string &zName, const integer id, const integer firstID, const integer lastID,
                 const cellZoneTypes::cellZoneType czT, const cellShapeType::shapeTypes cT);

        /**\brief Construct from components.*/
        cellZone(const integer id, const integer firstID, const integer lastID, const short czT,
                 const short cT);

        /**\brief Construct as copy.*/
        cellZone(const cellZone &);

        /**\brief Destructor.*/
        inline virtual ~cellZone() noexcept {}

        /**\brief Return the type of this cell zone, e.g., fluid, solid or other.*/
        hur_nodiscard inline cellZoneTypes::cellZoneType type() const noexcept { return type_; }

        hur_nodiscard inline cellShapeType::shapeTypes shapeType() const noexcept {
            return shapeType_;
        }

        hur_nodiscard inline bool isMixed() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::mixed);
        }

        hur_nodiscard inline bool isTriangular() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::triangular);
        }

        hur_nodiscard inline bool isTetrahedral() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::tetrahedral);
        }

        hur_nodiscard inline bool isQuadrilateral() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::quadrilateral);
        }

        hur_nodiscard inline bool isHexahedral() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::hexahedral);
        }

        hur_nodiscard inline bool isPyramid() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::pyramid);
        }

        hur_nodiscard inline bool isWedge() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::wedge);
        }

        hur_nodiscard inline bool isPolyhedral() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::polyhedral);
        }

        /*!\brief Reset cellZoneType.*/
        inline void setType(const short czT) noexcept {
            type_ = cellZoneTypes::getCellZoneType(czT);
        }

        /*!\brief Reset cellType.*/
        inline void setCellType(const short cT) noexcept {
            shapeType_ = cellShapeType::getshapeType(cT);
        }

        cellZone &operator=(const cellZone &cz);
    };

    /*!\brief Template specialization for feature<cellZone>.*/
    template <> class feature<cellZone> {
        cellZone p_;

    public:
        static void pack(const cellZone &, short *);

        static void objConstruct(List<cellZone> &, const integer *, const integer *, const short *,
                                 const char *, const integer &);

        declareClassNames;

        /*!\brief Construct from primitive.*/
        explicit feature(const cellZone &p) : p_(p) {}
    };

    /**
     * \brief The class of the face zone.
     */
    class faceZone : public zone {
    private:
        /**\brief The bc-bcType of faces in this face zone.*/
        faceBCType::bcTypes bcType_;

        /**\brief The shape of the faces in this face zone.*/
        faceShapeType::shapeTypes faceType_;

    public:
        /**\brief Construct from null.*/
        faceZone()
            : zone(), bcType_(faceBCType::bcTypes::INTERIOR),
              faceType_(faceShapeType::shapeTypes::quadrilateral) {}

        /**\brief Construct from components.*/
        faceZone(const string &zName, const integer id, const integer firstID, const integer lastID,
                 const faceBCType::bcTypes bcT, const faceShapeType::shapeTypes fT);

        /**\brief Construct from components.*/
        faceZone(const integer id, const integer firstID, const integer lastID, const short bcT,
                 const short fT);

        /**\brief Construct as copy.*/
        faceZone(const faceZone &);

        /**\brief Destructor.*/
        inline virtual ~faceZone() noexcept {}

        /**\brief Return the bc-bcType of faces in this face zone.*/
        hur_nodiscard inline faceBCType::bcTypes bcType() const noexcept { return bcType_; }

        /**\brief Return the faceType of faces in this face zone.*/
        hur_nodiscard inline faceShapeType::shapeTypes faceType() const noexcept {
            return faceType_;
        }

        /**\brief Return true if this face zone is boundary face zone.*/
        hur_nodiscard inline bool isBnd() const noexcept {
            return (!((bcType_ == faceBCType::bcTypes::INTERIOR) ||
                      (bcType_ == faceBCType::bcTypes::INTERFACE)));
        }

        /**\brief Return true if this face zone is parellel cut face zone.*/
        hur_nodiscard inline bool isCutFace() const noexcept {
            return (bcType_ == faceBCType::bcTypes::CUTFACE);
        }

        /**\brief Return true if this face zone is interior face zone or interface face zone.*/
        hur_nodiscard inline bool isInterior() const noexcept {
            return ((bcType_ == faceBCType::bcTypes::INTERIOR) ||
                    (bcType_ == faceBCType::bcTypes::INTERFACE));
        }

        /**\brief Return true if this face zone is interface face zone.*/
        hur_nodiscard inline bool isInterface() const noexcept {
            return (bcType_ == faceBCType::bcTypes::INTERFACE);
        }

        /**\brief Return true if this face zone is wall face zone.*/
        hur_nodiscard inline bool isWall() const noexcept {
            return (bcType_ == faceBCType::bcTypes::WALL);
        }

        /**\brief Return true if this face zone is periodic face zone.*/
        hur_nodiscard inline bool isPeriodic() const noexcept {
            return (bcType_ == faceBCType::bcTypes::PERIODIC);
        }

        /**\brief Return true if this face zone is periodic-shadow face zone.*/
        hur_nodiscard inline bool isPeriodicShadow() const noexcept {
            return (bcType_ == faceBCType::bcTypes::PERIODICSHADOW);
        }

        /**\brief Return true if this face zone is symmetric face zone.*/
        hur_nodiscard inline bool isSymmetric() const noexcept {
            return (bcType_ == faceBCType::bcTypes::SYMMETRY);
        }

        hur_nodiscard inline bool isMixed() const noexcept {
            return (faceType_ == faceShapeType::shapeTypes::mixed);
        }

        hur_nodiscard inline bool isLinear() const noexcept {
            return (faceType_ == faceShapeType::shapeTypes::linear);
        }

        hur_nodiscard inline bool isTriangular() const noexcept {
            return (faceType_ == faceShapeType::shapeTypes::triangular);
        }

        hur_nodiscard inline bool isQuadrilateral() const noexcept {
            return (faceType_ == faceShapeType::shapeTypes::quadrilateral);
        }

        hur_nodiscard inline bool isPolygonal() const noexcept {
            return (faceType_ == faceShapeType::shapeTypes::polygonal);
        }

        /*!\brief Reset bcType.*/
        inline void setBcType(const short bcT) noexcept { bcType_ = faceBCType::getBcTypes(bcT); }

        inline void setBcType(const faceBCType::bcTypes bcT) noexcept { bcType_ = bcT; }

        /*!\brief Reset faceType.*/
        inline void setFaceType(const short fT) noexcept {
            faceType_ = faceShapeType::getShapeTypes(fT);
        }

        faceZone &operator=(const faceZone &z);
    };

    /*!\brief Template specialization for feature<faceZone>.*/
    template <> class feature<faceZone> {
        faceZone p_;

    public:
        static void pack(const faceZone &, short *);

        static void objConstruct(List<faceZone> &, const integer *, const integer *, const short *,
                                 const char *, const integer &);

        declareClassNames;

        /*!\brief Construct from primitive.*/
        explicit feature(const faceZone &p) : p_(p) {}
    };

    /**
     * \brief The class of the cut zone.
     */
    class cutZone {
    private:
        // Private data

        /**\brief Index of the cut zone in the face zones.*/
        integer id_;

        /**\brief The sending processor number for current cutzone.*/
        integer sendProc_;

        /**\brief The receiving processor number for current cutzone.*/
        integer receivProc_;

        /**\brief the local cutbc number for sending  processor.*/
        integerList sor_;

        /**\brief The local cutbc number for receiving processor.*/
        integerList des_;

        // whether form two layer of ghost cell or not
        bool isTwoLayer_;

        mutable real *realSendBufPtr_;
        mutable real *realRecvBufPtr_;

    public:
        /**\brief Construct from null.*/
        inline cutZone()
            : id_(0), sendProc_(0), receivProc_(0), sor_(), des_(), isTwoLayer_(true),
              realSendBufPtr_(nullptr), realRecvBufPtr_(nullptr) {}

        /**\brief Construct from components.*/
        inline cutZone(const integer idx, const integer sP, const integer dP)
            : id_(idx), sendProc_(sP), receivProc_(dP), sor_(), des_(), isTwoLayer_(true),
              realSendBufPtr_(nullptr), realRecvBufPtr_(nullptr) {}

        /**\brief Construct from components.*/
        inline cutZone(const integer idx, const integer sP, const integer dP,
                       const integerList &sorBc, const integerList &desBc)
            : id_(idx), sendProc_(sP), receivProc_(dP), sor_(sorBc), des_(desBc), isTwoLayer_(true),
              realSendBufPtr_(nullptr), realRecvBufPtr_(nullptr) {}

        /**\brief Construct as copy.*/
        inline cutZone(const cutZone &cz)
            : id_(cz.id_), sendProc_(cz.sendProc_), receivProc_(cz.receivProc_), sor_(cz.sor_),
              des_(cz.des_), isTwoLayer_(cz.ghostLayer()), realSendBufPtr_(nullptr),
              realRecvBufPtr_(nullptr) {}

        /**\brief Destructor.*/
        inline ~cutZone() noexcept {
            HurDeleteDynArray(realSendBufPtr_);
            HurDeleteDynArray(realRecvBufPtr_);
            /*if (realSendBufPtr_ != nullptr)
            {
                    delete[] realSendBufPtr_;
            }*/
            /*if (realRecvBufPtr_ != nullptr)
            {
                    delete[] realRecvBufPtr_;
            }*/
        }

        // Access

        /**\brief Return the index of the cut zone in the face zones.*/
        hur_nodiscard inline integer id() const noexcept { return id_; }

        hur_nodiscard inline integer sendProc() const noexcept { return sendProc_; }

        inline integer receivProc() const noexcept { return receivProc_; }

        hur_nodiscard inline bool isSendFromThisProc() const noexcept {
            return sendProc_ == HurMPI::getProcRank();
        }

        hur_nodiscard inline bool isThisProcReceiv() const noexcept {
            return receivProc_ == HurMPI::getProcRank();
        }

        hur_nodiscard inline integerList &sor() noexcept { return sor_; }

        hur_nodiscard inline const integerList &sor() const noexcept { return sor_; }

        hur_nodiscard inline integerList &des() noexcept { return des_; }

        inline const hur_nodiscard integerList &des() const noexcept { return des_; }

        hur_nodiscard integer cutSize() const noexcept;

        hur_nodiscard inline bool ghostLayer() const noexcept { return isTwoLayer_; }

        inline void clearList() noexcept {
            sor_.clear();
            des_.clear();
        }

        // operator

        inline cutZone &operator=(const cutZone &cz) {
            if (this == &cz) {
                return *this;
            } else {
                sor_.resize(cz.sor_.size());
                des_.resize(cz.des_.size());
                sendProc_ = cz.sendProc_;
                receivProc_ = cz.receivProc_;
                sor_ = cz.sor_;
                des_ = cz.des_;
                isTwoLayer_ = cz.isTwoLayer_;
                return *this;
            }
        }

        /*!\brief Parallel transfer.*/

        /**\brief Transfer basic datatype.*/
        template <template <typename Type> class Form, typename Type>
        void transfer(Form<Type> &) const;

        integer sendNonBlocking(Array<real> &, HurMPI::Request *request) const;

        void recvNonBlocking(Array<real> &) const;

        /**\brief Transfer VertorSpace datatype.*/
        template <template <typename Type> class Form, typename Type>
        void transferVS(Form<Type> &) const;

        /**\brief Transfer complex datatype.*/
        template <template <typename Type> class Form, typename Type>
        void transferComp(Form<Type> &) const;

        /*!\brief Edit.*/

        /*!\brief set sendProc_.*/
        inline void setSendProc(const integer ip) noexcept { sendProc_ = ip; }

        /*!\brief set receivProc_.*/
        inline void setRecvProc(const integer ip) noexcept { receivProc_ = ip; }

        /*!\brief set isTwoLayer.*/
        inline void setGhostLayer(const bool layer) noexcept { isTwoLayer_ = layer; }
    };

    /*!\brief Transfer basic datatype.*/
    template <template <typename T> class Form, typename Type>
    void OpenHurricane::cutZone::transfer(Form<Type> &f) const {
        if (HurMPI::parRun()) {
            Type *send;
            Type *recv;
            if (HurMPI::getProcRank() == sendProc_) {

                send = new Type[sor_.size()];
                for (integer i = 0; i < sor_.size(); i++) {
                    send[i] = f[sor_[i]];
                }
                HurMPI::send(send, sor_.size(), feature<Type>::MPIType, receivProc_, sendProc_,
                             HurMPI::getComm());

                delete[] send;
            } else if (HurMPI::getProcRank() == receivProc_) {
                recv = new Type[des_.size()];
                HurMPI::Status status;
                HurMPI::recv(recv, des_.size(), feature<Type>::MPIType, sendProc_, sendProc_,
                             HurMPI::getComm(), &status);
                for (integer i = 0; i < sor_.size(); i++) {
                    f[des_[i]] = recv[i];
                }
                delete[] recv;
            }
            //HurMPI::barrier(HurMPI::getComm());
        }
    }

    /** \brief Transfer real array.*/
    template <> void OpenHurricane::cutZone::transfer(Array<real> &f) const;

    /**
     * \brief Transfer real array.
     *  The size of each f[i] must be the same.
     */
    template <> void OpenHurricane::cutZone::transfer(Array<realArray> &f) const;

    /*!\brief Transfer VertorSpace datatype.*/
    template <template <typename T> class Form, typename Type>
    void OpenHurricane::cutZone::transferVS(Form<Type> &f) const {
        if (HurMPI::parRun()) {
            typename feature<Type>::elementType *send;
            typename feature<Type>::elementType *recv;
            if (HurMPI::getProcRank() == sendProc_) {
                send = new
                    typename feature<Type>::elementType[sor_.size() * feature<Type>::nElements_];
                for (integer i = 0; i < sor_.size(); i++) {
                    for (int j = 0; j < feature<Type>::nElements_; j++) {
                        send[i * feature<Type>::nElements_ + j] = f[sor_[i]][j];
                    }
                }
                HurMPI::send(send, sor_.size() * feature<Type>::nElements_, feature<Type>::MPIType,
                             receivProc_, sendProc_, HurMPI::getComm());

                delete[] send;
            } else if (HurMPI::getProcRank() == receivProc_) {
                recv = new
                    typename feature<Type>::elementType[des_.size() * feature<Type>::nElements_];
                HurMPI::Status status;
                HurMPI::recv(recv, des_.size() * feature<Type>::nElements_, feature<Type>::MPIType,
                             sendProc_, sendProc_, HurMPI::getComm(), &status);
                for (integer i = 0; i < sor_.size(); i++) {
                    for (int j = 0; j < feature<Type>::nElements_; j++) {
                        f[des_[i]][j] = recv[i * feature<Type>::nElements_ + j];
                    }
                }
                delete[] recv;
            }
            //HurMPI::barrier(HurMPI::getComm());
        }
    }

    /*!\brief Transfer complex datatype.*/
    template <template <typename T> class Form, typename Type>
    void OpenHurricane::cutZone::transferComp(Form<Type> &f) const {
        if (HurMPI::parRun()) {
            typename feature<Type>::partType *send;
            typename feature<Type>::partType *recv;
            if (HurMPI::getProcRank() == sendProc_) {
                send = new typename feature<Type>::partType[sor_.size() * 2];
                for (integer i = 0; i < sor_.size(); i++) {
                    send[i * 2 + 0] = f[sor_[i]].Re();
                    send[i * 2 + 1] = f[sor_[i]].Im();
                }
                HurMPI::send(send, sor_.size() * 2, feature<Type>::MPIType, receivProc_, sendProc_,
                             HurMPI::getComm());

                delete[] send;
            } else if (HurMPI::getProcRank() == receivProc_) {
                recv = new typename feature<Type>::partType[des_.size() * 2];
                HurMPI::Status status;
                HurMPI::recv(recv, des_.size() * 2, feature<Type>::MPIType, sendProc_, sendProc_,
                             HurMPI::getComm(), &status);
                for (integer i = 0; i < sor_.size(); i++) {
                    f[des_[i]].Re() = recv[i * 2 + 0];
                    f[des_[i]].Im() = recv[i * 2 + 1];
                }
                delete[] recv;
            }
            //HurMPI::barrier(HurMPI::getComm());
        }
    }

    namespace periodicTypes {
        enum periodicType : short { TRANSLATIONAL = 0, ROTATIONAL = 1 };

        static const std::map<short, periodicType> periodicTypeMap =
            createMap<short, periodicType>(0, TRANSLATIONAL)(1, ROTATIONAL);

    } // namespace periodicTypes

    /**
     * \brief The class of periodic pair.
     */
    class periodicPair : public integerVector2DList {
    private:
        // Private data

        /**\brief The zone ID of the periodic face zone.*/
        integer periodicZone_;

        /**\brief The zone ID of the corresponding shadow face zone.*/
        integer shadowZone_;

        /**\brief The index of the first periodic face pair in the list.*/
        integer firstIndex_;

        /**\brief The index of the last periodic face pair in the list.*/
        integer lastIndex_;

        short type_;

    public:
        static void gatherPeriodicPair(const periodicPair &pp, periodicPair &gpp, const int root);

        // Constructors

        /**\brief Null constructor.*/
        periodicPair()
            : integerVector2DList(), periodicZone_(0), shadowZone_(0), firstIndex_(0),
              lastIndex_(0), type_(periodicTypes::TRANSLATIONAL) {}

        /**\brief Construct from components.*/
        inline periodicPair(const integer pZ, const integer sZ, const integer fI, const integer lI)
            : integerVector2DList(), periodicZone_(pZ), shadowZone_(sZ), firstIndex_(fI),
              lastIndex_(lI), type_(periodicTypes::TRANSLATIONAL) {
            integerVector2DList::resize(lI - fI + 1);
        }

        /**\brief Construct as copy.*/
        inline periodicPair(const periodicPair &pp)
            : integerVector2DList(pp), periodicZone_(pp.periodicZone_), shadowZone_(pp.shadowZone_),
              firstIndex_(pp.firstIndex_), lastIndex_(pp.lastIndex_), type_(pp.type_) {}

        /**\brief assignment.*/
        periodicPair &operator=(const periodicPair &);

        // Destructor
        inline ~periodicPair() noexcept {}

        hur_nodiscard inline integer periodicZone() const noexcept { return periodicZone_; }
        hur_nodiscard inline integer shadowZone() const noexcept { return shadowZone_; }
        hur_nodiscard inline integer firstIndex() const noexcept { return firstIndex_; }
        hur_nodiscard inline integer lastIndex() const noexcept { return lastIndex_; }

        inline void setFirstIndex(const integer fi) noexcept { firstIndex_ = fi; }

        inline void setLastIndex(const integer li) noexcept { lastIndex_ = li; }

        hur_nodiscard inline short type() const noexcept { return type_; }

        inline void setPeriodicZone(const integer pzid) noexcept { periodicZone_ = pzid; }

        inline void setShadowZone(const integer szid) noexcept { shadowZone_ = szid; }

        inline void setType(const short t) noexcept { type_ = t; }
        /**\brief Return one pair of faces. The index i must be greater than firstIndex_ and less than lastIndex_.*/
        hur_nodiscard inline const integerVector2D &operator()(const integer i) const noexcept {
#ifdef HUR_DEBUG
            if (i < firstIndex_ || i > lastIndex_) {
                LFatal("The index: %d is out of range: %d~%d", i, firstIndex_, lastIndex_);
            }
#endif // HUR_DEBUG

            return this->operator[](i - firstIndex_);
        }

        /**\brief Return one pair of faces. The index i must be greater than firstIndex_ and less than lastIndex_.*/
        hur_nodiscard inline integerVector2D &operator()(const integer i) noexcept {
#ifdef HUR_DEBUG
            if (i < firstIndex_ || i > lastIndex_) {
                LFatal("The index: %d is out of range: %d~%d", i, firstIndex_, lastIndex_);
            }
#endif // HUR_DEBUG

            return this->operator[](i - firstIndex_);
        }
    };

    /**
     * \brief The class of the periodic zone.
     */
    class perZone {
    public:
        enum axis : short { X_ROTATION = 0, Y_ROTATION, Z_ROTATION };

    private:
        /**\brief Index of the per zone in the face zones.*/
        integer id_;

        /**\brief Index of the pair zone in the face zones.*/
        integer pairId_;

        /**\brief Is this periodic shadow.*/
        bool isShadow_;

        /**\brief The sending processor number for current cutzone.*/
        integer sendProc_;

        /**\brief The receiving processor number for current cutzone.*/
        integer receivProc_;

        /**\brief the local perbc number of cell for sending  processor.*/
        integerList sor_;

        /**\brief The local perbc number of cell for receiving processor.*/
        integerList des_;

        /**\brief the local perbc number of face for sending  processor.*/
        integerList sorFace_;

        /**\brief The local perbc number of face for receiving processor.*/
        integerList desFace_;

        /*!\brief Rotation axis.*/
        short rotationAxis_;

        /**\brief Rotation matrix pointer.*/
        mutable tensor *RPtr_;

        /*!\brief Transform distance*/
        vector *tranDistancePtr_;

        /**\brief periodic bcType.*/
        short type_;

        /*!\brief Rotation radian.*/
        real phi_;

        // whether form two layer of ghost cell or not
        bool isTwoLayer_;

    public:

        inline perZone()
            : id_(0), pairId_(0), isShadow_(false), sendProc_(0), receivProc_(0), sor_(), des_(),
              sorFace_(), desFace_(), rotationAxis_(X_ROTATION), RPtr_(nullptr),
              tranDistancePtr_(nullptr), type_(0), phi_(0.0), isTwoLayer_(false) {}

        inline perZone(const integer idx, const integer sP, const integer dP)
            : id_(idx), pairId_(0), isShadow_(false), sendProc_(sP), receivProc_(dP), sor_(),
              des_(), sorFace_(), desFace_(), rotationAxis_(X_ROTATION), RPtr_(nullptr),
              tranDistancePtr_(nullptr), type_(0), phi_(0.0), isTwoLayer_(false) {}

        inline perZone(const integer idx, const integer sP, const integer dP,
                       const integerList &sorBc, const integerList &desBc)
            : id_(idx), pairId_(0), isShadow_(false), sendProc_(sP), receivProc_(dP), sor_(sorBc),
              des_(desBc), sorFace_(), desFace_(), rotationAxis_(X_ROTATION), RPtr_(nullptr),
              tranDistancePtr_(nullptr), type_(0), phi_(0.0), isTwoLayer_(false) {
            if (sor_.size() != des_.size()) {
                LFatal("The number of cut bc is not equal in periodic zone: %d sor.size() = %d   "
                       "des.size() = %d",
                       id_, sor_.size(), des_.size());
            }
        }

        perZone(const integer idx, const integer sP, const integer dP, const integerList &sorBc,
                const integerList &desBc, const short t, const tensor &R = tensor(Zero));

        inline perZone(const perZone &pz)
            : id_(pz.id_), pairId_(pz.pairId_), isShadow_(pz.isShadow_), sendProc_(pz.sendProc_),
              receivProc_(pz.receivProc_), sor_(pz.sor_), des_(pz.des_), sorFace_(pz.sorFace_),
              desFace_(pz.desFace_), rotationAxis_(X_ROTATION), RPtr_(nullptr),
              tranDistancePtr_(nullptr), type_(pz.type_), phi_(0.0), isTwoLayer_(pz.isTwoLayer_) {
            RPtr_ = pz.RPtr_;
            pz.RPtr_ = nullptr;
        }

        /**\brief Destructor.*/
        inline ~perZone() noexcept {
            HurDelete(RPtr_);
            HurDelete(tranDistancePtr_);
        }

        /**\brief Index of the per zone in the face zones.*/
        hur_nodiscard inline integer id() const noexcept { return id_; }

        /**\brief Index of the per zone in the face zones.*/
        hur_nodiscard inline integer &id() noexcept { return id_; }

        inline void setId(const integer id) noexcept { id_ = id; }

        inline void setPairId(const integer id) noexcept { pairId_ = id; }

        hur_nodiscard inline integer pairId() const noexcept { return pairId_; }

        hur_nodiscard inline integer &pairId() noexcept { return pairId_; }

        hur_nodiscard inline bool isShadow() const noexcept { return isShadow_; }

        inline void setShadow() noexcept { isShadow_ = true; }

        inline void setPeriodic() noexcept { isShadow_ = false; }

        hur_nodiscard inline integer sendProc() const noexcept { return sendProc_; }

        hur_nodiscard inline integer receivProc() const noexcept { return receivProc_; }

        hur_nodiscard inline integerList &sor() noexcept { return sor_; }

        hur_nodiscard inline const integerList &sor() const noexcept { return sor_; }

        hur_nodiscard inline integerList &des() noexcept { return des_; }

        hur_nodiscard inline const integerList &des() const noexcept { return des_; }

        hur_nodiscard inline integerList &sorFace() noexcept { return sorFace_; }

        hur_nodiscard inline const integerList &sorFace() const noexcept { return sorFace_; }

        hur_nodiscard inline integerList &desFace() noexcept { return desFace_; }

        hur_nodiscard inline const integerList &desFace() const noexcept { return desFace_; }

        hur_nodiscard inline short type() const noexcept { return type_; }

        inline void setType(const short t) noexcept { type_ = t; }

        hur_nodiscard inline bool isTranslational() const noexcept {
            return (type_ == periodicTypes::TRANSLATIONAL);
        }

        hur_nodiscard inline bool isRotational() const noexcept {
            return (type_ == periodicTypes::ROTATIONAL);
        }

        hur_nodiscard inline const tensor &RotionalMatrix() const {
            if (RPtr_ == nullptr) {
                LFatal("Attempt to access rotation tensor null pointer!");
            }

            return *RPtr_;
        }

        inline void setRotionalMatrix(const tensor _Rt) {
            HurDelete(RPtr_);
            type_ = periodicTypes::ROTATIONAL;
            RPtr_ = new tensor(_Rt);
        }

        hur_nodiscard inline const vector &tranDistance() const {
            if (tranDistancePtr_ == nullptr) {
                LFatal("Attempt to access translational distance null pointer!");
            }
            return *tranDistancePtr_;
        }

        inline void setTranDistance(const vector _tv) {
            HurDelete(tranDistancePtr_);
            type_ = periodicTypes::TRANSLATIONAL;
            tranDistancePtr_ = new vector(_tv);
        }

        hur_nodiscard inline real phi() const noexcept { return phi_; }

        hur_nodiscard inline integer perSize() {
#ifdef HUR_DEBUG
            if (sor_.size() != des_.size()) {
                LFatal("The number of cut bc is not equal in periodic zone: %d sor.size() = %d   "
                       "des.size() = %d",
                       id_, sor_.size(), des_.size());
            }
#endif // HUR_DEBUG
            return sor_.size();
        }

        hur_nodiscard inline bool isSameProc() const noexcept { return (sendProc_ == receivProc_); }

        hur_nodiscard inline bool ghostLayer() const noexcept { return isTwoLayer_; }

        hur_nodiscard inline bool isSendFromThisProc() const noexcept {
            return sendProc_ == HurMPI::getProcRank();
        }

        hur_nodiscard inline bool isThisProcReceiv() const noexcept {
            return receivProc_ == HurMPI::getProcRank();
        }

        inline void clearList() noexcept {
            sor_.clear();
            des_.clear();
        }

        perZone &operator=(const perZone &pz);

    public:
        /*!\brief Peridoic transfer.*/

        /**\brief Transfer basic datatype.*/
        template <template <typename Type> class Form, typename Type>
        void transfer(Form<Type> &) const;

        /**\brief Transfer (and rotate) VertorSpace datatype.*/
        template <template <typename Type> class Form, typename Type>
        void transferVS(Form<Type> &) const;

        /**\brief Transfer tensor datatype and (rotate the transposition of the tensor).*/
        void transferTensorTran(List<tensor> &) const;

        /**\brief Transfer complex datatype.*/
        template <template <typename Type> class Form, typename Type>
        void transferComp(Form<Type> &) const;

        void transferPoint(List<vector> &pointField, const vector &origin) const;

        /*!\brief Edit.*/

        /*!\brief set sendProc_.*/
        inline void setSendProc(const integer ip) noexcept { sendProc_ = ip; }

        /*!\brief set recvProc_.*/
        inline void setRecvProc(const integer ip) noexcept { receivProc_ = ip; }

        /*!\brief set isTwoLayer.*/
        inline void setGhostLayer(const bool layer) noexcept { isTwoLayer_ = layer; }

        template <typename Type> inline void rotate(Vector<Type> &v) const;

        template <typename Type> inline void rotate(Tensor<Type> &v) const;
    }; // namespace OpenHurricane

    /*!\brief Transfer basic datatype.*/
    template <template <typename T> class Form, typename Type>
    void OpenHurricane::perZone::transfer(Form<Type> &f) const {
        if (feature<Type>::nElements_ != 1) {
            std::string errMsg;
            errMsg = "Cannot transfer type: \"";
            errMsg += typeid(Type).name();
            errMsg += "\" (nElements = ";
            errMsg += std::to_string(feature<Type>::nElements_);
            errMsg += ") in periodic transfer fucntion.";
            errorAbortStr(errMsg);
        }
        Type *send;
        Type *recv;
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
                send = new Type[sor_.size()];
                for (integer i = 0; i < sor_.size(); i++) {
                    send[i] = f[sor_[i]];
                }
                HurMPI::send(send, sor_.size(), feature<Type>::MPIType, receivProc_, sendProc_,
                             HurMPI::getComm());
                delete[] send;
            } else if (HurMPI::getProcRank() == receivProc_) {
                recv = new Type[des_.size()];
                HurMPI::Status status;
                HurMPI::recv(recv, des_.size(), feature<Type>::MPIType, sendProc_, sendProc_,
                             HurMPI::getComm(), &status);
                for (integer i = 0; i < sor_.size(); i++) {
                    f[des_[i]] = recv[i];
                }
                delete[] recv;
            }
        }

        //HurMPI::barrier(HurMPI::getComm());
    } // End void transfer()

    /*!\brief Transfer VertorSpace datatype.*/
    template <template <typename T> class Form, typename Type>
    void OpenHurricane::perZone::transferVS(Form<Type> &f) const {
        typename feature<Type>::value_type *send;
        typename feature<Type>::value_type *recv;
        if (sendProc_ == receivProc_) //Periodic face pair on the same processor
        {
            if (HurMPI::getProcRank() == sendProc_) {
                for (integer i = 0; i < sor_.size(); i++) {
                    //std::cerr << "(";
                    for (int j = 0; j < feature<Type>::nElements_; j++) {
                        f[des_[i]][j] = f[sor_[i]][j];
                    }
                    if (RPtr_ != nullptr) {
                        f[des_[i]] = (*RPtr_) * f[des_[i]];                       
                    }
                }
            }
        } else // Periodic face pair on different processor
        {
            if (HurMPI::getProcRank() == sendProc_) {
                send = new
                    typename feature<Type>::elementType[sor_.size() * feature<Type>::nElements_];
                for (integer i = 0; i < sor_.size(); i++) {
                    for (int j = 0; j < feature<Type>::nElements_; j++) {
                        send[i * feature<Type>::nElements_ + j] = f[sor_[i]][j];
                    }
                }
                HurMPI::send(send, sor_.size() * feature<Type>::nElements_, feature<Type>::MPIType,
                             receivProc_, sendProc_, HurMPI::getComm());

                delete[] send;
            } else if (HurMPI::getProcRank() == receivProc_) {
                recv = new
                    typename feature<Type>::elementType[des_.size() * feature<Type>::nElements_];
                HurMPI::Status status;
                HurMPI::recv(recv, des_.size() * feature<Type>::nElements_, feature<Type>::MPIType,
                             sendProc_, sendProc_, HurMPI::getComm(), &status);
                for (integer i = 0; i < sor_.size(); i++) {
                    for (int j = 0; j < feature<Type>::nElements_; j++) {
                        f[des_[i]][j] = recv[i * feature<Type>::nElements_ + j];
                    }
                    if (RPtr_ != nullptr) {
                        f[des_[i]] = (*RPtr_) * f[des_[i]];
                    }
                }
                delete[] recv;
            }
        }
        //HurMPI::barrier(HurMPI::getComm());
    }

    /*!\brief Transfer complex datatype.*/
    template <template <typename T> class Form, typename Type>
    void OpenHurricane::perZone::transferComp(Form<Type> &f) const {

        typename feature<Type>::partType *send;
        typename feature<Type>::partType *recv;
        if (sendProc_ == receivProc_) //Periodic face pair on the same processor
        {
            if (HurMPI::getProcRank() == sendProc_) {
                for (integer i = 0; i < sor_.size(); i++) {
                    f[des_[i]].Re() = f[sor_[i]].Re();
                    f[des_[i]].Im() = f[sor_[i]].Im();
                }
            }
        } else // Periodic face pair on different processor
        {
            if (HurMPI::getProcRank() == sendProc_) {
                send = new typename feature<Type>::partType[sor_.size() * 2];
                for (integer i = 0; i < sor_.size(); i++) {
                    send[i * 2 + 0] = f[sor_[i]].Re();
                    send[i * 2 + 1] = f[sor_[i]].Im();
                }
                HurMPI::send(send, sor_.size() * 2, feature<Type>::MPIType, receivProc_, sendProc_,
                             HurMPI::getComm());

                delete[] send;
            } else if (HurMPI::getProcRank() == receivProc_) {
                recv = new typename feature<Type>::partType[des_.size() * 2];
                HurMPI::Status status;
                HurMPI::recv(recv, des_.size() * 2, feature<Type>::MPIType, sendProc_, sendProc_,
                             HurMPI::getComm(), &status);
                for (integer i = 0; i < sor_.size(); i++) {
                    f[des_[i]].Re() = recv[i * 2 + 0];
                    f[des_[i]].Im() = recv[i * 2 + 1];
                }
                delete[] recv;
            }
        }
        //HurMPI::barrier(HurMPI::getComm());
    }

    template <typename Type> inline void perZone::rotate(Vector<Type> &v) const {
        if (RPtr_ != nullptr) {
            v = (*RPtr_) * v;
        }
    }

    template <typename Type> inline void perZone::rotate(Tensor<Type> &v) const {
        if (RPtr_ != nullptr) {
            v = (*RPtr_) * v * (RPtr_)->transpose();
        }
    }

    /**
     * \brief Transfer basic datatype.
     *    The size of each f[i] must be the same.
     */
    template <> void OpenHurricane::perZone::transfer(Array<realArray> &f) const;
} // namespace OpenHurricane