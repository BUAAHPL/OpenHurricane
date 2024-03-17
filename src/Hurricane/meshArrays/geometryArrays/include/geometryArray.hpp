/*!
 * \file geometryArray.hpp
 * \brief Headers of the geometry array.
 *        The subroutines and functions are in the <i>geometryArray.cpp</i> file.
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
#include "Array.hpp"
#include "iteration.hpp"
#include "object.hpp"
#include "zoneMesh.hpp"

namespace OpenHurricane {
    /**\brief Forward declare*/
    template <class Type, class geoMesh> class boundary;

    /**
     * \brief The template class of array based on geometry.
     * \tparam Type - The component type.
     * \tparam GeometryMesh - The mesh type.
     */
    template <class Type, class GeometryMesh>
    class geometryArray : public object, public Array<Type> {
    public:
        using Base = Array<Type>;
        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;

        using Mesh = typename GeometryMesh::Mesh;

        using elementType = typename Base::elementType;
        using subArray = typename Base::subArray;
        using InternalArray = typename Base::subArray;
        using GhostArray = typename Base::subArray;

        using gradientType = geometryArray<typename outerProduct<vector, Type>::type, GeometryMesh>;

        template <class bType, class geoMesh>
        class boundaryCondition : public PtrList<boundary<bType, geoMesh>> {
        private:
            /**\brief Hold reference to face zone list.*/
            const faceZoneList &faceZones_;

        public:
            /**\brief Construct from components.*/
            boundaryCondition(const faceZoneList &fZL, geometryArray<bType, geoMesh> &gf,
                              const controller &cont)
                : faceZones_(fZL) {
                this->resize(fZL.size());
                for (integer i = 0; i < this->size(); i++) {
                    this->set(i, boundary<bType, geoMesh>::creator(fZL[i], gf, cont).release());
                }
            }

            /**\brief Destructor.*/
            ~boundaryCondition() noexcept {}

            /**\brief Return const reference to face zone list.*/
            inline const faceZoneList &faceZones() const { return faceZones_; }

            void calculate() {
                for (integer i = 0; i < this->size(); ++i) {
                    this->operator[](i).updateBoundary();
                }
            }
        };

    private:
        /*!\brief Const reference to mesh.*/
        const Mesh &mesh_;

        /*!\brief Is the geometry Array a entire internal Array.*/
        bool isOnlyInternalArray_;

        /**\brief Initial residual.*/
        Type rhs0_;

        /**\brief Current residual. */
        mutable Type curRhs_;

        mutable uniquePtr<Type> rhsAvePtr_;

        mutable uniquePtr<Type> rhsMaxPtr_;

        /*!\brief The internal Array pointer*/
        mutable uniquePtr<InternalArray> internalArrayPtr_;

        /*!\brief The ghost Array pointer*/
        mutable uniquePtr<GhostArray> ghostArrayPtr_;

        /**\brief The boundary condition pointer.*/
        mutable boundaryCondition<Type, GeometryMesh> *boundariesPtr_;

        /*!\brief The pointer to Array of the last time step .*/
        mutable uniquePtr<geometryArray<Type, GeometryMesh>> lastArrayPtr;

        /*!\brief The pointer to Array of the last last time step .*/
        mutable uniquePtr<geometryArray<Type, GeometryMesh>> lastLastArrayPtr;

        /*!\brief Residual term for primitive parameters.*/
        mutable uniquePtr<geometryArray<Type, GeometryMesh>> rhsPtr_;

        /*!\brief calculate time sum field for primitive parameters.*/
        mutable uniquePtr<geometryArray<Type, GeometryMesh>> timeSumPtr_;

        /*!\brief Limiters for primitive parameters.*/
        mutable uniquePtr<geometryArray<Type, GeometryMesh>> limiterPtr_;

        /*!\brief Diagnal source Jacobian for primitive parameters.*/
        mutable uniquePtr<geometryArray<Type, GeometryMesh>> diagSourcePtr_;

        /*!\brief Gradient for primitive parameters.*/
        mutable uniquePtr<gradientType> gradPtr_;

        /**\brief Only for internal Array.*/
        Array<Base> lastArrayArray_;

    protected:
        /**\brief Set current residual value. */
        virtual void setCurRhs(const Type &curRhs) const;

    public:
        declareClassName(geometryArray);

        hur_nodiscard inline static const geometryArray<Type, GeometryMesh> &nullObject() {
            return NullRefObj::nullRef<geometryArray<Type, GeometryMesh>>();
        }

        geometryArray(const object &ob, const Mesh &mesh);

        geometryArray(object &&ob, const Mesh &mesh);

        geometryArray(const object &ob, const Mesh &mesh, const zero);

        geometryArray(object &&ob, const Mesh &mesh, const zero);

        geometryArray(const object &ob, const Mesh &mesh, const Type &t);

        geometryArray(object &&ob, const Mesh &mesh, const Type &t);

        geometryArray(const object &ob, const Mesh &mesh, const Array<Type> &);

        geometryArray(const object &ob, const Mesh &mesh, Array<Type> &&);

        geometryArray(object &&ob, const Mesh &mesh, const Array<Type> &);

        geometryArray(object &&ob, const Mesh &mesh, Array<Type> &&);

        geometryArray(const geometryArray<Type, GeometryMesh> &);

        geometryArray(geometryArray<Type, GeometryMesh> &&) noexcept;

        geometryArray(const object &ob, const geometryArray<Type, GeometryMesh> &);

        hur_nodiscard uniquePtr<geometryArray<Type, GeometryMesh>> clone() const;

        /*!\brief Destructor.*/
        virtual ~geometryArray() noexcept;

        /*!\brief Return the const reference to the mesh.*/
        hur_nodiscard inline const Mesh &mesh() const noexcept;

        /*!\brief Return const access to the Array.*/
        hur_nodiscard inline const Array<Type> &array_ref() const noexcept;

        /*!\brief Return non-const access to the Array.*/
        hur_nodiscard inline Array<Type> &array_ref() noexcept;

        /*!\brief Return the size of the internal Array.*/
        hur_nodiscard inline integer internalArraySize() const noexcept;

        /*!\brief Is the geometry Array a entire internal Array.*/
        hur_nodiscard inline bool isOnlyInternalArray() const noexcept {
            return isOnlyInternalArray_;
        }

        /*!\brief Return access to the internal Array of the geometry Array.*/
        hur_nodiscard inline InternalArray &internalArray();

        /*!\brief Return const access to the internal Array of the geometry Array.*/
        hur_nodiscard inline const InternalArray &internalArray() const;

        /*!
         * \brief Return access to the ghost Array of the geometry Array.
         * If isOnlyInternalArray_==true, return nullptr.
         */
        hur_nodiscard inline GhostArray &ghostArray();

        /*!
         * \brief Return const access to the ghost Array of the geometry Array.
         * If isOnlyInternalArray_==true, return nullptr.
         */
        hur_nodiscard inline const GhostArray &ghostArray() const;

        hur_nodiscard inline virtual int nElements() const noexcept;

        void setOutputName();

        using Base::average;
        using Base::weightedAverage;

        /*!\brief The pointer to Array of the last time step .*/
        hur_nodiscard geometryArray<Type, GeometryMesh> &lastArray() const;

        void storeLastArray() const;

        hur_nodiscard geometryArray<Type, GeometryMesh> &rhs() const;

        virtual void calcTimeSumPtr(const real &dt) const;

        hur_nodiscard geometryArray<Type, GeometryMesh> &getTimeSumPtr() const;

    public:
        void getAveAndMaxRHS() const;

        /*!\brief Return the first step average residual.*/
        hur_nodiscard inline const Type &rhsAve0() const noexcept;

        /*!\brief Return the current step average residual.*/
        hur_nodiscard inline const Type &curRhsAve() const noexcept;

        /*!\brief Return the average residual.*/
        hur_nodiscard inline const Type &rhsAve() const noexcept;

        /*!\brief Return the maximum residual.*/
        hur_nodiscard inline const Type &rhsMax() const noexcept;

    public:
        hur_nodiscard geometryArray<Type, GeometryMesh> &limiter() const;

        hur_nodiscard geometryArray<Type, GeometryMesh> &diagSource() const;

        hur_nodiscard inline gradientType &grad() const {
            if (!gradPtr_) {
                string rn = object::name() + "Grad";
                gradPtr_.reset(
                    new gradientType(object(rn, object::tb(), object::NOT_WRITE), mesh_, Zero));
            }
            return *gradPtr_;
        }

        hur_nodiscard inline bool hasGradArray() const noexcept;

        inline void clearGradArray() noexcept;

        void setLastArrayArray(const integer lastI);

        void setLastArrayArray(const integer lastI, const zero);

        hur_nodiscard Array<Array<Type>> &lastArrayArray();
        hur_nodiscard const Array<Array<Type>> &lastArrayArray() const;

        void updateBoundary(integer layerI);

        void updateBoundary();

        Type initialize();

        hur_nodiscard virtual inline realArray realComponent(const int i) const { return realArray(); }

        void clear() noexcept;

        // Member operators

        geometryArray<Type, GeometryMesh> &operator=(const geometryArray<Type, GeometryMesh> &);
        geometryArray<Type, GeometryMesh> &operator=(const Array<Type> &);
        geometryArray<Type, GeometryMesh> &operator=(geometryArray<Type, GeometryMesh> &&) noexcept;
        geometryArray<Type, GeometryMesh> &operator=(Array<Type> &&) noexcept;
        geometryArray<Type, GeometryMesh> &operator=(const Type &);
        geometryArray<Type, GeometryMesh> &operator=(const zero);

        void operator+=(const geometryArray<Type, GeometryMesh> &);
        void operator+=(const Array<Type> &);
        void operator+=(const Type &);

        void operator-=(const geometryArray<Type, GeometryMesh> &);
        void operator-=(const Array<Type> &);
        void operator-=(const Type &);

        void operator*=(const geometryArray<real, GeometryMesh> &);
        void operator*=(const Array<real> &);
        void operator*=(const real &);

        void operator/=(const geometryArray<Type, GeometryMesh> &);
        void operator/=(const Array<Type> &);
        void operator/=(const Type &);

        // Write

        /*!
         * \brief Write object to output file.
         * Note: it is a virtual function, and should be rewritten in geometry Array class.
         */
        virtual void writeOutput(fileOsstream &fos) const;

        /*!
         * \brief Write object to output file by master.
         * Note: it is a virtual function, and should be rewritten in geometry Array class.
         */
        virtual void writeOutputByMaster(fileOsstream &fos) const;

        /*!
         * \brief Write object to output file.
         * Note: it is a virtual function, and should be rewritten in geometry Array class.
         */
        virtual void writeOutput(fileOsstream &fos, const integer fzid) const;

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in geometry Array class.
         */
        virtual void writeMinMaxOutput(fileOsstream &fos) const;

        /*!
         * \brief Write object minimum and maximum value to output file by master.
         * Note: it is a virtual function, and should be rewritten in geometry Array class.
         */
        virtual void writeMinMaxOutputByMaster(fileOsstream &fos) const;

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in geometry Array class.
         */
        virtual void writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;

        /*!
         * \brief Write object residuals to monitor.
         * Note: it is a virtual function, and should be rewritten in geometry Array class.
         */
        virtual void writeResiduals(fileOsstream &fos, integer intervalStep) const;

        /*!
         * \brief Write object to relay file.
         * Note: it is a virtual function, and should be rewritten in geometry Array class.
         */
        virtual void writeRelay(fileOsstream &fos) const;

        /*!
         * \brief Write object to relay file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeRelay(hdf5O &fos, const bool writeLast, const bool writeToGroup) const;

        virtual void readRelay(const hdf5I &fos, const bool readLast, const bool readFromGroup);

        virtual void interpolateRelay(const hdf5I &fos, const bool readLast,
                                      const bool readFromGroup);
    };

#ifndef checkMeshForArray
#define checkMeshForArray(f1, f2)       \
    if (&(f1).mesh() != &(f2).mesh()) { \
        LFatal("different mesh");       \
    }
#endif // !checkMeshForArray

    template <class Type, class GeometryMesh>
    geometryArray<Type, GeometryMesh> operator+(const geometryArray<Type, GeometryMesh> &f1,
                                                const geometryArray<Type, GeometryMesh> &f2) {
#ifdef HUR_DEBUG
        checkMeshForArray(f1, f2);
#endif // HUR_DEBUG
        geometryArray<Type, GeometryMesh> gf(
            object(f1.name() + "+" + f2.name(), f1.mesh(), object::NOT_WRITE, object::TEMPORARY),
            f1.mesh());
        for (integer i = 0; i < f1.size(); ++i) {
            gf[i] = f1[i] + f2[i];
        }
        return gf;
    }

    template <class Type, class GeometryMesh>
    geometryArray<Type, GeometryMesh> operator+(geometryArray<Type, GeometryMesh> &&f1,
                                                const geometryArray<Type, GeometryMesh> &f2) {
#ifdef HUR_DEBUG
        checkMeshForArray(f1, f2);
#endif // HUR_DEBUG
        geometryArray<Type, GeometryMesh> tf(std::move(f1));
        for (integer i = 0; i < tf.size(); ++i) {
            tf[i] += f2[i];
        }
        return tf;
    }

    template <class Type, class GeometryMesh>
    geometryArray<Type, GeometryMesh> operator+(const geometryArray<Type, GeometryMesh> &f1,
                                                geometryArray<Type, GeometryMesh> &&f2) {
#ifdef HUR_DEBUG
        checkMeshForArray(f1, f2);
#endif // HUR_DEBUG
        geometryArray<Type, GeometryMesh> tf(std::move(f2));
        for (integer i = 0; i < f1.size(); ++i) {
            tf[i] = f1[i] + tf[i];
        }
        return tf;
    }

    template <class Type, class GeometryMesh>
    geometryArray<Type, GeometryMesh> operator+(geometryArray<Type, GeometryMesh> &&f1,
                                                geometryArray<Type, GeometryMesh> &&f2) {
#ifdef HUR_DEBUG
        checkMeshForArray(f1, f2);
#endif // HUR_DEBUG
        geometryArray<Type, GeometryMesh> tf(std::move(f1));
        for (integer i = 0; i < tf.size(); ++i) {
            tf[i] += f2[i];
        }
        return tf;
    }

    template <class Type, class GeometryMesh>
    geometryArray<Type, GeometryMesh> operator-(const geometryArray<Type, GeometryMesh> &f1,
                                                const geometryArray<Type, GeometryMesh> &f2) {
#ifdef HUR_DEBUG
        checkMeshForArray(f1, f2);
#endif // HUR_DEBUG
        geometryArray<Type, GeometryMesh> gf(
            object(f1.name() + "-" + f2.name(), f1.mesh(), object::NOT_WRITE, object::TEMPORARY),
            f1.mesh());
        for (integer i = 0; i < f1.size(); ++i) {
            gf[i] = f1[i] - f2[i];
        }
        return gf;
    }

    template <class Type, class GeometryMesh>
    geometryArray<Type, GeometryMesh> operator-(const geometryArray<Type, GeometryMesh> &f1,
                                                geometryArray<Type, GeometryMesh> &&f2) {
#ifdef HUR_DEBUG
        checkMeshForArray(f1, f2);
#endif // HUR_DEBUG
        geometryArray<Type, GeometryMesh> tf(std::move(f2));
        for (integer i = 0; i < f1.size(); ++i) {
            tf[i] = f1[i] - tf[i];
        }
        return tf;
    }

    template <class Type, class GeometryMesh>
    geometryArray<Type, GeometryMesh> operator-(geometryArray<Type, GeometryMesh> &&f1,
                                                const geometryArray<Type, GeometryMesh> &f2) {
#ifdef HUR_DEBUG
        checkMeshForArray(f1, f2);
#endif // HUR_DEBUG
        geometryArray<Type, GeometryMesh> tf(std::move(f1));
        for (integer i = 0; i < tf.size(); ++i) {
            tf[i] -= f2[i];
        }
        return tf;
    }

    template <class Type, class GeometryMesh>
    geometryArray<Type, GeometryMesh> operator-(geometryArray<Type, GeometryMesh> &&f1,
                                                geometryArray<Type, GeometryMesh> &&f2) {
#ifdef HUR_DEBUG
        checkMeshForArray(f1, f2);
#endif // HUR_DEBUG
        geometryArray<Type, GeometryMesh> tf(std::move(f1));
        for (integer i = 0; i < tf.size(); ++i) {
            tf[i] -= f2[i];
        }
        return tf;
    }

    template <class Type, class GeometryMesh>
    inline void component(
        geometryArray<typename geometryArray<Type, GeometryMesh>::elementType, GeometryMesh> &comp,
        const geometryArray<Type, GeometryMesh> gf, const int d) {
        component(dynamic_cast<typename Array<Type>::elementType &>(comp),
                  dynamic_cast<const Array<Type> &>(gf), d);
    }

    template <class Type, class GeometryMesh>
    void transpose(geometryArray<Type, GeometryMesh> &tran,
                   const geometryArray<Type, GeometryMesh> gf) {
        transpose(dynamic_cast<typename Array<Type>::elementType &>(tran),
                  dynamic_cast<const Array<Type> &>(gf));
    }

#ifdef checkMeshForArray
#undef checkMeshForArray
#endif // checkMeshForArray

} //  namespace OpenHurricane
#include "geometryArray.inl"
