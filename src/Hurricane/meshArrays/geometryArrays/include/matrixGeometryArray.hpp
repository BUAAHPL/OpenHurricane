/*!
 * \file matrixGeometryArray.hpp
 * \brief Headers of the matrix geometryArray.
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

#include "geometryArray.hpp"

namespace OpenHurricane {

#define checkMeshForArray(f1, f2)       \
    if (&(f1).mesh() != &(f2).mesh()) { \
        LFatal("different mesh");       \
    }
    /**
     * \brief The template class of array based on geometry for matrix parameter.
     * \tparam matrices - The component type (matrix).
     * \tparam GeometryMesh - The mesh type.
     */
    template <class matrices, class GeometryMesh>
    class matrixGeometryArray : public Array<matrices> {
    public:
        using Base = Array<matrices>;
        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;
        using Mesh = typename GeometryMesh::Mesh;

        using elementType = typename Base::elementType;

    private:
        /*!\brief Const reference to mesh.*/
        const Mesh &mesh_;

    public:
        hur_nodiscard inline static const matrixGeometryArray &nullObject() {
            return NullRefObj::nullRef<matrixGeometryArray>();
        }

        matrixGeometryArray(const Mesh &mesh)
            : Base(GeometryMesh::internalArraySize(mesh)), mesh_(mesh) {}

        matrixGeometryArray(const Mesh &mesh, const zero)
            : Base(GeometryMesh::internalArraySize(mesh), Zero), mesh_(mesh) {}

        inline matrixGeometryArray(const Mesh &mesh, const matrices &m0)
            : Base(GeometryMesh::internalArraySize(mesh), m0), mesh_(mesh) {}

        inline matrixGeometryArray(const Mesh &mesh, const Base &M0) : Base(M0), mesh_(mesh) {
            if (M0.size() && (M0.size() != GeometryMesh::internalArraySize(mesh))) {
                LFatal("Cannot link a Array to mesh in different size");
            }
        }

        inline matrixGeometryArray(const Mesh &mesh, Base &&M0) noexcept
            : Base(std::move(M0)), mesh_(mesh) {
            if (Base::size() && (Base::size() != GeometryMesh::internalArraySize(mesh))) {
                LFatal("Cannot link a Array to mesh in different size");
            }
        }

        inline matrixGeometryArray(const matrixGeometryArray &ma) : Base(ma), mesh_(ma.mesh_) {}

        inline matrixGeometryArray(matrixGeometryArray &&ma) noexcept
            : Base(std::move(ma)), mesh_(ma.mesh_) {}

        hur_nodiscard inline uniquePtr<matrixGeometryArray> clone() const {
            return uniquePtr<matrixGeometryArray>(new matrixGeometryArray(*this));
        }

        /*!\brief Destructor.*/
        inline virtual ~matrixGeometryArray() noexcept {}

        /*!\brief Return the const reference to the mesh.*/
        hur_nodiscard inline const Mesh &mesh() const noexcept { return mesh_; }

        /*!\brief Return const access to the Array.*/
        hur_nodiscard inline const Base &array_ref() const noexcept { return *this; }

        /*!\brief Return non-const access to the Array.*/
        hur_nodiscard inline Base &array_ref() noexcept { return *this; }

        using Base::average;
        using Base::weightedAverage;

        inline void clear() noexcept { Base::clear(); }

        inline matrixGeometryArray &operator=(const matrixGeometryArray &ma) {
            if (this != std::addressof(ma)) {
#ifdef HUR_DEBUG
                checkMeshForArray(*this, ma);
#endif // HUR_DEBUG
                Base::operator=(ma);
            }
            return *this;
        }

        matrixGeometryArray &operator=(const Base &ma) {
            if (this != std::addressof(ma)) {
                if (Base::size() == ma.size()) {
                    Base::operator=(ma);
                } else {
                    const integer minSize = min(Base::size(), ma.size());
                    for (integer i = 0; i < minSize; ++i) {
                        Base::operator[](i) = ma[i];
                    }
                }
            }
            return *this;
        }

        inline matrixGeometryArray &operator=(matrixGeometryArray &&ma) noexcept {
#ifdef HUR_DEBUG
            checkMeshForArray(*this, ma);
#endif // HUR_DEBUG
            Base::operator=(std::move(ma));
            return *this;
        }

        matrixGeometryArray &operator=(Base &&ma) noexcept {
            if (Base::size() == ma.size()) {
                Base::operator=(std::move(ma));
            } else {
                const integer minSize = min(Base::size(), ma.size());
                for (integer i = 0; i < minSize; ++i) {
                    Base::operator[](i) = ma[i];
                }
                ma.clear();
            }
            return *this;
        }

        inline matrixGeometryArray &operator=(const matrices &m) {
            Base::operator=(m);
            return *this;
        }

        inline matrixGeometryArray &operator=(const zero) {
            Base::operator=(Zero);
            return *this;
        }

        inline void operator+=(const matrixGeometryArray &ma) {
#ifdef HUR_DEBUG
            checkMeshForArray(*this, ma);
#endif // HUR_DEBUG
            Base::operator+=(ma);
        }

        inline void operator+=(const Base &ma) { Base::operator+=(ma); }

        inline void operator+=(const matrices &m) { Base::operator+=(m); }

        inline void operator-=(const matrixGeometryArray &ma) {
#ifdef HUR_DEBUG
            checkMeshForArray(*this, ma);
#endif // HUR_DEBUG
            Base::operator-=(ma);
        }
        inline void operator-=(const Base &ma) { Base::operator-=(ma); }
        inline void operator-=(const matrices &m) { Base::operator-=(m); }
    };

} // namespace OpenHurricane

namespace OpenHurricane {
#undef checkMeshForArray
}