/*!
 * \file geometryArray.inl
 * \brief In-Line subroutines of the <i>geometryArray.hpp</i> file.
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

namespace OpenHurricane {
#define checkMeshForArray(f1, f2)       \
    if (&(f1).mesh() != &(f2).mesh()) { \
        LFatal("different mesh");       \
    }

} // namespace OpenHurricane

namespace OpenHurricane {
    template <class Type, class GeometryMesh>
    const std::string geometryArray<Type, GeometryMesh>::className_("geometryArray");
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::setCurRhs(const Type &curRhs) const {
    curRhs_ = curRhs;
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(const object &ob,
                                                                   const Mesh &mesh)
    : object(ob), Array<Type>(GeometryMesh::size(mesh)), mesh_(mesh), isOnlyInternalArray_(false),
      rhs0_(Type(0)), curRhs_(Type(0)), rhsAvePtr_(nullptr), rhsMaxPtr_(nullptr),
      internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr), boundariesPtr_(nullptr),
      lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr), timeSumPtr_(nullptr),
      limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr), lastArrayArray_() {
    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(object &&ob, const Mesh &mesh)
    : object(std::move(ob)), Array<Type>(GeometryMesh::size(mesh)), mesh_(mesh),
      isOnlyInternalArray_(false), rhs0_(Type(0)), curRhs_(Type(0)), rhsAvePtr_(nullptr),
      rhsMaxPtr_(nullptr), internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr),
      boundariesPtr_(nullptr), lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr),
      timeSumPtr_(nullptr), limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr),
      lastArrayArray_() {
    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(object &&ob, const Mesh &mesh,
                                                                   const zero)
    : object(std::move(ob)), Array<Type>(GeometryMesh::size(mesh), Zero), mesh_(mesh),
      isOnlyInternalArray_(false), rhs0_(Type(0)), curRhs_(Type(0)), rhsAvePtr_(nullptr),
      rhsMaxPtr_(nullptr), internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr),
      boundariesPtr_(nullptr), lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr),
      timeSumPtr_(nullptr), limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr),
      lastArrayArray_() {
    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(const object &ob,
                                                                   const Mesh &mesh, const zero)
    : object(ob), Array<Type>(GeometryMesh::size(mesh), Zero), mesh_(mesh),
      isOnlyInternalArray_(false), rhs0_(Type(0)), curRhs_(Type(0)), rhsAvePtr_(nullptr),
      rhsMaxPtr_(nullptr), internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr),
      boundariesPtr_(nullptr), lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr),
      timeSumPtr_(nullptr), limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr),
      lastArrayArray_() {
    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(const object &ob,
                                                                   const Mesh &mesh, const Type &t)
    : object(ob), Array<Type>(GeometryMesh::size(mesh), t), mesh_(mesh),
      isOnlyInternalArray_(false), rhs0_(Type(0)), curRhs_(Type(0)), rhsAvePtr_(nullptr),
      rhsMaxPtr_(nullptr), internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr),
      boundariesPtr_(nullptr), lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr),
      timeSumPtr_(nullptr), limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr),
      lastArrayArray_() {
    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(object &&ob, const Mesh &mesh,
                                                                   const Type &t)
    : object(std::move(ob)), Array<Type>(GeometryMesh::size(mesh), t), mesh_(mesh),
      isOnlyInternalArray_(false), rhs0_(Type(0)), curRhs_(Type(0)), rhsAvePtr_(nullptr),
      rhsMaxPtr_(nullptr), internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr),
      boundariesPtr_(nullptr), lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr),
      timeSumPtr_(nullptr), limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr),
      lastArrayArray_() {
    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(const object &ob,
                                                                   const Mesh &mesh,
                                                                   const Array<Type> &f)
    : object(ob), Array<Type>(f), mesh_(mesh), isOnlyInternalArray_(false), rhs0_(Type(0)),
      curRhs_(Type(0)), rhsAvePtr_(nullptr), rhsMaxPtr_(nullptr), internalArrayPtr_(nullptr),
      ghostArrayPtr_(nullptr), boundariesPtr_(nullptr), lastArrayPtr(nullptr),
      lastLastArrayPtr(nullptr), rhsPtr_(nullptr), timeSumPtr_(nullptr), limiterPtr_(nullptr),
      diagSourcePtr_(nullptr), gradPtr_(nullptr), lastArrayArray_() {
    if (f.size() && (f.size() != GeometryMesh::size(mesh))) {
        errorAbortStr(("Cannot link a Array to mesh in different size for \"" + name() + "\""));
    }

    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(const object &ob,
                                                                   const Mesh &mesh,
                                                                   Array<Type> &&f)
    : object(ob), Array<Type>(std::move(f)), mesh_(mesh), isOnlyInternalArray_(false),
      rhs0_(Type(0)), curRhs_(Type(0)), rhsAvePtr_(nullptr), rhsMaxPtr_(nullptr),
      internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr), boundariesPtr_(nullptr),
      lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr), timeSumPtr_(nullptr),
      limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr), lastArrayArray_() {
    if (Array<Type>::size() && (Array<Type>::size() != GeometryMesh::size(mesh))) {
        errorAbortStr(("Cannot link a Array to mesh in different size for \"" + name() + "\""));
    }

    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(object &&ob, const Mesh &mesh,
                                                                   const Array<Type> &f)
    : object(std::move(ob)), Array<Type>(f), mesh_(mesh), isOnlyInternalArray_(false),
      rhs0_(Type(0)), curRhs_(Type(0)), rhsAvePtr_(nullptr), rhsMaxPtr_(nullptr),
      internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr), boundariesPtr_(nullptr),
      lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr), timeSumPtr_(nullptr),
      limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr), lastArrayArray_() {
    if (f.size() && (f.size() != GeometryMesh::size(mesh))) {
        errorAbortStr(("Cannot link a Array to mesh in different size for \"" + name() + "\""));
    }

    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(object &&ob, const Mesh &mesh,
                                                                   Array<Type> &&f)
    : object(std::move(ob)), Array<Type>(std::move(f)), mesh_(mesh), isOnlyInternalArray_(false),
      rhs0_(Type(0)), curRhs_(Type(0)), rhsAvePtr_(nullptr), rhsMaxPtr_(nullptr),
      internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr), boundariesPtr_(nullptr),
      lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr), timeSumPtr_(nullptr),
      limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr), lastArrayArray_() {
    if (Array<Type>::size() && (Array<Type>::size() != GeometryMesh::size(mesh))) {
        errorAbortStr(("Cannot link a Array to mesh in different size for \"" + name() + "\""));
    }

    if (GeometryMesh::size(mesh) == GeometryMesh::internalArraySize(mesh)) {
        isOnlyInternalArray_ = true;
    }
    if (!object::hasSetOutputName()) {
        setOutputName();
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(
    const geometryArray<Type, GeometryMesh> &gf)
    : object(gf), Array<Type>(gf), mesh_(gf.mesh_), isOnlyInternalArray_(gf.isOnlyInternalArray_),
      rhs0_(gf.rhs0_), curRhs_(gf.curRhs_), rhsAvePtr_(nullptr), rhsMaxPtr_(nullptr),
      internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr), boundariesPtr_(nullptr),
      lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr), timeSumPtr_(nullptr),
      limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr), lastArrayArray_() {
    if (gf.rhsAvePtr_) {
        rhsAvePtr_.reset(new Type());
        *rhsAvePtr_ = *gf.rhsAvePtr_;
    }
    if (gf.rhsMaxPtr_) {
        rhsMaxPtr_.reset(new Type());
        *rhsMaxPtr_ = *gf.rhsMaxPtr_;
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(
    geometryArray<Type, GeometryMesh> &&gf) noexcept
    : object(std::move(gf)), Array<Type>(std::move(gf)), mesh_(gf.mesh_),
      isOnlyInternalArray_(std::move(gf.isOnlyInternalArray_)), rhs0_(std::move(gf.rhs0_)), curRhs_(std::move(gf.curRhs_)),
      rhsAvePtr_(nullptr), rhsMaxPtr_(nullptr), internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr),
      boundariesPtr_(nullptr), lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr),
      timeSumPtr_(nullptr), limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr),
      lastArrayArray_() {
    if (gf.rhsAvePtr_) {
        rhsAvePtr_ = std::move(gf.rhsAvePtr_);
    }
    if (gf.rhsMaxPtr_) {
        rhsMaxPtr_ = std::move(gf.rhsMaxPtr_);
    }
    if (gf.internalArrayPtr_) {
        internalArrayPtr_ = std::move(gf.internalArrayPtr_);
    }
    if (gf.ghostArrayPtr_) {
        ghostArrayPtr_ = std::move(gf.ghostArrayPtr_);
    }
    if (gf.lastArrayPtr) {
        lastArrayPtr = std::move(gf.lastArrayPtr);
    }
    if (gf.lastLastArrayPtr) {
        lastLastArrayPtr = std::move(gf.lastLastArrayPtr);
    }
    if (gf.rhsPtr_) {
        rhsPtr_ = std::move(gf.rhsPtr_);
    }

    if (gf.timeSumPtr_) {
        timeSumPtr_ = std::move(gf.timeSumPtr_);
    }
    if (gf.limiterPtr_) {
        limiterPtr_ = std::move(gf.limiterPtr_);
    }
    if (gf.diagSourcePtr_) {
        diagSourcePtr_ = std::move(gf.diagSourcePtr_);
    }
    if (gf.gradPtr_) {
        gradPtr_ = std::move(gf.gradPtr_);
    }
    if (gf.boundariesPtr_ != nullptr) {
        boundariesPtr_ = gf.boundariesPtr_;
        gf.boundariesPtr_ = nullptr;
    }
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::geometryArray(
    const object &ob, const geometryArray<Type, GeometryMesh> &gf)
    : object(ob), Array<Type>(gf), mesh_(gf.mesh_), isOnlyInternalArray_(gf.isOnlyInternalArray_),
      rhs0_(gf.rhs0_), curRhs_(gf.curRhs_), rhsAvePtr_(nullptr), rhsMaxPtr_(nullptr),
      internalArrayPtr_(nullptr), ghostArrayPtr_(nullptr), boundariesPtr_(nullptr),
      lastArrayPtr(nullptr), lastLastArrayPtr(nullptr), rhsPtr_(nullptr), timeSumPtr_(nullptr),
      limiterPtr_(nullptr), diagSourcePtr_(nullptr), gradPtr_(nullptr), lastArrayArray_() {
    if (gf.rhsAvePtr_) {
        rhsAvePtr_.reset(new Type());
        *rhsAvePtr_ = *gf.rhsAvePtr_;
    }
    if (gf.rhsMaxPtr_) {
        rhsMaxPtr_.reset(new Type());
        *rhsMaxPtr_ = *gf.rhsMaxPtr_;
    }
}

template <class Type, class GeometryMesh>
hur_nodiscard inline OpenHurricane::uniquePtr<OpenHurricane::geometryArray<Type, GeometryMesh>>
OpenHurricane::geometryArray<Type, GeometryMesh>::clone() const {
    return uniquePtr<geometryArray<Type, GeometryMesh>>(
        new geometryArray<Type, GeometryMesh>(*this));
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh>::~geometryArray() noexcept {
    HurDelete(boundariesPtr_);
}

template <class Type, class GeometryMesh>
hur_nodiscard inline const typename OpenHurricane::geometryArray<Type, GeometryMesh>::Mesh &
OpenHurricane::geometryArray<Type, GeometryMesh>::mesh() const noexcept {
    return mesh_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline const OpenHurricane::Array<Type> &
OpenHurricane::geometryArray<Type, GeometryMesh>::array_ref() const noexcept {
    return *this;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline OpenHurricane::Array<Type> &
OpenHurricane::geometryArray<Type, GeometryMesh>::array_ref() noexcept {
    return *this;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline OpenHurricane::integer
OpenHurricane::geometryArray<Type, GeometryMesh>::internalArraySize() const noexcept {
    return GeometryMesh::internalArraySize(mesh_);
}

template <class Type, class GeometryMesh>
hur_nodiscard inline typename OpenHurricane::geometryArray<Type, GeometryMesh>::InternalArray &
OpenHurricane::geometryArray<Type, GeometryMesh>::internalArray() {
    if (!internalArrayPtr_) {
        internalArrayPtr_.reset(new InternalArray(*this, this->internalArraySize()));
    }
    return *internalArrayPtr_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline const typename OpenHurricane::geometryArray<Type, GeometryMesh>::InternalArray &
OpenHurricane::geometryArray<Type, GeometryMesh>::internalArray() const {
    if (!internalArrayPtr_) {
        internalArrayPtr_.reset(new InternalArray(*this, this->internalArraySize()));
    }
    return *internalArrayPtr_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline typename OpenHurricane::geometryArray<Type, GeometryMesh>::GhostArray &
OpenHurricane::geometryArray<Type, GeometryMesh>::ghostArray() {
    if (!ghostArrayPtr_) {
        if (!isOnlyInternalArray_) {
            ghostArrayPtr_.reset(new GhostArray(*this, this->size() - this->internalArraySize(),
                                                this->internalArraySize()));
        }
    }
    return *ghostArrayPtr_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline const typename OpenHurricane::geometryArray<Type, GeometryMesh>::GhostArray &
OpenHurricane::geometryArray<Type, GeometryMesh>::ghostArray() const {
    if (!ghostArrayPtr_) {
        if (!isOnlyInternalArray_) {
            ghostArrayPtr_.reset(new GhostArray(*this, this->size() - this->internalArraySize(),
                                                this->internalArraySize()));
        }
    }
    return *ghostArrayPtr_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline int OpenHurricane::geometryArray<Type, GeometryMesh>::nElements() const noexcept {
    return feature<Type>::nElements_;
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::setOutputName() {
    if (object::isTemporary()) {
        return;
    }
    if (feature<Type>::nElements_ > 1) {
        object::outputVarNameL().resize(feature<Type>::nElements_);
        string nN;
        string nNl;
        for (int i = 0; i < feature<Type>::nElements_; ++i) {
            nN += "\"";
            nN += object::name();
            nN += "[";
            nN += std::to_string(i);
            nN += "]\"";
            if (i < feature<Type>::nElements_ - 1) {
                nN += ",";
            }
            nNl = "";
            nNl += object::name();
            nNl += "[";
            nNl += std::to_string(i);
            nNl += "]";
            object::outputVarNameL()[i] = nNl;
        }
        object::outputVarName() = nN;
    }
    object::hasSetOutputName() = true;
}

template <class Type, class GeometryMesh>
hur_nodiscard hur_forceinline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::lastArray() const {
    if (!lastArrayPtr) {
        string rn = object::name() + "lastArray";
        lastArrayPtr.reset(new geometryArray<Type, GeometryMesh>(
            object(rn, object::tb(), object::NOT_WRITE), mesh_));
    }
    return *lastArrayPtr;
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::storeLastArray() const {
    for (integer i = 0; i < this->size(); ++i) {
        lastArray()[i] = (*this)[i];
    }
}

template <class Type, class GeometryMesh>
hur_nodiscard hur_forceinline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::rhs() const {
    if (!rhsPtr_) {
        string rn = object::name() + "Rhs";
        rhsPtr_.reset(new geometryArray<Type, GeometryMesh>(
            object(rn, object::tb(), object::NOT_WRITE), mesh_, Zero));
    }
    return *rhsPtr_;
}

template <class Type, class GeometryMesh>
hur_nodiscard hur_forceinline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::getTimeSumPtr() const {
    if (!timeSumPtr_) {
        string rn = object::name() + "TimeSum";
        timeSumPtr_.reset(new geometryArray<Type, GeometryMesh>(
            object(rn, object::tb(), object::NOT_WRITE), mesh_, Zero));
    }
    return *timeSumPtr_;
}

template <class Type, class GeometryMesh>
hur_forceinline void
OpenHurricane::geometryArray<Type, GeometryMesh>::calcTimeSumPtr(const real &dt) const {
    if (!timeSumPtr_) {
        string rn = object::name() + "TimeSum";
        timeSumPtr_.reset(new geometryArray<Type, GeometryMesh>(
            object(rn, object::tb(), object::NOT_WRITE), mesh_, Zero));
    }
    (*timeSumPtr_) += dt * (*this);
}

template <class Type, class GeometryMesh>
void OpenHurricane::geometryArray<Type, GeometryMesh>::getAveAndMaxRHS() const {
    checkWarningStr(("The average and maximum residual for: " + name() + " is not defined"));
}

template <class Type, class GeometryMesh>
hur_nodiscard inline const Type &
OpenHurricane::geometryArray<Type, GeometryMesh>::rhsAve0() const noexcept {
    return rhs0_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline const Type &
OpenHurricane::geometryArray<Type, GeometryMesh>::curRhsAve() const noexcept {
    return curRhs_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline const Type &
OpenHurricane::geometryArray<Type, GeometryMesh>::rhsAve() const noexcept {
    if (!rhsAvePtr_) {
        rhsAvePtr_.reset(new Type());
    }
    return *rhsAvePtr_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline const Type &
OpenHurricane::geometryArray<Type, GeometryMesh>::rhsMax() const noexcept {
    if (!rhsMaxPtr_) {
        rhsMaxPtr_.reset(new Type());
    }
    return *rhsMaxPtr_;
}

template <class Type, class GeometryMesh>
hur_nodiscard hur_forceinline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::limiter() const {
    if (!limiterPtr_) {
        string rn = object::name() + "Limiter";
        limiterPtr_.reset(new geometryArray<Type, GeometryMesh>(
            object(rn, object::tb(), object::NOT_WRITE), mesh_, Zero));
    }
    return *limiterPtr_;
}

template <class Type, class GeometryMesh>
hur_nodiscard hur_forceinline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::diagSource() const {
    if (!diagSourcePtr_) {
        string rn = object::name() + "diagSource";
        diagSourcePtr_.reset(new geometryArray<Type, GeometryMesh>(
            object(rn, object::tb(), object::NOT_WRITE), mesh_, Zero));
    }
    return *diagSourcePtr_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline bool
OpenHurricane::geometryArray<Type, GeometryMesh>::hasGradArray() const noexcept {
    return !(gradPtr_.isNull());
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::clearGradArray() noexcept {
    gradPtr_.clear();
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::setLastArrayArray(const integer lastI) {
    lastArrayArray_.resize(lastI);
    for (integer i = 0; i < lastI; ++i) {
        lastArrayArray_[i].resize(internalArraySize());
    }
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::setLastArrayArray(const integer lastI,
                                                                            const zero) {
    lastArrayArray_.resize(lastI);
    for (integer i = 0; i < lastI; ++i) {
        lastArrayArray_[i].resize(internalArraySize(), Zero);
    }
}

template <class Type, class GeometryMesh>
hur_nodiscard hur_forceinline OpenHurricane::Array<OpenHurricane::Array<Type>> &
OpenHurricane::geometryArray<Type, GeometryMesh>::lastArrayArray() {
    return lastArrayArray_;
}

template <class Type, class GeometryMesh>
hur_nodiscard inline const OpenHurricane::Array<OpenHurricane::Array<Type>> &
OpenHurricane::geometryArray<Type, GeometryMesh>::lastArrayArray() const {
    return lastArrayArray_;
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::updateBoundary(integer layerI) {
    std::string errMsg;
    errMsg = "Ghost cell value calucation for this type of geometryArray named ";
    errMsg += this->name();
    errMsg += " has not been specified.";
    errorAbortStr(errMsg);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::updateBoundary() {
    for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
        this->updateBoundary(layerI);
    }
}

template <class Type, class GeometryMesh>
Type OpenHurricane::geometryArray<Type, GeometryMesh>::initialize() {
    return Zero;
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::clear() noexcept {
    Array<Type>::clear();
    object::clear();
    rhsAvePtr_.clear();
    rhsMaxPtr_.clear();
    internalArrayPtr_.clear();
    ghostArrayPtr_.clear();
    lastArrayPtr.clear();
    lastLastArrayPtr.clear();
    rhsPtr_.clear();
    timeSumPtr_.clear();
    limiterPtr_.clear();
    diagSourcePtr_.clear();
    gradPtr_.clear();
    HurDelete(boundariesPtr_);
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::operator=(
    const geometryArray<Type, GeometryMesh> &gf) {
    if (this != &gf) {
#ifdef HUR_DEBUG
        checkMeshForArray(*this, gf);
#endif // HUR_DEBUG

        Array<Type>::operator=(gf);
        isOnlyInternalArray_ = gf.isOnlyInternalArray_;

        if (gf.rhsPtr_) {
            rhs() = gf.rhs();
        }
    }
    return *this;
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::operator=(const Array<Type> &gf) {
    if (this != std::addressof(gf)) {
        if (Array<Type>::size() == gf.size()) {
            Array<Type>::operator=(gf);
        } else {
            const integer minSize = min(Array<Type>::size(), gf.size());
            for (integer i = 0; i < minSize; ++i) {
                Array<Type>::operator[](i) = gf[i];
            }
        }
    }
    return *this;
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::operator=(
    geometryArray<Type, GeometryMesh> &&gf) noexcept {
#ifdef HUR_DEBUG
    checkMeshForArray(*this, gf);
#endif // HUR_DEBUG
    object::operator=(std::move(gf));
    Array<Type>::operator=(std::move(gf));
    isOnlyInternalArray_ = std::move(gf.isOnlyInternalArray_);
    return *this;
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::operator=(Array<Type> &&gf) noexcept {
    if (Array<Type>::size() == gf.size()) {
        Array<Type>::operator=(std::move(gf));
    } else {
        const integer minSize = min(Array<Type>::size(), gf.size());
        for (integer i = 0; i < minSize; ++i) {
            Array<Type>::operator[](i) = gf[i];
        }
        gf.clear();
    }
    return *this;
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::operator=(const Type &t) {
    Array<Type>::operator=(t);
    return *this;
}

template <class Type, class GeometryMesh>
inline OpenHurricane::geometryArray<Type, GeometryMesh> &
OpenHurricane::geometryArray<Type, GeometryMesh>::operator=(const zero) {
    Array<Type>::operator=(Zero);
    return *this;
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator+=(
    const geometryArray<Type, GeometryMesh> &gf) {
#ifdef HUR_DEBUG
    checkMeshForArray(*this, gf);
#endif // HUR_DEBUG

    Array<Type>::operator+=(gf);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator+=(const Array<Type> &gf) {
    Array<Type>::operator+=(gf);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator+=(const Type &t) {
    Array<Type>::operator+=(t);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator-=(
    const geometryArray<Type, GeometryMesh> &gf) {
#ifdef HUR_DEBUG
    checkMeshForArray(*this, gf);
#endif // HUR_DEBUG

    Array<Type>::operator-=(gf);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator-=(const Array<Type> &gf) {
    Array<Type>::operator-=(gf);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator-=(const Type &t) {
    Array<Type>::operator-=(t);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator*=(
    const geometryArray<real, GeometryMesh> &gf) {
#ifdef HUR_DEBUG
    checkMeshForArray(*this, gf);
#endif // HUR_DEBUG

    Array<Type>::operator*=(gf);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator*=(const Array<real> &gf) {
    Array<Type>::operator*=(gf);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator*=(const real &t) {
    Array<Type>::operator*=(t);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator/=(
    const geometryArray<Type, GeometryMesh> &gf) {
#ifdef HUR_DEBUG
    checkMeshForArray(*this, gf);
#endif // HUR_DEBUG

    Array<Type>::operator/=(gf);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator/=(const Array<Type> &gf) {
    Array<Type>::operator/=(gf);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::operator/=(const Type &t) {
    Array<Type>::operator/=(t);
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::writeOutput(fileOsstream &fos) const {
    // Only write the internale filed's value.
    Array<Type>::writeToStream(fos, this->internalArraySize());
}

template <class Type, class GeometryMesh>
inline void
OpenHurricane::geometryArray<Type, GeometryMesh>::writeOutputByMaster(fileOsstream &fos) const {
    errorAbortStr(("Cannot write data : " + this->name() + " to file : " + fos.name()));
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::writeOutput(fileOsstream &fos,
                                                                      const integer fzid) const {}

template <class Type, class GeometryMesh>
inline void
OpenHurricane::geometryArray<Type, GeometryMesh>::writeMinMaxOutput(fileOsstream &fos) const {
    Array<Type>::writeMinMaxToStreamWithFactor(fos, this->internalArraySize());
}

template <class Type, class GeometryMesh>
inline void
OpenHurricane::geometryArray<Type, GeometryMesh>::writeMinMaxOutputByMaster(fileOsstream &fos) const {
    Array<Type>::writeMinMaxToStreamWithFactorByMaster(fos, this->internalArraySize());
}

template <class Type, class GeometryMesh>
inline void
OpenHurricane::geometryArray<Type, GeometryMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                                const integer fzid) const {}

template <class Type, class GeometryMesh>
inline void
OpenHurricane::geometryArray<Type, GeometryMesh>::writeResiduals(fileOsstream &fos,
                                                             integer intervalStep) const {
    if (Iteration().hasSubIteration()) {
        Array<Type>::writeAveToPout(fos, rhs(), mesh().cellVolume(), mesh().nCells(),
                                    mesh().allMeshCellNumber(), rhs0_,
                                    (Iteration().subIter().cSubStep() == intervalStep),
                                    (Iteration().subIter().cSubStep() == (2 * intervalStep)));
    } else {
        Array<Type>::writeAveToPout(fos, rhs(), mesh().cellVolume(), mesh().nCells(),
                                    mesh().allMeshCellNumber(), rhs0_,
                                    (Iteration().totalStep() == intervalStep),
                                    (Iteration().totalStep() == (2 * intervalStep)));
    }
}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::writeRelay(fileOsstream &fos) const {
    // Only write the internale filed's value.
    Array<Type>::writeToStream(fos, this->internalArraySize());
}

template <class Type, class GeometryMesh>
inline void
OpenHurricane::geometryArray<Type, GeometryMesh>::writeRelay(hdf5O &fos, const bool writeLast,
                                                         const bool writeToGroup) const {}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::readRelay(const hdf5I &fos,
                                                                    const bool readLast,
                                                                    const bool readFromGroup) {}

template <class Type, class GeometryMesh>
inline void OpenHurricane::geometryArray<Type, GeometryMesh>::interpolateRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {}

namespace OpenHurricane {
#undef checkMeshForArray
}