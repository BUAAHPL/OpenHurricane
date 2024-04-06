#include "timeMarching.hpp"
/*!
 * \file timeMarching.inl
 * \brief In-Line subroutines of the <i>timeMarching.hpp</i> file.
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

hur_nodiscard inline const OpenHurricane::runtimeMesh &
OpenHurricane::timeMarching::mesh() const noexcept {
    return mesh_;
}

hur_nodiscard inline OpenHurricane::realArray &OpenHurricane::timeMarching::dt() noexcept {
    return dt_;
}

hur_nodiscard inline const OpenHurricane::faceVector2DArray &
OpenHurricane::timeMarching::rai() const noexcept {
    return rai_;
}

hur_nodiscard inline const OpenHurricane::faceVector2DArray &
OpenHurricane::timeMarching::rav() const noexcept {
    return rav_;
}

hur_nodiscard inline bool OpenHurricane::timeMarching::isEuler() const noexcept {
    return isEuler_;
}

hur_nodiscard inline bool OpenHurricane::timeMarching::explicitSource() const noexcept {
    return sourceTermTreat_ == sourceTermTreating::EXPLICITSOURCE;
}

hur_nodiscard inline bool OpenHurricane::timeMarching::diagonalImpSource() const noexcept {
    return sourceTermTreat_ == sourceTermTreating::KIMDIAGIMPLICITSOURCE;
}

hur_nodiscard inline bool OpenHurricane::timeMarching::fullJacobianSource() const noexcept {
    return (sourceTermTreat_ == sourceTermTreating::FULLJACOBIAN);
}

hur_nodiscard inline bool OpenHurricane::timeMarching::fullJacobianSourceTable() const noexcept {
    return (sourceTermTreat_ == sourceTermTreating::FULLJACOBIANTABLE);
}

hur_nodiscard inline OpenHurricane::cellRealSquareMatrixArray &
OpenHurricane::timeMarching::Jacobian() {
    if (!JacPtr_) {
        LFatal("Attempt to access null pointer for Jacobian matrix array");
    }
    return *JacPtr_;
}

hur_nodiscard inline const OpenHurricane::cellRealSquareMatrixArray &
OpenHurricane::timeMarching::Jacobian() const {
    if (!JacPtr_) {
        LFatal("Attempt to access null pointer for Jacobian matrix array");
    }
    return *JacPtr_;
}

hur_nodiscard inline const OpenHurricane::cellRealArray &
OpenHurricane::timeMarching::shockFactor() const noexcept {
    return shockFactor_;
}

hur_nodiscard inline OpenHurricane::cellRealArray &
OpenHurricane::timeMarching::shockFactor() noexcept {
    return shockFactor_;
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::timeMarching::nPrimitives() const noexcept {
    return nPrimitives_;
}

hur_nodiscard inline OpenHurricane::integer
OpenHurricane::timeMarching::countParam() const noexcept {
    return countParam_;
}

hur_nodiscard inline OpenHurricane::realArrayArray &OpenHurricane::timeMarching::dq() noexcept {
    return dq_;
}

inline const OpenHurricane::realArrayArray &OpenHurricane::timeMarching::dq() const noexcept {
    return dq_;
}

inline void OpenHurricane::timeMarching::transferDq() {
    fv::transfer(mesh_, dq_, true);
}

template <class Type>
inline void OpenHurricane::timeMarching::addObject(geometryArray<Type, cellMesh> &gf) {
    objectList_.append(dynamic_cast<object *>(&gf));
    for (int i = 0; i < gf.nElements(); ++i) {
        integerVector2D tmp;
        tmp.x() = integer(objectList_.size() - 1);
        tmp.y() = i;
        paramMap_.append(tmp);
        countParam_++;
    }
}

hur_nodiscard inline bool OpenHurricane::timeMarching::updatedTemperature() const noexcept {
    return false;
}
