#include "reordering.hpp"
/*!
 * \file reordering.inl
 * \brief In-Line subroutines of the <i>reordering.hpp</i> file.
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
inline const OpenHurricane::integerList &
OpenHurricane::reordering::perm() const noexcept {
    return perm_;
}
inline OpenHurricane::integerList &OpenHurricane::reordering::perm() {
    return perm_;
}

inline const OpenHurricane::integerList &
OpenHurricane::reordering::iperm() const noexcept {
    return perm_;
}

inline OpenHurricane::integerList &OpenHurricane::reordering::iperm() {
    return iperm_;
}

inline OpenHurricane::reordering::reordering(decomposingMeshByMetis &dcpM)
    : decomposedMesh_(dcpM), nvtxs_(0), xadj_(), adjncy_(), perm_(), iperm_() {}