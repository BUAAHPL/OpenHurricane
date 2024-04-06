﻿#include "sourceTerm.hpp"
/*!
 * \file sourceTerm.inl
 * \brief The In-Line functions of the <i>sourceTerm.hpp</i> file.
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

hur_nodiscard inline const OpenHurricane::flowModel &OpenHurricane::sourceTerm::flows() const noexcept {
    return flows_;
}

hur_nodiscard inline const OpenHurricane::runtimeMesh &OpenHurricane::sourceTerm::mesh() const noexcept {
    return flows_.mesh();
}

hur_nodiscard inline const OpenHurricane::iteration &OpenHurricane::sourceTerm::iter() const noexcept {
    return iter_;
}

inline void OpenHurricane::sourceTerm::addSourceTerms(cellRealArray &rho) const {}

inline void OpenHurricane::sourceTerm::addSourceTerms(const cellRealArray &rho,
                                                  cellVectorArray &U) const {}

inline void OpenHurricane::sourceTerm::addSourceTerms(const cellRealArray &rho,
                                                  cellRealArray &phi) const {}
