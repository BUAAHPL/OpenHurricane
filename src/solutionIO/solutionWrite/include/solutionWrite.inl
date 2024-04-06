#include "solutionWrite.hpp"
/*!
 * \file solutionWrite.inl
 * \brief In-Line subroutines of the <i>solutionWrite.hpp</i> file.
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

hur_nodiscard inline const OpenHurricane::flowModel &
OpenHurricane::solutionWrite::flows() const noexcept {
    return flows_;
}

hur_nodiscard inline const OpenHurricane::iteration &
OpenHurricane::solutionWrite::iter() const noexcept {
    return iter_;
}

hur_nodiscard inline const OpenHurricane::runtimeMesh &
OpenHurricane::solutionWrite::mesh() const noexcept {
    return mesh_;
}

hur_nodiscard inline bool OpenHurricane::solutionWrite::found(const string &name) const {
    return sets_.find(name) != sets_.end();
}

hur_nodiscard inline const auto &OpenHurricane::solutionWrite::rho() const noexcept {
    return flows_.rho();
}

hur_nodiscard inline const auto &OpenHurricane::solutionWrite::p() const noexcept {
    return flows_.p();
}

hur_nodiscard inline const auto &OpenHurricane::solutionWrite::gama() const noexcept {
    return flows_.gama();
}

hur_nodiscard inline const auto &OpenHurricane::solutionWrite::T() const noexcept {
    return flows_.T();
}

hur_nodiscard inline const auto &OpenHurricane::solutionWrite::yi() const noexcept {
    return flows_.mixtures().Yi();
}

hur_nodiscard inline const auto &OpenHurricane::solutionWrite::v() const noexcept {
    return flows_.v();
}
