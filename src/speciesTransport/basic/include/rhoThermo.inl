/*!
 * \file rhoThermo.inl
 * \brief In-Line subroutines of the <i>rhoThermo.hpp</i> file.
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

inline const OpenHurricane::runtimeMesh &OpenHurricane::rhoThermo::mesh() const noexcept {
    return rho_.mesh();
}

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::rho() const noexcept {
    return rho_;
}

hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::rho() noexcept {
    return rho_;
}

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::p() const noexcept {
    return p_;
}

hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::p() noexcept {
    return p_;
}

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::T() const noexcept {
    return T_;
}

hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::T() noexcept {
    return T_;
}

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::E() const noexcept {
    return E_;
}

hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::E() noexcept {
    return E_;
}

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::gamma() const noexcept {
    return gamma_;
}

hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::gamma() noexcept {
    return gamma_;
}

hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::cp() noexcept {
    return cp_;
}

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::cp() const noexcept {
    return cp_;
}

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::mu() const noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::mu() noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::kappa() const noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::rhoThermo::kappa() noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard inline OpenHurricane::mixture &OpenHurricane::rhoThermo::mixtures() noexcept {
    return *mixPtr_;
}

hur_nodiscard inline bool OpenHurricane::rhoThermo::isMixturesPtrNUll() const noexcept {
    return mixPtr_.isNull();
}

hur_nodiscard inline bool OpenHurricane::rhoThermo::isSingleSpecie() const noexcept {
    return mixPtr_->isSingular();
}
