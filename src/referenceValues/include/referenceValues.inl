/*!
 * \file referenceValues.inl
 * \brief The In-Line functions of reference values
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

namespace OpenHurricane {
    inline referenceValues::~referenceValues() noexcept {}

    hur_nodiscard inline real referenceValues::Ma() const noexcept {
        return Ma_;
    }

    hur_nodiscard inline real referenceValues::vMag() const noexcept {
        return u_;
    }

    hur_nodiscard inline real referenceValues::p() const noexcept {
        return p_;
    }

    hur_nodiscard inline real referenceValues::T() const noexcept {
        return T_;
    }

    hur_nodiscard inline real referenceValues::Tt() const noexcept {
        return real();
    }

    hur_nodiscard inline real referenceValues::rho() const noexcept {
        return rho_;
    }

    hur_nodiscard inline real referenceValues::mu() const noexcept {
        return mu_;
    }

    hur_nodiscard inline real referenceValues::gama() const noexcept {
        return gamma_;
    }

    hur_nodiscard inline real referenceValues::area() const noexcept {
        return area_;
    }

    hur_nodiscard inline real referenceValues::length() const noexcept {
        return length_;
    }

} // namespace OpenHurricane