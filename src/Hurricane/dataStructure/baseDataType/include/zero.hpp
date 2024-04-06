/*!
 * \file zero.hpp
 * \brief Header of zero
 *       The subroutines and functions are in the <i>zero.inl</i> file.
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
#include "integer.hpp"

#include "string.hpp"

namespace OpenHurricane {
    /**
     * \brief The class of zero. You can set any basic data bcType to zero.
     */
    class zero {
    public:
        using value_type = zero;
        zero() = default;

        inline constexpr operator bool() const noexcept { return 0; }

        inline constexpr operator integer() const noexcept { return 0; }

        inline constexpr operator float() const noexcept { return 0.0f; }
        inline constexpr operator double() const noexcept { return 0.0; }
        inline constexpr operator long double() const noexcept { return 0.0l; }
    };
    template <> hur_nodiscard inline std::string toString(const zero &) {
        return "0";
    }
    static const zero Zero;

    inline zero operator+(const zero &, const zero &) noexcept {
        return Zero;
    }

    template <class Type> inline const Type &operator+(const Type &t, const zero &) noexcept {
        return t;
    }

    template <class Type> inline const Type &operator+(const zero &, const Type &t) noexcept {
        return t;
    }

    template <class Type> inline const Type &operator-(const Type &t, const zero &) noexcept {
        return t;
    }

    template <class Type> inline const Type &operator-(const zero &, const Type &t) noexcept {
        return -t;
    }

    inline zero operator*(const zero &, const zero &) noexcept {
        return Zero;
    }

    template <class Type> inline zero operator*(const Type &t, const zero &) noexcept {
        return Zero;
    }

    template <class Type> inline zero operator*(const zero &, const Type &t) noexcept {
        return Zero;
    }

    template <class Type> inline zero operator/(const zero &, const Type &t) noexcept {
        return Zero;
    }
} // namespace OpenHurricane