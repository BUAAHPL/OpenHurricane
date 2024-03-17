/*!
 * \file integer.hpp
 * \brief Header of lable
 *       The subroutines and functions are in the <i>integer.cpp</i> file.
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
#include "int32.hpp"
#include "int64.hpp"

namespace OpenHurricane {
#ifdef HURRICANE_32_INT
    using integer = int32_t;
    using uinteger = uint32_t;

#elif defined HURRICANE_64_INT
    using integer = int64_t;
    using uinteger = uint64_t;
#else
#error "integer.hpp: HURRICANE_32_INT or HURRICANE_64_INT must be define"
#endif

    static constexpr integer integerMin = feature<integer>::max;
    static constexpr integer integerMax = feature<integer>::max;
    static constexpr uinteger uintegerMax = feature<uinteger>::max;

    /**\brief Raise one integer to the power of another.*/
    hur_nodiscard constexpr integer pow(integer a, integer b) noexcept;

    inline constexpr integer &setComponent(integer &l) noexcept {
        return l;
    }
    inline constexpr integer component(const integer l) noexcept {
        return l;
    }

    /*!\brief Read a integer variable from character array.*/
    bool readInteger(const char *_str, integer &val);
    bool readInteger(const std::string &str, integer &val);

    hur_nodiscard constexpr inline integer sign(const integer &s) noexcept {
        return (s >= 0) ? 1 : -1;
    }

    hur_nodiscard constexpr inline integer max(const integer &l, const integer &k) noexcept {
        return (l > k) ? l : k;
    }

    hur_nodiscard constexpr inline integer min(const integer &l, const integer &k) noexcept {
        return (l < k) ? l : k;
    }

    hur_nodiscard constexpr inline uinteger max(const uinteger &l, const uinteger &k) noexcept {
        return (l > k) ? l : k;
    }

    hur_nodiscard constexpr inline uinteger max(const integer &l, const uinteger &k) noexcept {
        return ((uinteger)l > k) ? l : k;
    }

    hur_nodiscard constexpr inline uinteger max(const uinteger &l, const integer &k) noexcept {
        return (l > (uinteger)k) ? l : k;
    }

    hur_nodiscard constexpr inline uinteger min(const uinteger &l, const uinteger &k) noexcept {
        return (l < k) ? l : k;
    }

    hur_nodiscard constexpr inline uinteger min(const integer &l, const uinteger &k) noexcept {
        return ((uinteger)l < k) ? l : k;
    }

    hur_nodiscard constexpr inline uinteger min(const uinteger &l, const integer &k) noexcept {
        return (l < (uinteger)k) ? l : k;
    }
} // namespace OpenHurricane
