/*!
 * \file HurFormat.hpp
 * \brief Header of OpenHurricane formats
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
#include "string.hpp"
#include <functional>

namespace OpenHurricane {

    // For format is not supported in C++ 17

    template <typename... Args> constexpr std::string hurFormat(const char *fmt, Args... args) {
        if constexpr (sizeof...(args) > 0) {
            auto length = std::snprintf(nullptr, 0, fmt, args...);
            if (length <= 0) {
                return std::string(fmt);
            }
            std::string str;
            str.resize(static_cast<size_t>(length) + 1);
            std::snprintf(&str.front(), static_cast<size_t>(length) + 1, fmt, args...);
            str.resize(static_cast<size_t>(length));
            return str;
        } else {
            return std::string(fmt);
        }
    }
} // namespace OpenHurricane