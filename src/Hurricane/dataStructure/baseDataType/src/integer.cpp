/*!
 * \file integer.cpp
 * \brief The subroutines and functions of integer
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
#include "integer.hpp"
#include "errorAbort.hpp"

namespace OpenHurricane {

    hur_nodiscard constexpr integer pow(integer a, integer b) noexcept {
        integer result = 1;
        for (integer i = 0; i < b; i++) {
            result *= a;
        }
#ifdef HUR_DEBUG
        if (b < integer(0)) {
            LFatal("negative value for b is not supported");
        }
#endif // HUR_DEBUG
        return result;
    }

    bool readInteger(const char *_str, integer &val) {
        val = std::atoi(_str);
        if (errno == ERANGE || errno == EINVAL) {
            return false;
        }
        return true;
    }
    bool readInteger(const std::string &str, integer &val) {
        return readInteger(str.c_str(), val);
    }
} // namespace OpenHurricane