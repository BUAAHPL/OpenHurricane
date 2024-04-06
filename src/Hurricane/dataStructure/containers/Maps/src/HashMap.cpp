/*!
 * \file HashMap.cpp
 * \brief Main subroutines of HashMap.
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
#include "HashMap.hpp"
#include <cmath>

hur_nodiscard OpenHurricane::uinteger OpenHurricane::HashMapBase::mapSize(const uinteger n,
                                                                  const bool usePrime) {
    if (n <= 0) {
        return 0;
    }

    if (!usePrime) {
        uinteger size = n;
        if (size & (size - 1)) {
            size = 1;
            while (size < unsigned(n)) {
                size <<= 1;
            }
        }
        return size;
    } else {
        uinteger p = (n % 2) ? n + 2 : n + 1;
        uinteger i = 0;
        while (p <= maxSize_) {
            for (i = (uinteger)std::sqrt(p); i > 2; --i) {
                if ((p % i) == 0) {
                    break;
                }
            }
            if (i == 2) {
                break;
            } else {
                p += 2;
            }
        }

        return p;
    }
}
