/*!
 * \file basicFunctions.hpp
 * \brief Header of basic functions
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
#include "preset.hpp"
#include <typeinfo>
namespace OpenHurricane {
    /**\brief The global function of swap*/
    template <class T> inline void Swap(T &a, T &b) noexcept {
        T tmp = a;
        a = b;
        b = tmp;
    }

    template <class Type1, class Type2> hur_nodiscard inline bool isSameType(const Type2 &t) {
        return typeid(t) == typeid(Type1);
    }

    template <class Type1, class Type2> hur_nodiscard inline bool hasSameBase(const Type2 &t) {
        const Type2 *tPtr = &t;
        return dynamic_cast<const Type1 *>(tPtr);
    }
} // namespace OpenHurricane
