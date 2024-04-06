/*!
 * \file packs.cpp
 * \brief The subroutines and functions of packs
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

#include "packs.hpp"

void OpenHurricane::packageMap::unbind() noexcept {
    for (auto &e : objectList_) {
        e = nullptr;
    }
}

void OpenHurricane::packageMap::clear() noexcept {
    unbind();
    objectList_.clear();
    countPackage_ = 0;
    packageMap_.clear();
}

void OpenHurricane::packageMap::addObject(object &ob) {
    objectList_.append(std::addressof(ob));
    for (int i = 0; i < ob.nElements(); ++i) {
        integerVector2D tmp;
        tmp.x() = integer(objectList_.size() - 1);
        tmp.y() = i;
        packageMap_.append(tmp);
        countPackage_++;
    }
}
