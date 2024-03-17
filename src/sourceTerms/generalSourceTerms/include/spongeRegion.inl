﻿#include "spongeRegion.hpp"
/*!
 * \file spongeRegion.inl
 * \brief The In-Line functions of the <i>spongeRegion.hpp</i> file.
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

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::spongeRegion::d() const {
    if (dPtr_ == nullptr) {
        makeDist();
    }
    return *dPtr_;
}

inline const OpenHurricane::integerList &OpenHurricane::spongeRegion::spongeLayer() const noexcept {
    return regionSourceTerms::regionList()[0];
}

template <typename Type>
hur_nodiscard inline Type OpenHurricane::spongeRegion::getRefState(const string &UName) const {
    Type tmp(0);

    const auto &bcCont = iter().cont().subController("boundaryCondition");
    if (bcCont.found(faceZoneName_)) {
        const auto &fzCont = bcCont.subController(faceZoneName_);
        tmp = fzCont.findType<Type>(UName, tmp);
    } else {
        errorAbortStr(("Cannot find: " + faceZoneName_ + " in: " + bcCont.name()));
    }

    return tmp;
}