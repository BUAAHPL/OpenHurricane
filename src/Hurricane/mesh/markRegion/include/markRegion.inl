﻿#include "markRegion.hpp"
/*!
 * \file markRegion.inl
 * \brief In-Line subroutines of the <i>markRegion.hpp</i> file.
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

inline OpenHurricane::markRegion::markRegion() : id_(0), option_(regionOptions::NO_OPTIONs) {}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::markRegion::id() const noexcept {
    return id_;
}

hur_nodiscard inline bool OpenHurricane::markRegion::isOptionSet() const noexcept {
    return option_ != regionOptions::NO_OPTIONs;
}

hur_nodiscard inline OpenHurricane::markRegion::regionOptions
OpenHurricane::markRegion::option() const noexcept {
    return option_;
}

hur_nodiscard inline bool OpenHurricane::markRegion::isInside() const noexcept {
    return option_ == regionOptions::INSIDE;
}

hur_nodiscard inline bool OpenHurricane::markRegion::isOutside() const noexcept {
    return option_ == regionOptions::OUTSIDE;
}
