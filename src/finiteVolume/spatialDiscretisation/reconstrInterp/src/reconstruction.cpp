/*!
 * \file reconstruction.cpp
 * \brief Main subroutines for reconstruction.
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

#include "reconstruction.hpp"

namespace OpenHurricane {
    createClassNameStr(reconstruction, "reconstruction");
    createObjFty(reconstruction, controller);
} // namespace OpenHurricane

OpenHurricane::reconstruction::reconstruction() : gradPtr_(nullptr) {}

OpenHurricane::reconstruction::reconstruction(const controller &cont) : gradPtr_(nullptr) {
    gradPtr_ = gradient::creator(cont);
}

OpenHurricane::uniquePtr<OpenHurricane::reconstruction>
OpenHurricane::reconstruction::creator(const controller &cont) {
    string reconstructionType = cont.findWord(reconstruction::className_);
    defineInObjCreator(reconstruction, reconstructionType, controller, (cont));
}

OpenHurricane::reconstruction::~reconstruction() noexcept {}
