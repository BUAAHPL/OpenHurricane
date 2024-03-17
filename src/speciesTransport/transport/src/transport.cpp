/*!
 * \file transport.cpp
 * \brief Main subroutines for transport properties.
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
#include "transport.hpp"

namespace OpenHurricane {
    createClassNameStr(transport, "transport");
    createObjFty(transport, controller);
} // namespace OpenHurricane

OpenHurricane::transport::transport(const speciesList &sp, const integer index, const real Prl)
    : species_(sp), index_(index), Pr_(Prl) {}

OpenHurricane::transport::transport(const speciesList &sp, const integer index, const controller &cont)
    : species_(sp), index_(index), Pr_(cont.findOrDefault<real>("Prl", 0.72, true)) {}

OpenHurricane::transport::transport(const transport &tra)
    : species_(tra.species_), index_(tra.index_), Pr_(tra.Pr_) {}

OpenHurricane::transport::transport(const transport &tra, const speciesList &sp)
    : species_(sp), index_(tra.index_), Pr_(tra.Pr_) {}

OpenHurricane::uniquePtr<OpenHurricane::transport>
OpenHurricane::transport::creator(const speciesList &sp, const integer index, const controller &cont) {
    string transportType = cont.findWord("type");
    if (cont.found(transportType)) {
        defineInObjCreator(transport, transportType, controller,
                           (sp, index, cont.subController(transportType)));
    } else {
        defineInObjCreator(transport, transportType, controller, (sp, index, cont));
    }
}

OpenHurricane::transport &OpenHurricane::transport::operator=(const transport &t) {
    if (this != std::addressof(t)) {
        index_ = t.index_;
        Pr_ = t.Pr_;
    }
    return *this;
}