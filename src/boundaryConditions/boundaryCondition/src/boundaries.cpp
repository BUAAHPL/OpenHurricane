/*!
 * \file boundaries.cpp
 * \brief Main subroutines for boundaries.
 * \author Yang Hongzhen
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
#include "boundaries.hpp"
#include "parsingDirection.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

template <>
const std::string
    OpenHurricane::boundary<OpenHurricane::real, OpenHurricane::cellMesh>::className_ =
        "boundary";

template <>
const std::string
    OpenHurricane::boundary<OpenHurricane::vector, OpenHurricane::cellMesh>::className_ =
    "boundary";

namespace OpenHurricane {
    createObjFtyTmpl(realBoundary, controller);
    createObjFtyTmpl(vectorBoundary, controller);

    template <>
    void realBoundary::readValue(real &value, const std::string &name,
                                 const controller &cont) const {
        if (cont.found(name)) {
            value = cont.findType<real>(name, value);
        } else {
            const thermoList &thTable = varArray_.tb().findObject<thermoList>("thermoList");
            for (integer i = 0; i < thTable.species().size(); i++) {
                if (thTable.species().name(i) == name) {
                    value = Zero;
                    return;
                }
            }
            LFatal("Boundary %s does not specify %s value.\n Please check.",
                   boundaryZone_.name().c_str(), name.c_str());
        }
    }

    template <>
    void vectorBoundary::readValue(vector &value, const std::string &name,
                                   const controller &cont) const {
        if (cont.found(name)) {
            real valueMag;
            valueMag = cont.findType<real>(name, 0.0);
            vector direct;
            parsingDirection::getDirection(const_cast<controller &>(cont), this->varField().mesh(),
                                           this->boundaryfZ(), direct);
            value = direct * valueMag;

        } else {
            LFatal("Boundary %s does not specify %s value.\n Please check.",
                   boundaryZone_.name().c_str(), name.c_str());
        }
    }

    template <>
    void realBoundary::readVector(vector &value, const std::string &name,
                                  const controller &cont) const {
        if (cont.found(name)) {
            real valueMag;
            valueMag = cont.findType<real>(name, 0.0);
            vector direct;
            parsingDirection::getDirection(const_cast<controller &>(cont), mesh(),
                                           this->boundaryfZ(), direct);
            value = direct * valueMag;

        } else {
            LFatal("Boundary %s does not specify %s value.\n Please check.",
                   boundaryZone_.name().c_str(), name.c_str());
        }
    }
    template <>
    void vectorBoundary::readVector(vector &value, const std::string &name,
                                    const controller &cont) const {
        if (cont.found(name)) {
            real valueMag;
            valueMag = cont.findType<real>(name, 0.0);
            vector direct;
            parsingDirection::getDirection(const_cast<controller &>(cont), mesh(),
                                           this->boundaryfZ(), direct);
            value = direct * valueMag;

        } else {
            LFatal("Boundary %s does not specify %s value.\n Please check.",
                   boundaryZone_.name().c_str(), name.c_str());
        }
    }
} // namespace OpenHurricane