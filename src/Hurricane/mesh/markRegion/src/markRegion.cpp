/*!
 * \file markRegion.cpp
 * \brief Main subroutines for mark region of mesh.
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

#include "markRegion.hpp"
#include "runtimeMesh.hpp"

namespace OpenHurricane {
    createClassNameStr(markRegion, "markRegion");
}
namespace OpenHurricane {
    createObjFty(markRegion, controller);
}

OpenHurricane::markRegion::markRegion(const controller &cont)
    : id_(cont.findOrDefault<integer>("id", 0)), option_(regionOptions::NO_OPTIONs) {
    if (cont.found("option")) {
        const auto w = cont.findWord("option");
        if (w == "inside") {
            option_ = regionOptions::INSIDE;
        } else if (w == "outside") {
            option_ = regionOptions::OUTSIDE;
        } else {
            LFatal("Unknown option: %s", w.c_str());
        }
    }
}

OpenHurricane::uniquePtr<OpenHurricane::markRegion> OpenHurricane::markRegion::creator(const controller &cont) {
    string regineType = cont.findWord("shape");
#ifdef HUR_DEBUG
    PLInfo("    Info: creating %s shape regin\n", regineType.c_str());
#endif // HUR_DEBUG

    defineInObjCreator(markRegion, static_cast<std::string>(regineType), controller, (cont));
}

hur_nodiscard OpenHurricane::integerList
OpenHurricane::markRegion::regionCellId(const runtimeMesh &mesh) const {
    checkWarning("This function is not implemented");
    return integerList();
}