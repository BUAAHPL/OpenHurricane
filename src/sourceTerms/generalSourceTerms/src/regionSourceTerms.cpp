/*!
 * \file regionSourceTerms.cpp
 * \brief Main subroutines for regionSourceTerms.
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

#include "regionSourceTerms.hpp"
 namespace OpenHurricane{
	createClassNameStr(regionSourceTerms,"regionSourceTerms");
}

OpenHurricane::regionSourceTerms::regionSourceTerms(const flowModel &flows, const iteration &iter,
                                                const controller &cont)
    : sourceTerm(flows, iter, cont) {
    auto regionName = cont.findTextStr("regions");
    if (regionName.size() == 0) {
        errorAbortStr(("The size of \"regions\" is zero in " + cont.name()));
    }
    regionList_.resize(regionName.size());
    for (integer ire = 0; ire < regionName.size(); ++ire) {
        Pout << "        Adding mark region: " << regionName[ire] << std::endl;
        regionList_[ire] = mesh().region(regionName[ire]).regionCellId(mesh());
    }
}
