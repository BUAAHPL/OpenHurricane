/*!
 * \file sourceTerms.cpp
 * \brief Main subroutines for sourceTerms.
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

#include "sourceTerms.hpp"

void OpenHurricane::sourceTerms::setSourceTermList(const controller &cont) {
    auto sourceName = cont.findTextStr("sourceList");
    if (sourceName.size() == 0) {
#ifdef HUR_DEBUG
        checkWarningStr(("The size of \"sourceList\" is zero in " + cont.name()));
#endif // HUR_DEBUG
        return;
    }
    sorT_.resize(sourceName.size());
    for (integer i = 0; i < sourceName.size(); ++i) {
        Pout << "    Info: setting source term: " << sourceName[i] << std::endl;
        const auto &sCont = cont.subController(sourceName[i]);

        const auto sourceType = sCont.findWord("sourceType");

        sorT_.set(i, sourceTerm::creator(flows_, iter_, sCont, sourceType).release());
    }
}

OpenHurricane::sourceTerms::sourceTerms(const flowModel &flows, const iteration &iter,
                                    const controller &cont)
    : flows_(flows), iter_(iter), sorT_() {
    setSourceTermList(cont);
}