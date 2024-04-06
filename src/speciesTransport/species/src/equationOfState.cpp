/*!
 * \file equationOfState.cpp
 * \brief The subroutines and functions of base class of equation of state.
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

#include "equationOfState.hpp"

namespace OpenHurricane {
    createClassNameStr(equationOfState, "equationOfState");
    createObjFty(equationOfState, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::equationOfState>
OpenHurricane::equationOfState::creator(const speciesList &_species, const controller &cont) {
    string eosType = cont.findWord(equationOfState::className_);
#ifdef HUR_DEBUG
    PLInfo("    Info: creating %s equation of state\n", eosType.c_str());
#endif // HUR_DEBUG

    defineInObjCreator(equationOfState, static_cast<std::string>(eosType), controller,
                       (_species, cont));
}
