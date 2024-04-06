/*!
 * \file parameterContElement.cpp
 * \brief Main subroutines of parameter control elements.
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

#include "parameterContElement.hpp"

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::controlElement>
OpenHurricane::parameterContElement::clone(const controller &perentCont) const {
    return uniquePtr<controlElement>(new parameterContElement(*this));
}

hur_nodiscard OpenHurricane::IStringStream &
OpenHurricane::parameterContElement::ISStreamContEle() const {
    if (!hasSetIS_) {
        iss_.clear();
        iss_.str(parameterValue_);
    }
    const_cast<bool &>(hasSetIS_) = false;
    return iss_;
}

void OpenHurricane::parameterContElement::writeToXML(std::stringstream &sstr,
                                                 const integer ilayer) const {
    integer ispace = 4 * ilayer;
    std::string tmpSpac;
    for (integer i = 0; i < ispace; ++i) {
        tmpSpac += " ";
    }
    sstr << tmpSpac.c_str() << "<" << key_.c_str() << ">";
    sstr << parameterValue_.c_str();
    sstr << "</" << key_.c_str() << ">" << std::endl;
}
