/*!
 * \file textContElement.cpp
 * \brief Main subroutines of text control elements.
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

#include "textContElement.hpp"

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::controlElement>
OpenHurricane::textContElement::clone(const controller &perentCont) const {
    return uniquePtr<controlElement>(new textContElement(*this));
}

void OpenHurricane::textContElement::writeToXML(std::stringstream &sstr, const integer ilayer) const {
    integer ispace = 4 * ilayer;
    std::string tmpSpac;
    for (integer i = 0; i < ispace; ++i) {
        tmpSpac += " ";
    }
    sstr << tmpSpac.c_str() << "<" << key_.c_str() << ">";
    sstr << "<![CDATA[" << text_.c_str() << "]]>";
    sstr << "</" << key_.c_str() << ">" << std::endl;
}
