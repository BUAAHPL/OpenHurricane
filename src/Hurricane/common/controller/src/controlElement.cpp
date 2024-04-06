/*!
 * \file controlElement.cpp
 * \brief Main subroutines of control elements.
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

#include "controlElement.hpp"
#include "controller.hpp"
#include "realArray.hpp"

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::controlElement>
OpenHurricane::controlElement::clone(const controller &perentCont) const {
    return uniquePtr<controlElement>(new controlElement(*this));
}

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::controlElement>
OpenHurricane::controlElement::clone() const {
    return this->clone(controller::null);
}

hur_nodiscard const OpenHurricane::string &OpenHurricane::controlElement::wordContEle() const {
    return string::null;
}

hur_nodiscard OpenHurricane::string &OpenHurricane::controlElement::wordContEle() {
    return const_cast<string &>(string::null);
}

hur_nodiscard const std::string &OpenHurricane::controlElement::parameterContEle() const {
    return nullString;
}

hur_nodiscard std::string &OpenHurricane::controlElement::parameterContEle() {
    return const_cast<std::string &>(nullString);
}

hur_nodiscard OpenHurricane::IStringStream &OpenHurricane::controlElement::ISStreamContEle()const {
    return const_cast<IStringStream &>(nullIStringStream);
}

hur_nodiscard const std::string &OpenHurricane::controlElement::textContEle() const {
    return nullString;
}

hur_nodiscard std::string &OpenHurricane::controlElement::textContEle() {
    return const_cast<std::string &>(nullString);
}

hur_nodiscard const OpenHurricane::controller &OpenHurricane::controlElement::contContEle() const {
    return controller::null;
}

hur_nodiscard OpenHurricane::controller &OpenHurricane::controlElement::contContEle() {
    return const_cast<controller &>(controller::null);
}

hur_nodiscard const OpenHurricane::realArray &OpenHurricane::controlElement::rArrayContEle() const {
    return realArray::nullObject();
}

hur_nodiscard OpenHurricane::realArray &OpenHurricane::controlElement::rArrayContEle() {
    return const_cast<realArray &>(realArray::nullObject());
}
