
/*!
 * \file species.cpp
 * \brief The functions of the <i>species.hpp</i> file.
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

#include "species.hpp"

namespace OpenHurricane {
    species air("air", 28.97, (8.314e3));
}

OpenHurricane::species::species()
    : name_(), molecularWeight_(0.0), Ri_(constant::physicalConstant::Ru), elementList_() {}

OpenHurricane::species::species(const string &names, const real molWeight, const real Rr)
    : name_(names), molecularWeight_(molWeight), elementList_(0) {
    Ri_ = Rr / molecularWeight_;
}

OpenHurricane::species::species(const string &names, const real molWeight)
    : name_(names), molecularWeight_(molWeight), elementList_(0) {
    Ri_ = constant::physicalConstant::Ru / molecularWeight_;
}

OpenHurricane::species::species(const string &names, const real molWeight,
                            const speciesElementsList &eleList)
    : name_(names), molecularWeight_(molWeight), elementList_(eleList) {
    if (molecularWeight_ <= 0.0) {
        getMolWeight();
    }
    Ri_ = constant::physicalConstant::Ru / molecularWeight_;
}

OpenHurricane::species::species(const string &names, const speciesElementsList &eleList)
    : name_(names), molecularWeight_(0.0), elementList_(eleList) {
    getMolWeight();
    Ri_ = constant::physicalConstant::Ru / molecularWeight_;
}

OpenHurricane::species::species(const species &sp) {
    name_ = sp.name();
    molecularWeight_ = sp.W();
    elementList_.resize(sp.elementList().size());
    elementList_ = sp.elementList();
    Ri_ = sp.Ri_;
}

OpenHurricane::species::species(const species &sp, const string &name) {
    name_ = name;
    molecularWeight_ = sp.W();
    elementList_.resize(sp.elementList().size());
    elementList_ = sp.elementList();
    Ri_ = sp.Ri_;
}

OpenHurricane::species &OpenHurricane::species::operator=(const species &sp) {
    if (this != std::addressof(sp)) {
        name_ = sp.name();
        molecularWeight_ = sp.W();
        elementList_.resize(sp.elementList().size());
        elementList_ = sp.elementList();
        Ri_ = sp.Ri_;
    }
    return *this;
}

void OpenHurricane::species::getMolWeight() const {
    molecularWeight_ = 0.0;
    for (integer iEle = 0; iEle < elementList_.size(); iEle++) {
        molecularWeight_ += elementList_[iEle].nAtoms() * elementList_[iEle].atomicWeight();
    }
    Ri_ = constant::physicalConstant::Ru / molecularWeight_;
}

bool OpenHurricane::species::foundElement(const string &elementName, integer &nAtoms,
                                      real &elementWeight) const {
    for (integer iEle = 0; iEle < elementList_.size(); iEle++) {
        if (elementList_[iEle].name() == elementName) {
            nAtoms = elementList_[iEle].nAtoms();
            elementWeight = elementList_[iEle].atomicWeight();
            return true;
        }
    }
    return false;
}
