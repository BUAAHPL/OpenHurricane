/*!
 * \file speciesList.cpp
 * \brief Main subroutines for the specie list.
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

#include "speciesList.hpp"

hur_nodiscard OpenHurricane::string OpenHurricane::speciesList::nameList() const {
    // For example: O2,H2,N2

    string nl;
    for (integer i = 0; i < this->size(); i++) {
        nl += this->operator[](i).name();
        if (i != this->size() - 1) {
            nl += ",";
        }
    }
    return nl;
}

hur_nodiscard const OpenHurricane::stringList &
OpenHurricane::speciesList::elementsNameList() const noexcept {
    return elementsNameList_;
}

hur_nodiscard OpenHurricane::stringList &OpenHurricane::speciesList::elementsNameList() noexcept {
    return elementsNameList_;
}

hur_nodiscard bool OpenHurricane::speciesList::contains(const string &_spName) const {
    for (integer i = 0; i < this->size(); ++i) {
        if (this->operator[](i).name() == _spName) {
            return true;
        }
    }
    return false;
}

bool OpenHurricane::speciesList::contains(const string &_spName, integer &_id) const {
    for (integer i = 0; i < this->size(); ++i) {
        if (this->operator[](i).name() == _spName) {
            _id = i;
            return true;
        }
    }
    _id = -1;
    return false;
}

hur_nodiscard OpenHurricane::integer OpenHurricane::speciesList::index(const string &_spName) const {
    integer i;
    if (!contains(_spName, i)) {
        LFatal("Specie: \"%s\" does not exit in the table of species:\n {\n%s\n}", _spName.c_str(),
               nameList().c_str());
    }
    return i;
}

hur_nodiscard OpenHurricane::real OpenHurricane::speciesList::Rm(const PtrList<cellRealArray> &yi,
                                                         const integer cellI) const noexcept {
    return constant::physicalConstant::Ru / MWbyYi(yi, cellI);
}

hur_nodiscard OpenHurricane::real OpenHurricane::speciesList::MWbyXi(const realArray &xi) const noexcept {
#ifdef HUR_DEBUG

    if (xi.size() != (*this).size()) {
        LFatal("The size is not equal. Size of xi: %d != specie table size: %s", xi.size(), size());
    }

#endif // HUR_DEBUG
    if (size() == 1) {
        return this->operator[](0).W();
    }
    real Wmean = Zero;
    for (integer i = 0; i < size(); ++i) {
        Wmean += xi[i] * this->operator[](i).W();
    }
    return Wmean;
}

hur_nodiscard OpenHurricane::real OpenHurricane::speciesList::MWbyYi(const realArray &yi) const noexcept {
    if (size() == 1) {
        return this->operator[](0).W();
    }
    real Wmean = Zero;
    for (integer i = 0; i < size(); ++i) {
        Wmean += yi[i] / this->operator[](i).W();
    }
    Wmean = 1.0 / max(Wmean, veryTiny);
    return Wmean;
}

hur_nodiscard OpenHurricane::real OpenHurricane::speciesList::MWbyYi(const PtrList<cellRealArray> &yi,
                                                             const integer cellI) const noexcept {
    if (size() == 1) {
        return this->operator[](0).W();
    }
    real Wmean = Zero;
    for (integer i = 0; i < size(); ++i) {
        Wmean += yi[i][cellI] / this->operator[](i).W();
    }
    Wmean = 1 / max(Wmean, veryTiny);
    return Wmean;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::speciesList::Yi2Xi(const realArray &yi) const {
    if (size() == 1) {
        realArray xi(size(), Zero);
        xi = 1.0;
        return xi;
    }
    realArray xi(size(), Zero);
    Yi2Xi(yi, xi);
    return xi;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::speciesList::Yi2Xi(const PtrList<cellRealArray> &yi,
                                                                 const integer cellI) const {
    if (size() == 1) {
        realArray xi(size(), Zero);
        xi = 1.0;
        return xi;
    }
    realArray xi(size(), Zero);
    real Wmean = MWbyYi(yi, cellI);
    real xit = Zero;
    for (integer i = 0; i < size(); ++i) {
        xi[i] = yi[i][cellI] * Wmean / this->operator[](i).W();
        xi[i] = max(real(0.0), min(real(1.0), xi[i]));
        xit += xi[i];
    }
    if (xit != real(1.0)) {
        xi /= xit;
    }
    return xi;
}

void OpenHurricane::speciesList::Yi2Xi(const PtrList<cellRealArray> &yi, const integer cellI,
                                   realArray &xi) const {
    if (size() == 1) {
        xi = 1.0;
        return;
    }
    real Wmean = MWbyYi(yi, cellI);
    real xit = Zero;
    for (integer i = 0; i < size(); ++i) {
        xi[i] = yi[i][cellI] * Wmean / this->operator[](i).W();
        xi[i] = max(real(0.0), min(real(1.0), xi[i]));
        xit += xi[i];
    }
    if (xit != real(1.0)) {
        xi /= xit;
    }
}

void OpenHurricane::speciesList::Yi2Xi(const PtrList<cellRealArray> &yi, PtrList<cellRealArray> &xi,
                                   const bool isOnlyInternal) const {
    if (size() == 1) {
        xi[0] = 1.0;
        return;
    }

    integer fieldSize = yi[0].mesh().nCells();
    if (!isOnlyInternal) {
        fieldSize = yi[0].mesh().nTotalCells();
    }
    for (integer cellI = 0; cellI < fieldSize; ++cellI) {
        real Wmean = MWbyYi(yi, cellI);
        real xit = Zero;
        for (integer i = 0; i < size(); ++i) {
            xi[i][cellI] = yi[i][cellI] * Wmean / this->operator[](i).W();
            xi[i][cellI] = max(real(0.0), min(real(1.0), xi[i][cellI]));
            xit += xi[i][cellI];
        }
        if (xit != real(1.0)) {
            for (integer i = 0; i < size(); ++i) {
                xi[i][cellI] /= xit;
            }
        }
    }
}

void OpenHurricane::speciesList::Yi2Xi(const realArray &yi, realArray &xi) const {
   
    if (size() == 1) {
        xi = 1.0;
        return;
    }
    real Wmean = MWbyYi(yi);
    real xit = Zero;
    for (integer i = 0; i < size(); ++i) {
        xi[i] = yi[i] * Wmean / this->operator[](i).W();
        xi[i] = max(real(0.0), min(real(1.0), xi[i]));
        xit += xi[i];
    }
    if (xit != real(1.0)) {
        xi /= xit;
    }
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::speciesList::Xi2Yi(const realArray &xi) const {
    if (size() == 1) {
        realArray yi(size(), Zero);
        yi = 1.0;
        return yi;
    }
    realArray yi(xi.size(), Zero);
    Xi2Yi(xi, yi);
    return yi;
}

void OpenHurricane::speciesList::Xi2Yi(const realArray &xi, realArray &yi) const {
#ifdef HUR_DEBUG
    if (xi.size() != (*this).size()) {
        LFatal("The size is not equal. Size of xi: %d != specie table size: %s", xi.size(), size());
    }

#endif // HUR_DEBUG
    if (size() == 1) {
        yi = 1.0;
        return;
    }
    real Wmean = MWbyXi(xi);
    real yit = Zero;
    for (integer i = 0; i < size(); ++i) {
        yi[i] = xi[i] * this->operator[](i).W() / Wmean;
        yi[i] = max(real(0.0), min(real(1.0), yi[i]));
        yit += yi[i];
    }
    if (yit != real(1.0)) {
        yi /= yit;
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::speciesList::Zj(const string &elementName,
                                                         const realArray &yi) const {
    real zj = 0;
    real elementWeight = 0;
    for (integer isp = 0; isp < this->size(); ++isp) {
        integer nAtom = 0;
        if (this->operator[](isp).foundElement(elementName, nAtom, elementWeight)) {
            zj += nAtom * yi[isp] / this->operator[](isp).W();
        }
    }
    return zj * elementWeight;
}

hur_nodiscard OpenHurricane::real OpenHurricane::speciesList::Zj(const string &elementName,
                                                         const PtrList<cellRealArray> &yi,
                                                         const integer cellI) const {
    real zj = 0;
    real elementWeight = 0;
    for (integer isp = 0; isp < this->size(); ++isp) {
        integer nAtom = 0;
        if (this->operator[](isp).foundElement(elementName, nAtom, elementWeight)) {
            zj += nAtom * yi[isp][cellI] / this->operator[](isp).W();
        }
    }
    return zj * elementWeight;
}
