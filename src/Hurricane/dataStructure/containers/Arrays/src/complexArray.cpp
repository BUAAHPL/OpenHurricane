/*!
 * \file complexArray.cpp
 * \brief Main subroutines for complex Array.
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

#include "complexArray.hpp"

namespace OpenHurricane {
    template <>
    void OpenHurricane::Array<complex>::replace(const int i,
                                            const Array<feature<complex>::elementType> &l) {
        if (i == 0) {
            for (integer j = 0; j < this->size(); j++) {
                this->operator[](j).real(l[j]);
            }
        } else if (i == 1) {
            for (integer j = 0; j < this->size(); j++) {
                this->operator[](j).imag(l[j]);
            }
        } else {
            LFatal("Cannot access to the component of complex, for invalid component id: %d", i);
        }
    }

    template <>
    void OpenHurricane::Array<complex>::replace(const int i, const feature<complex>::elementType &c) {
        setComponent(i, c);
    }

    template <>
    void OpenHurricane::Array<complex>::setComponent(const int d,
                                                 const feature<complex>::elementType &c) {
        if (d == 0) {
            for (size_type j = 0; j < this->size(); j++) {
                this->operator[](j).real(c);
            }
        } else if (d == 1) {
            for (size_type j = 0; j < this->size(); j++) {
                this->operator[](j).imag(c);
            }
        } else {
            LFatal("Cannot access to the component of complex, for invalid component id: %d", d);
        }
    }

    template <>
    void OpenHurricane::Array<complex>::setComponent(const feature<complex>::elementType &c) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j).real(c);
            this->operator[](j).imag(c);
        }
    }

    template <>
    void component(Array<typename Array<complex>::elementType> &comp, const Array<complex> &lf,
                   const int d) {
        if (d == 0) {
            for (List<complex>::size_type i = 0; i < comp.size(); ++i) {
                comp[i] = lf[i].real();
            }
        } else if (d == 1) {
            for (List<complex>::size_type i = 0; i < comp.size(); ++i) {
                comp[i] = lf[i].imag();
            }
        } else {
            LFatal("Cannot access to the component of complex, for invalid component id: %d", d);
        }
    }

    template <>
    OpenHurricane::Array<typename OpenHurricane::Array<complex>::elementType>
    OpenHurricane::Array<complex>::component(const int d) const {
        Array<feature<complex>::elementType> cm(this->size());
        if (d == 0) {
            for (size_type i = 0; i < this->size(); ++i) {
                cm[i] = this->operator[](i).real();
            }
        } else if (d == 1) {
            for (size_type i = 0; i < this->size(); ++i) {
                cm[i] = this->operator[](i).imag();
            }
        } else {
            LFatal("Cannot access to the component of complex, for invalid component id: %d", d);
        }
        return cm;
    }

    template <>
    OpenHurricane::fileOsstream &OpenHurricane::Array<complex>::writeToStream(fileOsstream &fos) const {
        fos.setRealPrecision();

        if (fos.format() == IOsstream::ASCII_FORMAT) {
            for (size_type i = 0; i < this->size(); ++i) {
                fos.os() << this->operator[](i).real();
            }

            for (size_type i = 0; i < this->size(); ++i) {
                fos.os() << this->operator[](i).imag();
            }
        } else {
            if (this->size()) {
                Array<feature<complex>::elementType> cmRe = this->component(0);
                fos.write(reinterpret_cast<const char *>(&cmRe[0]), cmRe.byteSize());
                cmRe.clear();

                Array<feature<complex>::elementType> cmIm = this->component(1);
                fos.write(reinterpret_cast<const char *>(&cmIm[0]), cmIm.byteSize());
                cmIm.clear();
            }
        }

        fos.unsetRealPrecision();
        return fos;
    }
} // namespace OpenHurricane
