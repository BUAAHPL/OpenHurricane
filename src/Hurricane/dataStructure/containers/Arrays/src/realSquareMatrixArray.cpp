/*!
 * \file realSquareMatrixArray.cpp
 * \brief Main subroutines for realSquareMatrix Array.
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

#include "realSquareMatrixArray.hpp"

namespace OpenHurricane {
    template <>
    void Array<realSquareMatrix>::replace(const int i,
                                          const Array<realSquareMatrix::elementType> &l) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j).element(i) = l[j];
        }
    }

    template <>
    void Array<realSquareMatrix>::replace(const int i, const realSquareMatrix::elementType &c) {
        setComponent(i, c);
    }

    template <>
    void Array<realSquareMatrix>::setComponent(const int d,
                                               const realSquareMatrix::elementType &c) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j).element(d) = c;
        }
    }

    template <> void Array<realSquareMatrix>::setComponent(const realSquareMatrix::elementType &c) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c;
        }
    }

    hur_nodiscard realSquareMatrixArray transpose(const realSquareMatrixArray &lm) {
        realSquareMatrixArray TT;
        TT.resize(lm.size());
        for (realSquareMatrixArray::size_type j = 0; j < lm.size(); j++) {
            TT[j] = lm[j].transpose();
        }
        return TT;
    }

    template <>
    void component(Array<typename Array<realSquareMatrix>::elementType> &comp,
                   const Array<realSquareMatrix> &lf, const int d) {
        for (Array<realSquareMatrix>::size_type i = 0; i < comp.size(); ++i) {
            comp[i] = lf[i].element(d);
        }
    }

    template <>
    Array<typename Array<realSquareMatrix>::elementType>
    Array<realSquareMatrix>::component(const int d) const {
        Array<realSquareMatrix::elementType> cm(this->size());
        for (size_type i = 0; i < this->size(); ++i) {
            cm[i] = this->operator[](i).element(d);
        }
        return cm;
    }

    template <> fileOsstream &Array<realSquareMatrix>::writeToStream(fileOsstream &fos) const {
        LFatal("Cannot write \"realSquareMatrixArray\" to fileOsstream %s", fos.name().c_str());
        return fos;
    }

    template <>
    void Array<realSymmMatrix>::replace(const int i, const Array<realSymmMatrix::elementType> &l) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j).element(i) = l[j];
        }
    }

    template <>
    void Array<realSymmMatrix>::replace(const int i, const realSymmMatrix::elementType &c) {
        setComponent(i, c);
    }

    template <>
    void Array<realSymmMatrix>::setComponent(const int d, const realSymmMatrix::elementType &c) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j).element(d) = c;
        }
    }

    template <> void Array<realSymmMatrix>::setComponent(const realSymmMatrix::elementType &c) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c;
        }
    }

    template <>
    void component(Array<typename Array<realSymmMatrix>::elementType> &comp,
                   const Array<realSymmMatrix> &lf, const int d) {
        for (Array<realSymmMatrix>::size_type i = 0; i < comp.size(); ++i) {
            comp[i] = lf[i].element(d);
        }
    }

    template <>
    Array<typename Array<realSymmMatrix>::elementType>
    Array<realSymmMatrix>::component(const int d) const {
        Array<realSymmMatrix::elementType> cm(this->size());
        for (size_type i = 0; i < this->size(); ++i) {
            cm[i] = this->operator[](i).element(d);
        }
        return cm;
    }

    template <> fileOsstream &Array<realSymmMatrix>::writeToStream(fileOsstream &fos) const {
        LFatal("Cannot write \"realSymmMatrixArray\" to fileOsstream %s", fos.name().c_str());
        return fos;
    }
} // namespace OpenHurricane
