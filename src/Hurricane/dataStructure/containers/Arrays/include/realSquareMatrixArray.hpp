/*!
 * \file realSquareMatrixArray.hpp
 * \brief Headers of the realSquareMatrix Array.
 *        The subroutines and functions are in the <i>realSquareMatrixArray.cpp</i> file.
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

#pragma once

#include "Array.hpp"
#include "matrices.hpp"

namespace OpenHurricane {
    using realSquareMatrixArray = Array<realSquareMatrix>;
    using realSymmMatrixArray = Array<realSymmMatrix>;

    template <>
    void Array<realSquareMatrix>::replace(const int i,
                                          const Array<realSquareMatrix::elementType> &l);

    template <>
    void Array<realSquareMatrix>::replace(const int i, const realSquareMatrix::elementType &c);

    template <>
    void Array<realSquareMatrix>::setComponent(const int d, const realSquareMatrix::elementType &c);

    template <> void Array<realSquareMatrix>::setComponent(const realSquareMatrix::elementType &c);

    template <>
    void component(Array<typename Array<realSquareMatrix>::elementType> &comp,
                   const Array<realSquareMatrix> &lf, const int d);

    template <>
    Array<typename Array<realSquareMatrix>::elementType>
    Array<realSquareMatrix>::component(const int d) const;

    hur_nodiscard realSquareMatrixArray transpose(const realSquareMatrixArray &lm);

    template <> fileOsstream &Array<realSquareMatrix>::writeToStream(fileOsstream &fos) const;

    
    template <>
    void Array<realSymmMatrix>::replace(const int i, const Array<realSymmMatrix::elementType> &l);

    template <>
    void Array<realSymmMatrix>::replace(const int i, const realSymmMatrix::elementType &c);

    template <>
    void Array<realSymmMatrix>::setComponent(const int d, const realSymmMatrix::elementType &c);

    template <> void Array<realSymmMatrix>::setComponent(const realSymmMatrix::elementType &c);

    template <>
    void component(Array<typename Array<realSymmMatrix>::elementType> &comp,
                   const Array<realSymmMatrix> &lf, const int d);

    template <>
    Array<typename Array<realSymmMatrix>::elementType>
    Array<realSymmMatrix>::component(const int d) const;

    template <> fileOsstream &Array<realSymmMatrix>::writeToStream(fileOsstream &fos) const;

} // namespace OpenHurricane