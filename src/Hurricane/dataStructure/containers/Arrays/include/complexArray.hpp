/*!
 * \file complexArray.hpp
 * \brief Headers of the complex Array.
 *        The subroutines and functions are in the <i>complexArray.cpp</i> file.
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
#include "complex.hpp"

namespace OpenHurricane {
    using complexArray = Array<complex>;

    template <>
    void Array<complex>::replace(const int i, const Array<feature<complex>::elementType> &l);

    template <> void Array<complex>::replace(const int i, const feature<complex>::elementType &c);

    template <>
    void Array<complex>::setComponent(const int d, const feature<complex>::elementType &c);

    template <> void Array<complex>::setComponent(const feature<complex>::elementType &c);

    template <>
    void component(Array<typename Array<complex>::elementType> &comp, const Array<complex> &lf,
                   const int d);

    template <>
    Array<typename Array<complex>::elementType> Array<complex>::component(const int d) const;

    template <> fileOsstream &Array<complex>::writeToStream(fileOsstream &fos) const;

} // namespace OpenHurricane