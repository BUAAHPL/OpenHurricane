/*!
 * \file integerArray.hpp
 * \brief Headers of the integer Array.
 *        The subroutines and functions are in the <i>integerArray.cpp</i> file.
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

#pragma once

#include "Array.hpp"

namespace OpenHurricane {
    using integerArray = Array<integer>;
    using integerArrayArray = Array<integerArray>;
    using integerArrayArrayArray = Array<integerArrayArray>;

    template <> void Array<integer>::replace(const int i, const Array<integer> &l);

    template <> void Array<integer>::replace(const int i, const integer &c);

    template <> void Array<integer>::setComponent(const int d, const integer &c);

    template <> void Array<integer>::setComponent(const integer &c);

    template <>
    void component(Array<typename Array<integer>::value_type> &comp, const Array<integer> &lf,
                   const int d);

    template <>
    Array<typename Array<integer>::value_type> Array<integer>::component(const int d) const;

    template <> fileOsstream &Array<integer>::writeToStream(fileOsstream &fos) const;

    template <>
    fileOsstream &Array<integer>::writeToStream(fileOsstream &fos, const integer outSize) const;

    template <> fileOsstream &Array<integer>::writeToStreamWithFactor(fileOsstream &fos) const;

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &Array<integer>::writeToStreamWithFactor(fileOsstream &fos,
                                                          const size_type outSize) const;

    template <>
    void Array<integer>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                       const size_type outSize) const;

    template <>
    void Array<integer>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                               const size_type outSize) const;

    template <>
    void Array<integer>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const;

    hur_nodiscard inline integerArray operator*(const integerArray &f, const integer &t) {
        integerArray tf(f.size());
        for (integerArray::size_type i = 0; i < f.size(); ++i) {
            tf[i] = f[i] * t;
        }
        return tf;
    }

    hur_nodiscard inline integerArray operator*(integerArray &&f, const integer &t) {
        integerArray tf(std::move(f));
        for (List<integer>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] *= t;
        }
        return tf;
    }

    hur_nodiscard inline integerArray operator*(const integer t, const integerArray &f) {
        integerArray tf(f.size());
        for (integerArray::size_type i = 0; i < f.size(); ++i) {
            tf[i] = t * f[i];
        }
        return tf;
    }

    hur_nodiscard inline integerArray operator*(const integer t, integerArray &&f) noexcept {
        integerArray tf(std::move(f));
        for (List<integer>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = t * tf[i];
        }
        return tf;
    }

} // namespace OpenHurricane