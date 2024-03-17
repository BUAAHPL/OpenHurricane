/*!
 * \file realArray.hpp
 * \brief Headers of the real array.
 *        The subroutines and functions are in the <i>realArray.cpp</i> file.
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
#include "dataStructure.hpp"

namespace OpenHurricane {
    using realArray = Array<real>;
    using realArrayArray = Array<realArray>;

    template <> void Array<real>::replace(const int i, const realArray &l);

    template <> void realArrayArray::replace(const int i, const realArray &l);

    template <> void Array<real>::replace(const int i, const real &c);

    template <> void realArrayArray::replace(const int i, const real &c);

    template <> void Array<real>::setComponent(const int d, const real &c);

    template <> void Array<real>::setComponent(const int d, const Array<real> &c);

    template <> void realArrayArray::setComponent(const int d, const real &c);

    template <> void Array<real>::setComponent(const real &c);

    template <> void realArrayArray::setComponent(const real &c);

    template <>
    void component(Array<typename Array<real>::elementType> &comp, const realArray &lf,
                   const int d);

    template <>
    void component(Array<typename realArrayArray::elementType> &comp, const Array<realArray> &lf,
                   const int d);

    template <> Array<typename Array<real>::elementType> Array<real>::component(const int d) const;

    template <>
    Array<typename realArrayArray::elementType> realArrayArray::component(const int d) const;

    template <> fileOsstream &Array<real>::writeToStream(fileOsstream &fos) const;

    template <> fileOsstream &realArrayArray::writeToStream(fileOsstream &fos) const;

    template <>
    fileOsstream &Array<real>::writeToStream(fileOsstream &fos, const size_type outSize) const;

    template <>
    fileOsstream &realArrayArray::writeToStream(fileOsstream &fos, const size_type outSize) const;

    template <> fileOsstream &Array<real>::writeToStreamWithFactor(fileOsstream &fos) const;

    /*!\brief Start from v_[0].*/
    template <>
    void Array<real>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const;

    template <>
    void Array<real>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                    const size_type outSize) const;

    template <>
    void Array<real>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                            const size_type outSize) const;

    template <>
    fileOsstream &Array<real>::writeToStreamWithFactor(fileOsstream &fos,
                                                       const size_type outSize) const;

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &realArrayArray::writeToStreamWithFactor(fileOsstream &fos,
                                                          const size_type outSize) const;

    template <> fileOsstream &realArrayArray::writeToStreamWithFactor(fileOsstream &fos) const;

    template <>
    void Array<real>::writeAveToPout(fileOsstream &fos, const Array<real> &rhs,
                                     const Array<real> &cV, const size_type n, const size_type allN,
                                     const real &rhs0, const bool calRhs0,
                                     const bool modifyRhs0) const;

    hur_nodiscard realArray operator*(const realArray &f1, const realArray &f2);
    hur_nodiscard realArray operator*(Array<real> &&f1, Array<real> &&f2) noexcept;
    hur_nodiscard realArray operator*(Array<real> &&tft, const real t) noexcept;
    hur_nodiscard realArray operator*(const real t, const realArray &f);
    hur_nodiscard realArray operator*(const realArray &f, const real &t);
    hur_nodiscard realArray operator*(const real t, Array<real> &&tft) noexcept;

    hur_nodiscard realArray operator/(const realArray &f1, realArray &&f2) noexcept;
    hur_nodiscard realArray operator/(realArray &&f1, realArray &&f2) noexcept;
    hur_nodiscard realArray operator/(const real t, const realArray &f);
    hur_nodiscard realArray operator/(const real t, realArray &&f) noexcept;

    hur_nodiscard realArray exp(const realArray &f);
    hur_nodiscard realArray exp(realArray &&f);

    hur_nodiscard realArray sin(const realArray &f);
    hur_nodiscard realArray sin(realArray &&f);

    hur_nodiscard realArray sinh(const realArray &f);
    hur_nodiscard realArray sinh(realArray &&f);

    hur_nodiscard realArray cos(const realArray &f);
    hur_nodiscard realArray cos(realArray &&f);

    hur_nodiscard realArray cosh(const realArray &f);
    hur_nodiscard realArray cosh(realArray &&f);

    hur_nodiscard realArray tan(const realArray &f);
    hur_nodiscard realArray tan(realArray &&f);

    hur_nodiscard realArray tanh(const realArray &f);
    hur_nodiscard realArray tanh(realArray &&f);

    hur_nodiscard realArray asin(const realArray &f);
    hur_nodiscard realArray asin(realArray &&f);

    hur_nodiscard realArray asinh(const realArray &f);
    hur_nodiscard realArray asinh(realArray &&f);

    hur_nodiscard realArray acos(const realArray &f);
    hur_nodiscard realArray acos(realArray &&f);

    hur_nodiscard realArray acosh(const realArray &f);
    hur_nodiscard realArray acosh(realArray &&f);

    hur_nodiscard realArray atan(const realArray &f);
    hur_nodiscard realArray atan(realArray &&f);

    hur_nodiscard realArray atan2(const realArray &y, const realArray &x);
    hur_nodiscard realArray atan2(realArray &&y, const realArray &x);
    hur_nodiscard realArray atan2(const realArray &y, realArray &&x);
    hur_nodiscard realArray atan2(realArray &&y, realArray &&x);

    hur_nodiscard realArray atanh(const realArray &f);
    hur_nodiscard realArray atanh(realArray &&f);

    hur_nodiscard realArray log(const realArray &f);
    hur_nodiscard realArray log(realArray &&f);

    hur_nodiscard realArray log10(const realArray &f);
    hur_nodiscard realArray log10(realArray &&f);

    hur_nodiscard realArray pow(const realArray &f, const real p);
    hur_nodiscard realArray pow(realArray &&f, const real p);

    hur_nodiscard realArray pow(const realArray &f, const realArray &p);
    hur_nodiscard realArray pow(realArray &&f, const realArray &p);
    hur_nodiscard realArray pow(const realArray &f, realArray &&p);
    hur_nodiscard realArray pow(realArray &&f, realArray &&p);

    hur_nodiscard realArray sqr(const realArray &f);
    hur_nodiscard realArray sqr(realArray &&f);

    hur_nodiscard realArray sqrt(const realArray &f);
    hur_nodiscard realArray sqrt(realArray &&f);

} // namespace OpenHurricane
