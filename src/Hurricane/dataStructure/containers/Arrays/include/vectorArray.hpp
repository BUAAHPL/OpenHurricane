/*!
 * \file vectorArray.hpp
 * \brief Headers of the vector Array.
 *        The subroutines and functions are in the <i>vectorArray.cpp</i> file.
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
#include "real.hpp"
#include "realArray.hpp"

namespace OpenHurricane {
    using vectorArray = Array<vector>;
    using vectorArrayArray = Array<vectorArray>;
    using vector2DArray = Array<vector2D>;

    template <> fileOsstream &Array<vector>::writeToStream(fileOsstream &fos) const;

    template <>
    fileOsstream &Array<vector>::writeToStream(fileOsstream &fos, const size_type outSize) const;

    template <> fileOsstream &Array<vector>::writeToStreamWithFactor(fileOsstream &fos) const;

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &Array<vector>::writeToStreamWithFactor(fileOsstream &fos,
                                                         const size_type outSize) const;

    template <>
    void Array<vector>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                      const size_type outSize) const;

    template <>
    void Array<vector>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                              const size_type outSize) const;

    template <>
    void Array<vector>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const;

    template <>
    void Array<vector>::writeAveToPout(fileOsstream &fos, const Array<vector> &rhs,
                                       const Array<real> &cV, const size_type n,
                                       const size_type allN, const vector &rhs0, const bool calRhs0,
                                       const bool modifyRhs0) const;

    /*!\brief Inner product*/
    hur_nodiscard realArray operator*(const Array<vector> &f1, const Array<vector> &f2);

    /*!\brief Cross product*/
    hur_nodiscard vectorArray operator^(const Array<vector> &f1, const Array<vector> &f2);

    /*!\brief included angle of vectors
     * Return a real Array.
     */
    hur_nodiscard realArray cos(const Array<vector> &vf1, const Array<vector> &vf2);

    /*!\brief included angle of vectors
     * Return a real Array.
     */
    hur_nodiscard realArray cos(const Array<vector> &vf1, const vector &v2);

    /*!\brief included angle of vectors
     * Return a real Array.
     */
    hur_nodiscard realArray cos(const vector &v1, const Array<vector> &vf2);

    /*!\brief Distance between two points (vectors)
     * Return a real.
     */
    hur_nodiscard realArray dist(const Array<vector> &vf1, const Array<vector> &vf2);

    hur_nodiscard realArray div(const Array<vector> &vf1);

    template <> fileOsstream &Array<vector2D>::writeToStream(fileOsstream &fos) const;

    template <>
    fileOsstream &Array<vector2D>::writeToStream(fileOsstream &fos, const size_type outSize) const;

    template <> fileOsstream &Array<vector2D>::writeToStreamWithFactor(fileOsstream &fos) const;

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &Array<vector2D>::writeToStreamWithFactor(fileOsstream &fos,
                                                           const size_type outSize) const;

    template <>
    void Array<vector2D>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                        const size_type outSize) const;

    template <>
    void Array<vector2D>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                                const size_type outSize) const;

    template <>
    void Array<vector2D>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const;
} // namespace OpenHurricane