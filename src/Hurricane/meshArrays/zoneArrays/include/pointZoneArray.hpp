/*!
 * \file pointZoneArray.hpp
 * \brief Headers of the point zone Array.
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
#include "zoneMesh.hpp"
#include "zoneArray.hpp"

namespace OpenHurricane {
    template <class Type> using pointZoneArray = zoneArray<Type, pointZone>;

    // integer cell zone Array
    using integerPointZoneArray = pointZoneArray<integer>;
    using integerArrayPointZoneArray = pointZoneArray<integerArray>;
    // real cell zone Array
    using realPointZoneArray = pointZoneArray<real>;
    using realArrayPointZoneArray = pointZoneArray<realArray>;
    // symmTensor cell zone Array
    using symmTensorPointZoneArray = pointZoneArray<symmTensor>;
    using symmTensorArrayPointZoneArray = pointZoneArray<symmTensorArray>;
    // tensor cell zone Array
    using tensorPointZoneArray = pointZoneArray<tensor>;
    using tensorArrayPointZoneArray = pointZoneArray<tensorArray>;
    // vector cell zone Array
    using vectorPointZoneArray = pointZoneArray<vector>;
    using vectorArrayPointZoneArray = pointZoneArray<vectorArray>;
    // vector2D cell zone Array
    using vector2DPointZoneArray = pointZoneArray<vector2D>;
    using vector2DArrayPointZoneArray = pointZoneArray<vector2DArray>;

    template <> void pointZoneArray<integer>::writeOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<real>::writeOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<symmTensor>::writeOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<tensor>::writeOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<vector>::writeOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<vector2D>::writeOutput(fileOsstream &fos) const;

    template <>
    void pointZoneArray<integer>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <> void pointZoneArray<real>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void pointZoneArray<symmTensor>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void pointZoneArray<tensor>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void pointZoneArray<vector>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void pointZoneArray<vector2D>::writeOutput(fileOsstream &fos, const integer fzid) const;

    template <> void pointZoneArray<integer>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<real>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<symmTensor>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<tensor>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<vector>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void pointZoneArray<vector2D>::writeMinMaxOutput(fileOsstream &fos) const;

    template <>
    void pointZoneArray<integer>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void pointZoneArray<real>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void pointZoneArray<symmTensor>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void pointZoneArray<tensor>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void pointZoneArray<vector>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void pointZoneArray<vector2D>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;

} //  namespace OpenHurricane
