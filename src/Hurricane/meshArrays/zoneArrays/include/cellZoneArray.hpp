/*!
 * \file cellZoneArray.hpp
 * \brief Headers of the cell zone Array.
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
    template <class Type> using cellZoneArray = zoneArray<Type, cellZone>;

    // integer cell zone Array
    using integerCellZoneArray = cellZoneArray<integer>;
    using integerArrayCellZoneArray = cellZoneArray<integerArray>;
    // real cell zone Array
    using realCellZoneArray = cellZoneArray<real>;
    using realArrayCellZoneArray = cellZoneArray<realArray>;
    // symmTensor cell zone Array
    using symmTensorCellZoneArray = cellZoneArray<symmTensor>;
    using symmTensorArrayCellZoneArray = cellZoneArray<symmTensorArray>;
    // tensor cell zone Array
    using tensorCellZoneArray = cellZoneArray<tensor>;
    using tensorArrayCellZoneArray = cellZoneArray<tensorArray>;
    // vector cell zone Array
    using vectorCellZoneArray = cellZoneArray<vector>;
    using vectorArrayCellZoneArray = cellZoneArray<vectorArray>;
    // vector2D cell zone Array
    using vector2DCellZoneArray = cellZoneArray<vector2D>;
    using vector2DArrayCellZoneArray = cellZoneArray<vector2DArray>;

    template <> void cellZoneArray<integer>::writeOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<real>::writeOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<symmTensor>::writeOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<tensor>::writeOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<vector>::writeOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<vector2D>::writeOutput(fileOsstream &fos) const;

    template <>
    void cellZoneArray<integer>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <> void cellZoneArray<real>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void cellZoneArray<symmTensor>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void cellZoneArray<tensor>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void cellZoneArray<vector>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void cellZoneArray<vector2D>::writeOutput(fileOsstream &fos, const integer fzid) const;

    template <> void cellZoneArray<integer>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<real>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<symmTensor>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<tensor>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<vector>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void cellZoneArray<vector2D>::writeMinMaxOutput(fileOsstream &fos) const;

    template <>
    void cellZoneArray<integer>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void cellZoneArray<real>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void cellZoneArray<symmTensor>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void cellZoneArray<tensor>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void cellZoneArray<vector>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void cellZoneArray<vector2D>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;

} //  namespace OpenHurricane
