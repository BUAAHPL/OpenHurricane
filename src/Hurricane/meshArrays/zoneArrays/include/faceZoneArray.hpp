/*!
 * \file faceZoneArray.hpp
 * \brief Headers of the face zone Array.
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
    template <class Type> using faceZoneArray = zoneArray<Type, faceZone>;

    // integer face zone Array
    using integerFaceZoneArray = faceZoneArray<integer>;
    using integerArrayFaceZoneArray = faceZoneArray<integerArray>;
    // real face zone Array
    using realFaceZoneArray = faceZoneArray<real>;
    using realArrayFaceZoneArray = faceZoneArray<realArray>;
    // symmTensor face zone Array
    using symmTensorFaceZoneArray = faceZoneArray<symmTensor>;
    using symmTensorArrayFaceZoneArray = faceZoneArray<symmTensorArray>;
    // tensor face zone Array
    using tensorFaceZoneArray = faceZoneArray<tensor>;
    using tensorArrayFaceZoneArray = faceZoneArray<tensorArray>;
    // vector face zone Array
    using vectorFaceZoneArray = faceZoneArray<vector>;
    using vectorArrayFaceZoneArray = faceZoneArray<vectorArray>;
    // vector2D face zone Array
    using vector2DFaceZoneArray = faceZoneArray<vector2D>;
    using vector2DArrayFaceZoneArray = faceZoneArray<vector2DArray>;

    template <> void faceZoneArray<integer>::writeOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<real>::writeOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<symmTensor>::writeOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<tensor>::writeOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<vector>::writeOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<vector2D>::writeOutput(fileOsstream &fos) const;

    template <>
    void faceZoneArray<integer>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <> void faceZoneArray<real>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void faceZoneArray<symmTensor>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void faceZoneArray<tensor>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void faceZoneArray<vector>::writeOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void faceZoneArray<vector2D>::writeOutput(fileOsstream &fos, const integer fzid) const;

    template <> void faceZoneArray<integer>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<real>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<symmTensor>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<tensor>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<vector>::writeMinMaxOutput(fileOsstream &fos) const;
    template <> void faceZoneArray<vector2D>::writeMinMaxOutput(fileOsstream &fos) const;

    template <>
    void faceZoneArray<integer>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void faceZoneArray<real>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void faceZoneArray<symmTensor>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void faceZoneArray<tensor>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void faceZoneArray<vector>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;
    template <>
    void faceZoneArray<vector2D>::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;

} //  namespace OpenHurricane
