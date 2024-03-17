/*!
 * \file cellArrays.hpp
 * \brief Headers of the geometry arrays based on cell.
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

#include "cellMesh.hpp"
#include "geometryArrays.hpp"
#include "matrixGeometryArrays.hpp"

namespace OpenHurricane {
    using cellComplexArray = complexGeometryArray<cellMesh>;
    using cellIntegerArray = integerGeometryArray<cellMesh>;
    using cellRealArray = realGeometryArray<cellMesh>;
    using cellRealArrayArray = realArrayGeometryArray<cellMesh>;
    using cellVector2DArray = vector2DGeometryArray<cellMesh>;
    using cellVectorArray = vectorGeometryArray<cellMesh>;
    using cellTensorArray = tensorGeometryArray<cellMesh>;
    using cellSymmTensorArray = symmTensorGeometryArray<cellMesh>;
    using cellRealSquareMatrixArray = realSquareMatrixGeometryArray<cellMesh>;
    using cellRealSymmMatrixArray = realSymmMatrixGeometryArray<cellMesh>;

    // integer=============================================================
    template <>
    void geometryArray<integer, cellMesh>::calcTimeSumPtr(const real &dt) const;

    // real=============================================================
    template <> void geometryArray<real, cellMesh>::getAveAndMaxRHS() const;

    template <> void geometryArray<real, cellMesh>::updateBoundary(integer layerI);

    template <> real geometryArray<real, cellMesh>::initialize();

    template <> void geometryArray<real, cellMesh>::writeOutput(fileOsstream &fos) const;

    template <> void geometryArray<real, cellMesh>::writeOutputByMaster(fileOsstream &fos) const;

    template <>
    void geometryArray<real, cellMesh>::writeOutput(fileOsstream &fos, const integer fzid) const;

    template <>
    void geometryArray<real, cellMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                          const integer fzid) const;

    template <>
    void geometryArray<real, cellMesh>::writeRelay(hdf5O &fos, const bool writeLast,
                                                   const bool writeToGroup) const;

    template <>
    void geometryArray<real, cellMesh>::readRelay(const hdf5I &fos, const bool readLast,
                                                  const bool readFromGroup);

    template <>
    void geometryArray<real, cellMesh>::interpolateRelay(const hdf5I &fos, const bool readLast,
                                                         const bool readFromGroup);    
    
    template <>
    hur_nodiscard inline realArray geometryArray<real, cellMesh>::realComponent(const int i) const {
        return Base::component(i);
    }

    // vector2D=============================================================
    template <> void geometryArray<vector2D, cellMesh>::writeOutput(fileOsstream &fos) const;

    template <>
    void geometryArray<vector2D, cellMesh>::writeOutputByMaster(fileOsstream &fos) const;

    template <>
    void geometryArray<vector2D, cellMesh>::writeOutput(fileOsstream &fos,
                                                                   const integer fzid) const;

    template <>
    void geometryArray<vector2D, cellMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                                         const integer fzid) const;

    template <>
    void geometryArray<vector2D, cellMesh>::writeRelay(hdf5O &fos, const bool writeLast,
                                                                  const bool writeToGroup) const;

    template <>
    void geometryArray<vector2D, cellMesh>::readRelay(const hdf5I &fos,
                                                                 const bool readLast,
                                                                 const bool readFromGroup);

    template <>
    void geometryArray<vector2D, cellMesh>::interpolateRelay(const hdf5I &fos,
                                                                        const bool readLast,
                                                                        const bool readFromGroup);
    template <>
    hur_nodiscard inline realArray
    geometryArray<vector2D, cellMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    // vector=============================================================
    template <> void geometryArray<vector, cellMesh>::getAveAndMaxRHS() const;

    template <> void geometryArray<vector, cellMesh>::updateBoundary(integer layerI);

    template <> vector geometryArray<vector, cellMesh>::initialize();

    template <> void geometryArray<vector, cellMesh>::writeOutput(fileOsstream &fos) const;

    template <>
    void geometryArray<vector, cellMesh>::writeOutputByMaster(fileOsstream &fos) const;

    template <>
    void geometryArray<vector, cellMesh>::writeOutput(fileOsstream &fos,
                                                                 const integer fzid) const;

    template <>
    void geometryArray<vector, cellMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                                       const integer fzid) const;

    template <>
    void geometryArray<vector, cellMesh>::writeRelay(hdf5O &fos, const bool writeLast,
                                                                const bool writeToGroup) const;

    template <>
    void geometryArray<vector, cellMesh>::readRelay(const hdf5I &fos,
                                                               const bool readLast,
                                                               const bool readFromGroup);

    template <>
    void geometryArray<vector, cellMesh>::interpolateRelay(const hdf5I &fos,
                                                                      const bool readLast,
                                                                      const bool readFromGroup);
    template <>
    hur_nodiscard inline realArray
    geometryArray<vector, cellMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    // tensor=============================================================
    template <> void geometryArray<tensor, cellMesh>::writeOutput(fileOsstream &fos) const;

    template <>
    void geometryArray<tensor, cellMesh>::writeOutputByMaster(fileOsstream &fos) const;

    template <>
    void geometryArray<tensor, cellMesh>::writeOutput(fileOsstream &fos,
                                                                 const integer fzid) const;

    template <>
    void geometryArray<tensor, cellMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                                       const integer fzid) const;

    template <>
    void geometryArray<tensor, cellMesh>::writeRelay(hdf5O &fos, const bool writeLast,
                                                                const bool writeToGroup) const;

    template <>
    void geometryArray<tensor, cellMesh>::readRelay(const hdf5I &fos,
                                                               const bool readLast,
                                                               const bool readFromGroup);

    template <>
    void geometryArray<tensor, cellMesh>::interpolateRelay(const hdf5I &fos,
                                                                      const bool readLast,
                                                                      const bool readFromGroup);
    template <>
    hur_nodiscard inline realArray
    geometryArray<tensor, cellMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    // symmTensor=============================================================

    template <> void geometryArray<symmTensor, cellMesh>::writeOutput(fileOsstream &fos) const;

    template <>
    void
    geometryArray<symmTensor, cellMesh>::writeOutputByMaster(fileOsstream &fos) const;

    template <>
    void geometryArray<symmTensor, cellMesh>::writeOutput(fileOsstream &fos,
                                                                     const integer fzid) const;

    template <>
    void
    geometryArray<symmTensor, cellMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                                      const integer fzid) const;

    template <>
    void geometryArray<symmTensor, cellMesh>::writeRelay(hdf5O &fos,
                                                                    const bool writeLast,
                                                                    const bool writeToGroup) const;

    template <>
    void geometryArray<symmTensor, cellMesh>::readRelay(const hdf5I &fos,
                                                                   const bool readLast,
                                                                   const bool readFromGroup);

    template <>
    void geometryArray<symmTensor, cellMesh>::interpolateRelay(const hdf5I &fos,
                                                                          const bool readLast,
                                                                          const bool readFromGroup);
    template <>
    hur_nodiscard inline realArray
    geometryArray<symmTensor, cellMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
} // namespace OpenHurricane
