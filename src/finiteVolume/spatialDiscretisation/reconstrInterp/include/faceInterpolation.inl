/*!
 * \file faceInterpolation.inl
 * \brief In-Line subroutines of the <i>faceInterpolation.hpp</i> file.
 * \author Rao Sihang
 * \version V2.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of OpenHurricane project (Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation)
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
#include "faceInterpolation.hpp"

namespace OpenHurricane {
    namespace fv {
        /*!\brief Interpolation field onto face.*/
        template <class Type>
        geometryArray<Type, faceMesh> interpolate(const geometryArray<Type, cellMesh> &cf) {
            const runtimeMesh &mesh = cf.mesh();
            const realArray &faceWeight = mesh.faceWgt();
            const faceList &fl = mesh.faces();

            geometryArray<Type, faceMesh> ff(object("interpolation(" + cf.name() + ")", mesh,
                                                    object::NOT_WRITE, object::TEMPORARY),
                                             mesh);

            const integer nFaces = mesh.nFaces();
            for (integer fi = 0; fi < nFaces; ++fi) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();

                ff[fi] = faceWeight[fi] * cf[cl] + (real(1.0) - faceWeight[fi]) * cf[cr];
            }

            return ff;
        }

        template <class Type>
        static void interpolate(const geometryArray<Type, cellMesh> &cf,
                                geometryArray<Type, faceMesh> &ff) {
            const runtimeMesh &mesh = cf.mesh();
            const realArray &faceWeight = mesh.faceWgt();
            const faceList &fl = mesh.faces();

            const integer nFaces = mesh.nFaces();
            if (argParse::isSkewCorrect()) {
                const auto &fWSkew = mesh.faceSkewnessWeight();
                const auto &isFSkew = mesh.isSkewFace();
                const auto &fISFC = mesh.fIntSectFCentre();
                if (!cf.hasGradArray()) {
                    auto grdCf = fv::grad(cf, string("leastSquareGrad"));
                    for (integer fi = 0; fi < nFaces; ++fi) {
                        const auto &cl = fl[fi].leftCell();
                        const auto &cr = fl[fi].rightCell();

                        ff[fi] = fWSkew[fi] * cf[cl] + (real(1.0) - fWSkew[fi]) * cf[cr];
                        if (isFSkew[fi]) {
                            const auto grdff =
                                fWSkew[fi] * grdCf[cl] + (real(1.0) - fWSkew[fi]) * grdCf[cr];
                            ff[fi] += fISFC[fi] * grdff;
                        }
                    }
                } else {
                    const auto &grdCf = cf.grad();
                    for (integer fi = 0; fi < nFaces; ++fi) {
                        const auto &cl = fl[fi].leftCell();
                        const auto &cr = fl[fi].rightCell();

                        ff[fi] = fWSkew[fi] * cf[cl] + (real(1.0) - fWSkew[fi]) * cf[cr];
                        if (isFSkew[fi]) {
                            const auto grdff =
                                fWSkew[fi] * grdCf[cl] + (real(1.0) - fWSkew[fi]) * grdCf[cr];
                            ff[fi] += fISFC[fi] * grdff;
                        }
                    }
                }
            } else {
                for (integer fi = 0; fi < nFaces; ++fi) {
                    const auto &cl = fl[fi].leftCell();
                    const auto &cr = fl[fi].rightCell();

                    ff[fi] = faceWeight[fi] * cf[cl] + (real(1.0) - faceWeight[fi]) * cf[cr];
                }
            }
        }

        /*!\brief Interpolation field onto face zone by given zone id.*/
        template <class Type>
        Array<Type> interpolate(const geometryArray<Type, cellMesh> &cf, const integer faceZoneId) {
            const runtimeMesh &mesh = cf.mesh();
            const realArray &faceWeight = mesh.faceWgt();
            const faceList &fl = mesh.faces();
            const faceZone &fz = mesh.faceZones()[faceZoneId];
            Array<Type> ff(fz.size());
            integer ffi = 0;
            if (argParse::isSkewCorrect() && cf.hasGradArray()) {
                const auto &fWSkew = mesh.faceSkewnessWeight();
                const auto &isFSkew = mesh.isSkewFace();
                const auto &fISFC = mesh.fIntSectFCentre();
                const auto &grdCf = cf.grad();
                for (integer fi = fz.firstIndex(); fi <= fz.lastIndex(); fi++) {
                    const auto &cl = fl[fi].leftCell();
                    const auto &cr = fl[fi].rightCell();
                    ff[ffi++] = fWSkew[fi] * cf[cl] + (real(1.0) - fWSkew[fi]) * cf[cr];
                    if (isFSkew[fi]) {
                        const auto grdff =
                            fWSkew[fi] * grdCf[cl] + (real(1.0) - fWSkew[fi]) * grdCf[cr];
                        ff[fi] += fISFC[fi] * grdff;
                    }
                }
            } else {
                for (integer fi = fz.firstIndex(); fi <= fz.lastIndex(); fi++) {
                    const auto &cl = fl[fi].leftCell();
                    const auto &cr = fl[fi].rightCell();
                    ff[ffi++] = faceWeight[fi] * cf[cl] + (real(1.0) - faceWeight[fi]) * cf[cr];
                }
            }
            return ff;
        }

        template <class Type>
        Array<Type> interpolateNoCorr(const geometryArray<Type, cellMesh> &cf,
                                      const integer faceZoneId) {
            const runtimeMesh &mesh = cf.mesh();
            const realArray &faceWeight = mesh.faceWgt();
            const faceList &fl = mesh.faces();
            const faceZone &fz = mesh.faceZones()[faceZoneId];
            Array<Type> ff(fz.size());
            integer ffi = 0;
            for (integer fi = fz.firstIndex(); fi <= fz.lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                ff[ffi++] = faceWeight[fi] * cf[cl] + (real(1.0) - faceWeight[fi]) * cf[cr];
            }
            return ff;
        }

        /*!\brief Interpolation field onto ghost cell by given boundary face value.*/
        inline void extrapolate(const Array<Array<real>> &faceVar,
                                geometryArray<real, cellMesh> &var) {
            const auto &tb = var.mesh().Iteration();
            const auto &mesh = var.mesh();
            const auto &fL = mesh.faces();
            const auto &fZL = mesh.faceZones();

            for (integer fz = 0; fz < faceVar.size(); fz++) {
                for (integer i = 0; i < faceVar[fz].size(); i++) {
                    const integer index = fZL[fz].firstIndex();
                    const auto &inter = fL[i + index].leftCell();
                    const auto &ghost = fL[i + index].rightCell();
                    var[ghost] = real(2.0) * faceVar[fz][i] - var[inter];
                }
            }

            const cutZoneList &cZL = mesh.cutZones();
            integer ghostLayer = mesh.ghostCellsType();
            integer czs = cZL.size() / ghostLayer;

            for (integer layerI = 1; layerI < mesh.ghostCellsType(); layerI++) {               
                for (integer i = layerI * czs; i < (layerI + 1) * czs; i++) {
                    cZL[i].transfer<Array, real>(var);
                }
            }

            const perZoneList &pZL = mesh.perZones();
            for (integer i = 0; i < pZL.size(); i++) {
                pZL[i].transfer<Array, real>(var);
            }
        }

        /*!\brief Interpolation field onto ghost cell by given boundary face value.*/
        inline void extrapolate(const Array<Array<vector>> &faceVar,
                                geometryArray<vector, cellMesh> &var) {
            const auto &tb = var.mesh().Iteration();
            const auto &mesh = var.mesh();
            const auto &fL = mesh.faces();
            const auto &fZL = mesh.faceZones();
            //const integerListList& SNC = mesh.secondNeighbourCells();

            for (integer fz = 0; fz < faceVar.size(); fz++) {
                for (integer i = 0; i < faceVar[fz].size(); i++) {
                    const integer index = fZL[fz].firstIndex();
                    const auto &inter = fL[i + index].leftCell();
                    const auto &ghost = fL[i + index].rightCell();
                    var[ghost] = real(2.0) * faceVar[fz][i] - var[inter];
                }
            }

            const cutZoneList &cZL = mesh.cutZones();
            const integer ghostLayer = mesh.ghostCellsType();
            const integer czs = cZL.size() / ghostLayer;

            for (integer layerI = 1; layerI < mesh.ghostCellsType(); layerI++) {                
                for (integer i = layerI * czs; i < (layerI + 1) * czs; i++) {
                    cZL[i].transferVS<Array, vector>(var);
                }
            }

            const perZoneList &pZL = mesh.perZones();
            for (integer i = 0; i < pZL.size(); i++) {
                pZL[i].transferVS<Array, vector>(var);
            }
        }

        template <class Type>
        inline void extrapolate(const Array<Type> &faceVar, geometryArray<Type, cellMesh> &var,
                                const integer faceZoneId) {
            const auto &mesh = var.mesh();
            const auto &fs = mesh.faces();
            const auto &fZL = mesh.faceZones();
            integer count = 0;
            for (integer fi = fZL[faceZoneId].firstIndex(); fi <= fZL[faceZoneId].lastIndex();
                 ++fi) {
                const auto &cl = fs[fi].leftCell();
                const auto &cr = fs[fi].rightCell();
                var[cr] = real(2.0) * faceVar[count++] - var[cl];
            }
        }

    } // End namespace fv
} // namespace OpenHurricane
