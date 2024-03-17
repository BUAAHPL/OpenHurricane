/*!
 * \file faceInterpolation.cpp
 * \brief Main subroutines for face interpolation.
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

#include "faceInterpolation.hpp"

template <>
OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::faceMesh>
OpenHurricane::fv::interpolate(const geometryArray<real, cellMesh> &cf) {
    const runtimeMesh &mesh = cf.mesh();
    const realArray &faceWeight = mesh.faceWgt();
    const faceList &fl = mesh.faces();

    geometryArray<real, faceMesh> ff(
        object("interpolation(" + cf.name() + ")", mesh, object::NOT_WRITE, object::TEMPORARY),
        mesh);

    const integer nFaces = mesh.nFaces();

    if (argParse::isSkewCorrect()) {
        const auto &fWSkew = mesh.faceSkewnessWeight();
        const auto &isFSkew = mesh.isSkewFace();
        const auto &fISFC = mesh.fIntSectFCentre();
        if (!cf.hasGradArray()) {
            cellVectorArray grdCf = fv::grad(cf, string("leastSquareGrad"));
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

    return ff;
}

template <>
OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::faceMesh>
OpenHurricane::fv::interpolate(const geometryArray<vector, cellMesh> &cf) {
    const runtimeMesh &mesh = cf.mesh();
    const realArray &faceWeight = mesh.faceWgt();
    const faceList &fl = mesh.faces();

    geometryArray<vector, faceMesh> ff(
        object("interpolation(" + cf.name() + ")", mesh, object::NOT_WRITE, object::TEMPORARY),
        mesh);

    const integer nFaces = mesh.nFaces();

    if (argParse::isSkewCorrect()) {
        const auto &fWSkew = mesh.faceSkewnessWeight();
        const auto &isFSkew = mesh.isSkewFace();
        const auto &fISFC = mesh.fIntSectFCentre();
        if (!cf.hasGradArray()) {
            cellTensorArray grdCf = fv::grad(cf, string("leastSquareGrad"));
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

    return ff;
}