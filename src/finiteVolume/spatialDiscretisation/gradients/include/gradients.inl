#include "gradients.hpp"
/*!
 * \file gradients.inl
 * \brief In-Line subroutines of the <i>gradients.hpp</i> file.
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

namespace OpenHurricane {
    namespace fv {
        inline uniquePtr<gradient> gradScheme(const string &name) {
            controller cont("grad");
            cont.add(gradient::className_, name);
            return gradient::creator(cont);
        }

        template <class Type>
        inline geometryArray<typename outerProduct<vector, Type>::type, cellMesh>
        grad(const geometryArray<Type, cellMesh> &cellQ, const string &name) {
            string nn = "grad(" + cellQ.name() + ")";
            return gradScheme(name)->grad<Type>(cellQ, nn);
        }

        template <class Type>
        inline geometryArray<typename outerProduct<vector, Type>::type, faceMesh>
        gradf(const geometryArray<Type, cellMesh> &cellQ,
              const geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &cellGrad) {
            return faceGrad::grad(cellQ, cellGrad);
        }

        template <class Type>
        inline geometryArray<typename outerProduct<vector, Type>::type, faceMesh>
        gradf(const geometryArray<Type, cellMesh> &cellQ) {
            return faceGrad::grad(cellQ, cellQ.grad());
        }

        template <class Type>
        static void
        gradf(const geometryArray<Type, cellMesh> &cellQ,
              geometryArray<typename outerProduct<vector, Type>::type, faceMesh> &faceG) {
            faceGrad::grad(cellQ, cellQ.grad(), faceG);
        }

        template <class Type>
        Array<typename outerProduct<vector, Type>::type>
        gradf(const geometryArray<Type, cellMesh> &cellQ,
              const geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &cellGrad,
              const integer faceZoneId) {
            const runtimeMesh &mesh = cellQ.mesh();
            const faceList &fl = mesh.faces();
            const faceZone &fz = mesh.faceZones()[faceZoneId];
            Array<typename outerProduct<vector, Type>::type> fg(fz.size());
            integer fgi = 0;
            for (integer fi = fz.firstIndex(); fi <= fz.lastIndex(); ++fi) {
                fg[fgi++] = faceGrad::grad(cellQ, cellGrad, fi);
            }
            return fg;
        }

        template <class Type>
        Array<typename outerProduct<vector, Type>::type>
        gradf(const geometryArray<Type, cellMesh> &cellQ, const integer faceZoneId) {
            const runtimeMesh &mesh = cellQ.mesh();
            const faceList &fl = mesh.faces();
            const faceZone &fz = mesh.faceZones()[faceZoneId];
            Array<typename outerProduct<vector, Type>::type> fg(fz.size());
            integer fgi = 0;
            for (integer fi = fz.firstIndex(); fi <= fz.lastIndex(); ++fi) {
                fg[fgi++] = faceGrad::grad(cellQ, cellQ.grad(), fi);
            }
            return fg;
        }
    } // namespace fv
} // namespace OpenHurricane