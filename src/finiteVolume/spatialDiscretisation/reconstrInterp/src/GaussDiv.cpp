/*!
 * \file GaussDiv.cpp
 * \brief Main subroutines for computing divergence with Gauss theorem.
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

#include "GaussDiv.hpp"
#include "cellMesh.hpp"
#include "faceMesh.hpp"

namespace OpenHurricane {
    createClassNameStr(GaussDiv, "GaussDiv");
}

void OpenHurricane::GaussDiv::calcDiv(const geometryArray<vector, cellMesh> &cellQ,
                                      geometryArray<real, cellMesh> &div) {
    const auto &mesh = cellQ.mesh();
    const integer nCells = mesh.nCells();
    for (integer i = 0; i < nCells; ++i) {
        div[i] = Zero;
    }
    const integer nFaces = mesh.nFaces();

    const auto &fl = mesh.faces();
    const auto &fw = mesh.faceWeight();
    const auto &fA = mesh.faceArea();
    for (integer i = 0; i < nFaces; ++i) {
        const auto &cl = fl[i].leftCell();
        const auto &cr = fl[i].rightCell();

        const real wl = fw[i];
        const real wr = real(1.0) - wl;
        const vector qf = wl * cellQ[cl] + wr * cellQ[cr];

        div[cl] -= qf * fA[i];

        if (cr < nCells) {
            div[cr] += qf * fA[i];
        }
    }

    const auto &vol = mesh.cellVolume();
    for (integer i = 0; i < nCells; ++i) {
        div[i] /= vol[i];
    }
    updateBoundaryField<real>(div);
}

void OpenHurricane::GaussDiv::calcDiv(const geometryArray<tensor, cellMesh> &cellQ,
                                      geometryArray<vector, cellMesh> &div) {
    const auto &mesh = cellQ.mesh();
    const integer nCells = mesh.nCells();
    for (integer i = 0; i < nCells; ++i) {
        div[i] = Zero;
    }
    const integer nFaces = mesh.nFaces();

    const auto &fl = mesh.faces();
    const auto &fw = mesh.faceWeight();
    const auto &fA = mesh.faceArea();
    for (integer i = 0; i < nFaces; ++i) {
        const auto &cl = fl[i].leftCell();
        const auto &cr = fl[i].rightCell();

        const real wl = fw[i];
        const real wr = real(1.0) - wl;
        const tensor qf = wl * cellQ[cl] + wr * cellQ[cr];

        div[cl] -= qf * fA[i];

        if (cr < nCells) {
            div[cr] += qf * fA[i];
        }
    }

    const auto &vol = mesh.cellVolume();
    for (integer i = 0; i < nCells; ++i) {
        div[i] /= vol[i];
    }
    updateBoundaryField<vector>(div);
}

void OpenHurricane::GaussDiv::calcDiv(const geometryArray<vector, cellMesh> &cellQ,
                                      const geometryArray<vector, faceMesh> &faceQ,
                                      geometryArray<real, cellMesh> &div) {
    const auto &mesh = cellQ.mesh();
    const integer nCells = mesh.nCells();
    for (integer i = 0; i < nCells; ++i) {
        div[i] = Zero;
    }
    const integer nFaces = mesh.nFaces();

    const auto &fl = mesh.faces();
    const auto &fw = mesh.faceWeight();
    const auto &fA = mesh.faceArea();
    for (integer i = 0; i < nFaces; ++i) {
        const auto &cl = fl[i].leftCell();
        const auto &cr = fl[i].rightCell();

        div[cl] -= faceQ[i] * fA[i];

        if (cr < nCells) {
            div[cr] += faceQ[i] * fA[i];
        }
    }

    const auto &vol = mesh.cellVolume();
    for (integer i = 0; i < nCells; ++i) {
        div[i] /= vol[i];
    }
    updateBoundaryField<real>(div);
}

void OpenHurricane::GaussDiv::calcDiv(const geometryArray<tensor, cellMesh> &cellQ,
                                      const geometryArray<tensor, faceMesh> &faceQ,
                                      geometryArray<vector, cellMesh> &div) {
    const auto &mesh = cellQ.mesh();
    const integer nCells = mesh.nCells();
    for (integer i = 0; i < nCells; ++i) {
        div[i] = Zero;
    }
    const integer nFaces = mesh.nFaces();

    const auto &fl = mesh.faces();
    const auto &fw = mesh.faceWeight();
    const auto &fA = mesh.faceArea();
    for (integer i = 0; i < nFaces; ++i) {
        const auto &cl = fl[i].leftCell();
        const auto &cr = fl[i].rightCell();

        div[cl] -= faceQ[i] * fA[i];

        if (cr < nCells) {
            div[cr] += faceQ[i] * fA[i];
        }
    }

    const auto &vol = mesh.cellVolume();
    for (integer i = 0; i < nCells; ++i) {
        div[i] /= vol[i];
    }
    updateBoundaryField<vector>(div);
}

OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>
OpenHurricane::GaussDiv::div(const geometryArray<vector, cellMesh> &cellQ) {
    geometryArray<real, cellMesh> ff(
        object("div(" + cellQ.name() + ")", cellQ.mesh(), object::NOT_WRITE, object::TEMPORARY),
        cellQ.mesh());
    calcDiv(cellQ, ff);
    return ff;
}

OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::cellMesh>
OpenHurricane::GaussDiv::div(const geometryArray<tensor, cellMesh> &cellQ) {
    geometryArray<vector, cellMesh> ff(
        object("div(" + cellQ.name() + ")", cellQ.mesh(), object::NOT_WRITE, object::TEMPORARY),
        cellQ.mesh());
    calcDiv(cellQ, ff);
    return ff;
}

OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>
OpenHurricane::GaussDiv::div(const geometryArray<vector, cellMesh> &cellQ,
                             const geometryArray<vector, faceMesh> &faceQ) {
    geometryArray<real, cellMesh> ff(
        object("div(" + cellQ.name() + ")", cellQ.mesh(), object::NOT_WRITE, object::TEMPORARY),
        cellQ.mesh());
    calcDiv(cellQ, faceQ, ff);
    return ff;
}

OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::cellMesh>
OpenHurricane::GaussDiv::div(const geometryArray<tensor, cellMesh> &cellQ,
                             const geometryArray<tensor, faceMesh> &faceQ) {
    geometryArray<vector, cellMesh> ff(
        object("div(" + cellQ.name() + ")", cellQ.mesh(), object::NOT_WRITE, object::TEMPORARY),
        cellQ.mesh());
    calcDiv(cellQ, faceQ, ff);
    return ff;
}