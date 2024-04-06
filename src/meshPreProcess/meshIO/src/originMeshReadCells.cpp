/*!
 * \file originMeshReadCells.cpp
 * \brief Main subroutines for reading origin mesh cells.
 * \author Rao Sihang
 * \version V2.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of Hurricane project (Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation)
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

#include "logFile.hpp"
#include "originMeshRead.hpp"

void OpenHurricane::originMeshRead::assignCellType(cellList &cells) {
    for (integer celli = 0; celli < cells.size(); celli++) {
        if (cells[celli].isTriangular()) {
            cells[celli].nodesList().resize(3);
            cells[celli].facesList().resize(3);
        } else if (cells[celli].isTetrahedral()) {
            cells[celli].nodesList().resize(4);
            cells[celli].facesList().resize(4);
        } else if (cells[celli].isQuadrilateral()) {
            cells[celli].nodesList().resize(4);
            cells[celli].facesList().resize(4);
        } else if (cells[celli].isHexahedral()) {
            cells[celli].nodesList().resize(8);
            cells[celli].facesList().resize(6);
        } else if (cells[celli].isPyramid()) {
            cells[celli].nodesList().resize(5);
            cells[celli].facesList().resize(5);
        } else if (cells[celli].isWedge()) {
            cells[celli].nodesList().resize(6);
            cells[celli].facesList().resize(5);
        }
    }
}
void OpenHurricane::originMeshRead::formingNodes(const integer globalCells, const faceList &faces,
                                                 cellList &cells) {
    Pout << "    Info: computing forming node no. for each cell... " << std::endl;
    for (integer celli = 0; celli < globalCells; celli++) {
        if (cells[celli].isTriangular()) {
            triangularCellNodes(celli, faces, cells);
        } else if (cells[celli].isQuadrilateral()) {
            quadrilateralCellNodes(celli, faces, cells);
        } else if (cells[celli].isTetrahedral()) {
            tetrahedralCellNodes(celli, faces, cells);
        } else if (cells[celli].isPyramid()) {
            pyramidCellNodes(celli, faces, cells);
        } else if (cells[celli].isWedge()) {
            wedgeCellNodes(celli, faces, cells);
        } else if (cells[celli].isHexahedral()) {
            hexahedralCellNodes(celli, faces, cells);
        } else if (cells[celli].isPolyhedral()) {
            polyhedralCellNodes(celli, faces, cells);
        }
    }
}

void OpenHurricane::originMeshRead::triangularCellNodes(const integer celli, const faceList &faces,
                                                        cellList &cells) {
    LFatal("Triangular grid cell are not supported in current version!");
}

void OpenHurricane::originMeshRead::tetrahedralCellNodes(const integer celli, const faceList &faces,
                                                         cellList &cells) {
    integer i = cells[celli].facei(0);
    if (celli == faces[i].leftCell()) {
        for (integer j = 0; j < 3; j++) {
            cells[celli].nodesList()[j] = faces[i][j];
        }
    } else {
        for (integer j = 0; j < 3; j++) {
            cells[celli].nodesList()[j] = faces[i][2 - j];
        }
    }

    i = cells[celli].facei(1);
    for (integer j = 0; j < 3; j++) {
        integer k = faces[i][j];
        if (k != cells[celli].nodei(0) && k != cells[celli].nodei(1) &&
            k != cells[celli].nodei(2)) {
            cells[celli].nodesList()[3] = k;
            break;
        }
    }

    // Reordering faces
    integerList &fl = cells[celli].facesList();

    integerList nl = cells[celli].nodesList().sub(0, 2);
    integer iflag = 0;
    for (integer j = 1; j < 4; j++) {
        i = fl[j];
        if (isContained(nl, faces[i])) {
            iflag = j;
            break;
        }
    }

    if (iflag == 0) {
        LFatal("cannot forming nodes for tetrahedral");
    }
    if (iflag != 1) {
        integer temp = fl[1];
        fl[1] = fl[iflag];
        fl[iflag] = temp;
    }

    integerList nl1 = cells[celli].nodesList().sub(1, 2);
    i = fl[2];
    if (isContained(nl, faces[i])) {
        integer temp = fl[2];
        fl[2] = fl[3];
        fl[3] = temp;
    }
}

void OpenHurricane::originMeshRead::quadrilateralCellNodes(const integer celli,
                                                           const faceList &faces, cellList &cells) {
    LFatal("Quadrilateral grid cell are not supported in current version!");
}

void OpenHurricane::originMeshRead::pyramidCellNodes(const integer celli, const faceList &faces,
                                                     cellList &cells) {
    integer i;
    integer n = -1;
    for (integer m = 0; m < 5; m++) {
        i = cells[celli].facei(m);
        if (faces[i].isQuadrilateral()) {
            n = m;
            break;
        }
    }

    if (celli == faces[i].leftCell()) {
        for (integer j = 0; j < 4; j++) {
            cells[celli].nodesList()[j] = faces[i][j];
        }
    } else {
        for (integer j = 0; j < 4; j++) {
            cells[celli].nodesList()[j] = faces[i][3 - j];
        }
    }

    if (n == 1) {
        i = 2;
    } else {
        i = 1;
    }

    i = cells[celli].facei(i);

    for (integer j = 0; j < 3; j++) {
        integer k = faces[i][j];
        if (k != cells[celli].nodei(0) && k != cells[celli].nodei(1) &&
            k != cells[celli].nodei(2) && k != cells[celli].nodei(3)) {
            cells[celli].nodesList()[4] = k;
            break;
        }
    }

    /*!\brief Reordering faces.*/

    integerList &fl = cells[celli].facesList();

    if (n != 4) {
        integer temp = fl[n];
        fl[n] = fl[4];
        fl[4] = temp;
    }

    integerList nl(2);
    nl[0] = cells[celli].nodei(0);
    nl[1] = cells[celli].nodei(3);
    for (integer j = 0; j < 4; j++) {
        i = fl[j];
        if (isContained(nl, faces[i])) {
            if (j != 0) {
                integer temp = fl[j];
                fl[j] = fl[0];
                fl[0] = temp;
            }
            break;
        }
    }

    nl[0] = cells[celli].nodei(1);
    nl[1] = cells[celli].nodei(2);
    for (integer j = 1; j < 4; j++) {
        i = fl[j];
        if (isContained(nl, faces[i])) {
            if (j != 1) {
                integer temp = fl[j];
                fl[j] = fl[1];
                fl[1] = temp;
            }
            break;
        }
    }

    nl[0] = cells[celli].nodei(0);
    nl[1] = cells[celli].nodei(1);
    for (integer j = 2; j < 4; j++) {
        i = fl[j];
        if (isContained(nl, faces[i])) {
            if (j != 2) {
                integer temp = fl[j];
                fl[j] = fl[2];
                fl[2] = temp;
            }
            break;
        }
    }
}

void OpenHurricane::originMeshRead::wedgeCellNodes(const integer celli, const faceList &faces,
                                                   cellList &cells) {
    integer i1 = -1, i2 = -1;
    integer m1, m2;
    for (integer m = 0; m < 5; m++) {
        integer i = cells[celli].facei(m);
        if (faces[i].isTriangular()) {
            if (i1 == -1) {
                i1 = i;
                m1 = m;
            } else {
                i2 = i;
                m2 = m;
                break;
            }
        }
    }

    if (celli == faces[i1].leftCell()) {
        for (integer j = 0; j < 3; j++) {
            cells[celli].nodesList()[j] = faces[i1][j];
        }
    } else {
        for (integer j = 0; j < 3; j++) {
            cells[celli].nodesList()[j] = faces[i1][2 - j];
        }
    }

    if (celli == faces[i2].leftCell()) {
        for (integer j = 0; j < 3; j++) {
            cells[celli].nodesList()[j + 3] = faces[i2][2 - j];
        }
    } else {
        for (integer j = 0; j < 3; j++) {
            cells[celli].nodesList()[j + 3] = faces[i2][j];
        }
    }

    bool hasDone = false;
    while (true) {
        integerList nl1(2);
        nl1[0] = cells[celli].nodei(0);
        nl1[1] = cells[celli].nodei(1);
        integerList nl2(2);
        nl2[0] = cells[celli].nodei(3);
        nl2[1] = cells[celli].nodei(4);
        for (integer j = 0; j < 5; j++) {
            integer i = cells[celli].facei(j);
            if (i != i1 && i != i2) {
                if (isContained(nl1, faces[i]) && isContained(nl2, faces[i])) {
                    hasDone = true;
                    break;
                }
            }
        }
        if (hasDone) {
            break;
        } else {
            integer k = cells[celli].nodei(5);
            cells[celli].nodesList()[5] = cells[celli].nodesList()[4];
            cells[celli].nodesList()[4] = cells[celli].nodesList()[3];
            cells[celli].nodesList()[3] = k;
        }
    }

    /*!\brief Reordering faces.*/

    if (m2 != 4) {
        integer k = cells[celli].facei(4);
        cells[celli].facesList()[4] = cells[celli].facei(m2);
        cells[celli].facesList()[m2] = k;
    }
    if (m1 != 3) {
        integer k = cells[celli].facei(3);
        cells[celli].facesList()[3] = cells[celli].facei(m1);
        cells[celli].facesList()[m1] = k;
    }

    integerList nl(2);
    nl[0] = cells[celli].nodei(0);
    nl[1] = cells[celli].nodei(2);
    for (integer j = 0; j < 3; j++) {
        integer k = cells[celli].facei(j);
        if (isContained(nl, faces[k])) {
            if (j != 0) {
                cells[celli].facesList()[j] = cells[celli].facei(0);
                cells[celli].facesList()[0] = k;
                break;
            }
        }
    }

    nl[0] = cells[celli].nodei(1);
    nl[1] = cells[celli].nodei(2);
    for (integer j = 1; j < 3; j++) {
        integer k = cells[celli].facei(j);
        if (isContained(nl, faces[k])) {
            if (j != 1) {
                cells[celli].facesList()[j] = cells[celli].facei(1);
                cells[celli].facesList()[1] = k;
                break;
            }
        }
    }
}

void OpenHurricane::originMeshRead::hexahedralCellNodes(const integer celli, const faceList &faces,
                                                        cellList &cells) {
    integer i = cells[celli].facei(0);
    integer i1 = 0;
    integer i2;
    if (celli == faces[i].leftCell()) {
        for (integer j = 0; j < 4; j++) {
            cells[celli].nodesList()[j] = faces[i][j];
        }
    } else {
        for (integer j = 0; j < 4; j++) {
            cells[celli].nodesList()[j] = faces[i][3 - j];
        }
    }

    integer j;
    for (integer m = 1; m < 6; m++) {
        j = cells[celli].facei(m);
        i2 = m;
        integer equal = 0;
        for (integer n = 0; n < 4; n++) {
            for (integer t = 0; t < 4; t++) {
                if (faces[i][t] == faces[j][n]) {
                    equal += 1;
                    break;
                }
            }
            if (equal != 0) {
                break;
            }
        }
        if (equal == 0) {
            break;
        }
    }

    if (celli == faces[j].leftCell()) {
        for (integer m = 0; m < 4; m++) {
            cells[celli].nodesList()[m + 4] = faces[j][3 - m];
        }
    } else {
        for (integer m = 0; m < 4; m++) {
            cells[celli].nodesList()[m + 4] = faces[j][m];
        }
    }

    bool hasDone = false;
    while (true) {
        integerList nl1(2);
        nl1[0] = cells[celli].nodei(0);
        nl1[1] = cells[celli].nodei(1);
        integerList nl2(2);
        nl2[0] = cells[celli].nodei(4);
        nl2[1] = cells[celli].nodei(5);
        for (integer m = 1; m < 6; m++) {
            i = cells[celli].facei(m);

            if (isContained(nl1, faces[i]) && isContained(nl2, faces[i])) {
                hasDone = true;
                break;
            }
        }

        if (hasDone) {
            break;
        } else {
            integer k = cells[celli].nodei(7);
            cells[celli].nodesList()[7] = cells[celli].nodesList()[6];
            cells[celli].nodesList()[6] = cells[celli].nodesList()[5];
            cells[celli].nodesList()[5] = cells[celli].nodesList()[4];
            cells[celli].nodesList()[4] = k;
        }
    }

    // Reordering faces

    if (i2 != 5) {
        integer temp = cells[celli].facei(5);
        cells[celli].facesList()[5] = cells[celli].facei(i2);
        cells[celli].facesList()[i2] = temp;
    }

    if (i1 != 4) {
        integer temp = cells[celli].facei(4);
        cells[celli].facesList()[4] = cells[celli].facei(i1);
        cells[celli].facesList()[i1] = temp;
    }

    integerList nl(2);
    nl[0] = cells[celli].nodei(0);
    nl[1] = cells[celli].nodei(3);
    for (integer m = 0; m < 4; m++) {
        integer k = cells[celli].facei(m);
        if (isContained(nl, faces[k])) {
            if (m != 0) {
                cells[celli].facesList()[m] = cells[celli].facei(0);
                cells[celli].facesList()[0] = k;
                break;
            }
        }
    }

    nl[0] = cells[celli].nodei(1);
    nl[1] = cells[celli].nodei(2);
    for (integer m = 1; m < 4; m++) {
        integer k = cells[celli].facei(m);
        if (isContained(nl, faces[k])) {
            if (m != 1) {
                cells[celli].facesList()[m] = cells[celli].facei(1);
                cells[celli].facesList()[1] = k;
                break;
            }
        }
    }

    nl[0] = cells[celli].nodei(0);
    nl[1] = cells[celli].nodei(1);
    for (integer m = 2; m < 4; m++) {
        integer k = cells[celli].facei(m);
        if (isContained(nl, faces[k])) {
            if (m != 2) {
                cells[celli].facesList()[m] = cells[celli].facei(2);
                cells[celli].facesList()[2] = k;
                break;
            }
        }
    }
}

void OpenHurricane::originMeshRead::polyhedralCellNodes(const integer celli, const faceList &faces,
                                                        cellList &cells) {
    LFatal("Polyhedral grid cell are not supported in current version!");
}

/*!\brief Return true if list l1 is contained by l2.*/
bool OpenHurricane::originMeshRead::isContained(const integerList &l1, const integerList &l2) {
    integer j = 0;
    for (integer i = 0; i < l1.size(); i++) {
        for (j = 0; j < l2.size(); j++) {
            if (l2[j] == l1[i]) {
                break;
            }
        }

        if (j == l2.size()) {
            return false;
        }
    }
    return true;
}

void OpenHurricane::originMeshRead::formingFaces(const integer globalCells, const faceList &faces,
                                                 cellList &cells) {
    assignCellType(cells);

    Pout << "    Info: computing forming face no. for each cell... " << std::endl;

    integerList color(globalCells, OpenHurricane::Zero);

    for (integer facei = 0; facei < faces.size(); facei++) {
        integer cl = faces[facei].leftCell();
        integer cr = faces[facei].rightCell();

        integer k;
        if (cl != -1) {
            k = color[cl];
            if (cells[cl].isPolyhedral()) {
                cells[cl].facesList().append(facei);
            } else {
                cells[cl].facesList()[k] = facei;
            }
            color[cl] += 1;
        }
        if (cr != -1) {
            k = color[cr];
            if (cells[cr].isPolyhedral()) {
                cells[cr].facesList().append(facei);
            } else {
                cells[cr].facesList()[k] = facei;
            }
            color[cr] += 1;
        }
    }

    for (integer celli = 0; celli < cells.size(); celli++) {
        if (color[celli] < cells[celli].faceSize()) {
            errorAbortStr(("Cell is missing face: cell id " + toString(celli)));
        }
    }
}