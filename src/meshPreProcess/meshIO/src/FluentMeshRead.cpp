/*!
 * \file FluentMeshRead.cpp
 * \brief Main subroutines for reading fluent mesh.
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

#include "FluentMeshRead.hpp"
#include "logFile.hpp"
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <streambuf>
namespace OpenHurricane {
    createClassNameStr(FluentMeshRead, "fluent");
    registerObjFty(originMeshRead, FluentMeshRead, controller);
} // namespace OpenHurricane

void OpenHurricane::FluentMeshRead::createSectionIndexFluentMap() {
    sectionIndexFluentMap_ = createMap<short int, sectionIndexFluent>(0, XF_COMMENT)(1, XF_HEADER)(
        2, XF_DIMENSION)(0, XF_COMMENT)(10, XF_NODE)(18, XF_PERIODIC_FACE)(12, XF_CELL)(
        13, XF_FACE)(11, XF_EDGE)(59, XF_FACE_TREE)(58, XF_CELL_TREE)(61, XF_INTERFACE_PARENTS)(
        39, XF_RP_TV)(45, XF_RP_TV1);
}

void OpenHurricane::FluentMeshRead::removeBeginSpaceFluent(std::string &str,
                                                           size_t &firstSpacePos) const {
    while (true) {
        firstSpacePos = str.find_first_of(' ');
        if (firstSpacePos == 0) {
            str.erase(0, firstSpacePos + 1);
            continue;
        } else {
            break;
        }
    }
}

// Low efficiency
void OpenHurricane::FluentMeshRead::readFaceConnectionFluent(const faceZone &fZ,
                                                             const std::string &str,
                                                             const integer id, const short bcT) {
    std::stringstream sstr;
    sstr << std::hex << str;

    if (fZ.isLinear()) {
        faces_[id].resize(2);
        sstr >> faces_[id][0] >> faces_[id][1] >> faces_[id].leftCell() >> faces_[id].rightCell();

    } else if (fZ.isTriangular()) {
        faces_[id].resize(3);
        sstr >> faces_[id][0] >> faces_[id][1] >> faces_[id][2] >> faces_[id].leftCell() >>
            faces_[id].rightCell();
    } else if (fZ.isQuadrilateral()) {
        faces_[id].resize(4);
        sstr >> faces_[id][0] >> faces_[id][1] >> faces_[id][2] >> faces_[id][3] >>
            faces_[id].leftCell() >> faces_[id].rightCell();
    } else if (fZ.isMixed() || fZ.isPolygonal()) {
        short x;
        sstr >> x;
        faces_[id].resize(x);

        for (short i = 0; i < x; i++) {
            sstr >> faces_[id][i];
        }

        sstr >> faces_[id].leftCell() >> faces_[id].rightCell();
    }

    for (integer i = 0; i < faces_[id].size(); i++) {
        faces_[id][i] -= 1;
    }
    faces_[id].leftCell() -= 1;
    faces_[id].rightCell() -= 1;
    faces_[id].setBCType(bcT);
}

// High efficiency
void OpenHurricane::FluentMeshRead::readFaceConnectionNewFluent(const faceZone &fZ,
                                                                const std::string &str,
                                                                const integer fI, const integer lI,
                                                                const short bcT) {
    char *endptr;
    if (fZ.isLinear()) {
        for (integer i = fI - 1; i < lI; i++) {
            faces_[i].resize(2);
        }
        faces_[fI - 1][0] = (integer)strtol(str.c_str(), &endptr, 16);
        faces_[fI - 1][1] = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1].leftCell() = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1].rightCell() = (integer)strtol(endptr, &endptr, 16);
        for (integer i = fI; i < lI; i++) {
            faces_[i][0] = (integer)strtol(endptr, &endptr, 16);
            faces_[i][1] = (integer)strtol(endptr, &endptr, 16);
            faces_[i].leftCell() = (integer)strtol(endptr, &endptr, 16);
            faces_[i].rightCell() = (integer)strtol(endptr, &endptr, 16);
        }
    } else if (fZ.isTriangular()) {
        for (integer i = fI - 1; i < lI; i++) {
            faces_[i].resize(3);
        }
        faces_[fI - 1][0] = (integer)strtol(str.c_str(), &endptr, 16);
        faces_[fI - 1][1] = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1][2] = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1].leftCell() = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1].rightCell() = (integer)strtol(endptr, &endptr, 16);
        for (integer i = fI; i < lI; i++) {
            faces_[i][0] = (integer)strtol(endptr, &endptr, 16);
            faces_[i][1] = (integer)strtol(endptr, &endptr, 16);
            faces_[i][2] = (integer)strtol(endptr, &endptr, 16);
            faces_[i].leftCell() = (integer)strtol(endptr, &endptr, 16);
            faces_[i].rightCell() = (integer)strtol(endptr, &endptr, 16);
        }
    } else if (fZ.isQuadrilateral()) {
        for (integer i = fI - 1; i < lI; i++) {
            faces_[i].resize(4);
        }
        faces_[fI - 1][0] = (integer)strtol(str.c_str(), &endptr, 16);
        faces_[fI - 1][1] = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1][2] = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1][3] = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1].leftCell() = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1].rightCell() = (integer)strtol(endptr, &endptr, 16);
        for (integer i = fI; i < lI; i++) {
            faces_[i][0] = (integer)strtol(endptr, &endptr, 16);
            faces_[i][1] = (integer)strtol(endptr, &endptr, 16);
            faces_[i][2] = (integer)strtol(endptr, &endptr, 16);
            faces_[i][3] = (integer)strtol(endptr, &endptr, 16);
            faces_[i].leftCell() = (integer)strtol(endptr, &endptr, 16);
            faces_[i].rightCell() = (integer)strtol(endptr, &endptr, 16);
        }
    } else if (fZ.isMixed() || fZ.isPolygonal()) {
        integer x;
        x = (integer)strtol(str.c_str(), &endptr, 16);
        faces_[fI - 1].resize(x);
        for (integer j = 0; j < x; j++) {
            faces_[fI - 1][j] = (integer)strtol(endptr, &endptr, 16);
        }
        faces_[fI - 1].leftCell() = (integer)strtol(endptr, &endptr, 16);
        faces_[fI - 1].rightCell() = (integer)strtol(endptr, &endptr, 16);

        for (integer i = fI; i < lI; i++) {
            x = (integer)strtol(endptr, &endptr, 16);
            faces_[i].resize(x);
            for (integer j = 0; j < x; j++) {
                faces_[i][j] = (integer)strtol(endptr, &endptr, 16);
            }
            faces_[i].leftCell() = (integer)strtol(endptr, &endptr, 16);
            faces_[i].rightCell() = (integer)strtol(endptr, &endptr, 16);
        }
    }

    for (integer i = fI - 1; i < lI; i++) {
        for (integer j = 0; j < faces_[i].size(); j++) {
            faces_[i][j] -= 1;
        }
        faces_[i].leftCell() -= 1;
        faces_[i].rightCell() -= 1;
        faces_[i].setBCType(bcT);
    }
}

void OpenHurricane::FluentMeshRead::readMixedCellShapeFluent(const std::string &str,
                                                             const integer fI, const integer lI) {

    char *endptr;
    short cT;

    cT = (short)strtol(str.c_str(), &endptr, 16);
    cells_[fI - 1].setType(cT);
    for (integer i = fI; i <= lI - 1; i++) {
        cT = (short)strtol(endptr, &endptr, 16);
        cells_[i].setType(cT);
        cells_[i].setOrignIndex(i);
    }
}

void OpenHurricane::FluentMeshRead::readNodesCoordinatesFluent(const std::string &str,
                                                               const integer fI, const integer lI,
                                                               unsigned short nd) {

    char *endptr;
    if (nd == 2) {
        points_[fI - 1].x() = (real)strtod(str.c_str(), &endptr);
        endptr++;
        points_[fI - 1].y() = (real)strtod(endptr, &endptr);

    } else {
        points_[fI - 1].x() = (real)strtod(str.c_str(), &endptr);
        endptr++;
        points_[fI - 1].y() = (real)strtod(endptr, &endptr);
        endptr++;
        points_[fI - 1].z() = (real)strtod(endptr, &endptr);
    }

    for (integer i = fI; i <= lI - 1; i++) {
        if (nd == 2) {
            points_[i].x() = (real)strtod(endptr, &endptr);
            endptr++;
            points_[i].y() = (real)strtod(endptr, &endptr);
        } else {

            points_[i].x() = (real)strtod(endptr, &endptr);
            endptr++;
            points_[i].y() = (real)strtod(endptr, &endptr);
            endptr++;
            points_[i].z() = (real)strtod(endptr, &endptr);
        }
    }
}

bool OpenHurricane::FluentMeshRead::isFluentMeshHasBinaryFormat() const {
    std::ifstream myin(fileName_, std::ios::binary);
    if (!myin.good()) {
        LFatal("File: \"%s\" is not exist", fileName_.c_str());
    }
    std::string line;
    std::regex myReg("\\s*\\(\\s*([0-9]+).*");
    std::regex myReg1("(|20|30)10");
    std::regex myReg2("(|20|30)12");
    std::regex myReg3("(|20|30)13");
    bool isBinary = false;
    integer count = 0;
    while (std::getline(myin, line)) {
        replaceAllMarks(line, "\r", "");
        replaceAllMarks(line, "\n", "");
        std::smatch what;
        if (std::regex_match(line, what, myReg)) {
            std::string index = what[1];
            std::smatch what1;
            // Points
            // (10 (zone-id first-index last-index type ND)
            if (std::regex_match(index, what1, myReg1)) {
                if (what1[0] == "2010" || what1[0] == "3010") {
                    isBinary = true;
                    break;
                } else {
                    std::regex myRegP("\\s*\\(\\s*10\\s*\\(\\s*0.*");
                    if (!std::regex_match(line, myRegP)) {
                        isBinary = false;
                        break;
                    } else {
                        continue;
                    }
                }
            } else if (std::regex_match(index, what1, myReg2)) {
                if (what1[0] == "2012" || what1[0] == "3012") {
                    isBinary = true;
                    break;
                } else {
                    std::regex myRegP("\\s*\\(\\s*12\\s*\\(\\s*0");
                    if (!std::regex_match(line, myRegP)) {
                        isBinary = false;
                        break;
                    } else {
                        continue;
                    }
                }
            } else if (std::regex_match(index, what1, myReg3)) {
                if (what1[0] == "2013" || what1[0] == "3013") {
                    isBinary = true;
                    break;
                } else {
                    std::regex myRegP("\\s*\\(\\s*13\\s*\\(\\s*0");
                    if (!std::regex_match(line, myRegP)) {
                        isBinary = false;
                        break;
                    } else {
                        continue;
                    }
                }
            }
        } else {
            if (count++ > 100) {
                break;
            } else {
                continue;
            }
        }
    }
    myin.close();
    return isBinary;
}

void OpenHurricane::FluentMeshRead::countUnreferencedFaces(integer &n, integerList &unrefFL,
                                                           integerList &unrerfFinFZ) const {
    n = 0;
    unrerfFinFZ.resize(faceZones_.size(), Zero);
    unrefFL.resize(faces_.size(), 1);

    for (integer fzi = 0; fzi < faceZones_.size(); ++fzi) {
        for (integer fi = faceZones_[fzi].firstIndex(); fi <= faceZones_[fzi].lastIndex(); ++fi) {
            if (faces_[fi].leftCell() == -1 && faces_[fi].rightCell() == -1) {
                unrerfFinFZ[fzi]++;
                unrefFL[fi] = 0;
                n++;
            }
        }
    }
}

void OpenHurricane::FluentMeshRead::removeUnreferencedFaces() {
    integer n;
    integerList unrefFL;
    integerList unrerfFinFZ;

    countUnreferencedFaces(n, unrefFL, unrerfFinFZ);
    if (n != 0) {
        for (integer fzi = 0; fzi < unrerfFinFZ.size(); ++fzi) {
            if (unrerfFinFZ[fzi] != 0) {
                Pout << "    Removing " << unrerfFinFZ[fzi] << " unferenced faces in "
                     << faceZones_[fzi].name() << std::endl;
            }
        }

        faceList tmpFaces(unrefFL.size() - n);

        integer count = 0;
        for (integer fi = 0; fi < unrefFL.size(); ++fi) {
            if (unrefFL[fi] == 1) {
                tmpFaces[count].transfer(faces_[fi]);
                count++;
            }
        }
        faces_.transfer(tmpFaces);
        globalFaces_ -= n;

        for (integer fzi = 0; fzi < unrerfFinFZ.size(); ++fzi) {
            if (unrerfFinFZ[fzi] != 0) {
                // Firstly, renumber other face zones of which the first index is greater than
                // the last index of the zone including unreferenced faces.
                for (integer fzj = 0; fzj < faceZones_.size(); ++fzj) {
                    // If the same zones, continue.
                    if (fzi == fzj) {
                        continue;
                    }
                    if (faceZones_[fzj].firstIndex() > faceZones_[fzi].lastIndex()) {
                        faceZones_[fzj].setFirstIndex(faceZones_[fzj].firstIndex() -
                                                      unrerfFinFZ[fzi]);
                        faceZones_[fzj].setLastIndex(faceZones_[fzj].lastIndex() -
                                                     unrerfFinFZ[fzi]);
                    }
                }
                // Secondly, renumber the face zone that has unreferenced faces.
                faceZones_[fzi].setLastIndex(faceZones_[fzi].lastIndex() - unrerfFinFZ[fzi]);
            }
        }
    }
}

// Parsing section header for nodes, faces or cells section.
// (1) The format of nodes section header is as follows
//      (10 (zone-id first-index last-index bcType ND)(
//     or
//      (10 (zone-id first-index last-index bcType ND)
//      (
//   When zone-id = 0, first-index will be one, last-index will be
//   the total number of nodes in hexadecimal, bcType is zero, ND is omitted.
// (2) The format of faces section header is as follows
//      (13 (zone-id first-index last-index bc-bcType face-bcType))
//   When zone-id = 0, first-index will be one, last-index will be
//   the total number of faces in hexadecimal, bc-bcType and face-bcType are omitted.
// (3) The format of cells section header is as follows
//      (12 (zone-id first-index last-index bcType element-bcType))
//   When zone-id = 0, first-index will be one, last-index will be
//   the total number of cells in hexadecimal, bcType and element-bcType are omitted.
// Return 0 if it's the total number section.
// Return 1 if it's the single section.
short OpenHurricane::FluentMeshRead::parsingHeaderStrFluent(const std::string &str, integer &zoneId,
                                                            integer &fi, integer &li, short &type,
                                                            short &type2, short &beginListCount,
                                                            short &endListCount) const {
    std::string str1 = str;

    beginListCount = (short)std::count(str1.begin(), str1.end(), '(');
    endListCount = (short)std::count(str1.begin(), str1.end(), ')');

    size_t beginListPos = str1.find_first_of('(');
    str1.erase(0, beginListPos + 1);

    size_t firstSpacePos;
    removeBeginSpaceFluent(str1, firstSpacePos);

    std::string str2 = str1.substr(0, firstSpacePos);
    integer index = atoi(str2.c_str());

    /*!\brief zone id.*/
    beginListPos = str1.find_first_of('(');
    str1.erase(0, beginListPos + 1);
    removeBeginSpaceFluent(str1, firstSpacePos);
    str2 = str1.substr(0, firstSpacePos);
    str1.erase(0, firstSpacePos + 1);
    zoneId = (integer)strtol(str2.c_str(), NULL, 16);

    /*!\brief First id.*/
    removeBeginSpaceFluent(str1, firstSpacePos);
    str2 = str1.substr(0, firstSpacePos);
    str1.erase(0, firstSpacePos + 1);
    fi = (integer)strtol(str2.c_str(), NULL, 16);

    /*!\brief last id.*/
    removeBeginSpaceFluent(str1, firstSpacePos);
    str2 = str1.substr(0, firstSpacePos);
    str1.erase(0, firstSpacePos + 1);
    li = (integer)strtol(str2.c_str(), NULL, 16);
    if (zoneId == 0) {
        type = 0;
        type2 = 0;
        return 0;
    }

    /*!\brief bcType.*/
    removeBeginSpaceFluent(str1, firstSpacePos);
    str2 = str1.substr(0, firstSpacePos);
    str1.erase(0, firstSpacePos + 1);
    type = (integer)strtol(str2.c_str(), NULL, 16);

    /*!\brief type2.*/
    removeBeginSpaceFluent(str1, firstSpacePos);
    str2 = str1.substr(0, firstSpacePos);
    str1.erase(0, firstSpacePos + 1);
    type2 = (integer)strtol(str2.c_str(), NULL, 16);
    return 1;
}

void OpenHurricane::FluentMeshRead::readCFluent() {
    if (hasBeenRead_) {
        LFatal("Attempt to read the origin mesh file again!");
    }
    std::ifstream myf(fileName_, std::ios::binary);
    FILE *fp;
    fp = fopen(fileName_.c_str(), "rb");
    char *buffer;
    std::streampos length = myf.seekg(0, std::ios::end).tellg();
    buffer = new char[static_cast<std::streamsize>(length) + 1];

    myf.close();

    size_t len = fread(buffer, 1, length, fp);
    fclose(fp);

    buffer[static_cast<std::streamsize>(length)] = '\0';
    std::string str(buffer);
    delete[] buffer;

    parsingFluent(str);

    str.clear();

    hasBeenRead_ = true;
}

void OpenHurricane::FluentMeshRead::parsingFluent(std::string &str) {
    integer totalNodes = 0;
    integer totalFaces = 0;
    integer totalCells = 0;
    size_t NLPos = 0;
    size_t LBPos = 0;
    while (NLPos != std::string::npos) {

        NLPos = str.find_first_of('\n');
        std::string str1 = str.substr(0, NLPos);
        str.erase(0, NLPos + 1);
        LBPos = NLPos + 1;

        size_t beginListOps = str1.find_first_of('(');
        size_t endListOps = str1.find_last_of(')');
        if (beginListOps != std::string::npos) {
            size_t iter = beginListOps + 1;
            size_t str1Len = str1.size();
            while (iter < str1Len) {
                int Index;

                std::string str2 = str1.substr(iter, str1Len - 1);

                Index = (integer)strtol(str2.c_str(), nullptr, 10);
                str2.clear();

                if (Index == sectionIndexFluent::XF_COMMENT) {
                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_COMMENT

                else if (Index == sectionIndexFluent::XF_DIMENSION) {
                    if (endListOps != std::string::npos) {
                        if ((endListOps - 1) < (iter + 1)) {
                            Pout << "Error: "
                                 << "  The dimension of the grid: " << fileName_.name().c_str()
                                 << " is missing." << std::endl
                                 << "  The contents are \"" << str1.c_str() << "\". " << std::endl;
                            HurMPI::abort(HurMPI::getComm(), EXIT_FAILURE);
                        }
                        str2 = str1.substr(iter + 1, endListOps - 1);
                        ND_ = atoi(str2.c_str());
                        if (ND_ != DIMENSIONSET) {
                            Pout << "Error: "
                                 << "  The dimension of the grid: " << fileName_.name().c_str()
                                 << " is " << ND_ << "." << std::endl
                                 << "  But the dimension of this program is " << DIMENSIONSET
                                 << ". " << std::endl;
                            HurMPI::abort(HurMPI::getComm(), EXIT_FAILURE);
                        }
                    } else {
                        Pout << "Error: "
                             << "  The discription of the dimension of the "
                                "grid: "
                             << fileName_.name().c_str() << " is wrong!" << std::endl
                             << "  \"" << str1.c_str() << "\"" << std::endl;
                        HurMPI::abort(HurMPI::getComm(), EXIT_FAILURE);
                    }
                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_DIMENSION

                else if (Index == sectionIndexFluent::XF_HEADER) {
                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_HEADER

                else if (Index == sectionIndexFluent::XF_NODE) {
                    //size_t secondBeginListOps = str1.find_first_of('(');

                    integer zoneId;
                    integer fI = 0;
                    integer lI = 0;
                    short type;
                    short ND;
                    short bgListCounts = 0;
                    short endListCounts = 0;
                    parsingHeaderStrFluent(str1, zoneId, fI, lI, type, ND, bgListCounts,
                                           endListCounts);
                    if (zoneId == 0) {
                        globalNodes_ = lI;
                        //points_.clear();
                        points_.resize(globalNodes_);
                        break;
                    } else {
                        if (ND_ == 0) {
                            ND_ = ND;
                        }
                        if (ND_ != ND) {
                            std::string errMsg = "The dimension of grid: \"";
                            errMsg += fileName_.c_str();
                            errMsg += "\" is ambiguous. In the dimension section: " +
                                      std::to_string(ND_) + "\n";
                            errMsg += "In the nodes section: " + std::to_string(ND);
                            errorAbortStr(errMsg);
                        }
                        pointZone pZ(zoneId, fI - 1, lI - 1, ND, type);
                        pZ.resetName("points");
                        pointZones_.append(pZ);
                        totalNodes += (lI - fI + 1);

                        if (bgListCounts == 2) {
                            LBPos = str.find_first_of('(', 0);
                            NLPos = str.find_first_of(')', 0);
                            str1 = str.substr(LBPos + 1, NLPos - LBPos);
                        } else {
                            NLPos = str.find_first_of(')', 0);
                            str1 = str.substr(0, NLPos);
                        }
                        str.erase(0, NLPos + 1);
                        /*replaceAllMarks(str1, "\n", " ");
                        replaceAllMarks(str1, "\r", " ");*/
                        size_t firstS;
                        removeBeginSpaceFluent(str1, firstS);
                        readNodesCoordinatesFluent(str1, fI, lI, ND);
                    }
                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_NODE

                else if (Index == sectionIndexFluent::XF_FACE) {

                    integer zoneId;
                    integer fI = 0;
                    integer lI = 0;
                    short bcType;
                    short faceType;
                    short bgListCounts = 0;
                    short endListCounts = 0;
                    parsingHeaderStrFluent(str1, zoneId, fI, lI, bcType, faceType, bgListCounts,
                                           endListCounts);
                    if (zoneId == 0) {
                        globalFaces_ = lI;
                        //faces_.clear();
                        faces_.resize(globalFaces_);
                        break;
                    } else {
                        faceZone fZ(zoneId, fI - 1, lI - 1, bcType, faceType);
                        fZ.resetName("faces");
                        faceZones_.append(fZ);
                        totalFaces += (lI - fI + 1);

                        if (bcType == faceBCType::bcTypes::INTERIOR) {
                            globalInterior_ += (lI - fI + 1);
                        } else if (bcType == faceBCType::bcTypes::INTERFACE) {
                            if (!hasInterface_) {
                                hasInterface_ = true;
                            }
                        }

                        if (bgListCounts == 2) {
                            LBPos = str.find_first_of('(', 0);
                            NLPos = str.find_first_of(')', 0);
                            str1 = str.substr(LBPos + 1, NLPos - LBPos);
                        } else {
                            NLPos = str.find_first_of(')', 0);
                            str1 = str.substr(0, NLPos);
                        }
                        size_t firstS;
                        removeBeginSpaceFluent(str1, firstS);
                        readFaceConnectionNewFluent(fZ, str1, fI, lI, bcType);
                        str.erase(0, NLPos + 1);
                    }

                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_FACE

                else if ((Index == sectionIndexFluent::XF_RP_TV) ||
                         (Index == sectionIndexFluent::XF_RP_TV1)) {
                    int countBL = countMarks(str1, "(");
                    int countEL = countMarks(str1, ")");
                    replaceAllMarks(str1, "\r", " ");
                    replaceAllMarks(str1, "(", " ");
                    replaceAllMarks(str1, ")", " ");
                    char *endptr;
                    Index = (integer)strtol(str1.c_str(), &endptr, 10);
                    integer zID = (integer)strtol(endptr, &endptr, 10);
                    char *zName;
                    zName = strtok(endptr, " ");
                    zName = strtok(NULL, " ");
                    bool hasSet = false;
                    if (zName != NULL) {
                        for (integer j = 0; j < pointZones_.size(); j++) {
                            if (zID == pointZones_[j].index()) {
                                pointZones_[j].resetName(zName);
                                hasSet = true;
                            }
                        }

                        if (!hasSet) {
                            for (integer j = 0; j < faceZones_.size(); j++) {
                                if (zID == faceZones_[j].index()) {
                                    faceZones_[j].resetName(zName);
                                    hasSet = true;
                                }
                            }
                        }
                        if (!hasSet) {
                            for (integer j = 0; j < cellZones_.size(); j++) {
                                if (zID == cellZones_[j].index()) {
                                    cellZones_[j].resetName(zName);
                                }
                            }
                        }
                    }
                    if (countBL != countEL) {
                        size_t index0 = 0;
                        if (countBL != 3) {
                            index0 = str.find("(");
                            if (index0 == std::string::npos) {
                                LFatal("The format of mesh file in zone sections is wrong.");
                            }
                        }
                        size_t index1 = str.find("(", index0);
                        if (index1 == std::string::npos) {
                            LFatal("The format of mesh file in zone sections is wrong.");
                        }
                        size_t index2 = str.find(")", index1);
                        if (index2 == std::string::npos) {
                            LFatal("The format of mesh file in zone sections is wrong.");
                        }
                        std::string zstr = str.substr(index1, index2);
                        size_t indexfound = zstr.find("motion-spec .");
                        if (indexfound != std::string::npos) {
                            str.erase(0, index2 + 1);
                            index1 = str.find("(");
                            index2 = str.find(")", index1);
                            size_t indexXo = str.find("x-origin .");
                            char *endptrr;
                            real xo;
                            if (indexXo != std::string::npos) {
                                auto strx = str.substr(indexXo + 11, index2 - (indexXo + 11));
                                xo = (real)strtod(strx.c_str(), &endptrr);
                                str.erase(0, index2 + 1);
                            }

                            index1 = str.find("(");
                            index2 = str.find(")", index1);
                            size_t indexYo = str.find("y-origin .");
                            real yo;
                            if (indexYo != std::string::npos) {
                                auto strx = str.substr(indexYo + 11, index2 - (indexYo + 11));
                                yo = (real)strtod(strx.c_str(), &endptrr);
                                str.erase(0, index2 + 1);
                            }
                            index1 = str.find("(");
                            index2 = str.find(")", index1);
                            size_t indexZo = str.find("z-origin .");
                            real zo;
                            if (indexZo != std::string::npos) {
                                auto strx = str.substr(indexZo + 11, index2 - (indexZo + 11));
                                zo = (real)strtod(strx.c_str(), &endptrr);
                                str.erase(0, index2 + 1);
                            }
                            index1 = str.find("(");
                            index2 = str.find(")", index1);
                            size_t indexAi = str.find("ai .");
                            real xai;
                            if (indexAi != std::string::npos) {
                                auto strx = str.substr(indexAi + 5, index2 - (indexAi + 5));
                                xai = (real)strtod(strx.c_str(), &endptrr);
                                str.erase(0, index2 + 1);
                            }
                            index1 = str.find("(");
                            index2 = str.find(")", index1);
                            size_t indexAj = str.find("aj .");
                            real xaj;
                            if (indexAj != std::string::npos) {
                                auto strx = str.substr(indexAj + 5, index2 - (indexAj + 5));
                                xaj = (real)strtod(strx.c_str(), &endptrr);
                                str.erase(0, index2 + 1);
                            }
                            index1 = str.find("(");
                            index2 = str.find(")", index1);
                            size_t indexAk = str.find("ak .");
                            real xak;
                            if (indexAk != std::string::npos) {
                                auto strx = str.substr(indexAk + 5, index2 - (indexAk + 5));
                                xak = (real)strtod(strx.c_str(), &endptrr);
                                str.erase(0, index2 + 1);
                            }

                            origin_ = vector(xo, yo, zo);
                            axis_ = vector(xai, xaj, xak);
                        }
                        indexfound = zstr.find("angular? . #t");
                        if (indexfound != std::string::npos) {
                            for (integer ppi = 0; ppi < periodicPairZone_.size(); ++ppi) {
                                if (periodicPairZone_[ppi].periodicZone() == zID ||
                                    periodicPairZone_[ppi].shadowZone() == zID) {
                                    periodicPairZone_[ppi].setType(periodicTypes::ROTATIONAL);
                                    break;
                                }
                            }
                        }
                        indexfound = zstr.find("angular? . #f");
                        if (indexfound != std::string::npos) {
                            for (integer ppi = 0; ppi < periodicPairZone_.size(); ++ppi) {
                                if (periodicPairZone_[ppi].periodicZone() == zID ||
                                    periodicPairZone_[ppi].shadowZone() == zID) {
                                    periodicPairZone_[ppi].setType(periodicTypes::TRANSLATIONAL);
                                    break;
                                }
                            }
                        }
                    }
                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_RP_TV

                else if (Index == sectionIndexFluent::XF_PERIODIC_FACE) {
                    integer pZoneId, pSZoneId;
                    integer fI = 0;
                    integer lI = 0;
                    short bgListCounts = 0;
                    //short endListCounts = 0;
                    bgListCounts = (short)std::count(str1.begin(), str1.end(), '(');
                    //endListCounts = (short)std::count(str1.begin(), str1.end(), ')');

                    replaceAllMarks(str1, "(", " ");
                    replaceAllMarks(str1, ")", " ");

                    char *endptr;
                    Index = (integer)strtol(str1.c_str(), &endptr, 10);
                    fI = (integer)strtol(endptr, &endptr, 16);
                    lI = (integer)strtol(endptr, &endptr, 16);
                    pZoneId = (integer)strtol(endptr, &endptr, 16);
                    pSZoneId = (integer)strtol(endptr, &endptr, 16);

                    if (bgListCounts == 2) {
                        LBPos = str.find_first_of('(', 0);
                        NLPos = str.find_first_of(')', 0);
                        str1 = str.substr(LBPos + 1, NLPos - LBPos);
                    } else {
                        NLPos = str.find_first_of(')', 0);
                        str1 = str.substr(0, NLPos);
                    }
                    str.erase(0, NLPos + 1);

                    periodicPair pp(pZoneId, pSZoneId, fI, lI);
                    pp(fI).x() = (integer)strtol(str1.c_str(), &endptr, 16);
                    pp(fI).y() = (integer)strtol(endptr, &endptr, 16);

                    for (integer i = fI + 1; i <= lI; i++) {
                        pp(i).x() = (integer)strtol(endptr, &endptr, 16);
                        pp(i).y() = (integer)strtol(endptr, &endptr, 16);
                    }
                    for (integer i = fI; i <= lI; i++) {
                        pp(i).x() -= 1;
                        pp(i).y() -= 1;
                    }
                    periodicPairZone_.append(pp);
                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_PERIODIC_FACE

                else if (Index == sectionIndexFluent::XF_CELL) {
                    //size_t secondBeginListOps = str1.find_first_of('(');
                    integer zoneId;
                    integer fI = 0;
                    integer lI = 0;
                    short type;
                    short elementType;
                    short bgListCounts = 0;
                    short endListCounts = 0;
                    parsingHeaderStrFluent(str1, zoneId, fI, lI, type, elementType, bgListCounts,
                                           endListCounts);
                    if (zoneId == 0) {
                        globalCells_ = lI;
                        //cells_.clear();
                        cells_.resize(globalCells_);
                        break;
                    } else {
                        cellZone cZ(zoneId, fI - 1, lI - 1, type, elementType);
                        cZ.resetName("cells");
                        cellZones_.append(cZ);
                        totalCells += (lI - fI + 1);
                        if (cZ.isMixed()) {
                            if (bgListCounts == 2) {
                                LBPos = str.find_first_of('(', 0);
                                NLPos = str.find_first_of(')', 0);
                                str1 = str.substr(LBPos, NLPos - LBPos);
                            } else {
                                NLPos = str.find_first_of(')', 0);
                                str1 = str.substr(0, NLPos);
                            }
                            str.erase(0, NLPos + 1);
                            replaceAllMarks(str1, "(", " ");
                            replaceAllMarks(str1, ")", " ");
                            /*replaceAllMarks(str1, "\n", " ");
                            replaceAllMarks(str1, "\r", " ");*/

                            readMixedCellShapeFluent(str1, fI, lI);
                        } else {
                            for (integer id = fI - 1; id <= lI - 1; id++) {
                                cells_[id].setType(elementType);
                                cells_[id].setOrignIndex(id);
                            }
                        }
                    }
                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_CELL

                else if (Index == sectionIndexFluent::XF_FACE_TREE) {
                    LFatal("The face tree not supported yet");
                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_FACE_TREE

                else if (Index == sectionIndexFluent::XF_CELL_TREE) {
                    LFatal("The cell tree not supported yet");
                    break; //break current loop: while(iter < str1Len)
                }          // End if for XF_CELL_TREE

                if (sectionIndexFluentMap_.find(Index) == sectionIndexFluentMap_.end()) {
                    replaceAllMarks(str1, "\n", " ");
                    replaceAllMarks(str1, "\r", " ");
                    std::string errMsg = " The discriptiom: \"";
                    errMsg += str1.c_str();
                    errMsg += "\"is unknown.";
                    errorAbortStr(errMsg);
                }
            }
        } else {
            continue;
        }
    }

    /*!\brief To check whether the total number from the declaration section is equal to the sum of the number of all regular sections.*/

    /*!\brief Checking faces.*/
    if (totalFaces < globalFaces_) {
        globalFaces_ = totalFaces;
        faces_.resize(globalFaces_);
        if (report) {
            LWarning("Parsing mesh: The total number from the declaration section is not "
                     "equal to the sum of the number of all regular sections. nFaces_1 = %d, "
                     "nFaces_2 = %d",
                     globalFaces_, totalFaces);
        }
    }

    /*!\brief Checking nodes.*/
    if (totalNodes < globalNodes_) {
        LFatal("The total number from the declaration is not equal to that "
               "from the sum of the all regular sections. \n"
               " nPoints_from_decla = %d, nPoints_from_regu = %d",
               totalNodes, globalNodes_);
    }

    /*!\brief Checking cell.*/
    if (totalCells < globalCells_) {
        std::string errMsg;
        errMsg = "The total number from the declaration is not equal to that "
                 "from the sum of the all regular sections. \n";
        errMsg += " nCells_from_decla = " + std::to_string(totalCells);
        errMsg += ".  nCells_from_regu = " + std::to_string(globalCells_);
        errorAbortStr(errMsg);
    }
}

void OpenHurricane::FluentMeshRead::parsingFluentBinary() {
    std::ifstream myin(fileName_, std::ios::binary);
    if (!myin.good()) {
        LFatal("File: \"%s\" is not exist", fileName_.c_str());
    }
    std::string line;
    std::regex myReg("\\s*\\(\\s*([0-9]+).*");

    // node index
    std::regex myReg1("(|20|30)10");

    // cell index
    std::regex myReg2("(|20|30)12");

    // face index
    std::regex myReg3("(|20|30)13");

    // periodic face pair
    std::regex myReg4("(|20|30)18");

    // face tree
    std::regex myReg5("(|20|30)59");

    // cell tree
    std::regex myReg6("(|20|30)58");
    bool isBinary = false;
    integer count = 0;
    std::stringstream sstr;
    integer countNode = 0;
    integer countCell = 0;
    integer countFace = 0;
    while (std::getline(myin, line)) {
        replaceAllMarks(line, "\r", "");
        replaceAllMarks(line, "\n", "");
        std::smatch what;
        if (std::regex_match(line, what, myReg)) {
            std::string index = what[1];
            std::smatch what1;
            if (index == "0") {
                // Comment.
                continue;
            } else if (index == "1") {
                // Header.
                continue;
            } else if (index == "2") {
                // Dimensionality.
                replaceAllMarks(line, "(", "");
                replaceAllMarks(line, ")", "");
                sstr.str(line);
                integer ii;
                sstr >> ii >> ND_;
                if (ND_ != DIMENSIONSET) {
                    std::string errMsg = "The dimension of the grid: ";
                    errMsg += fileName_.name();
                    errMsg += " is ";
                    errMsg += toString(ND_);
                    errMsg += ". But the dimension of this program is";
                    errMsg += toString(integer(DIMENSIONSET));
                    errorAbortStr(errMsg);
                }
                continue;
            } else if (std::regex_match(index, what1, myReg1)) {
                parsingFluentPointBinary(myin, line, index, countNode);
            } else if (std::regex_match(index, what1, myReg2)) {
                // cells
                // (2012 (zone-id first-index last-index type element-type))
                parsingFluentCellBinary(myin, line, index, countCell);
            } else if (std::regex_match(index, what1, myReg3)) {
                parsingFluentFaceBinary(myin, line, index, countFace);
            } else if (index == "39" || index == "45") {
                parsingFluenZoneSectionBinary(myin, line);
            } else if (std::regex_match(index, what1, myReg4)) {
                parsingFluentPeriodicFaceBinary(myin, line, index);
            } else if (std::regex_match(index, what1, myReg5)) {
                LFatal("The face tree not supported yet");
            } else if (std::regex_match(index, what1, myReg6)) {
                LFatal("The cell tree not supported yet");
            }
        } else {
            continue;
        }
    }
    /*!\brief Checking faces.*/
    if (countFace < globalFaces_) {
        globalFaces_ = countFace;
        faces_.resize(globalFaces_);
        if (report) {
            LWarning("Parsing mesh: The total number from the declaration section is not "
                     "equal to the sum of the number of all regular sections. nFaces_1 = %d, "
                     "nFaces_2 = %d",
                     globalFaces_, countFace);
        }
    }

    /*!\brief Checking nodes.*/
    if (countNode < globalNodes_) {
        std::string errMsg;
        errMsg = "The total number from the declaration is not equal to that "
                 "from the sum of the all regular sections. \n";
        errMsg += " nPoints_from_decla = " + std::to_string(countNode);
        errMsg += ".  nPoints_from_regu = " + std::to_string(globalNodes_);
        errorAbortStr(errMsg);
    }

    /*!\brief Checking cell.*/
    if (countCell < globalCells_) {
        std::string errMsg;
        errMsg = "The total number from the declaration is not equal to that "
                 "from the sum of the all regular sections. \n";
        errMsg += " nCells_from_decla = " + std::to_string(countCell);
        errMsg += ".  nCells_from_regu = " + std::to_string(globalCells_);
        errorAbortStr(errMsg);
    }
    myin.close();
}

void OpenHurricane::FluentMeshRead::skipBrackets(std::ifstream &myin, const int nBrackets) const {
    integer nB = nBrackets;
    while (nB > 0) {
        char bChar[2];
        myin.read(bChar, 1);
        bChar[1] = '\0';
        if (bChar[0] == '(') {
            nB++;
        } else if (bChar[0] == ')') {
            nB--;
        }
    }
}

void OpenHurricane::FluentMeshRead::skipBrackets(std::ifstream &myin, std::string &str,
                                                 const int nBrackets) const {
    str = "";
    integer nB = nBrackets;
    while (nB > 0) {
        char bChar[2];
        myin.read(bChar, 1);
        bChar[1] = '\0';
        str += bChar[0];
        if (bChar[0] == '(') {
            nB++;
        } else if (bChar[0] == ')') {
            nB--;
        }
    }
}

void OpenHurricane::FluentMeshRead::skipTo(std::ifstream &myin, const char *c) const {
    char bChar[2];
    bChar[1] = '\0';
    while (bChar[0] != c[0]) {
        myin.read(bChar, 1);
    }
}

void OpenHurricane::FluentMeshRead::parsingFluentPointBinary(std::ifstream &myin, std::string &line,
                                                             std::string &index,
                                                             integer &countNode) {
    integer zoneId;
    integer fI = 0;
    integer lI = 0;
    short type;
    short ND;
    short bgListCounts = 0;
    short endListCounts = 0;
    parsingHeaderStrFluent(line, zoneId, fI, lI, type, ND, bgListCounts, endListCounts);

    if (zoneId == 0) {
        globalNodes_ = lI;
        //points_.clear();
        points_.resize(globalNodes_);
        return;
    } else {
        if (ND_ != ND) {
            std::string errMsg = "The dimension of grid: \"";
            errMsg += fileName_.c_str();
            errMsg += "\" is ambiguous. In the dimension section: " + std::to_string(ND_) + "\n";
            errMsg += "In the nodes section: " + std::to_string(ND);
            errorAbortStr(errMsg);
        }
        if (lI < fI) {
            LFatal("The point index is wrong in file: %s", fileName_.c_str());
        }
        if (lI > globalNodes_ + 1) {
            LFatal("The point index is wrong in file: %s", fileName_.c_str());
        }
        pointZone pZ(zoneId, fI - 1, lI - 1, ND, type);
        pZ.resetName("points");
        pointZones_.append(pZ);
        integer sizeP = (lI - fI + 1);
        countNode += sizeP;
        if (bgListCounts == 2) {
            char lastChar[2];
            myin.read(lastChar, 1);
            lastChar[1] = '\0';
        }
        if (index == "2010") {
            float *px = new float[ND * sizeP];
            myin.read((char *)px, ND * sizeP * sizeof(float));
            integer countPi = 0;
            for (integer pi = pZ.firstIndex(); pi <= pZ.lastIndex(); ++pi) {
                points_[pi].x() = (real)px[countPi * 3 + 0];
                points_[pi].y() = (real)px[countPi * 3 + 1];
                points_[pi].z() = (real)px[countPi * 3 + 2];
                countPi++;
            }
            delete[] px;
        } else if (index == "3010") {
            double *px = new double[ND * sizeP];
            myin.read((char *)px, ND * sizeP * sizeof(double));
            integer countPi = 0;
            for (integer pi = pZ.firstIndex(); pi <= pZ.lastIndex(); ++pi) {
                points_[pi].x() = (real)px[countPi * 3 + 0];
                points_[pi].y() = (real)px[countPi * 3 + 1];
                points_[pi].z() = (real)px[countPi * 3 + 2];
                countPi++;
            }
            delete[] px;
        }
    }
    skipBrackets(myin, 2);
}

void OpenHurricane::FluentMeshRead::parsingFluentCellBinary(std::ifstream &myin, std::string &line,
                                                            std::string &index,
                                                            integer &totalCells) {
    integer zoneId;
    integer fI = 0;
    integer lI = 0;
    short type;
    short elementType;
    short bgListCounts = 0;
    short endListCounts = 0;
    parsingHeaderStrFluent(line, zoneId, fI, lI, type, elementType, bgListCounts, endListCounts);
    if (zoneId == 0) {
        globalCells_ = lI;
        //cells_.clear();
        cells_.resize(globalCells_);
        return;
    } else {
        if (type == 0) {
            // dead zone
            return;
        }
        cellZone cZ(zoneId, fI - 1, lI - 1, type, elementType);
        cZ.resetName("cells");
        cellZones_.append(cZ);
        integer sizeC = (lI - fI + 1);
        totalCells += sizeC;

        if (cZ.isMixed()) {
            if (bgListCounts == 2) {
                skipTo(myin, "(");
            }

            if (index == "2012") {
                int32_t *cc = new int32_t[sizeC];
                myin.read((char *)cc, sizeC * sizeof(int32_t));
                integer countCi = 0;
                for (integer i = fI - 1; i <= lI - 1; i++) {
                    cells_[i].setType((integer)cc[countCi]);
                    cells_[i].setOrignIndex(i);
                    countCi++;
                }
                delete[] cc;
            } else if (index == "3012") {
                int32_t *cc = new int32_t[sizeC];
                myin.read((char *)cc, sizeC * sizeof(int32_t));
                integer countCi = 0;
                for (integer i = fI - 1; i <= lI - 1; i++) {
                    cells_[i].setType((integer)cc[countCi]);
                    cells_[i].setOrignIndex(i);
                    countCi++;
                }
                delete[] cc;
            }

            skipBrackets(myin, 2);
        } else {
            for (integer id = fI - 1; id <= lI - 1; id++) {
                cells_[id].setType(elementType);
                cells_[id].setOrignIndex(id);
            }
            if (bgListCounts == 2 && bgListCounts != endListCounts) {
                skipBrackets(myin, 1);
            }
        }
    }
}

void OpenHurricane::FluentMeshRead::parsingFluentFaceBinary(std::ifstream &myin, std::string &line,
                                                            std::string &index,
                                                            integer &totalFaces) {
    integer zoneId;
    integer fI = 0;
    integer lI = 0;
    short bcType;
    short faceType;
    short bgListCounts = 0;
    short endListCounts = 0;
    parsingHeaderStrFluent(line, zoneId, fI, lI, bcType, faceType, bgListCounts, endListCounts);
    if (zoneId == 0) {
        globalFaces_ = lI;
        //faces_.clear();
        faces_.resize(globalFaces_);
        return;
    } else {
        faceZone fZ(zoneId, fI - 1, lI - 1, bcType, faceType);
        fZ.resetName("faces");
        faceZones_.append(fZ);
        integer sizeF = (lI - fI + 1);
        totalFaces += sizeF;

        if (bcType == faceBCType::bcTypes::INTERIOR) {
            globalInterior_ += (lI - fI + 1);
        } else if (bcType == faceBCType::bcTypes::INTERFACE) {
            if (!hasInterface_) {
                hasInterface_ = true;
            }
        }

        if (bgListCounts == 2) {
            skipTo(myin, "(");
        }

        if (fZ.isLinear()) {
            for (integer i = fI - 1; i < lI; i++) {
                faces_[i].resize(2);
            }
            if (index == "2013") {
                int32_t *ff = new int32_t[4 * sizeF];
                myin.read((char *)ff, 4 * sizeF * sizeof(int32_t));
                integer countCi = 0;
                for (integer i = fI - 1; i <= lI - 1; i++) {
                    faces_[i][0] = (integer)ff[countCi * 4 + 0];
                    faces_[i][1] = (integer)ff[countCi * 4 + 1];
                    faces_[i].leftCell() = (integer)ff[countCi * 4 + 2];
                    faces_[i].rightCell() = (integer)ff[countCi * 4 + 3];
                    countCi++;
                }
                delete[] ff;
            } else if (index == "3013") {
                int32_t *ff = new int32_t[4 * sizeF];
                myin.read((char *)ff, 4 * sizeF * sizeof(int32_t));
                integer countCi = 0;
                for (integer i = fI - 1; i <= lI - 1; i++) {
                    faces_[i][0] = (integer)ff[countCi * 4 + 0];
                    faces_[i][1] = (integer)ff[countCi * 4 + 1];
                    faces_[i].leftCell() = (integer)ff[countCi * 4 + 2];
                    faces_[i].rightCell() = (integer)ff[countCi * 4 + 3];
                    countCi++;
                }
                delete[] ff;
            }
        } else if (fZ.isTriangular()) {
            for (integer i = fI - 1; i < lI; i++) {
                faces_[i].resize(3);
            }
            if (index == "2013") {
                int32_t *ff = new int32_t[5 * sizeF];
                myin.read((char *)ff, 5 * sizeF * sizeof(int32_t));
                integer countCi = 0;
                for (integer i = fI - 1; i <= lI - 1; i++) {
                    faces_[i][0] = (integer)ff[countCi * 5 + 0];
                    faces_[i][1] = (integer)ff[countCi * 5 + 1];
                    faces_[i][2] = (integer)ff[countCi * 5 + 2];
                    faces_[i].leftCell() = (integer)ff[countCi * 5 + 3];
                    faces_[i].rightCell() = (integer)ff[countCi * 5 + 4];
                    countCi++;
                }
                delete[] ff;
            } else if (index == "3013") {
                int32_t *ff = new int32_t[5 * sizeF];
                myin.read((char *)ff, 5 * sizeF * sizeof(int32_t));
                integer countCi = 0;
                for (integer i = fI - 1; i <= lI - 1; i++) {
                    faces_[i][0] = (integer)ff[countCi * 5 + 0];
                    faces_[i][1] = (integer)ff[countCi * 5 + 1];
                    faces_[i][2] = (integer)ff[countCi * 5 + 2];
                    faces_[i].leftCell() = (integer)ff[countCi * 5 + 3];
                    faces_[i].rightCell() = (integer)ff[countCi * 5 + 4];
                    countCi++;
                }
                delete[] ff;
            }
        } else if (fZ.isQuadrilateral()) {
            for (integer i = fI - 1; i < lI; i++) {
                faces_[i].resize(4);
            }
            if (index == "2013") {
                int32_t *ff = new int32_t[6 * sizeF];
                myin.read((char *)ff, 6 * sizeF * sizeof(int32_t));
                integer countCi = 0;
                for (integer i = fI - 1; i <= lI - 1; i++) {
                    faces_[i][0] = (integer)ff[countCi * 6 + 0];
                    faces_[i][1] = (integer)ff[countCi * 6 + 1];
                    faces_[i][2] = (integer)ff[countCi * 6 + 2];
                    faces_[i][3] = (integer)ff[countCi * 6 + 3];
                    faces_[i].leftCell() = (integer)ff[countCi * 6 + 4];
                    faces_[i].rightCell() = (integer)ff[countCi * 6 + 5];
                    countCi++;
                }
                delete[] ff;
            } else if (index == "3013") {
                int32_t *ff = new int32_t[6 * sizeF];
                myin.read((char *)ff, 6 * sizeF * sizeof(int32_t));
                integer countCi = 0;
                for (integer i = fI - 1; i <= lI - 1; i++) {
                    faces_[i][0] = (integer)ff[countCi * 6 + 0];
                    faces_[i][1] = (integer)ff[countCi * 6 + 1];
                    faces_[i][2] = (integer)ff[countCi * 6 + 2];
                    faces_[i][3] = (integer)ff[countCi * 6 + 3];
                    faces_[i].leftCell() = (integer)ff[countCi * 6 + 4];
                    faces_[i].rightCell() = (integer)ff[countCi * 6 + 5];
                    countCi++;
                }
                delete[] ff;
            }
        } else if (fZ.isMixed() || fZ.isPolygonal()) {
            for (integer i = fI - 1; i <= lI - 1; i++) {
                integer x;
                if (index == "2013") {
                    int32_t ff;
                    myin.read((char *)(&ff), sizeof(int32_t));
                    x = (integer)ff;

                } else if (index == "3013") {
                    int32_t ff;
                    myin.read((char *)(&ff), sizeof(int32_t));
                    x = (integer)ff;
                }
                faces_[i].resize(x);
                if (index == "2013") {
                    int32_t *ff = new int32_t[x + 2];
                    myin.read((char *)(ff), (x + 2) * sizeof(int32_t));
                    for (integer j = 0; j < x; j++) {
                        faces_[i][j] = (integer)ff[j];
                    }
                    faces_[i].leftCell() = (integer)ff[x];
                    faces_[i].rightCell() = (integer)ff[x + 1];
                    delete[] ff;
                } else if (index == "3013") {
                    int32_t *ff = new int32_t[x + 2];
                    myin.read((char *)(ff), (x + 2) * sizeof(int32_t));
                    for (integer j = 0; j < x; j++) {
                        faces_[i][j] = (integer)ff[j];
                    }
                    faces_[i].leftCell() = (integer)ff[x];
                    faces_[i].rightCell() = (integer)ff[x + 1];
                    delete[] ff;
                }
            }
        }

        // 因为在C/C++中数组坐标是从0开始的，所以face对应的points和cells的索引需要左移一位
        for (integer i = fI - 1; i < lI; i++) {
            for (integer j = 0; j < faces_[i].size(); j++) {
                faces_[i][j] -= 1;
            }
            faces_[i].leftCell() -= 1;
            faces_[i].rightCell() -= 1;
            faces_[i].setBCType(bcType);
        }
        skipBrackets(myin, 2);
    }
}

void OpenHurricane::FluentMeshRead::parsingFluenZoneSectionBinary(std::ifstream &myin,
                                                                  std::string &line) {
    const auto beginListCount = (short)std::count(line.begin(), line.end(), '(');
    const auto endListCount = (short)std::count(line.begin(), line.end(), ')');

    replaceAllMarks(line, "(", " ");
    replaceAllMarks(line, ")", " ");
    stdStringList tk;
    split(line, tk, " ");
    char *endptr = nullptr;
    auto zID = (integer)strtol(tk[1].c_str(), &endptr, 10);

    std::string zName;
    if (tk.size() == 3) {
        zName = tk[2];
    } else if (tk.size() >= 4) {
        zName = tk[3];
    } else {
        LFatal("The zone sections are wrong in file: %s", fileName_.c_str());
    }

    bool hasSet = false;

    for (integer j = 0; j < pointZones_.size(); j++) {
        if (zID == pointZones_[j].index()) {
            pointZones_[j].resetName(zName);
            hasSet = true;
        }
    }

    if (!hasSet) {
        for (integer j = 0; j < faceZones_.size(); j++) {
            if (zID == faceZones_[j].index()) {
                faceZones_[j].resetName(zName);
                hasSet = true;
            }
        }
    }
    if (!hasSet) {
        for (integer j = 0; j < cellZones_.size(); j++) {
            if (zID == cellZones_[j].index()) {
                cellZones_[j].resetName(zName);
                hasSet = true;
            }
        }
    }

    if (!hasSet) {
        checkWarningStr(("Cannot find zone: " + zName + " in zone list"));
    }

    if (beginListCount == endListCount && beginListCount == 3) {
        return;
    }

    std::string zoneInfoStr;
    skipBrackets(myin, zoneInfoStr, 2);
    auto indexfound = zoneInfoStr.find("motion-spec .");
    if (indexfound != std::string::npos) {
        auto index2 = zoneInfoStr.find(")", indexfound + 14 - 1);
        auto msstr = zoneInfoStr.substr(indexfound + 14 - 1, index2 - (indexfound + 14) + 1);
        std::stringstream ssstr;
        ssstr.str(msstr);
        real motionSpec;
        ssstr >> motionSpec;
    }
    auto indexXo = zoneInfoStr.find("x-origin .");
    if (indexXo != std::string::npos) {
        auto index2 = zoneInfoStr.find(")", indexXo + 11 - 1);
        auto msstr = zoneInfoStr.substr(indexXo + 11 - 1, index2 - (indexXo + 11) + 1);
        std::stringstream ssstr;
        ssstr.str(msstr);
        real xo;
        ssstr >> xo;
        origin_.x() = xo;
    }
    auto indexYo = zoneInfoStr.find("y-origin .");
    if (indexYo != std::string::npos) {
        auto index2 = zoneInfoStr.find(")", indexYo + 11 - 1);
        auto msstr = zoneInfoStr.substr(indexYo + 11 - 1, index2 - (indexYo + 11) + 1);
        std::stringstream ssstr;
        ssstr.str(msstr);
        real yo;
        ssstr >> yo;
        origin_.y() = yo;
    }
    auto indexZo = zoneInfoStr.find("z-origin .");
    if (indexZo != std::string::npos) {
        auto index2 = zoneInfoStr.find(")", indexZo + 11 - 1);
        auto msstr = zoneInfoStr.substr(indexZo + 11 - 1, index2 - (indexZo + 11) + 1);
        std::stringstream ssstr;
        ssstr.str(msstr);
        real zo;
        ssstr >> zo;
        origin_.z() = zo;
    }
    auto indexAi = zoneInfoStr.find("ai .");
    if (indexAi != std::string::npos) {
        auto index2 = zoneInfoStr.find(")", indexAi + 5 - 1);
        auto msstr = zoneInfoStr.substr(indexAi + 5 - 1, index2 - (indexAi + 5) + 1);
        std::stringstream ssstr;
        ssstr.str(msstr);
        real ai;
        ssstr >> ai;
        axis_.x() = ai;
    }
    auto indexAj = zoneInfoStr.find("aj .");
    if (indexAj != std::string::npos) {
        auto index2 = zoneInfoStr.find(")", indexAj + 5 - 1);
        auto msstr = zoneInfoStr.substr(indexAj + 5 - 1, index2 - (indexAj + 5) + 1);
        std::stringstream ssstr;
        ssstr.str(msstr);
        real aj;
        ssstr >> aj;
        axis_.y() = aj;
    }
    auto indexAk = zoneInfoStr.find("ak .");
    if (indexAk != std::string::npos) {
        auto index2 = zoneInfoStr.find(")", indexAk + 5 - 1);
        auto msstr = zoneInfoStr.substr(indexAk + 5 - 1, index2 - (indexAk + 5) + 1);
        std::stringstream ssstr;
        ssstr.str(msstr);
        real ak;
        ssstr >> ak;
        axis_.z() = ak;
    }
    auto indexZt = zoneInfoStr.find("angular? . #t");
    if (indexZt != std::string::npos) {
        for (integer ppi = 0; ppi < periodicPairZone_.size(); ++ppi) {
            if (periodicPairZone_[ppi].periodicZone() == zID ||
                periodicPairZone_[ppi].shadowZone() == zID) {
                periodicPairZone_[ppi].setType(periodicTypes::ROTATIONAL);
                break;
            }
        }
    }
    auto indexZf = zoneInfoStr.find("angular? . #f");
    if (indexZf != std::string::npos) {
        for (integer ppi = 0; ppi < periodicPairZone_.size(); ++ppi) {
            if (periodicPairZone_[ppi].periodicZone() == zID ||
                periodicPairZone_[ppi].shadowZone() == zID) {
                periodicPairZone_[ppi].setType(periodicTypes::TRANSLATIONAL);
                break;
            }
        }
    }
}

void OpenHurricane::FluentMeshRead::parsingFluentPeriodicFaceBinary(std::ifstream &myin,
                                                                    std::string &line,
                                                                    std::string &index) {
    integer pZoneId, pSZoneId;
    integer fI = 0;
    integer lI = 0;
    short bgListCounts = 0;
    //short endListCounts = 0;
    bgListCounts = (short)std::count(line.begin(), line.end(), '(');
    //endListCounts = (short)std::count(line.begin(), line.end(), ')');

    replaceAllMarks(line, "(", " ");
    replaceAllMarks(line, ")", " ");

    char *endptr;
    auto Index = (integer)strtol(line.c_str(), &endptr, 10);
    fI = (integer)strtol(endptr, &endptr, 16);
    lI = (integer)strtol(endptr, &endptr, 16);
    pZoneId = (integer)strtol(endptr, &endptr, 16);
    pSZoneId = (integer)strtol(endptr, &endptr, 16);

    if (bgListCounts == 2) {
        skipTo(myin, "(");
    }

    integer sizeF = (lI - fI + 1);
    periodicPair pp(pZoneId, pSZoneId, fI, lI);
    if (index == "2018") {
        int32_t *ff = new int32_t[2 * sizeF];
        myin.read((char *)ff, 2 * sizeF * sizeof(int32_t));
        integer countCi = 0;
        for (integer i = fI; i <= lI; i++) {
            pp(i).x() = (integer)ff[countCi * 2 + 0];
            pp(i).y() = (integer)ff[countCi * 2 + 1];
            countCi++;
        }
        delete[] ff;
    } else if (index == "3018") {
        int32_t *ff = new int32_t[2 * sizeF];
        myin.read((char *)ff, 2 * sizeF * sizeof(int32_t));
        integer countCi = 0;
        for (integer i = fI; i <= lI; i++) {
            pp(i).x() = (integer)ff[countCi * 2 + 0];
            pp(i).y() = (integer)ff[countCi * 2 + 1];
            countCi++;
        }
        delete[] ff;
    }
    for (integer i = fI; i <= lI; i++) {
        pp(i).x() -= 1;
        pp(i).y() -= 1;
    }
    periodicPairZone_.append(pp);
    skipBrackets(myin, 2);
}

void OpenHurricane::FluentMeshRead::readFluent() {
    if (HurMPI::master()) {
        bool isBinary = false;
        if (hasBeenRead_) {
            std::string errMsg = "Attempt to read the origin mesh file again!";
            errorAbortStr(errMsg);
        }
        if (!fileName_.empty()) {
            isBinary = isFluentMeshHasBinaryFormat();
            if (!isBinary) {
                meshStr_ = readFileToString(fileName_);
            }
        } else if (meshStr_.empty()) {
            LFatal("Attempt to read the null mesh string");
        }
        if (isBinary) {
            parsingFluentBinary();
        } else {
            parsingFluent(meshStr_);
        }

        removeUnreferencedFaces();

        FluentMeshRead::formingFaces(cells_.size(), faces_, cells_);
        FluentMeshRead::formingNodes(cells_.size(), faces_, cells_);
        //computingCellWall();

        /*decomposingMeshByMetis dMM(cells_, faces_, globalInterior_, 2);

        dMM.decomposing();*/

        meshStr_.clear();

        hasBeenRead_ = true;

        printMeshInfo();

        originCellIndex_.resize(cellZones_.size());
        for (integer czi = 0; czi < cellZones_.size(); ++czi) {
            originCellIndex_[czi].resize(cellZones_[czi].size());

            integer countCi = 0;
            for (integer ci = cellZones_[czi].firstIndex(); ci <= cellZones_[czi].lastIndex();
                 ++ci) {
                originCellIndex_[czi][countCi] = ci;
                countCi++;
            }
        }
    }
}

OpenHurricane::FluentMeshRead::FluentMeshRead() : originMeshRead() {
    createSectionIndexFluentMap();
}

OpenHurricane::FluentMeshRead::FluentMeshRead(const fileName &fN, const int nP)
    : originMeshRead(fN, nP) {
    createSectionIndexFluentMap();
}

OpenHurricane::FluentMeshRead::FluentMeshRead(const std::string &str, const int nP)
    : originMeshRead(str, nP) {
    createSectionIndexFluentMap();
}

void OpenHurricane::FluentMeshRead::clear() noexcept {
    originMeshRead::clear();
    sectionIndexFluentMap_.clear();
}

void OpenHurricane::FluentMeshRead::reading(string &gridUnit) {
    readFluent();
}