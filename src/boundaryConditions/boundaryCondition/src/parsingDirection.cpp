/*!
 * \file parsingDirection.cpp
 * \brief Main subroutines for parsing boundary velocity direction.
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
#include "parsingDirection.hpp"
#include "constants.hpp"

OpenHurricane::string OpenHurricane::parsingDirection::getDirection(controller &bcCont,
                                                            const runtimeMesh &mesh,
                                                            const faceZone &fz, vector &direct,
                                                            const bool addToCont) {
    string dtw;
    if (bcCont.found("directionType")) {
        string str = bcCont.findWord("directionType");
        if (str == string("normalToBoundary")) {
            integer flag = Zero;
            integer rightproc;
            integerList buf(HurMPI::getProcSize(), Zero);
            realList direction(3);

            for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
                if (HurMPI::isThisProc(ip)) {
                    if (fz.size() != 0) {
                        flag = 1;
                        vector normal = mesh.faceArea()[fz.firstIndex()].normalized();
                        direction[0] = normal.x();
                        direction[1] = normal.y();
                        direction[2] = normal.z();
                        break;
                    }
                }
            }
            HurMPI::gather(&flag, 1, feature<integer>::MPIType, buf.data(), 1,
                           feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

            if (HurMPI::master()) {
                for (integer ip = 0; ip < buf.size(); ip++) {
                    if (buf[ip] != 0) {
                        rightproc = ip;
                        break;
                    }
                }
            }
            HurMPI::bcast(&rightproc, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                          HurMPI::getComm());
            HurMPI::bcast(direction.data(), 3, feature<real>::MPIType, rightproc,
                          HurMPI::getComm());

            direct.x() = direction[0];
            direct.y() = direction[1];
            direct.z() = direction[2];
            if (addToCont) {
                string directWord = string("car") + toString(direct);
                if (bcCont.found("direction")) {
                    bcCont.remove("direction");
                }
                bcCont.addWord(string("direction"), directWord);
            }
        } else if (str == string("directionVector")) {
            if (!bcCont.found("direction")) {
                errorAbortStr(("The direction must be specified in face zone: " + fz.name()));
            }
            const auto directWord = bcCont.findWord("direction");
            const std::regex contFlag("(.*?)\\s*\\((.*?)(?:,\\s*?)(.*?)(?:,\\s*?)(.*?)\\s*?\\)");
            std::smatch what;
            std::regex_search(directWord, what, contFlag);
            std::string type = what[1];
            std::istringstream cmpt1(what[2]);
            std::istringstream cmpt2(what[3]);
            std::istringstream cmpt3(what[4]);
            real x = static_cast<real>(feature<real>(cmpt1));
            real y = static_cast<real>(feature<real>(cmpt2));
            real z = static_cast<real>(feature<real>(cmpt3));
            vector direction(x, y, z);

            if (type == "car") /*!Cartesian coordinates.*/
            {
                direct = direction.normalized();
            } else {
                std::string errMsg;
                errMsg = "Current program do not support coordinate "
                         "specification in ";
                errMsg += type;
                errMsg += " way in faceZone: ";
                errMsg += fz.name();
                errorAbortStr(errMsg);
            }
        } else if (str == "angle") {
            if (!bcCont.found("alpha")) {
                errorAbortStr(("The algle of attack must be specified in face zone: " + fz.name()));
            }
            if (!bcCont.found("beta")) {
                errorAbortStr(
                    ("The algle of sideslip must be specified in face zone: " + fz.name()));
            }
            real alpha = bcCont.findType<real>("alpha", real(0));
            real beta = bcCont.findType<real>("beta", real(0));
            const real pi = constant::mathConstants::pi;
            alpha *= (pi / 180.0);
            beta *= (pi / 180.0);
            direct.x() = cos(beta) * cos(alpha);
            direct.y() = cos(beta) * sin(alpha);
            direct.z() = sin(beta);
            if (addToCont) {
                string directWord = string("car") + toString(direct);
                if (bcCont.found("direction")) {
                    bcCont.remove("direction");
                }
                bcCont.addWord(string("direction"), directWord);
            }
        } else {
            errorAbortStr(("Unknown direction type: " + str + " in face zone:" + fz.name()));
        }
        dtw = str;
    } else {
        if (!bcCont.found("direction")) {
            errorAbortStr(("The direction must be specified in face zone: " + fz.name()));
        }
        const auto directWord = bcCont.findWord("direction");
        const std::regex contFlag("(.*?)\\s*\\((.*?)(?:,\\s*?)(.*?)(?:,\\s*?)(.*?)\\s*?\\)");
        std::smatch what;
        std::regex_search(directWord, what, contFlag);
        std::string type = what[1];
        std::istringstream cmpt1(what[2]);
        std::istringstream cmpt2(what[3]);
        std::istringstream cmpt3(what[4]);
        real x = static_cast<real>(feature<real>(cmpt1));
        real y = static_cast<real>(feature<real>(cmpt2));
        real z = static_cast<real>(feature<real>(cmpt3));
        vector direction(x, y, z);

        if (type == "car") /*!Cartesian coordinates.*/
        {
            direct = direction.normalized();
        } else {
            std::string errMsg;
            errMsg = "Current program do not support coordinate specification in ";
            errMsg += type;
            errMsg += " way in faceZone: ";
            errMsg += fz.name();
            errorAbortStr(errMsg);
        }
        dtw = string("directionVector");
    }
    return dtw;
}

OpenHurricane::string OpenHurricane::parsingDirection::getDirectionType(const controller &bcCont,
                                                                const faceZone &fz) {
    if (bcCont.found("directionType")) {
        string str = bcCont.findWord("directionType");
        if (str == string("normalToBoundary")) {
            // Do nothing
        } else if (str == string("directionVector")) {
            // Do nothing
        } else if (str == "angle") {
            // Do nothing
        } else {
            errorAbortStr(("Unknown direction type: " + str + " in face zone:" + fz.name()));
        }
        return str;
    } else {
        return string("directionVector");
    }
}

void OpenHurricane::parsingDirection::getFaceNormal(const runtimeMesh &mesh, const faceZone &fz,
                                                vector &direct) {
    integer flag = Zero;
    integer rightproc;
    integerList buf(HurMPI::getProcSize(), Zero);
    realList direction(3);

    for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
        if (HurMPI::isThisProc(ip)) {
            if (fz.size() != 0) {
                flag = 1;
                vector normal = mesh.faceArea()[fz.firstIndex()].normalized();
                direction[0] = normal.x();
                direction[1] = normal.y();
                direction[2] = normal.z();
                break;
            }
        }
    }
    HurMPI::gather(&flag, 1, feature<integer>::MPIType, buf.data(), 1, feature<integer>::MPIType,
                   HurMPI::masterNo(), HurMPI::getComm());

    if (HurMPI::master()) {
        for (integer ip = 0; ip < buf.size(); ip++) {
            if (buf[ip] != 0) {
                rightproc = ip;
                break;
            }
        }
    }
    HurMPI::bcast(&rightproc, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
    HurMPI::bcast(direction.data(), 3, feature<real>::MPIType, rightproc, HurMPI::getComm());

    direct.x() = direction[0];
    direct.y() = direction[1];
    direct.z() = direction[2];
}
