/*!
 * \file FluentProfiles.cpp
 * \brief Main subroutines for Ansys Fluent Profiles.
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

#include "FluentProfiles.hpp"
#include "fileOsstream.hpp"

namespace OpenHurricane {
    createClassNameStr(FluentProfiles, "Fluent");
    registerObjFty(profiles, FluentProfiles, controller);
} // namespace OpenHurricane

OpenHurricane::FluentProfiles::FluentProfiles(const fileName &file) : profiles(file) {}

void OpenHurricane::FluentProfiles::readProfiles(controller &cont) const {
    Pout << "    Reading Ansys Fluent profile file..." << std::endl;
    auto fstr = readFileToString(files_);
    bool valid = true;
    integer blCount = 0;
    integerList bPosCount;
    integerList ePosCount;
    for (size_t i = 0; i < fstr.size(); ++i) {
        const auto &ich = fstr[i];
        if (blCount == 0 && ich == ')') {
            valid = false;
            break;
        } else if (blCount == 0 && ich == '(') {
            blCount++;
            bPosCount.append(static_cast<integer>(i));
        } else if (blCount == 1 && ich == ')') {
            blCount--;
            ePosCount.append(static_cast<integer>(i));
        } else if (blCount == 1 && ich == '(') {
            blCount++;
        } else if (blCount == 2 && ich == ')') {
            blCount--;
        } else if (blCount == 2 && ich == '(') {
            valid = false;
            break;
        }
    }

    if (valid) {
        if (blCount != 0) {
            valid = false;
        }
        if (bPosCount.size() == 0 || ePosCount.size() == 0) {
            valid = false;
        }
        if (bPosCount.size() != ePosCount.size()) {
            valid = false;
        }
    }

    if (!valid) {
        LFatal("The file %s is not valid Ansys Fluent standard profile file, "
               "please check!",
               files_.c_str());
    }
    std::string profileList;
    for (integer i = 0; i < bPosCount.size(); ++i) {
        auto proName =
            readFluentProfile(fstr.substr(bPosCount[i] + 1, ePosCount[i] - bPosCount[i] - 1), cont);
        if (profileList.empty()) {
            profileList = proName;
        } else {
            profileList += ",";
            profileList += proName;
        }
    }
    if (cont.found("profileList")) {
        auto oldProList = cont.findText("profileList");
        profileList = oldProList + "," + profileList;
    }
    cont.setText("profileList", profileList);
}

void OpenHurricane::FluentProfiles::writeProfiles(const controller &cont) const {
    if (!HurMPI::master()) {
        return;
    }
    if (cont.found("profileList")) {
        fileOsstream myOut(files_);
        const auto pfl = cont.findTextStr("profileList");
        myOut.setRealPrecision();
        for (integer i = 0; i < pfl.size(); ++i) {
            const auto &pflCont = cont.subController(pfl[i]);
            if (!pflCont.found("fieldNameList")) {
                LFatal("Cannot found field name list in %s for writting profile to %s",
                       pflCont.name().c_str(), files_.c_str());
            }
            const auto pflNL = pflCont.findTextStr("fieldNameList");

            if (pflNL.size() == 0) {
                LFatal("Cannot found field name list in %s for writting profile to %s",
                       pflCont.name().c_str(), files_.c_str());
            }

            const auto &pfl0 = pflCont.findRealArray(pflNL[0]);
            myOut.os() << "((" << pfl[i].c_str() << " point " << pfl0.size() << ")" << std::endl;

            myOut.os() << " (" << pflNL[0].c_str() << std::endl << " \t";
            for (integer j = 0; j < pfl0.size(); ++j) {
                myOut.os() << pfl0[j] << " ";
                if (j == pfl0.size() - 1) {
                    myOut.os() << ")" << std::endl;
                } else if ((j != 0) && (j % 4 == 0)) {
                    myOut.os() << std::endl << " \t";
                }
            }

            for (integer k = 1; k < pflNL.size(); ++k) {
                const auto &pflk = pflCont.findRealArray(pflNL[k]);
                myOut.os() << " (" << pflNL[k].c_str() << std::endl << " \t";
                for (integer j = 0; j < pflk.size(); ++j) {
                    myOut.os() << pflk[j] << " ";
                    if (j == pflk.size() - 1) {
                        myOut.os() << ")" << std::endl;
                    } else if ((j != 0) && (j % 4 == 0)) {
                        myOut.os() << std::endl << " \t";
                    }
                }
            }

            myOut.os() << ")" << std::endl;
        }
        myOut.close();
    } else {
        LFatal("The profile list is not found in %s", cont.name().c_str());
    }
}

std::string OpenHurricane::FluentProfiles::readFluentProfile(const std::string &str,
                                                             controller &cont) const {
    auto blPos = str.find_first_of('(');
    auto elPos = str.find_first_of(')');

    std::string headerStr = str.substr(blPos + 1, elPos - blPos - 1);
    replaceAllMarks(headerStr, "\n", " ");
    replaceAllMarks(headerStr, "\t", " ");
    replaceAllMarks(headerStr, "\r", " ");
    stdStringList hsl;
    split(headerStr, hsl, " ");
    if (hsl.size() >= 3) {
        char *endptr;
        integer size = 0;

        if (hsl[1] == "point" || hsl[1] == "line") {
            size = (integer)strtol(hsl[2].c_str(), &endptr, 10);
            Pout << "        \"" << hsl[0] << "\" " << hsl[1] << "-profile with n = " << size
                 << ": ";
        } else if (hsl[1] == "mesh") {
            auto mm = (integer)strtol(hsl[2].c_str(), &endptr, 10);
            auto nn = (integer)strtol(hsl[3].c_str(), &endptr, 10);
            size = mm * nn;
            Pout << "        \"" << hsl[0] << "\" " << hsl[1] << "-profile with m = " << mm
                 << ", n = " << nn << ": ";
        } else if (hsl[1] == "transient") {
            size = (integer)strtol(hsl[2].c_str(), &endptr, 10);
            Pout << "        \"" << hsl[0] << "\" " << hsl[1] << "-profile with n = " << size
                 << ": ";
        } else {
            errorAbortStr(("The type of profile: " + hsl[1] + " is not supported yet"));
        }
        controller proCont(hsl[0], cont);
        std::string nameList;
        std::string outnameList;
        while (true) {
            blPos = str.find_first_of('(', elPos + 1);
            elPos = str.find_first_of(')', elPos + 1);
            bool valid = true;
            if (blPos == std::string::npos && elPos == std::string::npos) {
                break;
            } else if (blPos != std::string::npos && elPos == std::string::npos) {
                valid = false;
            } else if (blPos == std::string::npos && elPos != std::string::npos) {
                valid = false;
            }
            if (!valid) {
                LFatal("The file is not valid Ansys Fluent standard "
                       "profile file, please check!");
            }

            std::string varStr = str.substr(blPos + 1, elPos - blPos - 1);
            replaceAllMarks(varStr, "\n", " ");
            replaceAllMarks(varStr, "\t", " ");
            replaceAllMarks(varStr, "\r", " ");
            trim(varStr);
            stdStringList vsl;
            split(varStr, vsl, " ");
            if (vsl.size() != size + 1) {
                errorAbortStr(("The size of filed variable " + vsl[0] +
                               " in Ansys Fluent profile file: " + files_ + " is not equal to " +
                               toString(size)));
            }
            realArray v(size, Zero);
            for (integer i = 1; i <= size; ++i) {
                v[i - 1] = (real)strtod(vsl[i].c_str(), &endptr);
            }
            proCont.addRealArray(vsl[0], v);
            if (nameList.empty()) {
                nameList = vsl[0];
                outnameList = vsl[0];
            } else {
                nameList += ",";
                nameList += vsl[0];
                outnameList += ", ";
                outnameList += vsl[0];
            }
        }
        proCont.addText("fieldNameList", nameList);
        Pout << outnameList << std::endl;
        if (hsl[1] == "point" || hsl[1] == "line" || hsl[1] == "mesh") {
            if (!proCont.found("x") || !proCont.found("y") || !proCont.found("z")) {
                LFatal("The Ansys Fluent profile of type \"point\", "
                       "\"line\" and \"mesh\" must contain fields with "
                       "names \"x\", \"y\", and \"z\"");
            }
        } else if (hsl[1] == "transient") {
            if (!proCont.found("t")) {
                LFatal("The Ansys Fluent profile of type \"transient\" must contain field with "
                       "name \"t\"");
            }
        }
        cont.addController(hsl[0], proCont);
    }
    return hsl[0];
}
