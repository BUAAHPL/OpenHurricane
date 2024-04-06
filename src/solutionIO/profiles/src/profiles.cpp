/*!
 * \file profiles.cpp
 * \brief Main subroutines for profiles.
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

#include "profiles.hpp"
namespace OpenHurricane {
    createClassNameStr(profiles, "profiles");
    createObjFty(profiles, controller);
} // namespace OpenHurricane

OpenHurricane::profiles::profiles(const fileName &file) : files_(file) {}

OpenHurricane::uniquePtr<OpenHurricane::profiles>
OpenHurricane::profiles::creator(const fileName &file, const string &profileType) {

    defineInObjCreator(profiles, static_cast<std::string>(profileType), controller, (file));
}

void OpenHurricane::profiles::read(controller &cont) const {
    Pout << "    Reading " << files_ << std::endl;
    readProfiles(cont);
}

void OpenHurricane::profiles::write(const flowModel &flows, const controller &cont) const {
    controller profCont("profCont", cont);
    const auto pfl = cont.findTextStr("profileList");
    const auto pflStr = cont.findText("profileList");
    profCont.addText("profileList", pflStr);

    for (integer i = 0; i < pfl.size(); ++i) {
        controller profConti(pfl[i], profCont);
        getProfiles(flows, cont.subController(pfl[i]), profConti);
        profCont.addController(pfl[i], profConti);
    }

    writeProfiles(profCont);
}

void OpenHurricane::profiles::getProfiles(const flowModel &flows, const controller &cont,
                                          controller &profCont) const {
    const auto &mesh = flows.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fs = mesh.faces();
    const auto &fC = mesh.faceCentre();
    if (cont.found("profileFace")) {
        const auto pfw = cont.findWord("profileFace");
        integer fzi = 0;
        for (; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].name() == pfw) {
                break;
            }
        }
        const auto vsize = fzl[fzi].size();

        realArray xx(vsize, Zero);
        realArray yy(vsize, Zero);
        realArray zz(vsize, Zero);

        integer cnt = 0;
        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            xx[cnt] = fC[fi].x();
            yy[cnt] = fC[fi].y();
            zz[cnt] = fC[fi].z();
            cnt++;
        }

        profCont.addRealArray("x", xx);
        profCont.addRealArray("y", yy);
        profCont.addRealArray("z", zz);

        const auto fnl = cont.findTextStr("fieldNameList");
        std::string fnlStr;
        fnlStr = "x,y,z,";
        for (integer inl = 0; inl < fnl.size(); ++inl) {
            fnlStr += getFieldVar(flows, fzi, fnl[inl], profCont);
            if (inl != fnl.size() - 1) {
                fnlStr += ",";
            }
        }
        profCont.addText("fieldNameList", fnlStr);
    } else {
        LFatal("The profileFace is not set in %s", cont.name().c_str());
    }
}

std::string OpenHurricane::profiles::getFieldVar(const flowModel &flows, const integer fzi,
                                                 const string &fieldName,
                                                 controller &profCont) const {
    const auto &mesh = flows.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fs = mesh.faces();
    std::string nam;
    if (mesh.foundOnlyObject(fieldName)) {
        const auto &fob = mesh.findOnlyObject(fieldName);
        if (fob.nElements() == 1) {
            const auto &f = mesh.findObject<cellRealArray>(fieldName);
            if (HurMPI::parRun()) {
                realArray allVFC;
                allGatherVList(f, allVFC);
                profCont.addRealArray(fieldName, allVFC);
            } else {
                profCont.addRealArray(fieldName, fv::interpolateNoCorr(f, fzi));
            }
            nam = fieldName;
        } else if (fob.nElements() == 2) {
            const auto &f = mesh.findObject<cellVector2DArray>(fieldName);
            auto vf = fv::interpolateNoCorr(f, fzi);
            for (integer ii = 0; ii < fob.nElements(); ++ii) {
                if (HurMPI::parRun()) {
                    realArray allVFC;
                    allGatherVList(vf.component(ii), allVFC);
                    profCont.addRealArray(fob.outputVarNameL()[ii], allVFC);
                } else {
                    profCont.addRealArray(fob.outputVarNameL()[ii], vf.component(ii));
                }
                nam += fob.outputVarNameL()[ii];
                if (ii != fob.nElements() - 1) {
                    nam += ",";
                }
            }
        } else if (fob.nElements() == 3) {
            const auto &f = mesh.findObject<cellVectorArray>(fieldName);
            auto vf = fv::interpolateNoCorr(f, fzi);
            for (integer ii = 0; ii < fob.nElements(); ++ii) {
                if (HurMPI::parRun()) {
                    realArray allVFC;
                    allGatherVList(vf.component(ii), allVFC);
                    profCont.addRealArray(fob.outputVarNameL()[ii], allVFC);
                } else {
                    profCont.addRealArray(fob.outputVarNameL()[ii], vf.component(ii));
                }
                nam += fob.outputVarNameL()[ii];
                if (ii != fob.nElements() - 1) {
                    nam += ",";
                }
            }
        }

        else if (fob.nElements() == 6) {
            const auto &f = mesh.findObject<cellSymmTensorArray>(fieldName);
            auto vf = fv::interpolateNoCorr(f, fzi);
            for (integer ii = 0; ii < fob.nElements(); ++ii) {
                if (HurMPI::parRun()) {
                    realArray allVFC;
                    allGatherVList(vf.component(ii), allVFC);
                    profCont.addRealArray(fob.outputVarNameL()[ii], allVFC);
                } else {
                    profCont.addRealArray(fob.outputVarNameL()[ii], vf.component(ii));
                }
                nam += fob.outputVarNameL()[ii];
                if (ii != fob.nElements() - 1) {
                    nam += ",";
                }
            }
        } else if (fob.nElements() == 9) {
            const auto &f = mesh.findObject<cellTensorArray>(fieldName);
            auto vf = fv::interpolateNoCorr(f, fzi);
            for (integer ii = 0; ii < fob.nElements(); ++ii) {
                if (HurMPI::parRun()) {
                    realArray allVFC;
                    allGatherVList(vf.component(ii), allVFC);
                    profCont.addRealArray(fob.outputVarNameL()[ii], allVFC);
                } else {
                    profCont.addRealArray(fob.outputVarNameL()[ii], vf.component(ii));
                }
                nam += fob.outputVarNameL()[ii];
                if (ii != fob.nElements() - 1) {
                    nam += ",";
                }
            }
        } else {
            LFatal("The type of array is not supported yet");
        }
    }

    return nam;
}
