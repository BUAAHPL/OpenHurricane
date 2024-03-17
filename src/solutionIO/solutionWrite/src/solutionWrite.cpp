/*!
 * \file solutionWrite.cpp
 * \brief Main subroutines for solution write-out.
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

#include "solutionWrite.hpp"
#include "Lapacians.hpp"
#include "calculateFieldVar.hpp"
#include "profiles.hpp"

namespace OpenHurricane {
    createClassNameStr(solutionWrite, "solutionWrite");
}

namespace OpenHurricane {
    createObjFty(solutionWrite, controller);
}

void OpenHurricane::solutionWrite::getOutSetMap(const stringList &outFieldVarList) {
    for (integer i = 0; i < outFieldVarList.size(); ++i) {
        if (!mesh().foundOnlyObject(outFieldVarList[i])) {
            if (outFieldVarMap_.find(outFieldVarList[i]) != outFieldVarMap_.end()) {
#ifdef HUR_DEBUG
                checkWarningStr(
                    ("The variable " + outFieldVarList[i] + "has already been set in outFieldMap"));
#endif // HUR_DEBUG
            } else {
                outFieldVarMap_.emplace(outFieldVarList[i], nullptr);
            }
        }
        if (sets_.find(outFieldVarList[i]) != sets_.end()) {
#ifdef HUR_DEBUG
            checkWarningStr(
                ("The variable " + outFieldVarList[i] + "has already been set in outFieldMap"));
#endif // HUR_DEBUG
        } else {
            sets_.emplace(outFieldVarList[i]);
        }

        if (solType_ != solutionType::onlyInstantaneous) {
            if (outAveFieldVarMap_.find(outFieldVarList[i]) != outAveFieldVarMap_.end()) {
#ifdef HUR_DEBUG
                checkWarningStr(
                    ("The variable " + outFieldVarList[i] + "has already been set in outFieldMap"));
#endif // HUR_DEBUG
            } else {
                outAveFieldVarMap_.emplace(outFieldVarList[i], nullptr);
            }
            if (solType_ == solutionType::onlyPulse ||
                solType_ == solutionType::InstantAndAveAndPulse) {
                if (outPulseFieldVarMap_.find(outFieldVarList[i]) != outPulseFieldVarMap_.end()) {
#ifdef HUR_DEBUG
                    checkWarningStr(("The variable " + outFieldVarList[i] +
                                     "has already been set in outFieldMap"));
#endif // HUR_DEBUG
                } else {
                    outPulseFieldVarMap_.emplace(outFieldVarList[i], nullptr);
                }
            }
        }
    }
}

hur_nodiscard OpenHurricane::stringList OpenHurricane::solutionWrite::outVarNameList() const {
    stringList nl;
    if (solType_ == solutionType::onlyInstantaneous) {
        return outVarNameList_;
    } else if (solType_ == solutionType::onlyAverage) {
        nl.resize(outVarNameList_.size());
        for (integer i = 0; i < nl.size(); ++i) {
            nl[i] = outVarNameList_[i] + "_Ave";
        }
    } else if (solType_ == solutionType::onlyPulse) {
        nl.resize(outVarNameList_.size());
        for (integer i = 0; i < nl.size(); ++i) {
            nl[i] = outVarNameList_[i] + "_Pulse";
        }
    } else if (solType_ == solutionType::InatsntAndAve) {
        nl.resize(outVarNameList_.size());
        for (integer i = 0; i < nl.size(); ++i) {
            nl[i] = outVarNameList_[i];
        }
        stringList nlAve(outVarNameList_.size());
        for (integer i = 0; i < nlAve.size(); ++i) {
            nlAve[i] = outVarNameList_[i] + "_Ave";
        }
        nl.append(nlAve);
    } else if (solType_ == solutionType::InstantAndAveAndPulse) {
        nl.resize(outVarNameList_.size());
        for (integer i = 0; i < nl.size(); ++i) {
            nl[i] = outVarNameList_[i];
        }
        stringList nlAve(outVarNameList_.size());
        for (integer i = 0; i < nlAve.size(); ++i) {
            nlAve[i] = outVarNameList_[i] + "_Ave";
        }
        nl.append(nlAve);

        stringList nlPul(outVarNameList_.size());
        for (integer i = 0; i < nlPul.size(); ++i) {
            nlPul[i] = outVarNameList_[i] + "_Pulse";
        }
        nl.append(nlPul);
    }
    return nl;
}

void OpenHurricane::solutionWrite::removeDuplicate(stdStringList &writeList) const {
    if (writeList.size() == 0) {
        return;
    }
    integer m = 0;
    stdStringList tmpList(writeList.size());

    for (integer i = 0; i < writeList.size(); ++i) {
        integer j = 0;
        for (; j < m; ++j) {
            if (writeList[i] == tmpList[j]) {
                break;
            }
        }

        if (j == m) {
            tmpList[m] = writeList[i];
            m++;
        }
    }
    tmpList.resize(m);
    writeList.transfer(tmpList);
}

void OpenHurricane::solutionWrite::getVarList(const controller &cont) {
    if (cont.found("varList")) {
        std::string rnl = cont.findText("varList");
        replaceAllMarks(rnl, "\n", " ");
        if (!rnl.empty()) {
            stdStringList rll;
            split(rnl, rll, ",");

            for (integer i = 0; i < rll.size(); ++i) {
                trim(rll[i]);
            }

            removeDuplicate(rll);
            stdStringList spN;
            integer iAllSpc = -1;
            for (integer i = 0; i < rll.size(); ++i) {
                if (rll[i] == "allSpecies") {
                    iAllSpc = i;
                    spN.resize(flows_.mixtures().species().size());
                    for (integer isp = 0; isp < spN.size(); ++isp) {
                        spN[isp] = flows_.mixtures().species().name(isp);
                    }
                    break;
                }
            }

            if (iAllSpc != -1) {
                if (iAllSpc < rll.size() - 1) {
                    rll.insert(&rll[iAllSpc], spN);
                } else {
                    rll.remove(rll[iAllSpc]);
                    rll.append(spN);
                }
            }

            iAllSpc = -1;
            for (integer i = 0; i < rll.size(); ++i) {
                if (rll[i] == "allSpecies-mole") {
                    iAllSpc = i;
                    spN.resize(flows_.mixtures().species().size());
                    for (integer isp = 0; isp < spN.size(); ++isp) {
                        spN[isp] = flows_.mixtures().species().name(isp) + "-mole";
                    }
                    break;
                }
            }

            if (iAllSpc != -1) {
                if (iAllSpc < rll.size() - 1) {
                    rll.insert(&rll[iAllSpc], spN);
                } else {
                    rll.remove(rll[iAllSpc]);
                    rll.append(spN);
                }
            }

            outVarNameList_.resize(rll.size());
            for (integer i = 0; i < rll.size(); ++i) {
                outVarNameList_[i] = trimCopy(rll[i]);
            }
        }
    } else {
        outVarNameList_ = mesh().outputTitleNameDocList();
    }
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::solutionWrite::outFile(const char *c) const {
    fileName outN = iter_.outputName();
    auto pathOut = outN.parentPath();
    auto fname = outN.name(true);
    fname += "-";
    fname += toString(iter_.totalStep());
    string fext = c;
    outN = fname + fext;
    outN = pathOut / outN;
    return outN;
}

void OpenHurricane::solutionWrite::clearOutVarMap(
    std::map<std::string, object *> &outFieldVarMap) const {
    for (auto &iter : outFieldVarMap) {
        HurDelete(iter.second);
    }
}

void OpenHurricane::solutionWrite::calcAveVar() const {
    if (!iter().hasPhysicalTimeStep()) {
        return;
    }

    if (solType_ == solutionType::onlyInstantaneous) {
        return;
    }

    const auto cTT = iter().pTStep().totalTime();
    if (cTT <= startTime_) {
        timeElasped_ = 0;
        lastTime_ = cTT;
        return;
    }
    const auto dt = cTT - lastTime_;

    lastTime_ = cTT;

    timeElasped_ += dt;
    for (auto &e : outAveFieldVarMap_) {
        if (!mesh().foundOnlyObject(e.first)) {
            HurDelete(e.second);
            continue;
        }
        const auto &ob = mesh().findOnlyObject(e.first);
        const auto &obfonl = ob.outputVarNameL();
        if (ob.nElements() == 1) {
            const auto &obf = mesh().findObject<cellRealArray>(e.first);
            if (e.second == nullptr) {
                e.second =
                    new cellRealArray(object(e.first + "_Ave", string("\"" + obfonl[0] + "_Ave\""),
                                             mesh(), object::WRITE_OUTPUT),
                                      mesh(), obf);
                (*e.second).calcTimeSumPtr(dt);
            } else {
                auto &obfav = static_cast<cellRealArray &>(*e.second);
                obfav = obf;
                obfav.calcTimeSumPtr(dt);
                obfav = obfav.getTimeSumPtr() / timeElasped_;
            }
        } else if (ob.nElements() == 3) {
            const auto &obf = mesh().findObject<cellVectorArray>(e.first);

            if (e.second == nullptr) {
                e.second = new cellVectorArray(
                    object(e.first + "_Ave",
                           string("\"" + obfonl[0] + "_Ave\"," + "\"" + obfonl[1] + "_Ave\"," +
                                  "\"" + obfonl[2] + "_Ave\""),
                           mesh(), object::WRITE_OUTPUT),
                    mesh(), obf);
                (*e.second).calcTimeSumPtr(dt);
            } else {
                auto &obfav = static_cast<cellVectorArray &>(*e.second);
                obfav = obf;
                obfav.calcTimeSumPtr(dt);
                obfav = obfav.getTimeSumPtr() / timeElasped_;
            }
        } else if (ob.nElements() == 2) {
            const auto &obf = mesh().findObject<cellVector2DArray>(e.first);

            if (e.second == nullptr) {
                e.second = new cellVector2DArray(
                    object(e.first + "_Ave",
                           string("\"" + obfonl[0] + "_Ave\"," + "\"" + obfonl[1] + "_Ave\""),
                           mesh(), object::WRITE_OUTPUT),
                    mesh(), obf);
                (*e.second).calcTimeSumPtr(dt);
            } else {
                auto &obfav = static_cast<cellVector2DArray &>(*e.second);
                obfav = obf;
                obfav.calcTimeSumPtr(dt);
                obfav = obfav.getTimeSumPtr() / timeElasped_;
            }
        } else if (ob.nElements() == 9) {
            const auto &obf = mesh().findObject<cellTensorArray>(e.first);

            if (e.second == nullptr) {
                e.second = new cellTensorArray(
                    object(e.first + "_Ave",
                           string("\"" + obfonl[0] + "_Ave\"," + "\"" + obfonl[1] + "_Ave\"," +
                                  "\"" + obfonl[2] + "_Ave\"," + "\"" + obfonl[3] + "_Ave\"," +
                                  "\"" + obfonl[4] + "_Ave\"," + "\"" + obfonl[5] + "_Ave\"," +
                                  "\"" + obfonl[6] + "_Ave\"," + "\"" + obfonl[7] + "_Ave\"," +
                                  "\"" + obfonl[8] + "_Ave\""),
                           mesh(), object::WRITE_OUTPUT),
                    mesh(), obf);
                (*e.second).calcTimeSumPtr(dt);
            } else {
                auto &obfav = static_cast<cellTensorArray &>(*e.second);
                obfav = obf;
                obfav.calcTimeSumPtr(dt);
                obfav = obfav.getTimeSumPtr() / timeElasped_;
            }
        } else {
            errorAbortStr(("The type of cell array, which has " + toString(ob.nElements()) +
                           " components, is not support yet."));
        }
    }
}

void OpenHurricane::solutionWrite::calcPulseVar() const {
    if (!iter().hasPhysicalTimeStep()) {
        return;
    }

    if (solType_ != solutionType::onlyPulse && solType_ != solutionType::InstantAndAveAndPulse) {
        return;
    }
    const auto cTT = iter().pTStep().totalTime();
    if (cTT <= startTime_) {
        return;
    }

    for (auto &e : outPulseFieldVarMap_) {
        if (!mesh().foundOnlyObject(e.first)) {
            HurDelete(e.second);
            continue;
        }
        const auto &ob = mesh().findOnlyObject(e.first);
        const auto &obfonl = ob.outputVarNameL();
        if (ob.nElements() == 1) {
            const auto &obf = mesh().findObject<cellRealArray>(e.first);
            const auto &obfav = mesh().findObject<cellRealArray>(e.first + "_Ave");
            HurDelete(e.second);
            e.second =
                new cellRealArray(object(e.first + "_Pulse", string("\"" + obfonl[0] + "_Pulse\""),
                                         mesh(), object::WRITE_OUTPUT),
                                  mesh(), obf - obfav);
        } else if (ob.nElements() == 3) {
            const auto &obf = mesh().findObject<cellVectorArray>(e.first);
            const auto &obfav = mesh().findObject<cellVectorArray>(e.first + "_Ave");
            HurDelete(e.second);
            e.second = new cellVectorArray(
                object(e.first + "_Pulse",
                       string("\"" + obfonl[0] + "_Pulse\"," + "\"" + obfonl[1] + "_Pulse\"," +
                              "\"" + obfonl[2] + "_Pulse\""),
                       mesh(), object::WRITE_OUTPUT),
                mesh(), obf - obfav);
        } else if (ob.nElements() == 2) {
            const auto &obf = mesh().findObject<cellVector2DArray>(e.first);
            const auto &obfav = mesh().findObject<cellVector2DArray>(e.first + "_Ave");
            HurDelete(e.second);
            e.second = new cellVector2DArray(
                object(e.first + "_Pulse",
                       string("\"" + obfonl[0] + "_Pulse\"," + "\"" + obfonl[1] + "_Pulse\""),
                       mesh(), object::WRITE_OUTPUT),
                mesh(), obf - obfav);
        } else if (ob.nElements() == 9) {
            const auto &obf = mesh().findObject<cellTensorArray>(e.first);
            const auto &obfav = mesh().findObject<cellTensorArray>(e.first + "_Ave");
            HurDelete(e.second);
            e.second = new cellTensorArray(
                object(e.first + "_Pulse",
                       string("\"" + obfonl[0] + "_Pulse\"," + "\"" + obfonl[1] + "_Pulse\"," +
                              "\"" + obfonl[2] + "_Pulse\"," + "\"" + obfonl[3] + "_Pulse\"," +
                              "\"" + obfonl[4] + "_Pulse\"," + "\"" + obfonl[5] + "_Pulse\"," +
                              "\"" + obfonl[6] + "_Pulse\"," + "\"" + obfonl[7] + "_Pulse\"," +
                              "\"" + obfonl[8] + "_Pulse\""),
                       mesh(), object::WRITE_OUTPUT),
                mesh(), obf - obfav);
        } else {
            errorAbortStr(("The type of cell array, which has " + toString(ob.nElements()) +
                           " components, is not support yet."));
        }
    }
}

void OpenHurricane::solutionWrite::updating() const {
    calcStagnationParameters();
    calcViscousRatio();
    calcVorticity();
    calcMoleSpecies();

    // These two subroutines must be called at the end of this funtion
    calcAveVar();
    calcPulseVar();
}

OpenHurricane::solutionWrite::solutionWrite(const flowModel &flows, const iteration &iter,
                                            const runtimeMesh &mesh, const controller &cont)
    : flows_(flows), iter_(iter), mesh_(mesh), sets_(), outFieldVarMap_(), outAveFieldVarMap_(),
      outPulseFieldVarMap_(), solType_(solutionType::onlyInstantaneous),
      startTime_(cont.findOrDefault<real>("startTime", 0)), timeElasped_(0), lastTime_(0),
      firstCalling_(true), outVarNameList_(), isWriteProfile_(false), isSampling_(false),
      samplingStep_(10) {
    getVarList(cont);
    if (iter.isUnsteadyFlow()) {
        if (cont.found("solutionType")) {
            auto wtw = cont.findWord("solutionType");
            trim(wtw);
            if (wtw == "onlyInstantaneous") {
                solType_ = solutionType::onlyInstantaneous;
            } else if (wtw == "onlyAverage") {
                solType_ = solutionType::onlyAverage;
            } else if (wtw == "onlyPulse") {
                solType_ = solutionType::onlyPulse;
            } else if (wtw == "InatsntAndAve") {
                solType_ = solutionType::InatsntAndAve;
            } else if (wtw == "InstantAndAveAndPulse") {
                solType_ = solutionType::InstantAndAveAndPulse;
            } else {
                errorAbortStr(("Unknown writeType: " + wtw + " in " + cont.name()));
            }

            if (solType_ != solutionType::onlyInstantaneous) {
                if (cont.found("sampling")) {
                    const auto &samCont = cont.subController("sampling");

                    controllerSwitch issamCont(samCont);
                    isSampling_ = issamCont("isSampling", isSampling_);
                    samplingStep_ = samCont.findOrDefault<integer>("samplingStep", samplingStep_);
                }
            }
        }
    }
    getOutSetMap(outVarNameList_);

    if (cont.found("writeProfile")) {
        isWriteProfile_ = true;
    }
}

OpenHurricane::uniquePtr<OpenHurricane::solutionWrite>
OpenHurricane::solutionWrite::creator(const flowModel &flows, const iteration &iter,
                                      const runtimeMesh &mesh, const controller &cont) {
    string monitorType = cont.findWord("writeType");

    defineInObjCreator(solutionWrite, static_cast<std::string>(monitorType), controller,
                       (flows, iter, mesh, cont));
}

void OpenHurricane::solutionWrite::write() const {
    if (iter().hasPhysicalTimeStep() && firstCalling_) {
        const auto cTT = iter().pTStep().totalTime();
        timeElasped_ = 0;
        lastTime_ = cTT - iter().pTStep().pTimeStep();
        firstCalling_ = false;
    }
    if (writeNow() || samplingNow()) {
        updating();
    }

    if (writeNow()) {
        writeToFile();
        writeProfiles();
        clearOutVarMap(outFieldVarMap_);
        clearOutVarMap(outPulseFieldVarMap_);
    }
}

void OpenHurricane::solutionWrite::writeProfiles() const {
    if (!isWriteProfile_) {
        return;
    }

    const auto &wCont = iter().cont().subController("writeControl");
    const auto &solCont = wCont.subController("solutionWrite");
    const auto &prfCont = solCont.subController("writeProfile");

    fileName outPrfName;
    if (prfCont.found("fileName")) {
        outPrfName = prfCont.findWord("fileName");
    } else {
        outPrfName = "HurricaneProfile.dat";
    }

    fileName outPath;
    if (!outPrfName.isAbsolute()) {
        outPath = iter().outputPath();
    } else {
        outPath = outPrfName.parentPath();
    }
    auto fname = outPrfName.name(true);
    fname += "-";
    fname += toString(iter_.totalStep());
    string fext = outPrfName.ext();
    outPrfName = fname + fext;
    outPrfName = outPath / outPrfName;
    auto prfTw = prfCont.findWordOrDefault("profileType", "OpenHurricane");
    auto myProfPtr = profiles::creator(outPrfName, prfTw);
    myProfPtr->write(flows(), prfCont);
}

hur_nodiscard bool OpenHurricane::solutionWrite::samplingNow() const noexcept {
    if (outVarNameList_.size() == 0) {
        return false;
    }
    if (isSampling_ && iter_.hasPhysicalTimeStep()) {
        if (iter_.cStep() % samplingStep_ == 0) {
            return true;
        }
    }
    return false;
}

hur_nodiscard bool OpenHurricane::solutionWrite::writeNow() const noexcept {
    if (outVarNameList_.size() == 0) {
        return false;
    } else {
        if (iter_.cStep() % iter_.writeOutputStep() == 0 ||
            (iter_.end() && !argParse::noWriteAtEnd())) {
            return true;
        }
    }
    return false;
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName,
                                                  const realArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellRealArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, realArray &&value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellRealArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(),
                                          std::move(value));
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  const realArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellRealArray(object(varName, outName, mesh(), object::WRITE_OUTPUT), mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  realArray &&value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellRealArray(object(varName, outName, mesh(), object::WRITE_OUTPUT),
                                          mesh(), std::move(value));
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName,
                                                  const vectorArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellVectorArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, vectorArray &&val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellVectorArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(),
                                            std::move(val));
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  const vectorArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellVectorArray(object(varName, outName, mesh(), object::WRITE_OUTPUT),
                                            mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  vectorArray &&val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellVectorArray(object(varName, outName, mesh(), object::WRITE_OUTPUT),
                                            mesh(), std::move(val));
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName,
                                                  const vector2DArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellVector2DArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName,
                                                  vector2DArray &&val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellVector2DArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(),
                                              std::move(val));
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  const vector2DArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellVector2DArray(
            object(varName, outName, mesh(), object::WRITE_OUTPUT), mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  vector2DArray &&val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellVector2DArray(
            object(varName, outName, mesh(), object::WRITE_OUTPUT), mesh(), std::move(val));
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName,
                                                  const tensorArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellTensorArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, tensorArray &&val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellTensorArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(),
                                            std::move(val));
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  const tensorArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellTensorArray(object(varName, outName, mesh(), object::WRITE_OUTPUT),
                                            mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  tensorArray &&val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellTensorArray(object(varName, outName, mesh(), object::WRITE_OUTPUT),
                                            mesh(), std::move(val));
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName,
                                                  const symmTensorArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellSymmTensorArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName,
                                                  symmTensorArray &&val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellSymmTensorArray(object(varName, mesh(), object::WRITE_OUTPUT),
                                                mesh(), std::move(val));
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  const symmTensorArray &val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellSymmTensorArray(
            object(varName, outName, mesh(), object::WRITE_OUTPUT), mesh(), val);
    }
}

void OpenHurricane::solutionWrite::setOutputField(const string &varName, const string &outName,
                                                  symmTensorArray &&val) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellSymmTensorArray(
            object(varName, outName, mesh(), object::WRITE_OUTPUT), mesh(), std::move(val));
    }
}

void OpenHurricane::solutionWrite::calcStagnationParameters() const {
    const auto iter = sets_.find("Ma");
    if (iter != sets_.end()) {
        setOutputField(*iter, calculateFieldVar::calcMachNumber(v(), gama(), p(), rho()));
    }
    const auto iter2 = sets_.find("pt");
    if (iter2 != sets_.end()) {
        setOutputField(*iter2, calculateFieldVar::calcTotalPressure(flows_));
    }

    const auto iterPtCp = sets_.find("pt_constCp");
    if (iterPtCp != sets_.end()) {
        setOutputField(*iterPtCp,
                       calculateFieldVar::calcTotalPressure(p(), gama(), v(), rho(), flows_));
    }

    const auto iter3 = sets_.find("Tt");
    if (iter3 != sets_.end()) {
        if (flows_.mixtures().isSingular()) {
            setOutputField(*iter3, calculateFieldVar::calcTotalTemperature(T(), gama(), v(), p(),
                                                                           rho(), flows_));
        } else {
            setOutputField(*iter3,
                           calculateFieldVar::totalTemperature(flows_, p(), T(), v(), yi()));
        }
    }
    const auto iterTtCp = sets_.find("Tt_constCp");
    if (iterTtCp != sets_.end()) {
        setOutputField(*iterTtCp, calculateFieldVar::calcTotalTemperature(T(), gama(), v(), p(),
                                                                          rho(), flows_));
    }

    const auto iterHat = sets_.find("hat");
    if (iterHat != sets_.end()) {
        setOutputField(*iterHat, calculateFieldVar::calcTotalEnthalpy(flows_));
    }

    const auto iterCV = sets_.find("cellVolume");
    if (iterCV != sets_.end()) {
        setOutputField(*iterCV, calculateFieldVar::cellVolume(flows_));
    }

    const auto iterShadowG = sets_.find("Shadowgraph");
    if (iterShadowG != sets_.end()) {
        setOutputField(*iterShadowG, fv::lapacian(rho()));
    }
}

void OpenHurricane::solutionWrite::calcViscousRatio() const {
    const auto iter = sets_.find("viscousRatio");
    if (iter != sets_.end()) {
        setOutputField(*iter, calculateFieldVar::calcViscousRatio(flows_));
    }
}

void OpenHurricane::solutionWrite::calcVorticity() const {
    const auto iter = sets_.find("vorticity");
    if (iter != sets_.end()) {
        setOutputField(*iter, calculateFieldVar::calcVorticity(flows_));
    }

    const auto iterQ = sets_.find("QCriterion");
    if (iterQ != sets_.end()) {
        setOutputField(*iterQ, calculateFieldVar::calcQCriterion(flows_));
    }

    const auto iterO = sets_.find("OmegaCriterion");
    if (iterO != sets_.end()) {
        setOutputField(*iterO, calculateFieldVar::calcOmegaCriterion(flows_));
    }

    const auto iterDel = sets_.find("DeltaCriterion");
    if (iterDel != sets_.end()) {
        setOutputField(*iterDel, calculateFieldVar::calcDeltaCriterion(flows_));
    }
}

void OpenHurricane::solutionWrite::calcMoleSpecies() const {
    if (flows_.mixtures().isSingular()) {
        return;
    }
    bool isOutXi = false;
    for (integer i = 0; i < flows_.mixtures().species().size(); ++i) {
        auto iter = sets_.find(flows_.mixtures().species()[i].name() + "-mole");
        if (iter != sets_.end()) {
            isOutXi = true;
            break;
        }
    }
    if (!isOutXi) {
        return;
    }
    PtrList<cellRealArray> tmpXi;
    for (integer i = 0; i < flows_.mixtures().species().size(); ++i) {
        tmpXi.append(new cellRealArray(
            object(flows_.mixtures().species()[i].name() + "-mole", mesh(), object::WRITE_OUTPUT),
            mesh()));
    }
    flows_.mixtures().species().Yi2Xi(yi(), tmpXi);
    for (integer i = 0; i < flows_.mixtures().species().size(); ++i) {
        auto iter = outFieldVarMap_.find(flows_.mixtures().species()[i].name() + "-mole");
        if (iter != outFieldVarMap_.end()) {
            iter->second = tmpXi.set(i, nullptr);
        }
    }
}
