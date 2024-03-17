/*!
 * \file writeFaceZone.cpp
 * \brief Main subroutines for writeFaceZone.
 * \author Chen Zhenyi
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

#include "writeFaceZone.hpp"
#include "calculateFaceZoneVar.hpp"
#include "calculateFieldVar.hpp"
#include "faceZoneArray.hpp"
#include "fieldParameterFuncMap.hpp"
#include "tecplotWriter.hpp"

void OpenHurricane::writeFaceZone::checkDuplicate(const faceZone &fz,
                                                  const stringList &outV) const {
    const auto &mesh = flows_.mesh();
    stringList checkN(1);
    for (integer i = 0; i < outV.size(); ++i) {
        checkN[0] = outV[i];
        auto nlist = mesh.outputTitleNameDocList(checkN);
        checkN[0] += fz.name();

        auto nlist2 = mesh.outputTitleNameDocList(checkN);

        if (nlist.size() != 0 && nlist2.size() != 0) {
            LFatal("Duplicate variable name: %s written for face zone ", outV[i].c_str(),
                   fz.name().c_str());
        }
    }
}

void OpenHurricane::writeFaceZone::setFaceZoneList() {
    if (cont_.found("writeControl")) {
        auto &interCont = cont_.subController("writeControl");
        if (interCont.found("faceZones")) {
            std::string rnl = interCont.findText("faceZones");
            replaceAllMarks(rnl, "\n", " ");
            if (!rnl.empty()) {
                size_t pos = 0;
                stdStringList rll;
                split(rnl, rll, ",");
                writeFaceZoneNameList_.resize(rll.size());
                for (integer i = 0; i < rll.size(); ++i) {
                    writeFaceZoneNameList_[i] = trimCopy(rll[i]);
                }
            }
            writeFaceZoneVarList_.resize(writeFaceZoneNameList_.size());
            for (integer i = 0; i < writeFaceZoneNameList_.size(); ++i) {
                std::string nstr = interCont.findText(writeFaceZoneNameList_[i]);
                replaceAllMarks(nstr, "\n", " ");
                if (!rnl.empty()) {
                    size_t pos = 0;
                    stdStringList rll;
                    split(nstr, rll, ",");
                    writeFaceZoneVarList_[i].resize(rll.size());
                    for (integer j = 0; j < rll.size(); ++j) {
                        writeFaceZoneVarList_[i][j] = trimCopy(rll[j]);
                    }
                } else {
                    writeFaceZoneVarList_[i].resize(0);
                }
            }
        }
    }
}

void OpenHurricane::writeFaceZone::setFaceZoneId() {
    const auto &mesh = flows_.mesh();
    const auto &fz = mesh.faceZones();
    writeFaceZoneIdList_.resize(writeFaceZoneNameList_.size(), -1);
    for (integer i = 0; i < writeFaceZoneNameList_.size(); ++i) {
        for (integer j = 0; j < fz.size(); ++j) {
            if (fz[j].name() == writeFaceZoneNameList_[i]) {
                writeFaceZoneIdList_[i] = fz[j].index();
                break;
            }
        }
    }
}

void OpenHurricane::writeFaceZone::writeToFile() const {
    for (integer i = 0; i < writeFaceZoneIdList_.size(); ++i) {
        if (writeFaceZoneIdList_[i] == -1) {
            Pout << "   Warning: cannot find face zone: \"" << writeFaceZoneNameList_[i]
                 << "\", and the output skipping" << std::endl;
            continue;
        }
        const auto &mesh = flows_.mesh();
        fileName outNWall = getFileName(writeFaceZoneNameList_[i]);
        Pout << "    Writting result to file: " << outNWall << std::endl;

        fileOsstream foswall(outNWall, IOsstream::BINARY_FORMAT);

        const auto &iter = mesh.Iteration();
        if (isSampling_ && iter.hasPhysicalTimeStep()) {
            stringList outVarNL(writeFaceZoneVarList_[i].size() * 2);
            for (integer j = 0; j < writeFaceZoneVarList_[i].size(); ++j) {
                outVarNL[j] = writeFaceZoneVarList_[i][j];
                outVarNL[j + writeFaceZoneVarList_[i].size()] =
                    writeFaceZoneVarList_[i][j] + "_sampling";
            }
            writeToTecplot(foswall, writeFaceZoneIdList_[i], outVarNL);
        } else {
            writeToTecplot(foswall, writeFaceZoneIdList_[i], writeFaceZoneVarList_[i]);
        }
        foswall.close();
    }
}

void OpenHurricane::writeFaceZone::writeToTecplot(fileOsstream &fos, const integer fzid,
                                                  const stringList &outVarName) const {
    const auto &mesh = flows_.mesh();
    const auto &iter = flows_.mesh().Iteration();

    integer id = -1;
    for (integer i = 0; i < mesh.faceZones().size(); ++i) {
        if (fzid == mesh.faceZones()[i].index()) {
            id = i;
            break;
        }
    }
    if (id == -1) {
        LFatal("Cannot find face zone with index: %d", fzid);
    }

    checkDuplicate(mesh.faceZones()[id], outVarName);

    int nVar = 0;
    nVar = mesh.outputTitleSize(outVarName);
    stringList outFzN(outVarName.size());
    for (integer i = 0; i < outVarName.size(); ++i) {
        outFzN[i] = outVarName[i] + mesh.faceZones()[id].name();
    }
    nVar += mesh.outputTitleSize(outFzN);
    nVar += 3;

    if (HurMPI::master()) {
        const auto &name = mesh.globalFaceZoneInfo(fzid).fZ().name();
        string title = iter.configName().name(true);
        title += "-";
        title += name;
        stdStringList varName(nVar);
        varName[0] = "x";
        varName[1] = "y";
        varName[2] = "z";
        auto nlist = mesh.outputTitleNameDocList(outVarName);
        for (integer i = 0; i < nlist.size(); ++i) {
            integer j = i + 3;
            varName[j] = nlist[i];
        }
        auto nlist2 = mesh.outputTitleNameDocList(outFzN);
        for (integer i = 0; i < nlist2.size(); ++i) {
            integer j = i + 3 + nlist.size();
            replaceAllMarks(nlist2[i], mesh.faceZones()[id].name(), "");
            varName[j] = nlist2[i];
        }
        // Tecplot file header
        tecplotWriter::writeFileHeader(fos, title, nVar, varName, tecplotFormat::FULL);

        real st = 0.0;
        if (iter.hasPhysicalTimeStep()) {
            st = iter.pTStep().totalTime();
        }

        const auto nPts = mesh.globalFaceZoneInfo(fzid).totalNodes();
        const auto nEles = mesh.globalFaceZoneInfo(fzid).totalFaces();
        /*
        *	+-----------+
                | INT32		|	 ZoneType 0=ORDERED,		1=FELINESEG,
                +-----------+			  2=FETRIANGLE,		3=FEQUADRILATERAL,
                                                                  4=FETETRAHEDRON,	5=FEBRICK,
                                                                  6=FEPOLYGON,		7=FEPOLYHEDRON
        */
        int ZoneType = 3;
        int strandid = -1;
        if (iter.hasPhysicalTimeStep()) {
            strandid = fzid;
        }
        tecplotWriter::writeZoneHeader(fos, name, st, nVar, HurMPI::masterNo(), nPts, nEles, 0,
                                       ZoneType, strandid);

        tecplotWriter::writeEOHMARKER(fos);
    }
    tecplotWriter::writeDataSection(fos, mesh, nVar, fzid);
    mesh.writeMinMaxOutput(fos, fzid, outVarName);
    mesh.writeMinMaxOutput(fos, fzid, outFzN);

    const auto &nodel = mesh.globalFaceZoneInfo(fzid).facePoints();

    if (HurMPI::master()) {
        // Points
        tecplotWriter::writeData(fos, nodel, nodel.size());
    }

    mesh.writeOutput(fos, fzid, outVarName);
    mesh.writeOutput(fos, fzid, outFzN);

    integerList faceMap;
    tecplotWriter::faceConnectList(faceMap, mesh, fzid);

    if (HurMPI::master()) {
        fos.write(reinterpret_cast<const char *>(faceMap.data()), faceMap.byteSize());
    }
}

void OpenHurricane::writeFaceZone::updating() const {
    const auto &mesh = flows_.mesh();
    for (integer i = 0; i < writeFaceZoneIdList_.size(); ++i) {
        if (writeFaceZoneIdList_[i] == -1) {
            Pout << "   Warning: cannot find face zone: \"" << writeFaceZoneNameList_[i]
                 << "\", and the output skipping" << std::endl;
            continue;
        }
        auto fzid = writeFaceZoneIdList_[i];
        integer fzi = -1;
        for (integer i = 0; i < mesh.faceZones().size(); ++i) {
            if (fzid == mesh.faceZones()[i].index()) {
                fzi = i;
                break;
            }
        }
        if (fzi == -1) {
            LFatal("Cannot find face zone with index: %d", fzid);
        }

        const auto &varNL = writeFaceZoneVarList_[i];
        for (integer jj = 0; jj < varNL.size(); ++jj) {
            const auto keyN = calculateFaceZoneVar::setKeyName(varNL[jj]);
            if (faceZoneVarFuncMap_.found(keyN)) {
                const auto &fz = mesh.faceZones()[fzi];
                setOutputField(fz, varNL[jj], faceZoneVarFuncMap_.find(keyN)(flows_, fzi));
            }
        }
    }
    /*calcWallFaceZoneHeatFlux();
    calcWallFaceZoneFrictionCoefficient();
    calcWallFaceZoneAbsPressCoefficient();
    calcWallFaceZoneRelPressCoefficient();
    calcWallFaceZoneHeatCoefficient();
    calcWallYPlus();
    calcWallUPlus();*/

    if (samplingNow()) {
        sampling();
    }
}

OpenHurricane::fileName OpenHurricane::writeFaceZone::getFileName(const string &zone) const {
    const auto &iter = flows_.mesh().Iteration();
    fileName outN = iter.outputName();

    auto pathOut = outN.parentPath();
    auto fname = outN.name(true);
    string fext;
    fext = ".plt";

    fname += "-";
    fname += zone;
    fname += "-";
    fname += toString(iter.totalStep());
    outN = fname + fext;
    outN = pathOut / outN;
    return outN;
}

void OpenHurricane::writeFaceZone::setFaceZoneVarMap() {
    const auto &mesh = flows_.mesh();
    const auto &fz = mesh.faceZones();
    for (integer i = 0; i < writeFaceZoneIdList_.size(); ++i) {
        const auto fzi = writeFaceZoneIdList_[i];
        if (fzi != -1) {
            std::set<std::string> fzSet;
            for (integer j = 0; j < writeFaceZoneVarList_[i].size(); ++j) {
                fzSet.emplace(writeFaceZoneVarList_[i][j]);
            }
            integer id = -1;
            for (integer k = 0; k < fz.size(); ++k) {
                if (fzi == fz[k].index()) {
                    id = k;
                    break;
                }
            }
            if (id == -1) {
                LFatal("Cannot find face zone with index: %d", fzi);
            }
            faceZoneVarMap_.emplace(fz[id].name(), fzSet);
        }
    }
}

void OpenHurricane::writeFaceZone::clearOutFieldVarMap() const {
    for (auto &iter : outFieldVarMap_) {
        HurDelete(iter.second);
    }

    for (auto &iter : outSamplingFieldVarMap_) {
        HurDelete(iter.second);
    }
}

void OpenHurricane::writeFaceZone::resetOutFieldVarMap() const {
    clearOutFieldVarMap();
    outFieldVarMap_.clear();
}

OpenHurricane::writeFaceZone::writeFaceZone(const flowModel &flows, const controller &cont)
    : flows_(flows), cont_(cont), writeFaceZoneNameList_(), writeFaceZoneVarList_(),
      writeFaceZoneIdList_(), outFieldVarMap_(), outSamplingFieldVarMap_(), faceZoneVarMap_(),
      faceZoneVarFuncMap_(), isSampling_(false), samplingStep_(10), startTime_(0), timeElasped_(0),
      lastTime_(0), firstCalling_(true) {
    setFaceZoneList();
    setFaceZoneId();
    setFaceZoneVarMap();
    calculateFaceZoneVar myffc(faceZoneVarFuncMap_);

    if (flows_.mesh().Iteration().hasPhysicalTimeStep()) {
        if (cont.found("writeControl")) {
            const auto &writeControlCont = cont.subController("writeControl");
            startTime_ = writeControlCont.findOrDefault<real>("startTime", startTime_);

            if (writeControlCont.found("sampling")) {
                const auto &samCont = writeControlCont.subController("sampling");
                controllerSwitch issamCont(samCont);
                isSampling_ = issamCont("isSampling", isSampling_);
                samplingStep_ = samCont.findOrDefault<integer>("samplingStep", samplingStep_);
            }
        }
    }
}

void OpenHurricane::writeFaceZone::sampling() const {
    const auto &iter = flows_.mesh().Iteration();
    const auto cTT = iter.pTStep().totalTime();
    if (cTT <= startTime_) {
        timeElasped_ = 0;
        lastTime_ = cTT;
        return;
    }
    const auto dt = cTT - lastTime_;

    lastTime_ = cTT;

    timeElasped_ += dt;

    for (integer i = 0; i < writeFaceZoneIdList_.size(); ++i) {
        if (writeFaceZoneIdList_[i] == -1) {
            continue;
        }
        const auto &mesh = flows_.mesh();

        auto &outVarNL = writeFaceZoneVarList_[i];
        auto fzid = writeFaceZoneIdList_[i];

        integer id = -1;
        for (integer j = 0; j < mesh.faceZones().size(); ++j) {
            if (fzid == mesh.faceZones()[j].index()) {
                id = j;
                break;
            }
        }
        if (id == -1) {
            LFatal("Cannot find face zone with index: %d", fzid);
        }

        checkDuplicate(mesh.faceZones()[id], outVarNL);
        for (integer j = 0; j < outVarNL.size(); ++j) {
            if (mesh.foundOnlyObject(outVarNL[j])) {
                samplingFromField(mesh.faceZones()[id], outVarNL[j], dt, id);
            } else {
                string outFzN = outVarNL[j] + mesh.faceZones()[id].name();
                if (mesh.foundOnlyObject(outFzN)) {
                    samplingFromMap(mesh.faceZones()[id], outVarNL[j], dt, id);
                }
            }
        }
    }
}

void OpenHurricane::writeFaceZone::samplingFromField(const faceZone &fz, const string &varN,
                                                     const real dt, const integer fzid) const {
    const auto &mesh = flows_.mesh();
    const auto &ob = mesh.findOnlyObject(varN);
    if (ob.nElements() == 1) {
        const auto &obf = mesh.findObject<cellRealArray>(varN);
        auto inter = outSamplingFieldVarMap_.find(varN + "_sampling" + fz.name());
        if (inter != outSamplingFieldVarMap_.end()) {
            auto &obfav = static_cast<realFaceZoneArray &>(*inter->second);
            obfav = fv::interpolate(obf, fzid);
            obfav.calcTimeSumPtr(dt);
            obfav = obfav.timeSum() / timeElasped_;
        } else {
            outSamplingFieldVarMap_.emplace(
                varN + "_sampling" + fz.name(),
                new realFaceZoneArray(fz, varN + "_sampling", mesh, fv::interpolate(obf, fzid)));
            auto inter2 = outSamplingFieldVarMap_.find(varN + "_sampling" + fz.name());
            (*inter2->second).calcTimeSumPtr(dt);
        }
    } else if (ob.nElements() == 3) {
        const auto &obf = mesh.findObject<cellVectorArray>(varN);
        auto inter = outSamplingFieldVarMap_.find(varN + "_sampling" + fz.name());
        if (inter != outSamplingFieldVarMap_.end()) {
            auto &obfav = static_cast<vectorFaceZoneArray &>(*inter->second);
            obfav = fv::interpolate(obf, fzid);
            obfav.calcTimeSumPtr(dt);
            obfav = obfav.timeSum() / timeElasped_;
        } else {
            outSamplingFieldVarMap_.emplace(
                varN + "_sampling" + fz.name(),
                new vectorFaceZoneArray(fz, varN + "_sampling", mesh, fv::interpolate(obf, fzid)));
            auto inter2 = outSamplingFieldVarMap_.find(varN + "_sampling" + fz.name());
            (*inter2->second).calcTimeSumPtr(dt);
        }
    } else {
        LFatal("The type of face zone array, which has %d components, is not support yet.",
               ob.nElements());
    }
}

void OpenHurricane::writeFaceZone::samplingFromMap(const faceZone &fz, const string &varN,
                                                   const real dt, const integer fzid) const {
    const auto &mesh = flows_.mesh();
    const auto &ob = mesh.findOnlyObject(varN + fz.name());
    if (ob.nElements() == 1) {
        const auto &obf = mesh.findObject<realFaceZoneArray>(varN + fz.name());
        auto inter = outSamplingFieldVarMap_.find(varN + "_sampling" + fz.name());
        if (inter != outSamplingFieldVarMap_.end()) {
            auto &obfav = static_cast<realFaceZoneArray &>(*inter->second);
            obfav = obf;
            obfav.calcTimeSumPtr(dt);
            obfav = obfav.timeSum() / timeElasped_;
        } else {
            outSamplingFieldVarMap_.emplace(
                varN + "_sampling" + fz.name(),
                new realFaceZoneArray(fz, varN + "_sampling", mesh, obf));
            auto inter2 = outSamplingFieldVarMap_.find(varN + "_sampling" + fz.name());
            (*inter2->second).calcTimeSumPtr(dt);
        }
    } else if (ob.nElements() == 3) {
        const auto &obf = mesh.findObject<vectorFaceZoneArray>(varN + fz.name());
        auto inter = outSamplingFieldVarMap_.find(varN + "_sampling" + fz.name());
        if (inter != outSamplingFieldVarMap_.end()) {
            auto &obfav = static_cast<vectorFaceZoneArray &>(*inter->second);
            obfav = obf;
            obfav.calcTimeSumPtr(dt);
            obfav = obfav.timeSum() / timeElasped_;
        } else {
            outSamplingFieldVarMap_.emplace(
                varN + "_sampling" + fz.name(),
                new vectorFaceZoneArray(fz, varN + "_sampling", mesh, obf));
            auto inter2 = outSamplingFieldVarMap_.find(varN + "_sampling" + fz.name());
            (*inter2->second).calcTimeSumPtr(dt);
        }
    } else {
        LFatal("The type of face zone array, which has %d components, is not support yet.",
               ob.nElements());
    }
}

bool OpenHurricane::writeFaceZone::findInMap(const string &fzName, const string &varName) const {
    const auto iter = faceZoneVarMap_.find(fzName);
    if (iter != faceZoneVarMap_.end()) {
        const auto iter2 = iter->second.find(varName);
        if (iter2 != iter->second.end()) {
            return true;
        }
    }
    return false;
}

void OpenHurricane::writeFaceZone::setOutputField(const faceZone &fz, const string &varName,
                                                  const realArray &value) const {
    const auto &mesh = flows_.mesh();

    integer id = -1;
    for (integer i = 0; i < writeFaceZoneIdList_.size(); ++i) {
        if (fz.index() == writeFaceZoneIdList_[i]) {
            id = i;
            break;
        }
    }

    if (id == -1) {
        LFatal("Cannot find face zone: %s in writing zone list", fz.name().c_str());
    }
    auto inter = outFieldVarMap_.find(varName + fz.name());

    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new realFaceZoneArray(fz, varName, mesh, value);
    } else {
        outFieldVarMap_.emplace(varName + fz.name(),
                                new realFaceZoneArray(fz, varName, mesh, value));
    }
}

void OpenHurricane::writeFaceZone::setOutputField(const faceZone &fz, const string &varName,
                                                  realArray &&value) const {
    const auto &mesh = flows_.mesh();
    integer id = -1;
    for (integer i = 0; i < writeFaceZoneIdList_.size(); ++i) {
        if (fz.index() == writeFaceZoneIdList_[i]) {
            id = i;
            break;
        }
    }

    if (id == -1) {
        LFatal("Cannot find face zone: %s in writing zone list", fz.name().c_str());
    }
    auto inter = outFieldVarMap_.find(varName + fz.name());

    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new realFaceZoneArray(fz, varName, mesh, std::move(value));
    } else {
        outFieldVarMap_.emplace(varName + fz.name(),
                                new realFaceZoneArray(fz, varName, mesh, std::move(value)));
    }
}

void OpenHurricane::writeFaceZone::setOutputField(const faceZone &fz, const string &varName,
                                                  const vectorArray &value) const {
    const auto &mesh = flows_.mesh();
    integer id = -1;
    for (integer i = 0; i < writeFaceZoneIdList_.size(); ++i) {
        if (fz.index() == writeFaceZoneIdList_[i]) {
            id = i;
            break;
        }
    }

    if (id == -1) {
        LFatal("Cannot find face zone: %s in writing zone list", fz.name().c_str());
    }

    auto inter = outFieldVarMap_.find(varName + fz.name());
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new vectorFaceZoneArray(fz, varName, mesh, value);
    } else {
        outFieldVarMap_.emplace(varName + fz.name(),
                                new vectorFaceZoneArray(fz, varName, mesh, value));
    }
}

void OpenHurricane::writeFaceZone::setOutputField(const faceZone &fz, const string &varName,
                                                  vectorArray &&value) const {
    const auto &mesh = flows_.mesh();
    integer id = -1;
    for (integer i = 0; i < writeFaceZoneIdList_.size(); ++i) {
        if (fz.index() == writeFaceZoneIdList_[i]) {
            id = i;
            break;
        }
    }

    if (id == -1) {
        LFatal("Cannot find face zone: %s in writing zone list", fz.name().c_str());
    }

    auto inter = outFieldVarMap_.find(varName + fz.name());
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new vectorFaceZoneArray(fz, varName, mesh, std::move(value));
    } else {
        outFieldVarMap_.emplace(varName + fz.name(),
                                new vectorFaceZoneArray(fz, varName, mesh, std::move(value)));
    }
}
//
//void OpenHurricane::writeFaceZone::calcWallFaceZoneHeatFlux() const
//{
//	const auto& mesh = flows_.mesh();
//	for (integer fzi = 0; fzi < mesh.faceZones().size(); ++fzi)
//	{
//		if (mesh.faceZones()[fzi].isWall())
//		{
//			const auto& fz = mesh.faceZones()[fzi];
//			if (findInMap(fz.name(), "Qw"))
//			{
//				setOutputField(fz, "Qw", calculateFieldVar::wallFaceZoneHeatFlux(flows_, fzi));
//			}
//		}
//	}
//}
//
//void OpenHurricane::writeFaceZone::calcWallFaceZoneFrictionCoefficient() const
//{
//	const auto& mesh = flows_.mesh();
//	for (integer fzi = 0; fzi < mesh.faceZones().size(); ++fzi)
//	{
//		if (mesh.faceZones()[fzi].isWall())
//		{
//			const auto& fz = mesh.faceZones()[fzi];
//			if (findInMap(fz.name(), "cf"))
//			{
//				setOutputField(fz, "cf", calculateFieldVar::wallFaceZoneFrictionCoefficient(flows_, fzi));
//			}
//		}
//	}
//}
//
//void OpenHurricane::writeFaceZone::calcWallFaceZoneAbsPressCoefficient() const
//{
//	const auto& mesh = flows_.mesh();
//	for (integer fzi = 0; fzi < mesh.faceZones().size(); ++fzi)
//	{
//		if (mesh.faceZones()[fzi].isWall())
//		{
//			const auto& fz = mesh.faceZones()[fzi];
//			if (findInMap(fz.name(), "cpAbs"))
//			{
//				setOutputField(fz, "cpAbs", calculateFieldVar::wallFaceZoneAbsPressCoefficient(flows_, fzi));
//			}
//		}
//	}
//}
//
//void OpenHurricane::writeFaceZone::calcWallFaceZoneRelPressCoefficient() const
//{
//	const auto& mesh = flows_.mesh();
//	for (integer fzi = 0; fzi < mesh.faceZones().size(); ++fzi)
//	{
//		if (mesh.faceZones()[fzi].isWall())
//		{
//			const auto& fz = mesh.faceZones()[fzi];
//			if (findInMap(fz.name(), "cpRel"))
//			{
//				setOutputField(fz, "cpRel", calculateFieldVar::wallFaceZoneRelPressCoefficient(flows_, fzi));
//			}
//		}
//	}
//}
//
//void OpenHurricane::writeFaceZone::calcWallFaceZoneHeatCoefficient() const
//{
//	const auto& mesh = flows_.mesh();
//	for (integer fzi = 0; fzi < mesh.faceZones().size(); ++fzi)
//	{
//		if (mesh.faceZones()[fzi].isWall())
//		{
//			const auto& fz = mesh.faceZones()[fzi];
//			if (findInMap(fz.name(), "ch"))
//			{
//				setOutputField(fz, "ch", calculateFieldVar::wallFaceZoneHeatCoefficient(flows_, fzi));
//			}
//		}
//	}
//}
//
//void OpenHurricane::writeFaceZone::calcWallYPlus() const
//{
//	const auto& mesh = flows_.mesh();
//	for (integer fzi = 0; fzi < mesh.faceZones().size(); ++fzi)
//	{
//		if (mesh.faceZones()[fzi].isWall())
//		{
//			const auto& fz = mesh.faceZones()[fzi];
//			if (findInMap(fz.name(), "YPlus"))
//			{
//				setOutputField(fz, "YPlus", calculateFieldVar::wallFaceZoneYPlus(flows_, fzi));
//			}
//		}
//	}
//}
//
//void OpenHurricane::writeFaceZone::calcWallUPlus() const
//{
//	const auto& mesh = flows_.mesh();
//	for (integer fzi = 0; fzi < mesh.faceZones().size(); ++fzi)
//	{
//		if (mesh.faceZones()[fzi].isWall())
//		{
//			const auto& fz = mesh.faceZones()[fzi];
//			if (findInMap(fz.name(), "uPlus"))
//			{
//				setOutputField(fz, "uPlus", calculateFieldVar::wallFaceZoneUPlus(flows_, fzi));
//			}
//		}
//	}
//}
