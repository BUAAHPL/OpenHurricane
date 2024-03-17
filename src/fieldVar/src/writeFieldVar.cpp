/*!
 * \file writeFieldVar.cpp
 * \brief Main subroutines for writeFieldVar.
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

#include "writeFieldVar.hpp"
#include "calculateFieldVar.hpp"
#include "tecplotWriter.hpp"
namespace OpenHurricane {
    createClassNameStr(writeFieldVar, "writeFieldVar");
    createObjFty(writeFieldVar, controller);
} // namespace OpenHurricane

void OpenHurricane::writeFieldVar::removeDuplicate(stdStringList &writeList) const {
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

void OpenHurricane::writeFieldVar::getVarList(const controller &cont) {
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

            writeFieldVarList_.resize(rll.size());
            for (integer i = 0; i < rll.size(); ++i) {
                writeFieldVarList_[i] = trimCopy(rll[i]);
            }
        }
    } else {
        writeFieldVarList_ = mesh().outputTitleNameDocList();
    }
}

void OpenHurricane::writeFieldVar::setoutVarMap(const stringList &outFieldVarList) {
    for (integer i = 0; i < outFieldVarList.size(); ++i) {
        if (outFieldVarMap_.find(outFieldVarList[i]) != outFieldVarMap_.end()) {
            checkWarningStr(
                ("The variable " + outFieldVarList[i] + "has already been set in outFieldMap"));
        } else {
            outFieldVarMap_.emplace(outFieldVarList[i], nullptr);
        }
        if (sets_.find(outFieldVarList[i]) != sets_.end()) {
            checkWarningStr(
                ("The variable " + outFieldVarList[i] + "has already been set in outFieldMap"));
        } else {
            sets_.emplace(outFieldVarList[i]);
        }
    }
}

//void OpenHurricane::writeFieldVar::clearoutVarMap() const
//{
//    for (auto iter = outFieldVarMap_.begin(); iter != outFieldVarMap_.end(); iter++)
//    {
//        HurDelete(iter->second);
//    }
//}

OpenHurricane::writeFieldVar::writeFieldVar(const flowModel &flows, const iteration &iter,
                                            const controller &cont, const string &writeId,
                                            std::map<std::string, object *> &outFieldVarMap)
    : flows_(flows), iter_(iter), writeFieldVarList_(), sets_(), isWritting_(false),
      writeId_(writeId), yiAvePtr_(nullptr), outFieldVarMap_(outFieldVarMap) {
    getVarList(cont);
    isWritting_ = controllerSwitch(cont)("write", isWritting_);
}

OpenHurricane::uniquePtr<OpenHurricane::writeFieldVar>
OpenHurricane::writeFieldVar::creator(const flowModel &flows, const iteration &iter,
                                      const controller &cont, const string &writeId,
                                      std::map<std::string, object *> &outFieldVarMap) {
    string writeType = cont.findWord("writeType");

    defineInObjCreator(writeFieldVar, writeType, controller,
                       (flows, iter, cont, writeId, outFieldVarMap));
}

void OpenHurricane::writeFieldVar::fieldOtherVarsWriting(const string &type,
                                                         const combustionModel *chemtryPtr) {}

void OpenHurricane::writeFieldVar::fieldOtherVarsWriting(const string &type) {}

void OpenHurricane::writeFieldVar::calcTotalPressure(const string &Ma, const string &pt,
                                                     const realArray &p, const realArray &gama,
                                                     const vectorArray &v,
                                                     const realArray &rho) const {
    for (auto iter = outFieldVarMap_.begin(); iter != outFieldVarMap_.end(); iter++) {
        if (iter->first == Ma) {
            if (iter->second == nullptr) {
                setOutputField(Ma, calculateFieldVar::calcMachNumber(v, gama, p, rho));
                auto iter2 = outFieldVarMap_.find(Ma);
                cellRealArray *MaValue = static_cast<cellRealArray *>(iter2->second);
                setOutputField(pt, calculateFieldVar::calcTotalPressure(p, gama, *MaValue, flows_));
                return;
            } else {
                cellRealArray *MaValue = static_cast<cellRealArray *>(iter->second);
                setOutputField(pt, calculateFieldVar::calcTotalPressure(p, gama, *MaValue, flows_));
                return;
            }
        }
    }
    setOutputField(pt, calculateFieldVar::calcTotalPressure(p, gama, v, rho, flows_));
}

void OpenHurricane::writeFieldVar::calcTotalTemperature(const string &Ma, const string &Tt,
                                                        const realArray &T, const realArray &gama,
                                                        const vectorArray &v, const realArray &p,
                                                        const realArray &rho) const {
    for (auto iter = outFieldVarMap_.begin(); iter != outFieldVarMap_.end(); iter++) {
        if (iter->first == Ma) {
            if (iter->second == nullptr) {
                setOutputField(Ma, calculateFieldVar::calcMachNumber(v, gama, p, rho));
                auto iter2 = outFieldVarMap_.find(Ma);
                cellRealArray *MaValue = static_cast<cellRealArray *>(iter2->second);
                setOutputField(Tt,
                               calculateFieldVar::calcTotalTemperature(T, gama, *MaValue, flows_));
                return;
            } else {
                cellRealArray *MaValue = static_cast<cellRealArray *>(iter->second);
                setOutputField(Tt,
                               calculateFieldVar::calcTotalTemperature(T, gama, *MaValue, flows_));
                return;
            }
        }
    }
    setOutputField(Tt, calculateFieldVar::calcTotalTemperature(T, gama, v, p, rho, flows_));
}

//void OpenHurricane::writeFieldVar::setCurretTimeGap()
//{
//    curretTimeGap_ = iter_.pTStep().totalTime() - startTime_;
//}

//void OpenHurricane::writeFieldVar::calcTimeAveSpecies()
//{
//    yiAve_.resize(yi().size());
//    for (integer i = 0; i < yi().size(); ++i)
//    {
//        yiAve_[i] = yi()[i].getTimeSumPtr()/ curretTimeGap_;
//    }
//    integer fieldSize = yi()[0].mesh().nCells();
//    for (integer cellI = 0; cellI < fieldSize; ++cellI)
//    {
//        real sum = 0;
//        for (integer i = 0; i < yi().size(); ++i)
//        {
//            sum += yiAve_[i][cellI];
//        }
//        for (integer i = 0; i < yi().size(); ++i)
//        {
//            yiAve_[i][cellI] = yiAve_[i][cellI] / sum;
//        }
//    }
//    for (integer i = 0; i < yi().size(); ++i)
//    {
//        yiAve_[i] = const_cast<cellRealArray&>(yiAve_[i]);
//    }
//}

void OpenHurricane::writeFieldVar::correct() const {}

void OpenHurricane::writeFieldVar::reseting() const {
    HurDelete(yiAvePtr_);
}

void OpenHurricane::writeFieldVar::writeToFile() const {
    if (writeFieldVarList_.size() != 0) {
        const auto outN = getFileName();
        Pout << "    Writting result to file: " << outN << std::endl;
        fileOsstream fos(outN, IOsstream::BINARY_FORMAT);
        if (HurMPI::multiNodes() || argParse::tecplotFileWriteOnMaster()) {
            tecplotWriter::writeResultToTecplotByMaster(fos, mesh(), writeFieldVarList_);
        } else {
            tecplotWriter::writeResultToTecplot(fos, mesh(), writeFieldVarList_);
            //writeToTecplot(fos);
        }
        fos.close();
    }
}
