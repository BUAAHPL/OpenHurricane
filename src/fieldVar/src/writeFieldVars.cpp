/*!
 * \file writeFieldVars.cpp
 * \brief Main subroutines for writeFieldVars.
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

#include "writeFieldVars.hpp"

OpenHurricane::writeFieldVars::writeFieldVars(const flowModel &flows_, const controller &cont_,
                                              const iteration &iter_)
    : flows_(flows_), cont_(cont_), iter_(iter_), timeSumVarList_(), endTimeAveCal_(0.0),
      writeFieldPtrList_(), outFieldVarMap_() {
    if (cont_.found("writeControl")) {
        const auto &iterCont = cont_.subController("writeControl");
        setWriteFieldVarList(iterCont);
    }
}

void OpenHurricane::writeFieldVars::setWriteFieldVarList(const controller &iterCont) {
    if (iterCont.found("writeList")) {
        std::string rnl = iterCont.findText("writeList");
        replaceAllMarks(rnl, "\n", " ");
        stringList writeList;
        if (!rnl.empty()) {
            size_t pos = 0;
            stdStringList rll;
            split(rnl, rll, ",");
            writeList.resize(rll.size());
            for (integer i = 0; i < rll.size(); ++i) {
                writeList[i] = trimCopy(rll[i]);
            }
        }
        removeDuplicate(writeList);
        if (writeFieldPtrList_.size() == 0) {
            for (integer i = 0; i < writeList.size(); ++i) {
                if (iterCont.found(writeList[i])) {
                    const auto &icont = iterCont.subController(writeList[i]);
                    writeFieldPtrList_.append(
                        writeFieldVar::creator(flows_, iter_, icont, writeList[i], outFieldVarMap_)
                            .release());
                }
            }
        } else {
            for (integer i = 0; i < writeFieldPtrList_.size(); ++i) {
                if (iterCont.found(writeList[i])) {
                    const auto &icont = iterCont.subController(writeList[i]);
                    integer m = 0;
                    bool found = false;
                    for (integer j = 0; j < writeList.size(); ++j) {
                        if (writeFieldPtrList_[j].writeId() == writeList[i]) {
                            found = true;
                            m = j;
                            break;
                        }
                    }
                    if (found) {
                        auto oldWrite = writeFieldPtrList_.set(
                            m, writeFieldVar::creator(flows_, iter_, icont, writeList[i],
                                                      outFieldVarMap_)
                                   .release());
                        HurDelete(oldWrite);
                    } else {
                        writeFieldPtrList_.append(writeFieldVar::creator(flows_, iter_, icont,
                                                                         writeList[i],
                                                                         outFieldVarMap_)
                                                      .release());
                    }
                }
            }
        }
    }
}

void OpenHurricane::writeFieldVars::write() const {
    for (integer i = 0; i < writeFieldPtrList_.size(); ++i) {
        writeFieldPtrList_[i].correct();
    }
    if (iter_.cStep() % iter_.writeOutputStep() == 0 ||
        (iter_.end() && !argParse::noWriteAtEnd())) {
        updating();

        for (integer i = 0; i < writeFieldPtrList_.size(); ++i) {
            writeFieldPtrList_[i].write();
        }

        for (integer i = 0; i < writeFieldPtrList_.size(); ++i) {
            writeFieldPtrList_[i].reseting();
        }

        clearoutVarMap();
    }
}

void OpenHurricane::writeFieldVars::setOutputField(const string &varName,
                                                   const realArray &value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellRealArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(), value);
    }
}

void OpenHurricane::writeFieldVars::setOutputField(const string &varName, realArray &&value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellRealArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(),
                                          std::move(value));
    }
}

void OpenHurricane::writeFieldVars::setOutputField(const string &varName,
                                                   const vectorArray &value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellVectorArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(), value);
    }
}

void OpenHurricane::writeFieldVars::setOutputField(const string &varName,
                                                   vectorArray &&value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellVectorArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(),
                                            std::move(value));
    }
}

void OpenHurricane::writeFieldVars::removeDuplicate(stringList &writeList) const {
    if (writeList.size() == 0) {
        return;
    }
    integer m = 0;
    stringList tmpList(writeList.size());

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

//void OpenHurricane::writeFieldVars::setTimeSumVarList(const stringList& varList)
//{
//	for (integer i = 0; i < varList.size(); ++i)
//	{
//		if (mesh().foundOnlyObject(getPrefix(varList[i])))
//		{
//			integer j = 0;
//			for (j; j < timeSumVarList_.size(); ++j)
//			{
//				if (getPrefix(varList[i]) == timeSumVarList_[j])
//				{
//					break;
//				}
//			}
//			if (j == timeSumVarList_.size())
//			{
//				timeSumVarList_.append(getPrefix(varList[i]));
//			}
//		}
//		else if (getPrefix(varList[i]) == "Ma")
//		{
//			stringList varnameList = { "v","gama","p","rho" };
//			addTimeSumPriVar(varnameList);
//		}
//		else if (getPrefix(varList[i]) == "Pt")
//		{
//			stringList varnameList = { "p","v","gama","p","rho" };
//			addTimeSumPriVar(varnameList);
//		}
//		else if (getPrefix(varList[i]) == "Tt")
//		{
//			if (flows_.mixtures().isSingular())
//			{
//				stringList varnameList = { "T","v","gama","p","rho" };
//				addTimeSumPriVar(varnameList);
//			}
//			else
//			{
//				stringList varnameList;
//				for (integer j = 0; j < flows_.mixtures().species().size(); ++j)
//				{
//					varnameList.append(flows_.mixtures().species()[i].name());
//				}
//				addTimeSumPriVar(varnameList);
//			}
//		}
//	}
//}
//
//void OpenHurricane::writeFieldVars::addTimeSumPriVar(const stringList& varnameList)
//{
//	for (integer i = 0; i < varnameList.size(); ++i)
//	{
//		integer j = 0;
//		for (j; j < timeSumVarList_.size(); ++j)
//		{
//			if (varnameList[i] == timeSumVarList_[j])
//			{
//				break;
//			}
//		}
//		if (j == timeSumVarList_.size())
//		{
//			timeSumVarList_.append(varnameList[i]);
//		}
//	}
//}
//
//void OpenHurricane::writeFieldVars::calcTimeSumVar() const
//{
//	bool calcSpecies = false;
//	bool isSpecies = false;
//	for (integer i = 0; i < timeSumVarList_.size(); ++i)
//	{
//		if (!calcSpecies)
//		{
//			for (integer j = 0; j < flows_.mixtures().species().size(); ++j)
//			{
//				if (timeSumVarList_[i] == flows_.mixtures().species()[j].name())
//				{
//					calcSpecies = true;
//					isSpecies = true;
//					for (integer k = 0; k < flows_.mixtures().Yi().size(); ++k)
//					{
//						flows_.mixtures().Yi()[k].calcTimeSumPtr(iter_.pTStep().pTimeStep());
//					}
//					break;
//				}
//			}
//		}
//		if (isSpecies)
//		{
//			isSpecies = false;
//			return;
//		}
//		for (auto iter = mesh().table().begin(); iter != mesh().table().end(); ++iter)
//		{
//			if (iter->first == timeSumVarList_[i])
//			{
//				if (iter->second->nElements() == 1)
//				{
//					cellRealArray* f = static_cast<cellRealArray*>(iter->second);
//					(*f).calcTimeSumPtr(iter_.pTStep().pTimeStep());
//				}
//				if (iter->second->nElements() == 3)
//				{
//					cellVectorArray* f = static_cast<cellVectorArray*>(iter->second);
//					(*f).calcTimeSumPtr(iter_.pTStep().pTimeStep());
//				}
//				break;
//			}
//		}
//	}
//}

void OpenHurricane::writeFieldVars::updating() const {
    if (writeFieldPtrList_.size() != 0) {
        for (integer i = 0; i < writeFieldPtrList_.size(); ++i) {
            writeFieldPtrList_[i].updating();
        }
    }
}

void OpenHurricane::writeFieldVars::clearoutVarMap() const {
    for (auto iter = outFieldVarMap_.begin(); iter != outFieldVarMap_.end(); iter++) {
        HurDelete(iter->second);
    }
}
