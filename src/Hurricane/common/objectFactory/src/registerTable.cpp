/*!
 * \file registerTable.cpp
 * \brief The subroutines and functions of register table
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

#include "registerTable.hpp"
#include "cellArrays.hpp"
#include "clock.hpp"
#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "iteration.hpp"
#include <iomanip>

hur_nodiscard bool OpenHurricane::registerTable::parentNotIteration() const {
    return (&parent_ != dynamic_cast<const registerTable *>(&iteration_));
}

hur_nodiscard const OpenHurricane::registerTable &OpenHurricane::registerTable::nullObject() {
    return NullRefObj::nullRef<registerTable>();
}

OpenHurricane::registerTable::registerTable()
    : object(*this, table_.get(), true), iteration_(iteration::nullObject()),
      parent_(registerTable::nullObject()) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating register table\n");
    }
#endif // HUR_LESS_LOGGER || HUR_FULL_LOGGER
}

OpenHurricane::registerTable::registerTable(const iteration &_iteration)
    : object(_iteration.pName(), _iteration, table_.get(), true), iteration_(_iteration),
      parent_(_iteration) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating register table: %s\n", _iteration.pName().c_str());
    }
#endif // HUR_LESS_LOGGER || HUR_FULL_LOGGER
}

OpenHurricane::registerTable::registerTable(const object &ob)
    : object(ob), iteration_(ob.Iteration()), parent_(ob.Iteration()), table_(nullptr) {}

OpenHurricane::registerTable::registerTable(object &&ob)
    : object(std::move(ob)), iteration_(object::Iteration()), parent_(object::Iteration()),
      table_(nullptr) {}

OpenHurricane::registerTable::~registerTable() noexcept {
    if (!table_) {
        return;
    }
    List<object *> removeList(integer(table_->size()));
    integer removeSize = 0;
    for (typename tableType::iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->ownedByTable()) {
            removeList[removeSize++] = iter->second;
        }
    }

    for (integer i = 0; i < removeSize; ++i) {
        removeFromTable(*removeList[i]);
    }
    residuals_.unbind();
    residuals_.clear();
    table_.clear();
}

hur_nodiscard const OpenHurricane::iteration &OpenHurricane::registerTable::Iteration() const {
    return iteration_;
}

hur_nodiscard std::string OpenHurricane::registerTable::outputTitleNameDoc() const {
    std::string onStr;

    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->writeOption() == WRITE_OUTPUT ||
            iter->second->writeOption() == WRITE_RELAY_OUTPUT) {
            onStr += ",";
            onStr += iter->second->outputVarName();
        }
    }
    return onStr;
}

hur_nodiscard std::string
OpenHurricane::registerTable::outputTitleNameDoc(const stringList &outVarName) const {
    std::string onStr;
    for (integer i = 0; i < outVarName.size(); ++i) {
        auto iter = table_->find(outVarName[i]);
        if (iter != table_->end()) {
            onStr += ",";
            onStr += iter->second->outputVarName();
        }
    }

    return onStr;
}

hur_nodiscard OpenHurricane::stringList OpenHurricane::registerTable::outputTitleNameDocList() const {
    stringList onStr;

    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->writeOption() == WRITE_OUTPUT ||
            iter->second->writeOption() == WRITE_RELAY_OUTPUT) {
            onStr.append(iter->second->outputVarNameL());
        }
    }
    return onStr;
}

hur_nodiscard OpenHurricane::stringList
OpenHurricane::registerTable::outputTitleNameDocList(const stringList &outVarName) const {
    stringList onStr;

    for (integer i = 0; i < outVarName.size(); ++i) {
        auto iter = table_->find(outVarName[i]);
        if (iter != table_->end()) {
            onStr.append(iter->second->outputVarNameL());
        }
    }
    return onStr;
}

hur_nodiscard int OpenHurricane::registerTable::outputTitleSize() const {
    int is = 0;

    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->writeOption() == WRITE_OUTPUT ||
            iter->second->writeOption() == WRITE_RELAY_OUTPUT) {
            is += iter->second->nElements();
        }
    }
    return is;
}

hur_nodiscard int OpenHurricane::registerTable::outputTitleSize(const stringList &outVarName) const {
    int is = 0;
    for (integer i = 0; i < outVarName.size(); ++i) {
        auto iter = table_->find(outVarName[i]);
        if (iter != table_->end()) {
            is += iter->second->nElements();
        }
    }
    return is;
}

hur_nodiscard OpenHurricane::List<OpenHurricane::object *> OpenHurricane::registerTable::outRelayList() const {
    List<object *> outList;
    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->writeOption() == WRITE_RELAY ||
            iter->second->writeOption() == WRITE_RELAY_OUTPUT) {
            outList.append(iter->second);
        }
    }
    return outList;
}

hur_nodiscard OpenHurricane::List<OpenHurricane::object *>
OpenHurricane::registerTable::readRelayList(const stringList &varName) const {
    List<object *> rList;
    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        for (integer i = 0; i < varName.size(); ++i) {
            if (iter->second->name() == varName[i]) {
                rList.append(iter->second);
                break;
            }
        }
    }
    return rList;
}

hur_nodiscard int OpenHurricane::registerTable::outputRelaySize() const {
    int is = 0;
    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->writeOption() == WRITE_RELAY ||
            iter->second->writeOption() == WRITE_RELAY_OUTPUT) {
            is++;
        }
    }
    return is;
}

hur_nodiscard OpenHurricane::stringList OpenHurricane::registerTable::outputRelayNameDocList() const {
    stringList onStr;

    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->writeOption() == WRITE_RELAY ||
            iter->second->writeOption() == WRITE_RELAY_OUTPUT) {
            onStr.append(iter->second->name());
        }
    }
    return onStr;
}

hur_nodiscard OpenHurricane::List<OpenHurricane::object *> OpenHurricane::registerTable::outResultMap() const {
    List<object *> outList;
    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->writeOption() == WRITE_OUTPUT ||
            iter->second->writeOption() == WRITE_RELAY_OUTPUT) {
            outList.append(iter->second);
        }
    }
    return outList;
}

hur_nodiscard OpenHurricane::List<OpenHurricane::object *>
OpenHurricane::registerTable::outResultMap(const stringList &outVarName) const {
    List<object *> outList;
    for (integer i = 0; i < outVarName.size(); ++i) {
        auto iter = table_->find(outVarName[i]);
        if (iter != table_->end()) {
            outList.append(iter->second);
        }
    }

    return outList;
}

hur_nodiscard bool OpenHurricane::registerTable::checkOutVarName(const stringList &outVarName) const {
    for (integer i = 0; i < outVarName.size(); ++i) {
        auto iter = table_->find(outVarName[i]);
        if (iter == table_->end()) {
            return false;
        }
    }
    return true;
}

hur_nodiscard bool
OpenHurricane::registerTable::checkAndReportOutVarName(const stringList &outVarName) const {
    bool fl = true;
    for (integer i = 0; i < outVarName.size(); ++i) {
        auto iter = table_->find(outVarName[i]);
        if (iter == table_->end()) {
            checkWarningStr((" Variable: " + outVarName[i] + " is not exist"));
            fl = false;
        }
    }
    return fl;
}

hur_nodiscard typename OpenHurricane::registerTable::tableType
OpenHurricane::registerTable::primitiveParamMap() const {
    tableType pM;
    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->objectType() == object::PRIMITIVE) {
            pM.emplace(iter->first, iter->second);
        }
    }
    return pM;
}

hur_nodiscard OpenHurricane::integer OpenHurricane::registerTable::primitiveParamSize() const {
    integer n = 0;
    for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end(); ++iter) {
        if (iter->second->objectType() == object::PRIMITIVE) {
            n += iter->second->nElements();
        }
    }
    return n;
}

void OpenHurricane::registerTable::initResidualsList(const stringList &wl) const {
    for (integer i = 0; i < wl.size(); i++) {
        typename tableType::const_iterator iter = table_->find(wl[i]);
        if (iter != table_->end()) {
            residuals_.append(iter->second);
        }
    }
}

OpenHurricane::integerList
OpenHurricane::registerTable::initResidualsListAndCheck(stringList &wl, stringList &cmptName) const {
    stringList tmpWL(wl.size());
    integerList cmptNum(tmpWL.size(), 0);
    cmptName.resize(tmpWL.size());

    integer m = 0;
    for (integer i = 0; i < wl.size(); i++) {
        typename tableType::const_iterator iter = table_->find(wl[i]);
        if (iter != table_->end()) {
            cmptNum[m] = iter->second->nElements();
            cmptName[m] = iter->second->outputVarName();
            residuals_.append(iter->second);
            tmpWL[m] = wl[i];
            m++;
        }
    }
    if (m < wl.size()) {
        tmpWL.resize(m);
        cmptNum.resize(m);
        cmptName.resize(m);
        wl.transfer(tmpWL);
    }
    return cmptNum;
}

hur_nodiscard const OpenHurricane::registerTable &
OpenHurricane::registerTable::subRegisterTable(const string &name) const {
    return findObject<registerTable>(name);
}

bool OpenHurricane::registerTable::addToTable(object &ob) const {
    if (!table_) {
        table_.reset(new tableType());
    }
#ifdef HUR_DEBUG
    if (report) {
        LInfo("Adding a new object: \"%s\" to the register table: %s", ob.name().c_str(),
              name().c_str());
    }
#endif // HUR_DEBUG

    auto inserResult = table_->emplace(static_cast<std::string>(ob.name()), &ob);
    return inserResult.second;
}

bool OpenHurricane::registerTable::removeFromTable(object &ob) const {
    if (!table_) {
        return false;
    }
    typename tableType::iterator iter = const_cast<tableType &>(*table_).find(ob.name());

    if (iter != table_->end()) {
        if (iter->second != &ob) {
#ifdef HUR_DEBUG
            PLWarning(" Attempted to erase a copy object: \"%s\" in table: %s", ob.name().c_str(),
                      name().c_str());
#else
            if (report) {
                LWarning(" Attempted to erase a copy object: \"%s\" in table: %s",
                         ob.name().c_str(), name().c_str());
            }
#endif // HUR_DEBUG

            return false;
        } else {
            object *obj = iter->second;
            bool erased = const_cast<tableType &>(*table_).erase(iter->first);
            if (ob.ownedByTable()) {
                delete obj;
            }
            return erased;
        }
    }

    return false;
}

hur_nodiscard bool OpenHurricane::registerTable::foundOnlyObject(const string &name) const {
    typename tableType::const_iterator iter = table_->find(name);

    if (iter != table_->end()) {
        const object *Ptr_ = iter->second;
        if (Ptr_ != nullptr) {
            return true;
        }
    } else if (parentNotIteration()) {
        return parent_.foundOnlyObject(name);
    }
    return false;
}

hur_nodiscard const OpenHurricane::object &
OpenHurricane::registerTable::findOnlyObject(const string &name) const {
    typename tableType::const_iterator iter = table_->find(name);

    if (iter != table_->end()) {
        const object *Ptr_ = iter->second;
        if (Ptr_ != nullptr) {
            return *Ptr_;
        } else {
            LFatal("Null pointer of object: %s in %s", name.c_str(), this->name().c_str());
        }
    } else {
        if (parentNotIteration()) {
            return parent_.findOnlyObject(name);
        } else {
            LFatal("Cannot find object: \"%s\" in register table: %s.", name.c_str(),
                   this->name().c_str());
        }
    }
    return NullRefObj::nullRef<object>();
}

void OpenHurricane::registerTable::writeOutput(fileOsstream &fos) const {
    List<object *> outList = outResultMap();
    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeOutput(fos);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeOutputByMaster(fileOsstream &fos) const {
    List<object *> outList = outResultMap();

    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeOutputByMaster(fos);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeOutput(fileOsstream &fos, const stringList &outVarName) const {
    List<object *> outList = outResultMap(outVarName);

    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeOutput(fos);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeOutputByMaster(fileOsstream &fos,
                                                   const stringList &outVarName) const {
    List<object *> outList = outResultMap(outVarName);

    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeOutputByMaster(fos);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeOutput(fileOsstream &fos, const integer fzid) const {
    List<object *> outList = outResultMap();
    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeOutput(fos, fzid);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeOutput(fileOsstream &fos, const integer fzid,
                                           const stringList &outVarName) const {
    List<object *> outList = outResultMap(outVarName);
    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeOutput(fos, fzid);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeMinMaxOutput(fileOsstream &fos) const {
    List<object *> outList = outResultMap();
    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeMinMaxOutput(fos);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeMinMaxOutputByMaster(fileOsstream &fos) const {
    List<object *> outList = outResultMap();

    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeMinMaxOutputByMaster(fos);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeMinMaxOutput(fileOsstream &fos,
                                                 const stringList &outVarName) const {
    List<object *> outList = outResultMap(outVarName);

    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeMinMaxOutput(fos);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeMinMaxOutputByMaster(fileOsstream &fos,
                                                         const stringList &outVarName) const {
    List<object *> outList = outResultMap(outVarName);

    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeMinMaxOutputByMaster(fos);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeMinMaxOutput(fileOsstream &fos, const integer fzid,
                                                 const stringList &outVarName) const {
    List<object *> outList = outResultMap(outVarName);

    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeMinMaxOutput(fos, fzid);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeResiduals(fileOsstream &fos, integer intervalStep) const {
    for (integer i = 0; i < residuals_.size(); ++i) {
        residuals_(i)->writeResiduals(fos, intervalStep);
    }
}

void OpenHurricane::registerTable::writeResidualsName() const {
    for (integer i = 0; i < residuals_.size(); ++i) {
        Pout << residuals_[i].outputVarName().c_str() << "  ";
    }
}

void OpenHurricane::registerTable::writeRelay(fileOsstream &fos) const {
    auto outList = outRelayList();
    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeRelay(fos);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::writeRelay(hdf5O &fos, const bool writeLast,
                                          const bool writeToGroup) const {
    if (HurMPI::master()) {
        if (!writeToGroup) {
            fos.writeIntegerAttributeToFile(0, "nTimeGroups");
        } else if (!writeLast) {
            fos.writeIntegerAttributeToFile(1, "nTimeGroups");
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.createGroup(groupName);
        } else {
            const auto iter = table_->find("rho");
            integer ng = 0;
            if (iter != table_->end()) {
                auto ob = iter->second;
                ng = static_cast<const cellRealArray &>(*ob).lastArrayArray().size();
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.createGroup(groupName);

                for (integer kk = 0; kk < ng; ++kk) {
                    groupName = "timeGroupm";
                    groupName += toString(kk + 1);
                    fos.createGroup(groupName);
                }

                ob = nullptr;
            }
            fos.writeIntegerAttributeToFile(ng, "nTimeGroups");
        }
    }
    auto outList = outRelayList();
    for (integer i = 0; i < outList.size(); ++i) {
        outList[i]->writeRelay(fos, writeLast, writeToGroup);
        outList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::readRelay(const hdf5I &fos, const stringList &varN,
                                         const bool readLast, const bool readFromGroup) {
    auto rList = readRelayList(varN);
    for (integer i = 0; i < rList.size(); ++i) {
        rList[i]->readRelay(fos, readLast, readFromGroup);
        rList[i] = nullptr;
    }
}

void OpenHurricane::registerTable::interpolateRelay(const hdf5I &fos, const stringList &varN,
                                                const bool readLast, const bool readFromGroup) {
    auto rList = readRelayList(varN);
    if (HurMPI::master()) {
        Pout("    Interpolating solution from other calculation:        ");
    }
    for (integer i = 0; i < rList.size(); ++i) {
        rList[i]->interpolateRelay(fos, readLast, readFromGroup);
        rList[i] = nullptr;
        if (HurMPI::master()) {
            std::cout << "\b\b\b\b\b";
            std::cout << std::setw(4) << std::fixed << std::setprecision(1)
                      << (float)(i + 1) / rList.size() * 100 << "%";
        }
    }
    if (HurMPI::master()) {
        std::cout << std::endl << std::endl;
        std::cout.unsetf(std::ios::fixed);
        Pout.unsetReal();
    }
}
