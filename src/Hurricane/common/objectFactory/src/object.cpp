#include "object.hpp"
/*!
 * \file object.cpp
 * \brief The subroutines and functions of object
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

#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "iteration.hpp"
#include "object.hpp"
#include "realArray.hpp"
#include "registerTable.hpp"

OpenHurricane::object::object(const char *_c, const registerTable &_tb)
    : name_(_c), outputVarName_(_c), outputVarNameL_(1, _c), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(OTHERS), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", _c,
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const std::string &nam, const registerTable &_tb)
    : name_(nam), outputVarName_(nam), outputVarNameL_(1, nam), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(OTHERS), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(std::string &&nam, const registerTable &_tb)
    : name_(std::move(nam)), outputVarNameL_(1), tb_(_tb), registered_(false), ownedByTable_(false),
      writeOption_(NOT_WRITE), objectType_(OTHERS), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    outputVarName_ = name_;
    outputVarNameL_[0] = name_;
    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const char *_c, const registerTable &_tb, tableType *tab,
                          const bool setNullTablePtr)
    : name_(_c), outputVarName_(_c), outputVarNameL_(1, _c), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(OTHERS), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", _c,
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    if (setNullTablePtr) {
        const_cast<registerTable &>(tb_).unsafeSetNUllTable();
    }
    addToTable();
    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const std::string &nam, const registerTable &_tb, tableType *tab,
                          const bool setNullTablePtr)
    : name_(nam), outputVarName_(nam), outputVarNameL_(1, nam), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(OTHERS), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    if (setNullTablePtr) {
#if defined(HUR_FULL_LOGGER)
        if (report) {
            PLDebug("    Info: setting null table pointer for \"%s\".", name_.c_str());
        }
#endif // HUR_FULL_LOGGER
        const_cast<registerTable &>(tb_).unsafeSetNUllTable();

#if defined(HUR_FULL_LOGGER)
        if (report) {
            PLDebug("    Info: setting null table pointer for \"%s\". Done.", name_.c_str());
        }
#endif // HUR_FULL_LOGGER
    }

#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLDebug("    Info: adding \"%s\" to register table \"%s\".", name_.c_str(),
                _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLDebug("    Info: adding \"%s\" to register table \"%s\". Done.", name_.c_str(),
                _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const char *_c, const registerTable &_tb, const ObjectType ot)
    : name_(_c), outputVarName_(_c), outputVarNameL_(1, _c), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(ot), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", _c,
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const std::string &nam, const registerTable &_tb, const ObjectType ot)
    : name_(nam), outputVarName_(nam), outputVarNameL_(1, nam), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(ot), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const char *_c, const char *_oc, const registerTable &_tb)
    : name_(_c), outputVarName_(_oc), outputVarNameL_(1, _oc), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(OTHERS), hasSetOutputName_(true) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", _c,
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();

    std::string nl(_oc);
    stdStringList taken;
    split(nl, taken, ",");
    outputVarNameL_.resize(taken.size());
    for (integer i = 0; i < taken.size(); ++i) {
        replaceAllMarks(taken[i], "\"", "");
        outputVarNameL_[i] = taken[i];
    }
}

OpenHurricane::object::object(const std::string &nam, const std::string &outnam,
                          const registerTable &_tb)
    : name_(nam), outputVarName_(outnam), outputVarNameL_(1, outnam), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(OTHERS), hasSetOutputName_(true) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();

    std::string nl(outnam);
    stdStringList taken;
    split(nl, taken, ",");
    outputVarNameL_.resize(taken.size());
    for (integer i = 0; i < taken.size(); ++i) {
        replaceAllMarks(taken[i], "\"", "");
        outputVarNameL_[i] = taken[i];
    }
}

OpenHurricane::object::object(const char *_c, const char *_oc, const registerTable &_tb,
                          const ObjectType ot)
    : name_(_c), outputVarName_(_oc), outputVarNameL_(1, _oc), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(ot), hasSetOutputName_(true) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n",_c,
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    std::string nl(_oc);
    stdStringList taken;
    split(nl, taken, ",");
    outputVarNameL_.resize(taken.size());
    for (integer i = 0; i < taken.size(); ++i) {
        replaceAllMarks(taken[i], "\"", "");
        outputVarNameL_[i] = taken[i];
    }
}

OpenHurricane::object::object(const std::string &nam, const std::string &outnam,
                          const registerTable &_tb, const ObjectType ot)
    : name_(nam), outputVarName_(outnam), outputVarNameL_(1, outnam), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(NOT_WRITE), objectType_(ot), hasSetOutputName_(true) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    std::string nl(outnam);
    stdStringList taken;
    split(nl, taken, ",");
    outputVarNameL_.resize(taken.size());
    for (integer i = 0; i < taken.size(); ++i) {
        replaceAllMarks(taken[i], "\"", "");
        outputVarNameL_[i] = taken[i];
    }
}

OpenHurricane::object::object(const char *_c, const registerTable &_tb, const writeOptions wo)
    : name_(_c), outputVarName_(_c), outputVarNameL_(1, _c), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(wo), objectType_(OTHERS), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", _c,
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();

    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const std::string &nam, const registerTable &_tb, const writeOptions wo)
    : name_(nam), outputVarName_(nam), outputVarNameL_(1, nam), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(wo), objectType_(OTHERS), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();

    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const char *_c, const registerTable &_tb, const writeOptions wo,
                          const ObjectType ot)
    : name_(_c), outputVarName_(_c), outputVarNameL_(1, _c), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(wo), objectType_(ot), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", _c,
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const std::string &nam, const registerTable &_tb, const writeOptions wo,
                          const ObjectType ot)
    : name_(nam), outputVarName_(nam), outputVarNameL_(1, nam), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(wo), objectType_(ot), hasSetOutputName_(false) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    outputVarName_.insert(0, "\"");
    outputVarName_ += "\"";
}

OpenHurricane::object::object(const char *_c, const char *_oc, const registerTable &_tb,
                          const writeOptions wo)
    : name_(_c), outputVarName_(_oc), outputVarNameL_(1, _oc), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(wo), objectType_(OTHERS), hasSetOutputName_(true) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", _c,
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    std::string nl(_oc);
    stdStringList taken;
    split(nl, taken, ",");
    outputVarNameL_.resize(taken.size());
    for (integer i = 0; i < taken.size(); ++i) {
        replaceAllMarks(taken[i], "\"", "");
        outputVarNameL_[i] = taken[i];
    }
}

OpenHurricane::object::object(const std::string &nam, const std::string &outnam,
                          const registerTable &_tb, const writeOptions wo)
    : name_(nam), outputVarName_(outnam), outputVarNameL_(1, outnam), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(wo), objectType_(OTHERS), hasSetOutputName_(true) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();
    std::string nl(outnam);
    stdStringList taken;
    split(nl, taken, ",");
    outputVarNameL_.resize(taken.size());
    for (integer i = 0; i < taken.size(); ++i) {
        replaceAllMarks(taken[i], "\"", "");
        outputVarNameL_[i] = taken[i];
    }
}

OpenHurricane::object::object(const char *_c, const char *_oc, const registerTable &_tb,
                          const writeOptions wo, const ObjectType ot)
    : name_(_c), outputVarName_(_oc), outputVarNameL_(1, _oc), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(wo), objectType_(ot), hasSetOutputName_(true) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", _c,
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();

    std::string nl(_oc);
    stdStringList taken;
    split(nl, taken, ",");
    outputVarNameL_.resize(taken.size());
    for (integer i = 0; i < taken.size(); ++i) {
        replaceAllMarks(taken[i], "\"", "");
        outputVarNameL_[i] = taken[i];
    }
}

OpenHurricane::object::object(const std::string &nam, const std::string &outnam,
                          const registerTable &_tb, const writeOptions wo, const ObjectType ot)
    : name_(nam), outputVarName_(outnam), outputVarNameL_(1, outnam), tb_(_tb), registered_(false),
      ownedByTable_(false), writeOption_(wo), objectType_(ot), hasSetOutputName_(true) {
#if defined(HUR_FULL_LOGGER)
    if (report) {
        PLInfo("    Info: creating object \"%s\" in register table \"%s\"\n", name_.c_str(),
               _tb.name().c_str());
    }
#endif // HUR_FULL_LOGGER
    addToTable();

    std::string nl(outnam);
    stdStringList taken;
    split(nl, taken, ",");
    outputVarNameL_.resize(taken.size());
    for (integer i = 0; i < taken.size(); ++i) {
        replaceAllMarks(taken[i], "\"", "");
        outputVarNameL_[i] = taken[i];
    }
}

OpenHurricane::object::object(const object &ob)
    : name_(ob.name_), outputVarName_(ob.outputVarName_), outputVarNameL_(ob.outputVarNameL_),
      tb_(ob.tb_), registered_(false), ownedByTable_(false), writeOption_(ob.writeOption_),
      objectType_(ob.objectType_), hasSetOutputName_(ob.hasSetOutputName_) {}

OpenHurricane::object::object(const object &ob, tableType *tab, const bool setNullTablePtr)
    : name_(ob.name_), outputVarName_(ob.outputVarName_), outputVarNameL_(ob.outputVarNameL_),
      tb_(ob.tb_), registered_(false), ownedByTable_(false), writeOption_(ob.writeOption_),
      objectType_(ob.objectType_), hasSetOutputName_(ob.hasSetOutputName_) {
    if (setNullTablePtr) {
        const_cast<registerTable &>(tb_).unsafeSetNUllTable();
    }
}

OpenHurricane::object &OpenHurricane::object::operator=(const object &ob) {
    if (this != std::addressof(ob)) {
        removeFromTable();
        name_ = ob.name_;
        outputVarName_ = ob.outputVarName_;
        writeOption_ = ob.writeOption_;
        ownedByTable_ = false;
        addToTable();
    }
    return *this;
}

OpenHurricane::object::object(object &&ob) noexcept
    : name_(ob.name_), outputVarName_(std::move(ob.outputVarName_)),
      outputVarNameL_(std::move(ob.outputVarNameL_)), tb_(ob.tb_), registered_(false),
      ownedByTable_(false), writeOption_(std::move(ob.writeOption_)),
      objectType_(std::move(ob.objectType_)), hasSetOutputName_(std::move(ob.hasSetOutputName_)) {
    if (ob.registered_) {
        ob.removeFromTable();
        addToTable();
    }
}

OpenHurricane::object::object(object &&ob, tableType *tab, const bool setNullTablePtr) noexcept
    : name_(ob.name_), outputVarName_(std::move(ob.outputVarName_)),
      outputVarNameL_(std::move(ob.outputVarNameL_)), tb_(ob.tb_), registered_(false),
      ownedByTable_(false), writeOption_(ob.writeOption_), objectType_(ob.objectType_),
      hasSetOutputName_(ob.hasSetOutputName_) {
    if (setNullTablePtr) {
        const_cast<registerTable &>(tb_).unsafeSetNUllTable();
    }
    if (ob.registered_) {
        const_cast<object &>(ob).removeFromTable();
        addToTable();
    }
}

OpenHurricane::object &OpenHurricane::object::operator=(object &&ob) noexcept {
    removeFromTable();
    name_ = std::move(ob.name_);
    outputVarName_ = std::move(ob.outputVarName_);
    writeOption_ = std::move(ob.writeOption_);

    if (ob.registered_) {
        ob.removeFromTable();
    }

    addToTable();

    return *this;
}

OpenHurricane::object::object(const object &ob, bool registerCopy)
    : name_(ob.name_), outputVarName_(ob.outputVarName_), outputVarNameL_(ob.outputVarNameL_),
      tb_(ob.tb_), registered_(false), ownedByTable_(false), writeOption_(ob.writeOption_),
      objectType_(ob.objectType_), hasSetOutputName_(ob.hasSetOutputName_) {
    if (registerCopy && ob.registered_) {
        const_cast<object &>(ob).removeFromTable();
        addToTable();
    }
}

OpenHurricane::object::~object() noexcept {
    if (!ownedByTable_) {
        removeFromTable();
    }
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::object::realComponent(const int i) const {
    return realArray();
}

hur_nodiscard const OpenHurricane::registerTable &OpenHurricane::object::tb() const {
    return tb_;
}

hur_nodiscard const OpenHurricane::iteration &OpenHurricane::object::Iteration() const {
    return tb_.Iteration();
}

bool OpenHurricane::object::addToTable() {
    if (!registered_) {
        registered_ = tb().addToTable(*this);
        if (!registered_) {
#ifdef HUR_DEBUG
            std::string errMsg;
            errMsg = " Cannot add the object: \"";
            errMsg += name_;
            errMsg += "\" to the table: ";
            errMsg += tb().name();
            LFatal(errMsg.c_str());
#else
            if (report) {
                LError("Cannot add the object: \"%s\" to the table: %s", name_.c_str(),
                       tb().name().c_str());
            }

#endif // HUR_DEBUG
        }
    }
    return registered_;
}

bool OpenHurricane::object::removeFromTable() {
    if (registered_) {
        registered_ = false;
        return tb().removeFromTable(*this);
    }
    return false;
}

void OpenHurricane::object::clear() noexcept {
    if (!ownedByTable_) {
        removeFromTable();
    }
}

void OpenHurricane::object::writeOutput(fileOsstream &fos) const {
    if (writeOption_ == WRITE_OUTPUT || writeOption_ == WRITE_RELAY_OUTPUT) {
        fos.os() << name_.c_str() << std::endl;
    }
}

void OpenHurricane::object::writeOutputByMaster(fileOsstream &fos) const {
    if (HurMPIBase::master()) {
        if (writeOption_ == WRITE_OUTPUT || writeOption_ == WRITE_RELAY_OUTPUT) {
            fos.os() << name_.c_str() << std::endl;
        }
    }
}

void OpenHurricane::object::writeOutput(fileOsstream &fos, const integer fzid) const {
    if (writeOption_ == WRITE_OUTPUT || writeOption_ == WRITE_RELAY_OUTPUT) {
        fos.os() << name_.c_str() << std::endl;
    }
}

void OpenHurricane::object::writeMinMaxOutput(fileOsstream &fos) const {}

void OpenHurricane::object::writeMinMaxOutputByMaster(fileOsstream &fos) const {}

void OpenHurricane::object::writeMinMaxOutput(fileOsstream &fos, const integer fzid) const {}

void OpenHurricane::object::writeResiduals(fileOsstream &fos, integer intervalStep) const {}

void OpenHurricane::object::writeRelay(fileOsstream &fos) const {
    if (writeOption_ == WRITE_RELAY || writeOption_ == WRITE_RELAY_OUTPUT) {
        fos.os() << name_.c_str() << std::endl;
    }
}

void OpenHurricane::object::writeRelay(hdf5O &fos, const bool writeLast,
                                   const bool writeToGroup) const {}

void OpenHurricane::object::readRelay(const hdf5I &fos, const bool readLast, const bool readFromGroup) {
}

void OpenHurricane::object::calcTimeSumPtr(const real &dt) const {
    LFatal("This function must be implemented in derived class");
}

void OpenHurricane::object::interpolateRelay(const hdf5I &fos, const bool readLast,
                                         const bool readFromGroup) {}
