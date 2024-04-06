#include "controller.hpp"
/*!
 * \file controller.cpp
 * \brief Main subroutines of controller.
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

#include "controller.hpp"
#include "controllerContElement.hpp"
#include "fileOsstream.hpp"
#include "parameterContElement.hpp"
#include "realArrayContElement.hpp"
#include "Lists.hpp"
#include "textContElement.hpp"
#include "wordContElement.hpp"

const OpenHurricane::controller OpenHurricane::controller::null;

OpenHurricane::controller::controller(const fileName &name, const controller &parent)
    : parent_(parent), mapEntries_(), name_(name) {}

OpenHurricane::controller::controller(const controller &other)
    : parent_(other.parent_), mapEntries_(), name_(other.name_) {
    for (auto &e : other.mapEntries_) {
        mapEntries_.emplace(e.first, sharedPtr<controlElement>(e.second->clone(*this).release()));
    }
}
OpenHurricane::controller::controller(const controller &newParent, const controller &other)
    : parent_(newParent), mapEntries_(), name_(other.name_) {
    for (auto &e : other.mapEntries_) {
        mapEntries_.emplace(e.first, sharedPtr<controlElement>(e.second->clone(*this).release()));
    }
}

OpenHurricane::controller &OpenHurricane::controller::operator=(const controller &other) {
    if (this != std::addressof(other)) {
        name_ = other.name_;
        for (auto &e : other.mapEntries_) {
            mapEntries_.emplace(e.first,
                                sharedPtr<controlElement>(e.second->clone(*this).release()));
        }
    }
    return *this;
}

hur_nodiscard const OpenHurricane::controller &OpenHurricane::controller::topCont() const {
    const auto &pc = this->parent_;
    if (&pc != this && !pc.name().empty()) {
        return pc.topCont();
    } else {
        return *this;
    }
}

bool OpenHurricane::controller::merge(const controller &other) {
    if (this == std::addressof(other)) {
        LFatal("Attempt to merge a controller to itself: %s", name().c_str());
        return false;
    }

    bool isChanged = false;

    for (auto &e : other.mapEntries_) {
        auto mapIter = mapEntries_.find(e.first);

        if (mapIter != mapEntries_.end()) {
            if ((mapIter->second)->isControllerCE() && e.second->isControllerCE()) {
                if ((mapIter->second)->contContEle().merge(e.second->contContEle())) {
                    isChanged = true;
                }
            } else {
                add(e.second->clone(*this), true);
                isChanged = true;
            }
        } else {
            add(e.second->clone(*this));
            isChanged = true;
        }
    }

    return isChanged;
}

void OpenHurricane::controller::transfer(controller &other) {
    name_.swap(other.name());

    mapEntries_.swap(other.mapEntries_);
    other.clear();
}

bool OpenHurricane::controller::add(sharedPtr<controlElement> conElePtr, bool isMerge) {
    if (!conElePtr) {
        return false;
    }
    auto iter = mapEntries_.find(conElePtr->key());
    if (isMerge && iter != mapEntries_.end()) {
        if ((iter->second)->isControllerCE() && conElePtr->isControllerCE()) {
            (iter->second)->contContEle().merge(conElePtr->contContEle());
            return true;
        } else {
            mapEntries_.erase((iter->first));

            std::pair<typename map_type::iterator, bool> inserted =
                mapEntries_.insert(std::pair<std::string, sharedPtr<controlElement>>(
                    static_cast<std::string>(conElePtr->key()), conElePtr));
            if (inserted.second) {
                if (conElePtr->isControllerCE()) {
                    conElePtr->contContEle().name() = name() + '.' + conElePtr->key();
                }
                return true;
            } else {
                checkWarningStr(("Cannot replace control entry: " + conElePtr->key() +
                                 " in controller: " + name()));
                return false;
            }
        }
    }

    std::pair<typename map_type::iterator, bool> inserted =
        mapEntries_.insert(std::pair<std::string, sharedPtr<controlElement>>(
            static_cast<std::string>(conElePtr->key()), conElePtr));

    if (inserted.second) {
        if (conElePtr->isControllerCE()) {
            conElePtr->contContEle().name() = name() + '.' + conElePtr->key();
        }
        return true;
    } else {
        checkWarningStr(
            ("Cannot replace control entry: " + conElePtr->key() + " in controller: " + name()));
        return false;
    }

    return false;
}

bool OpenHurricane::controller::add(uniquePtr<controlElement> conElePtr, bool isMerge) {
    return add(sharedPtr<controlElement>(conElePtr.release()), isMerge);
}

bool OpenHurricane::controller::add(const controlElement &conEle, bool isMerge) {
    return add(conEle.clone(*this), isMerge);
}

bool OpenHurricane::controller::add(const std::string &key, const string &wValue, bool overwrite) {
    return add(sharedPtr<controlElement>(new wordContElement(key, wValue)), overwrite);
}

bool OpenHurricane::controller::addWord(const std::string &key, const string &wValue, bool overwrite) {
    return add(sharedPtr<controlElement>(new wordContElement(key, wValue)), overwrite);
}

bool OpenHurricane::controller::addWord(const std::string &key, string &&wValue, bool overwrite) {
    return add(sharedPtr<controlElement>(new wordContElement(key, std::move(wValue))), overwrite);
}

bool OpenHurricane::controller::addText(const std::string &key, const std::string &wValue,
                                    bool overwrite) {
    return add(sharedPtr<controlElement>(new textContElement(key, wValue)), overwrite);
}

bool OpenHurricane::controller::addText(const std::string &key, std::string &&wValue, bool overwrite) {
    return add(sharedPtr<controlElement>(new textContElement(key, std::move(wValue))), overwrite);
}

bool OpenHurricane::controller::addRealArray(const std::string &key, const realArray &raValue,
                                         bool overwrite) {
    return add(sharedPtr<controlElement>(new realArrayContElement(key, raValue)), overwrite);
}

bool OpenHurricane::controller::addRealArray(const std::string &key, realArray &&raValue,
                                         bool overwrite) {
    return add(sharedPtr<controlElement>(new realArrayContElement(key, std::move(raValue))),
               overwrite);
}

bool OpenHurricane::controller::add(const std::string &key, const std::string &pValue, bool overwrite) {
    return add(sharedPtr<controlElement>(new parameterContElement(key, pValue)), overwrite);
}

bool OpenHurricane::controller::addParameter(const std::string &key, const std::string &pValue,
                                         bool overwrite) {
    return add(sharedPtr<controlElement>(new parameterContElement(key, pValue)), overwrite);
}

bool OpenHurricane::controller::add(const std::string &key, const integer pValue, bool overwrite) {
    return add(sharedPtr<controlElement>(new parameterContElement(key, pValue)), overwrite);
}

bool OpenHurricane::controller::addParameter(const std::string &key, const integer pValue,
                                         bool overwrite) {
    return add(sharedPtr<controlElement>(new parameterContElement(key, pValue)), overwrite);
}

bool OpenHurricane::controller::add(const std::string &key, const real pValue, bool overwrite) {
    return add(sharedPtr<controlElement>(new parameterContElement(key, pValue)), overwrite);
}

bool OpenHurricane::controller::addParameter(const std::string &key, const real pValue,
                                         bool overwrite) {
    return add(sharedPtr<controlElement>(new parameterContElement(key, pValue)), overwrite);
}

bool OpenHurricane::controller::add(const std::string &key, const controller &cont, bool overwrite) {
    return add(sharedPtr<controlElement>(new controllerContElement(key, *this, cont)), overwrite);
}

bool OpenHurricane::controller::addController(const std::string &key, const controller &cont,
                                          bool overwrite) {
    return add(sharedPtr<controlElement>(new controllerContElement(key, *this, cont)), overwrite);
}

void OpenHurricane::controller::set(sharedPtr<controlElement> conElePtr) {
    auto existingCEPtr = findControlElePtr(conElePtr->key(), false);
    if (existingCEPtr && existingCEPtr->isControllerCE()) {
        existingCEPtr->contContEle().clear();
    }
    add(conElePtr, true);
}

void OpenHurricane::controller::set(const controlElement &conEle) {
    set(sharedPtr<controlElement>(conEle.clone(*this).release()));
}

void OpenHurricane::controller::set(const std::string &key, const controller &cont) {
    set(sharedPtr<controlElement>(new controllerContElement(key, *this, cont)));
}

void OpenHurricane::controller::setController(const std::string &key, const controller &cont) {
    set(sharedPtr<controlElement>(new controllerContElement(key, *this, cont)));
}

void OpenHurricane::controller::set(const std::string &key, const string &w) {
    set(sharedPtr<controlElement>(new wordContElement(key, w)));
}

void OpenHurricane::controller::setWord(const std::string &key, const string &w) {
    set(sharedPtr<controlElement>(new wordContElement(key, w)));
}

void OpenHurricane::controller::setText(const std::string &key, const std::string &w) {
    set(sharedPtr<controlElement>(new textContElement(key, w)));
}

void OpenHurricane::controller::setRealArray(const std::string &key, const realArray &w) {
    set(sharedPtr<controlElement>(new realArrayContElement(key, w)));
}

void OpenHurricane::controller::setRealArray(const std::string &key, realArray &&w) {
    set(sharedPtr<controlElement>(new realArrayContElement(key, std::move(w))));
}

bool OpenHurricane::controller::remove(const std::string &key) {
    auto iter = mapEntries_.find(key);
    if (iter != mapEntries_.end()) {
        iter->second.reset(nullptr);
        mapEntries_.erase(iter->first);
        return true;
    }
    return false;
}

bool OpenHurricane::controller::changeKey(const std::string &oldKey, const std::string &newKey,
                                      bool forceOverwrite) {
    if (oldKey == newKey) {
        return false;
    }

    auto iter = mapEntries_.find(oldKey);

    if (iter == mapEntries_.end()) {
        return false;
    }

    auto iter2 = mapEntries_.find(newKey);

    if (iter2 != mapEntries_.end()) {
        if (forceOverwrite) {
            iter2->second.reset(nullptr);
            mapEntries_.erase(iter2->first);
        } else {
            LFatal("Cannot change key:  %s with the new key %s in the controller: %s",
                   oldKey.c_str(), newKey.c_str(), name().c_str());
            return false;
        }
    }

    iter->second->key() = newKey;
    if (iter->second->isControllerCE()) {
        iter->second->contContEle().name() = name() + '.' + newKey;
    }
    mapEntries_.erase(oldKey);
    mapEntries_.emplace(newKey, iter->second);

    return true;
}

hur_nodiscard bool OpenHurricane::controller::found(const std::string &key, bool recursive) const {
    if (mapEntries_.find(key) != mapEntries_.end()) {
        return true;
    } else {
        if (recursive && &parent_ != &controller::null) {
            return parent_.found(key, recursive);
        } else {
            return false;
        }
    }
    return false;
}

hur_nodiscard const OpenHurricane::controlElement &
OpenHurricane::controller::findControlEle(const std::string &key, bool recursive) const {
    const auto cPtr = findControlElePtr(key, recursive);
    if (!cPtr) {
        LFatal("Cannot find option:  %s in the controller: %s", key.c_str(), name().c_str());
    }
    return *cPtr;
}

hur_nodiscard const OpenHurricane::sharedPtr<OpenHurricane::controlElement>
OpenHurricane::controller::findControlElePtr(const std::string &key, bool recursive) const {
    auto iter = mapEntries_.find(key);
    if (iter == mapEntries_.end()) {
        if (recursive && &parent_ != &controller::null) {
            return parent_.findControlElePtr(key, recursive);
        } else {
            return nullptr;
        }
    }
    return iter->second;
}

hur_nodiscard OpenHurricane::sharedPtr<OpenHurricane::controlElement>
OpenHurricane::controller::findControlElePtr(const std::string &key, bool recursive) {
    auto iter = mapEntries_.find(key);
    if (iter == mapEntries_.end()) {
        if (recursive && &parent_ != &controller::null) {
            return const_cast<controller &>(parent_).findControlElePtr(key, recursive);
        } else {
            return nullptr;
        }
    }
    return iter->second;
}

hur_nodiscard OpenHurricane::string OpenHurricane::controller::findWord(const std::string &key,
                                                                bool recursive) const {
    const auto &cE = findControlEle(key, recursive);
    if (!cE.isWordCE()) {
        LFatal("The value of option:  %s in the controller: %s is not a word element", key.c_str(),
               name().c_str());
    }
    return cE.wordContEle();
}

hur_nodiscard OpenHurricane::string OpenHurricane::controller::findWord(const std::string &key,
                                                                const string &defaultWord,
                                                                bool recursive) const {
    const auto cPtr = findControlElePtr(key, recursive);
    if (cPtr) {
        if (!cPtr->isWordCE()) {
            LFatal("The value of option:  %s in the controller: %s is not a word element",
                   key.c_str(), name().c_str());
        }
        return cPtr->wordContEle();
    }
    return defaultWord;
}

hur_nodiscard std::string OpenHurricane::controller::findText(const std::string &key,
                                                          bool recursive) const {
    const auto &cE = findControlEle(key, recursive);
    if (!cE.isTextCE()) {
        LFatal("The value of option:  %s in the controller: %s is not a text element", key.c_str(),
               name().c_str());
    }
    return cE.textContEle();
}

hur_nodiscard std::string OpenHurricane::controller::findText(const std::string &key,
                                                          const std::string &defaultText,
                                                          bool recursive) const {
    const auto cPtr = findControlElePtr(key, recursive);
    if (cPtr) {
        if (!cPtr->isTextCE()) {
            LFatal("The value of option:  %s in the controller: %s is not a text element",
                   key.c_str(), name().c_str());
        }
        return cPtr->textContEle();
    }
    return defaultText;
}

hur_nodiscard OpenHurricane::List<OpenHurricane::string>
OpenHurricane::controller::findTextStr(const std::string &key, bool recursive) const {
    std::string rnl = findText(key, recursive);
    List<string> rll;
    replaceAllMarks(rnl, "\n", " ");
    if (!rnl.empty()) {
        size_t pos = 0;
        split(rnl, rll, ",");
        for (integer i = 0; i < rll.size(); ++i) {
            trim(rll[i]);
        }
    }
    return rll;
}

hur_nodiscard std::string OpenHurricane::controller::findParameter(const std::string &key,
                                                               bool recursive) const {
    const auto &cE = findControlEle(key, recursive);
    if (!cE.isParameterCE()) {
        LFatal("The value of option:  %s in the controller: %s is not a parameter element",
               key.c_str(), name().c_str());
    }
    return cE.parameterContEle();
}

hur_nodiscard std::string OpenHurricane::controller::findParameter(const std::string &key,
                                                               const std::string &defaultParameter,
                                                               bool recursive) const {
    const auto cPtr = findControlElePtr(key, recursive);
    if (cPtr) {
        if (!cPtr->isParameterCE()) {
            LFatal("The value of option:  %s in the controller: %s is not a parameter element",
                   key.c_str(), name().c_str());
        }
        return cPtr->parameterContEle();
    }
    return defaultParameter;
}

hur_nodiscard const OpenHurricane::realArray &
OpenHurricane::controller::findRealArray(const std::string &key, bool recursive) const {
    const auto &cE = findControlEle(key, recursive);
    if (!cE.isRealArrayCE()) {
        LFatal("The value of option:  %s in the controller: %s is not a realArray element",
               key.c_str(), name().c_str());
    }
    return cE.rArrayContEle();
}

hur_nodiscard bool OpenHurricane::controller::isController(const std::string &key) const {
    const auto cEPtr = findControlElePtr(key, false);
    if (cEPtr) {
        return cEPtr->isControllerCE();
    }
    return false;
}

hur_nodiscard const OpenHurricane::controller &
OpenHurricane::controller::subController(const std::string &key) const {
    const auto cEPtr = findControlElePtr(key, false);
    if (!cEPtr) {
        errorAbortStr(
            ("The value of option: " + key + " can not be found in the controller: " + name()));
    } else if (!cEPtr->isControllerCE()) {
        errorAbortStr(("The value of option: " + key + " in the controller: " + name() +
                       " is not a controller"));
    }
    return cEPtr->contContEle();
}

hur_nodiscard OpenHurricane::controller &OpenHurricane::controller::subController(const std::string &key) {
    auto cEPtr = findControlElePtr(key, false);
    if (!cEPtr) {
        errorAbortStr(
            ("The value of option: " + key + " can not be found in the controller: " + name()));
    } else if (!cEPtr->isControllerCE()) {
        errorAbortStr(("The value of option: " + key + " in the controller: " + name() +
                       " is not a controller"));
    }
    return cEPtr->contContEle();
}

hur_nodiscard OpenHurricane::stringList OpenHurricane::controller::keysOfMap() const {
    stringList keys(static_cast<integer>(mapEntries_.size()));
    integer nKey = 0;

    for (const auto &e : mapEntries_) {
        keys[nKey++] = e.first;
    }
    return keys;
}

OpenHurricane::controller& OpenHurricane::controller::operator+=(const controller &other) {
    if (this != std::addressof(other)) {
        for (const auto &e : other.mapEntries_) {
            add(e.second->clone(*this));
        }
    }
    return *this;
}

void OpenHurricane::controller::replaceXMLComments(std::string &str) const {
    replaceAllMarks(str, "\r", " ");
    replaceAllMarks(str, "\t", " ");

    std::string::size_type pos = str.find("<!--", 0);
    std::string::size_type lastPos = str.find("-->", pos);
    while (pos != std::string::npos || lastPos != std::string::npos) {
        str.replace(pos, lastPos - pos + 3, " ");
        pos = str.find("<!--", pos + 1);
        lastPos = str.find("-->", pos);
    }
}

void OpenHurricane::controller::readXMLStatement(std::string &str) const {
    std::string::size_type pos = str.find("<?xml", 0);
    std::string::size_type lastPos = str.find("?>", pos);
    if (pos != std::string::npos || lastPos != std::string::npos) {
#ifdef HUR_DEBUG
        std::string stateStr = str.substr(pos, lastPos - pos + 2);
        Pout << stateStr << std::endl;
#endif // HUR_DEBUG
        str.replace(pos, lastPos - pos + 2, " ");
    }
}

void OpenHurricane::controller::parsingXML(std::string &str) {
    const std::regex eleFlag(
        "<((?!!\\[CDATA\\[).*?[^<>])(?:\\s+(?:(.*?[^<>])\\s*=\\s*(?:\"|\')(.*?)"
        "(?:\"|\')))*>");
    const std::regex nullEleFlag("<(\\w+)/>");

    const std::regex CDATA("<!\\[CDATA\\[");
    std::smatch what;
    auto beginIter = str.cbegin();
    auto endIter = str.cend();

    while (std::regex_search(beginIter, endIter, what, eleFlag)) {
        std::string newSubContName = what[1];
        std::string::const_iterator iter1 = what[0].second;

        if (std::regex_search(beginIter, iter1, what, nullEleFlag)) {
            beginIter = what[0].second;
            continue;
        }

        std::string eleBeginFlagStr = "<" + contNameForRegInXML(newSubContName) + "(?:\\s+|>)";
        const std::regex eleBeginFlag(eleBeginFlagStr.c_str());
        std::string eleEndFlagStr = "</" + contNameForRegInXML(newSubContName) + ">";
        const std::regex eleEndFlag(eleEndFlagStr.c_str());
        if (!std::regex_search(iter1, endIter, what, eleEndFlag)) {
            errorAbortStr(("The format of the control section is wrong: " + newSubContName +
                           " not balance in " + this->name()));
        }
        std::string::const_iterator iter2 = what[0].first;
        std::string::const_iterator iter2Erase = what[0].second;

        if (std::regex_search(iter1, iter2, eleBeginFlag)) {
            if (!std::regex_search(iter2Erase, endIter, what, eleEndFlag)) {
                errorAbortStr(("The format of the control section is wrong: " + newSubContName +
                               " not balance in " + this->name()));
            }
            iter2 = what[0].first;
            iter2Erase = what[0].second;
        }

        std::string subStr = std::string(iter1, iter2);
        if (newSubContName == "Config") {
            this->parsingXML(subStr);
        } else if (std::regex_search(subStr, what, eleFlag)) {
            controller newSubCont(newSubContName, *this);
            newSubCont.parsingXML(subStr);
            this->addController(newSubContName, newSubCont, true);
        } else if (std::regex_search(subStr, what, CDATA)) {
            auto pos = subStr.find("<![CDATA[", 0);
            auto lastPos = subStr.find("]]>", pos);

            std::string cdataStr = subStr.substr(pos + 9, lastPos - pos - 9);
            this->addText(newSubContName, cdataStr, true);
        } else {
            const std::regex equExp(
                "(?:(^(?:-?[0-9]+|\\.)\\.*[0-9]*(?:[eEgGdD][-+]?[0-9]*)*"
                "$)" // To match numbers (numerical digit), e.g. "10" or "1e1"
                "|(^\\s*\\S+?\\(\\s*(?:.*?)\\s*\\)\\s*$)" // To match the string in the form "string(...)", e.g. "car(1,0,0)"
                "|(?:^(\\w+)$))"                          // To match string, e.g. "controler"
            );
            replaceAllMarks(subStr, "\n", "");
            trim(subStr);
            std::string::const_iterator start = subStr.begin();
            std::string::const_iterator end = subStr.end();
            if (std::regex_search(start, end, what, equExp)) {
                if (what[1].matched) {
                    // To match numbers (numerical digit), e.g. "10" or "1e1"
                    std::string newP = what[1]; // The numerical digit must be added as string.
                    this->addParameter(newSubContName, newP, true);
                } else if (what[2].matched) {
                    // To match the string in the form "string(...)", e.g. "car(1,0,0)"
                    string newW = std::string(what[2]); // String must be added as string.
                    replaceAllMarks(newW, ";", " ");
                    trim(newW);
                    this->addWord(newSubContName, newW, true);
                } else if (what[3].matched) { // To match string, e.g. "controler"
                    string newW = std::string(what[3]);
                    replaceAllMarks(newW, ";", " ");
                    trim(newW);
                    this->addWord(newSubContName, newW, true);
                }
            } else {
                // To match the vector, e.g. (1,0,0)
                const std::regex vectorExp("(?:^\\s*?(\\([\\s\\S]*?\\))\\s*?$)");
                std::string::const_iterator start = subStr.begin();
                std::string::const_iterator end = subStr.end();
                if (std::regex_search(start, end, what, vectorExp)) {
                    std::string newW = subStr; // The vector must added as string.
                    replaceAllMarks(newW, "\n", " ");
                    this->addParameter(newSubContName, newW, true);
                } else {
                    string newW = subStr;
                    this->addWord(newSubContName, newW, true);
                }
            }
        }

        beginIter = iter2Erase;
    }
}

hur_nodiscard std::string OpenHurricane::controller::contNameForRegInXML(const std::string &na) const {
    std::string nam = na;
    replaceAllMarks(nam, "(", "\\(");
    replaceAllMarks(nam, ")", "\\)");
    replaceAllMarks(nam, "*", "\\*");
    replaceAllMarks(nam, ".", "\\.");
    replaceAllMarks(nam, "?", "\\?");
    return nam;
}

void OpenHurricane::controller::readXML(std::string &str) {
    HurMPIBase::bcastString(str, HurMPIBase::masterNo(), HurMPIBase::getComm());
    replaceXMLComments(str);
    readXMLStatement(str);
    parsingXML(str);
}

void OpenHurricane::controller::writeToXML(std::stringstream &sstr, const integer ilayer) const {
    integer ispace = 4 * ilayer;
    std::string tmpSpac;
    for (integer i = 0; i < ispace; ++i) {
        tmpSpac += " ";
    }

    for (auto iter = mapEntries_.begin(); iter != mapEntries_.end(); ++iter) {
        sstr << tmpSpac.c_str() << "<" << iter->first.c_str() << ">" << std::endl;
        iter->second->writeToXML(sstr, ilayer + 1);
        sstr << tmpSpac.c_str() << "</" << iter->first.c_str() << ">" << std::endl;
    }
}

void OpenHurricane::controller::write(fileOsstream &fos, bool onlyMasterNode) const {
    if (onlyMasterNode && !HurMPIBase::master()) {
        return;
    }
    std::stringstream sstr;
    write(sstr);
    fos.os() << sstr.str().c_str();
}

void OpenHurricane::controller::write(std::stringstream &sstr) const {
    sstr << "<?xml version=\"1.0\" encoding=\"UTF - 8\" standalone=\"no\"?>" << std::endl;
    sstr << "<Config>" << std::endl;
    writeToXML(sstr);
    sstr << "</Config>" << std::endl;
}
