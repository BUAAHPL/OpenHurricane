/*!
 * \file controller.hpp
 * \brief Headers of controller.
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

#pragma once
#include "controlElement.hpp"
#include "parameterContElement.hpp"
#include "Lists.hpp"

namespace OpenHurricane {
    class fileOsstream;

    class controller {
    public:
        using map_type = std::map<std::string, sharedPtr<controlElement>>;

    private:
        const controller &parent_;

        map_type mapEntries_;
        fileName name_;

    public:
        static const controller null;

        inline controller() : parent_(null), mapEntries_() {}
        inline controller(const fileName &name) : parent_(null), mapEntries_(), name_(name) {}
        controller(const fileName &name, const controller &parent);

        controller(const controller &other);

        inline controller(const controller *otherPtr) : parent_(null), mapEntries_() {
            if (otherPtr != nullptr) {
                operator=(*otherPtr);
            }
        }

        inline controller(const uniquePtr<controller> &otherPtr) : parent_(null), mapEntries_() {
            if (otherPtr) {
                operator=(*otherPtr);
            }
        }

        inline controller(const sharedPtr<controller> &otherPtr) : parent_(null), mapEntries_() {
            if (otherPtr) {
                operator=(*otherPtr);
            }
        }

        controller(const controller &newParent, const controller &other);

        controller &operator=(const controller &other);

        inline virtual ~controller() noexcept { clear(); }

        hur_nodiscard inline const controller &parent() const noexcept { return parent_; }
        hur_nodiscard const controller &topCont() const;
        hur_nodiscard inline const fileName &name() const noexcept { return name_; }
        hur_nodiscard inline fileName &name() noexcept { return name_; }
        hur_nodiscard inline const map_type &mapEntries() const noexcept { return mapEntries_; }

        inline void clear() noexcept {
            mapEntries_.clear();
            name_.clear();
        }

        bool merge(const controller &other);
        void transfer(controller &other);

    private:
        bool add(sharedPtr<controlElement> conElePtr, bool isMerge = false);
        bool add(uniquePtr<controlElement> conElePtr, bool isMerge = false);
        bool add(const controlElement &conEle, bool isMerge = false);

    public:
        bool add(const std::string &key, const string &wValue, bool overwrite = false);
        bool addWord(const std::string &key, const string &wValue, bool overwrite = false);
        bool addWord(const std::string &key, string &&wValue, bool overwrite = false);

        bool addText(const std::string &key, const std::string &wValue, bool overwrite = false);
        bool addText(const std::string &key, std::string &&wValue, bool overwrite = false);

        bool addRealArray(const std::string &key, const realArray &raValue, bool overwrite = false);
        bool addRealArray(const std::string &key, realArray &&raValue, bool overwrite = false);

        bool add(const std::string &key, const std::string &pValue, bool overwrite = false);
        bool addParameter(const std::string &key, const std::string &pValue,
                          bool overwrite = false);

        bool add(const std::string &key, const integer pValue, bool overwrite = false);
        bool addParameter(const std::string &key, const integer pValue, bool overwrite = false);

        bool add(const std::string &key, const real pValue, bool overwrite = false);
        bool addParameter(const std::string &key, const real pValue, bool overwrite = false);

        bool add(const std::string &key, const controller &cont, bool overwrite = false);
        bool addController(const std::string &key, const controller &cont, bool overwrite = false);

        template <class Type>
        bool add(const std::string &key, const Type &t, bool overwrite = false) {
            return add(sharedPtr<controlElement>(new parameterContElement(key, t)), overwrite);
        }

    private:
        void set(sharedPtr<controlElement> conElePtr);
        void set(const controlElement &conEle);

    public:
        void set(const std::string &key, const controller &cont);
        void setController(const std::string &key, const controller &cont);

        void set(const std::string &key, const string &w);
        void setWord(const std::string &key, const string &w);

        void setText(const std::string &key, const std::string &w);

        void setRealArray(const std::string &key, const realArray &w);
        void setRealArray(const std::string &_key, realArray &&w);

        template <class Type> void set(const std::string &key, const Type &t) {
            set(sharedPtr<controlElement>(new parameterContElement(key, t)));
        }
        template <class Type> void setParameter(const std::string &key, const Type &t) {
            set(sharedPtr<controlElement>(new parameterContElement(key, t)));
        }

        bool remove(const std::string &key);

        bool changeKey(const std::string &oldKey, const std::string &newKey,
                       bool forceOverwrite = false);

        hur_nodiscard bool found(const std::string &key, bool recursive = false) const;

        hur_nodiscard const controlElement &findControlEle(const std::string &key,
                                                           bool recursive = false) const;

        hur_nodiscard const sharedPtr<controlElement>
        findControlElePtr(const std::string &key, bool recursive = false) const;

        hur_nodiscard sharedPtr<controlElement> findControlElePtr(const std::string &key,
                                                                  bool recursive = false);

        hur_nodiscard string findWord(const std::string &key, bool recursive = false) const;
        hur_nodiscard string findWord(const std::string &key, const string &defaultWord,
                                      bool recursive = false) const;
        hur_nodiscard inline string findWordOrDefault(const std::string &key,
                                                      const string &defaultWord,
                                                      bool recursive = false) const {
            return findWord(key, defaultWord, recursive);
        }

        hur_nodiscard std::string findText(const std::string &key, bool recursive = false) const;
        hur_nodiscard std::string findText(const std::string &key, const std::string &defaultText,
                                           bool recursive = false) const;
        hur_nodiscard List<string> findTextStr(const std::string &key,
                                               bool recursive = false) const;

        hur_nodiscard std::string findParameter(const std::string &key,
                                                bool recursive = false) const;
        hur_nodiscard std::string findParameter(const std::string &key,
                                                const std::string &defaultParameter,
                                                bool recursive = false) const;

        hur_nodiscard const realArray &findRealArray(const std::string &key,
                                                     bool recursive = false) const;

        template <class Type>
        hur_nodiscard Type find(const std::string &key, const Type &whatType,
                                bool recursive = false) const {
            const auto cPtr = findControlElePtr(key, recursive);
            if (!cPtr) {
                LFatal(" Cannot find the option: %s in the controller: %s", key.c_str(),
                       name().c_str());
            }
            return static_cast<Type>(feature<Type>(cPtr->ISStreamContEle()));
        }

        template <class Type>
        hur_nodiscard Type findType(const std::string &key, const Type &whatType,
                                    bool recursive = false) const {
            const auto cPtr = findControlElePtr(key, recursive);
            if (!cPtr) {
                LFatal(" Cannot find the option: %s in the controller: %s", key.c_str(),
                       name().c_str());
            }
            return static_cast<Type>(feature<Type>(cPtr->ISStreamContEle()));
        }

        template <class Type>
        hur_nodiscard Type findOrDefault(const std::string &key, const Type &defatValue,
                                         bool recursive = false) const {
            const auto cPtr = findControlElePtr(key, recursive);
            if (!cPtr) {
                return defatValue;
            }
            return static_cast<Type>(feature<Type>(cPtr->ISStreamContEle()));
        }

        template <class Type>
        Type findOrAddDefault(const std::string &key, const Type &defatValue,
                              bool recursive = false) const {
            const auto cPtr = findControlElePtr(key, recursive);
            if (!cPtr) {

                return defatValue;
            }
            return static_cast<Type>(feature<Type>(cPtr->ISStreamContEle()));
        }

        hur_nodiscard bool isController(const std::string &key) const;

        hur_nodiscard const controller &subController(const std::string &key) const;

        hur_nodiscard controller &subController(const std::string &key);

        hur_nodiscard const controller &subController(const char *const c) const {
            return subController(std::string(c));
        }

        hur_nodiscard controller &subController(const char *const c) {
            return subController(std::string(c));
        }

        hur_nodiscard stringList keysOfMap() const;

        OpenHurricane::controller &operator+=(const controller &other);

    private:
        void replaceXMLComments(std::string &str) const;
        void readXMLStatement(std::string &str) const;
        void parsingXML(std::string &str);
        hur_nodiscard std::string contNameForRegInXML(const std::string &nam) const;

    public:
        void readXML(std::string &str);

    protected:
        void writeToXML(std::stringstream &sstr, const integer ilayer = 1) const;

    public:
        void write(fileOsstream &fos, bool onlyMasterNode = true) const;
        void write(std::stringstream &sstr) const;
    };
} // namespace OpenHurricane
