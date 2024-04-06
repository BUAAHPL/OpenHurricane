/*!
 * \file registerTable.hpp
 * \brief Header of register table.
 *       The subroutines and functions are in the <i>registerTable.cpp</i> file.
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
#pragma once

#include "object.hpp"
#include "smartPointerList.hpp"
#include "stdMaps.hpp"
namespace OpenHurricane {
    class iteration;
    class hdf5O;
    class hdf5I;

    /*!\brief The class of registerTable.*/
    class registerTable : public object {
    public:
        using tableType = std::unordered_map<std::string, object *>;

    private:
        // Private data

        /*!\brief Const reference to the iteration.*/
        const iteration &iteration_;

        /*!\brief Parent registerTable.*/
        const registerTable &parent_;

        /*!\brief If the parent is not an iteration.*/
        hur_nodiscard bool parentNotIteration() const;

        mutable uniquePtr<tableType> table_;

        mutable sharedPtrList<object> residuals_;

    public:
        hur_nodiscard static const registerTable &nullObject();

        registerTable();

        /*!\brief Construct from object.*/
        explicit registerTable(const iteration &_iteration);

        /*!\brief Construct from object.*/
        explicit registerTable(const object &ob);

        /*!\brief Construct from object.*/
        explicit registerTable(object &&ob);

        /*!\brief Destructor.*/
        virtual ~registerTable() noexcept;

        /*!\brief Iteration.*/
        hur_nodiscard const iteration &Iteration() const;

        hur_nodiscard inline tableType &table() noexcept { return *table_; }

        hur_nodiscard inline const tableType &table() const noexcept { return *table_; }

        /**
         * \brief Use with care.
         */
        inline auto unsafeSetNUllTable () noexcept-> tableType * { return table_.unsafeRelease(); }

        /*!\brief Parent registerTable.*/
        hur_nodiscard inline const registerTable &parent() const noexcept { return parent_; }

        /*!\brief Return the object name doc in the map.*/
        hur_nodiscard inline std::string nameDoc() const { return stringMapDoc(*table_); }

        /*!\brief Return the object output title name doc in the map.*/
        hur_nodiscard std::string outputTitleNameDoc() const;

        /*!\brief Return the object output title name doc in the map.*/
        hur_nodiscard std::string outputTitleNameDoc(const stringList &outVarName) const;

        /*!\brief Return the object output title name doc in the map.*/
        hur_nodiscard stringList outputTitleNameDocList() const;

        /*!\brief Return the object output title name doc in the map.*/
        hur_nodiscard stringList outputTitleNameDocList(const stringList &outVarName) const;

        /*!\brief Return the size of the object to be output in the map.*/
        hur_nodiscard int outputTitleSize() const;

        /*!\brief Return the size of the object to be output in the map.*/
        hur_nodiscard int outputTitleSize(const stringList &outVarName) const;

        /*!\brief Output relay object map.*/
        hur_nodiscard List<object *> outRelayList() const;
        hur_nodiscard List<object *> readRelayList(const stringList &varName) const;

        /*!\brief Return the size of the object to be output to relay file in the map.*/
        hur_nodiscard int outputRelaySize() const;

        /*!\brief Return the object output relay name doc in the map.*/
        hur_nodiscard stringList outputRelayNameDocList() const;

        /*!\brief Output result object map.*/
        hur_nodiscard List<object *> outResultMap() const;

        /*!\brief Output result object map.*/
        hur_nodiscard List<object *> outResultMap(const stringList &outVarName) const;

        /*!\brief Output result object map.*/
        hur_nodiscard bool checkOutVarName(const stringList &outVarName) const;

        /*!\brief Output result object map.*/
        hur_nodiscard bool checkAndReportOutVarName(const stringList &outVarName) const;

        /*!\brief Primitive paramaters object map.*/
        hur_nodiscard tableType primitiveParamMap() const;

        /*!\brief Primitive paramaters object map.*/
        hur_nodiscard integer primitiveParamSize() const;

        void initResidualsList(const stringList &wl) const;
        integerList initResidualsListAndCheck(stringList &wl, stringList &cmptName) const;

        hur_nodiscard inline bool isResidualsListSet() const noexcept {
            return residuals_.size() != 0;
        }

        /*!\brief Find and return a sub register table.*/
        hur_nodiscard const registerTable &subRegisterTable(const string &name) const;

        /*!\breif Add the object to the register table, and return true if success.*/
        bool addToTable(object &ob) const;

        /*!\brief Remove the object from the register table, and return true if success.*/
        bool removeFromTable(object &ob) const;

        /*!\brief Change the register table name.*/
        inline virtual void changeName(const string &nName) { object::changeName(nName); }

        /*!\brief Find and return the given Type.*/
        template <class Type>
        hur_nodiscard std::unordered_map<std::string, Type *>
        findClass(const bool strict = false) const {
            std::unordered_map<std::string, Type *> classMap;
            for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end();
                 ++iter) {
                if ((strict && isSameType<Type>(*iter->second)) ||
                    (!strict && hasSameBase<Type>(*iter->second))) {
                    classMap.emplace(iter->first, iter->second);
                }
            }
            return classMap;
        }

        /*!\brief Find and return the given Type.*/
        template <class Type>
        hur_nodiscard std::unordered_map<std::string, Type *> findClass(const bool strict = false) {
            std::unordered_map<std::string, Type *> classMap;
            for (typename tableType::const_iterator iter = table_->begin(); iter != table_->end();
                 ++iter) {
                if ((strict && isSameType<Type>(*iter->second)) ||
                    (!strict && hasSameBase<Type>(*iter->second))) {
                    classMap.emplace(iter->first, iter->second);
                }
            }
            return classMap;
        }

        /*!\brief Return true if the name has been found in the map.*/
        template <class Type> hur_nodiscard bool foundObject(const string &name) const {
            typename tableType::const_iterator iter = table_->find(name);
            if (iter != table_->end()) {
                const Type *Ptr_ = dynamic_cast<const Type *>(iter->second);
                if (Ptr_ != nullptr) {
                    return true;
                }
            } else if (parentNotIteration()) {
                return parent_.foundObject<Type>(name);
            }
            return false;
        }

        /*!\brief Find and return the object of given bcType*/
        template <class Type> hur_nodiscard const Type &findObject(const string &name) const {
            typename tableType::const_iterator iter = table_->find(name);
            if (iter != table_->end()) {
                const Type *Ptr_ = dynamic_cast<const Type *>(iter->second);
                if (Ptr_ != nullptr) {
                    return *Ptr_;
                } else {
                    LFatal("Finding object: \"%s\" in register table: %s successful. But it is not "
                           "a Type: %s.",
                           name.c_str(), this->name().c_str(), typeid(Type).name());
                }
            } else {
                if (parentNotIteration()) {
                    return parent_.findObject<Type>(name);
                } else {
                    LFatal("Cannot find object: \"%s\" in register table: %s.", name.c_str(),
                           this->name().c_str());
                }
            }
            return NullRefObj::nullRef<Type>();
        }

        /*!\brief Find and return the object of given bcType*/
        template <class Type> hur_nodiscard Type &findObjectRef(const string &name) const {
            return const_cast<Type &>(findObject<Type>(name));
        }

        /*!\brief Return true if the name has been found in the map.*/
        hur_nodiscard bool foundOnlyObject(const string &name) const;

        /*!\brief Find and return the object of given name*/
        hur_nodiscard const object &findOnlyObject(const string &name) const;

        // Wtite

        /*!
         * \brief Write object to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeOutput(fileOsstream &fos) const;

        /*!
         * \brief Write object to output file by master.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeOutputByMaster(fileOsstream &fos) const;

        /*!
         * \brief Write object to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeOutput(fileOsstream &fos, const stringList &outVarName) const;

        /*!
         * \brief Write object to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeOutputByMaster(fileOsstream &fos, const stringList &outVarName) const;

        /*!
         * \brief Write object to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeOutput(fileOsstream &fos, const integer fzid) const;

        /*!
         * \brief Write object to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeOutput(fileOsstream &fos, const integer fzid,
                                 const stringList &outVarName) const;

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeMinMaxOutput(fileOsstream &fos) const;

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeMinMaxOutputByMaster(fileOsstream &fos) const;

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeMinMaxOutput(fileOsstream &fos, const stringList &outVarName) const;

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeMinMaxOutputByMaster(fileOsstream &fos,
                                               const stringList &outVarName) const;

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeMinMaxOutput(fileOsstream &fos, const integer fzid,
                                       const stringList &outVarName) const;

        /*!
         * \brief Write object residuals to monitor.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeResiduals(fileOsstream &fos, integer intervalStep) const;

        /*!
         * \brief Write object residuals to monitor.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeResidualsName() const;

        /*!
         * \brief Write object to relay file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeRelay(fileOsstream &fos) const;

        /*!
         * \brief Write object to relay file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeRelay(hdf5O &fos, const bool writeLast, const bool writeToGroup) const;

        /*!
         * \brief Write object to relay file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void readRelay(const hdf5I &fos, const stringList &varN, const bool readLast,
                               const bool readFromGroup);

        /*!\brief Interpolate relay file to current geometryArray.*/
        virtual void interpolateRelay(const hdf5I &fos, const stringList &varN, const bool readLast,
                                      const bool readFromGroup);
    };
} // namespace OpenHurricane
