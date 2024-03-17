/*!
 * \file object.hpp
 * \brief Header of object.
 *       The subroutines and functions are in the <i>object.cpp</i> file.
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

#include "dataStructure.hpp"
#include "fileOsstream.hpp"
#include "logFile.hpp"
#include "objectFactory.hpp"

namespace OpenHurricane {
    class registerTable;
    class iteration;
    class hdf5O;
    class hdf5I;

    template<class Type>
    class Array;

    /*!\brief The class of object.*/
    class object {
    public:
        /*!\brief Write to file options enum.*/
        enum writeOptions : short {
            WRITE_RELAY,        /**< Write object to relay file.*/
            WRITE_OUTPUT,       /**< Write object to result output file.*/
            WRITE_RELAY_OUTPUT, /**< Write object to relay and result output file.*/
            NOT_WRITE           /**< Do not write object.*/
        };

        /*!\brief Write to file options enum.*/
        enum ObjectType : short {
            PRIMITIVE,     /**< primitive parameters.*/
            NOT_PRIMITIVE, /**< not primitive parameters.*/
            TEMPORARY,     /**< Temporary parameters.*/
            OTHERS,        /**< not specified paramaters.*/
        };

        using string_type = std::string;
        using obj_table = registerTable;
        using tableType = std::unordered_map<std::string, object *>;

    private:
        /*!\brief The object name.*/
        string_type name_;

        /*!\brief The object output var name.*/
        string_type outputVarName_;

        /*!\brief The object output var name list.*/
        List<string_type> outputVarNameL_;

        /*!\brief Const reference to the register table.*/
        const obj_table &tb_;

        /*!\brief If the object registered in the register table.*/
        bool registered_;

        /*!\brief If the object owned by register table.*/
        bool ownedByTable_;

        /*!\brief The option for writting to file*/
        writeOptions writeOption_;

        /*!\brief The bcType of the object*/
        ObjectType objectType_;

        /*!\brief If the output name of the object has been specified by manually.*/
        bool hasSetOutputName_;

    protected:
        /*!\brief Return true if the output name of the object has been specified by manually.*/
        hur_nodiscard inline bool &hasSetOutputName() noexcept { return hasSetOutputName_; }

    public:

        /*!\brief Construct from components.*/
        object(const char *_c, const registerTable &_tb);
        object(const std::string &nam, const registerTable &_tb);
        object(std::string &&nam, const registerTable &_tb);
        object(const char *_c, const registerTable &_tb, tableType *tab,
               const bool setNullTablePtr);
        object(const std::string &nam, const registerTable &_tb, tableType *tab,
               const bool setNullTablePtr);

        /*!\brief Construct from components.*/
        object(const char *_c, const registerTable &_tb, const ObjectType ot);
        object(const std::string &nam, const registerTable &_tb, const ObjectType ot);

        /*!\brief Construct from components.*/
        object(const char *_c, const char *_oc, const registerTable &_tb);
        object(const std::string &nam, const std::string &outnam, const registerTable &_tb);

        /*!\brief Construct from components.*/
        object(const char *_c, const char *_oc, const registerTable &_tb, const ObjectType ot);
        object(const std::string &nam, const std::string &outnam, const registerTable &_tb,
               const ObjectType ot);

        /*!\brief Construct from components.*/
        object(const char *_c, const registerTable &_tb, const writeOptions wo);
        object(const std::string &nam, const registerTable &_tb, const writeOptions wo);

        /*!\brief Construct from components.*/
        object(const char *_c, const registerTable &_tb, const writeOptions wo,
               const ObjectType ot);
        object(const std::string &nam, const registerTable &_tb, const writeOptions wo,
               const ObjectType ot);

        /*!\brief Construct from components.*/
        object(const char *_c, const char *_oc, const registerTable &_tb, const writeOptions wo);
        object(const std::string &nam, const std::string &outnam, const registerTable &_tb,
               const writeOptions wo);

        /*!\brief Construct from components.*/
        object(const char *_c, const char *_oc, const registerTable &_tb, const writeOptions wo,
               const ObjectType ot);
        object(const std::string &nam, const std::string &outnam, const registerTable &_tb,
               const writeOptions wo, const ObjectType ot);

        /*!\brief Construct as copy.*/
        object(const object &ob);
        object(const object &ob, tableType *tab, const bool setNullTablePtr);
        object &operator=(const object &ob);

        /*!\brief Construct as copy.*/
        object(object &&ob) noexcept;
        object(object &&ob, tableType *tab, const bool setNullTablePtr) noexcept;
        object &operator=(object &&ob) noexcept;

        /*!\brief Construct as copy, transferring register to copy if registerCopy is true.*/
        object(const object &ob, bool registerCopy);

        /*!\brief Destructor.*/
        virtual ~object() noexcept;

        /*!\brief Return const access to the object name.*/
        hur_nodiscard inline const string_type &name() const noexcept { return name_; }

        /*!\brief Return const access to the object output title name.*/
        hur_nodiscard inline virtual const string_type &outputVarName() const noexcept {
            return outputVarName_;
        }

        /*!\brief Return const access to the object output title name.*/
        hur_nodiscard inline string_type &outputVarName() noexcept { return outputVarName_; }

        /*!\brief Return access to the object output title name.*/
        hur_nodiscard inline List<string_type> &outputVarNameL() noexcept {
            return outputVarNameL_;
        }

        /*!\brief Return access to the object output title name.*/
        hur_nodiscard inline const List<string_type> &outputVarNameL() const noexcept {
            return outputVarNameL_;
        }

        /*!
         * \brief Return the number of components of object.
         *   For example, the number of components of vector is 3.
         */
        hur_nodiscard inline virtual int nElements() const noexcept { return 1; }

        hur_nodiscard virtual Array<real> realComponent(const int i) const;

        /*!\brief Return const reference to the register table.*/
        hur_nodiscard const registerTable &tb() const;

        /*!\brief Return const access to the iteration.*/
        hur_nodiscard const iteration &Iteration() const;

        /*!\brief If the object owned by register table.*/
        hur_nodiscard inline bool ownedByTable() const noexcept { return ownedByTable_; }

        /*!\brief The option for writting to file*/
        hur_nodiscard inline writeOptions writeOption() const noexcept { return writeOption_; }

        /*!\brief Set NO_WRITE option.*/
        inline void setNoWrite() noexcept { writeOption_ = NOT_WRITE; }

        /*!\brief Set WRITE_OUTPUT option.*/
        inline void setWriteResult() noexcept { writeOption_ = WRITE_OUTPUT; }

        /*!\brief The option for writting to file*/
        hur_nodiscard inline ObjectType objectType() const noexcept { return objectType_; }

        hur_nodiscard inline bool isTemporary() const noexcept { return objectType_ == TEMPORARY; }

        /*!\brief Change the object name.*/
        inline virtual void changeName(const string &nName) { name_ = nName; }

        /*!\brief Add the object to the register table.*/
        bool addToTable();

        /*!\brief Remove the object from the register table.*/
        bool removeFromTable();

        /*!\brief Change the ownership of thie object to its register table.*/
        inline void store() noexcept { ownedByTable_ = true; }

        /*!
         * \brief Change the ownership of thie object to its register table.
         * And return the reference to the object.
         */
        template <class Type> inline static Type &store(Type *tPtr) {
            if (!tPtr) {
                LFatal(" Attempted to get access to a deallocated object");
            }
            tPtr->object::ownedByTable_ = true;

            return *tPtr;
        }

        /*!
         * \brief Change the ownership of thie object to its register table.
         * And return the reference to the object.
         */
        template <class Type> inline static Type &store(uniquePtr<Type> &atPtr) {
            Type *tPtr = atPtr.release();

            if (!tPtr) {
                LFatal(" Attempted to get access to a deallocated object");
            }
            tPtr->object::ownedByTable_ = true;

            return *tPtr;
        }

        /*!\brief Release ownership of this object from register table.*/
        inline void release() noexcept { ownedByTable_ = false; }

        void clear() noexcept;

        // Write object

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
        virtual void writeOutput(fileOsstream &fos, const integer fzid) const;

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeMinMaxOutput(fileOsstream &fos) const;

        /*!
         * \brief Write object minimum and maximum value to output file by master.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeMinMaxOutputByMaster(fileOsstream &fos) const;

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeMinMaxOutput(fileOsstream &fos, const integer fzid) const;

        /*!
         * \brief Write object residuals to monitor.
         * Note: it is a virtual function, and should be rewritten in geometry field class.
         */
        virtual void writeResiduals(fileOsstream &fos, integer intervalStep) const;

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

        virtual void readRelay(const hdf5I &fos, const bool readLast, const bool readFromGroup);

        virtual void calcTimeSumPtr(const real &dt) const;

        virtual void interpolateRelay(const hdf5I &fos, const bool readLast,
                                      const bool readFromGroup);
    };
} // namespace OpenHurricane
