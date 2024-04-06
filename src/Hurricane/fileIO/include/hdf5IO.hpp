/*!
 * \file hdf5IO.hpp
 * \brief Headers of input and output hdf5 file.
 *        The subroutines and functions are in the <i>hdf5IO.cpp</i> file.
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
#include <iostream>

#include "H5Cpp.h"

#include "List.hpp"
#include "dataStructure.hpp"

namespace OpenHurricane {

    /*!
     * \brief The base class of input and output hdf5 file.
     */
    class hdf5IO {
    public:
        enum accessState { OPENED, CLOSED };

    protected:
        /**\brief File name.*/
        fileName filename_;

        /**\brief H5 file pointer.*/
        uniquePtr<H5::H5File> filePtr_;

        /**\brief H5 file flag.*/
        unsigned int flag_;

        /*! \brief The access state of the file.*/
        accessState openClosed_;

        /*! \brief Set file opened.*/
        inline void setOpened() noexcept;

        /*! \brief Set file opened.*/
        inline void setClosed() noexcept;

    public:

        /**\brief NUll constructor*/
        inline hdf5IO();

        /**\brief Construct from file name.*/
        inline hdf5IO(const fileName &fN);

        /**\brief Construct from file name.*/
        inline hdf5IO(const fileName &fN, const unsigned int flg);

        /**\brief Destructor.*/
        inline ~hdf5IO() noexcept;

        inline void open(const fileName &fN);
        inline void open();

        inline void open(const unsigned int flg);

        inline void close();

        /*!\brief Return true if the file opened.*/
        hur_nodiscard inline bool opened() const noexcept;

        /*!\brief Return true if the file closed.*/
        hur_nodiscard inline bool closed() const noexcept;

        inline void createGroup(const string &groupName);

        // Write

        template <class Type>
        void write(const Type *data, const integer size, const string &dataName);

        template <class Type>
        void write(const Type *data, const integer size, const string &groupName,
                   const string &dataName);

        template <class Type>
        void writeSingleCmpt(const List<Type> &l, const Type factor, const integer size,
                             const string &dataName);

        template <class Type>
        void writeSingleCmpt(const List<Type> &l, const Type factor, const integer size,
                             const string &groupName, const string &dataName);

        template <class Type>
        void writeSingleCmpt(const List<Type> &l, const integer size, const string &dataName);

        template <class Type>
        void writeSingleCmpt(const List<Type> &l, const integer size, const string &groupName,
                             const string &dataName);

        template <class Type> void writeSingleCmpt(const List<Type> &l, const string &dataName);

        template <class Type>
        void writeSingleCmpt(const List<Type> &l, const string &groupName, const string &dataName);

        template <class Type>
        void writeMultipleCmpt(const List<Type> &l,
                               const typename feature<Type>::elementType factor, const integer size,
                               const string &dataName);

        template <class Type>
        void writeMultipleCmpt(const List<Type> &l,
                               const typename feature<Type>::elementType factor, const integer size,
                               const string &groupName, const string &dataName);

        template <class Type>
        void writeMultipleCmpt(const List<Type> &l, const integer size, const string &dataName);

        template <class Type>
        void writeMultipleCmpt(const List<Type> &l, const integer size, const string &groupName,
                               const string &dataName);

        template <class Type> void writeMultipleCmpt(const List<Type> &l, const string &dataName);

        template <class Type>
        void writeMultipleCmpt(const List<Type> &l, const string &groupName,
                               const string &dataName);

        void writeString(const std::string &str, const string &dataName);

        void writeString(const std::string &str, const string &groupName, const string &dataName);

        template <class Type>
        void write(const List<Type> &l, const typename feature<Type>::elementType factor,
                   const integer size, const string &dataName);

        template <class Type>
        void write(const List<Type> &l, const typename feature<Type>::elementType factor,
                   const integer size, const string &groupName, const string &dataName);

        template <class Type>
        void write(const List<Type> &l, const integer size, const string &dataName);

        template <class Type>
        void write(const List<Type> &l, const integer size, const string &groupName,
                   const string &dataName);

        template <class Type> void write(const List<Type> &l, const string &dataName);

        template <class Type>
        void write(const List<Type> &l, const string &groupName, const string &dataName);

        template <template <typename T> class Form, typename Type>
        void writeArrayArray(const Form<Type> &l, const string &dataName);

        template <template <typename T> class Form, typename Type>
        void writeArrayArray(const Form<Type> &l, const string &groupName, const string &dataName);

        void writeIntegerAttributeToFile(const integer l, const string &attrName);

        void writeRealAttributeToFile(const real l, const string &attrName);
        void writeVectorAttributeToFile(const vector &l, const string &attrName);

        void writeStringAttributeToFile(const std::string &l, const string &attrName);

        void writeIntegerAttributeToDataset(const integer l, const string &attrName,
                                            const string &dataName);
        void writeIntegerAttributeToDataset(const integer l, const string &attrName,
                                            const string &groupName, const string &dataName);

        void writeRealAttributeToDataset(const real l, const string &attrName,
                                         const string &dataName);
        void writeRealAttributeToDataset(const real l, const string &attrName,
                                         const string &groupName, const string &dataName);

        void writeVectorAttributeToDataset(const vector &l, const string &attrName,
                                           const string &dataName);
        void writeVectorAttributeToDataset(const vector &l, const string &attrName,
                                           const string &groupName, const string &dataName);

        void writeStringAttributeToDataset(const std::string &l, const string &attrName,
                                           const string &dataName);
        void writeStringAttributeToDataset(const std::string &l, const string &attrName,
                                           const string &groupName, const string &dataName);

        void writeIntegerAttributeToGroup(const integer l, const string &attrName,
                                          const string &groupName);

        void writeRealAttributeToGroup(const real l, const string &attrName,
                                       const string &groupName);
        void writeVectorAttributeToGroup(const vector &l, const string &attrName,
                                         const string &groupName);

        void writeStringAttributeToGroup(const std::string &l, const string &attrName,
                                         const string &groupName);

        // Read

        void readIntegerAttributeFromFile(integer &l, const string &attrName) const;

        void readRealAttributeFromFile(real &l, const string &attrName) const;
        void readVectorAttributeFromFile(vector &l, const string &attrName) const;

        void readStringAttributeFromFile(std::string &l, const string &attrName) const;

        void readIntegerAttributeFromDataset(integer &l, const string &attrName,
                                             const string &dataName) const;
        void readIntegerAttributeFromDataset(integer &l, const string &attrName,
                                             const string &groupName, const string &dataName) const;

        void readRealAttributeFromDataset(real &l, const string &attrName,
                                          const string &dataName) const;
        void readRealAttributeFromDataset(real &l, const string &attrName, const string &groupName,
                                          const string &dataName) const;

        void readVectorAttributeFromDataset(vector &l, const string &attrName,
                                            const string &dataName) const;
        void readVectorAttributeFromDataset(vector &l, const string &attrName,
                                            const string &groupName, const string &dataName) const;

        void readStringAttributeFromDataset(std::string &l, const string &attrName,
                                            const string &dataName) const;
        void readStringAttributeFromDataset(std::string &l, const string &attrName,
                                            const string &groupName, const string &dataName) const;

        void readIntegerAttributeFromGroup(integer &l, const string &attrName,
                                           const string &groupName) const;

        void readRealAttributeFromGroup(real &l, const string &attrName,
                                        const string &groupName) const;
        void readVectorAttributeFromGroup(vector &l, const string &attrName,
                                          const string &groupName) const;

        void readStringAttributeFromGroup(std::string &l, const string &attrName,
                                          const string &groupName) const;

        void readTypeFromAttribute(const H5::DataSet &dataset, int &dataTypr) const;

        template <class Type> void readSingleCmpt(List<Type> &l, const string &dataName) const;
        template <class Type>
        void readSingleCmpt(List<Type> &l, const string &groupName, const string &dataName) const;

        template <class Type> void readMultipleCmpt(List<Type> &l, const string &dataName) const;

        template <class Type>
        void readMultipleCmpt(List<Type> &l, const string &groupName, const string &dataName) const;

        template <template <typename T> class Form, typename Type>
        void readArrayArray(Form<Type> &l, const string &dataName) const;

        template <template <typename T> class Form, typename Type>
        void readArrayArray(Form<Type> &l, const string &groupName, const string &dataName) const;

        void readString(std::string &str, const string &dataName) const;

        void readString(std::string &str, const string &groupName, const string &dataName) const;

        template <class Type> void read(List<Type> &l, const string &dataName) const;

        template <class Type>
        void read(List<Type> &l, const string &groupName, const string &dataName) const;

        bool exist(const string &name) const;
        bool exist(const string &groupName, const string &name) const;

        bool isHDF5File() const;
    };

} // namespace OpenHurricane

template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const typename feature<real>::elementType factor,
                              const integer size, const string &dataName);
template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const typename feature<real>::elementType factor,
                              const integer size, const string &groupName, const string &dataName);

template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const integer size, const string &dataName);
template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const integer size, const string &groupName,
                              const string &dataName);

template <> void OpenHurricane::hdf5IO::write(const List<real> &l, const string &dataName);

template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const string &groupName, const string &dataName);

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l,
                              const typename feature<integer>::elementType factor,
                              const integer size, const string &dataName);

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l,
                              const typename feature<integer>::elementType factor,
                              const integer size, const string &groupName, const string &dataName);

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l, const integer size, const string &dataName);

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l, const integer size, const string &groupName,
                              const string &dataName);

template <> void OpenHurricane::hdf5IO::write(const List<integer> &l, const string &dataName);

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l, const string &groupName,
                              const string &dataName);

template <> void OpenHurricane::hdf5IO::write(const List<std::string> &l, const string &dataName);

template <>
void OpenHurricane::hdf5IO::write(const List<std::string> &l, const string &groupName,
                              const string &dataName);

template <> void OpenHurricane::hdf5IO::write(const List<string> &l, const string &dataName);

template <>
void OpenHurricane::hdf5IO::write(const List<string> &l, const string &groupName,
                              const string &dataName);

template <> void OpenHurricane::hdf5IO::read(List<real> &l, const string &dataName) const;

template <>
void OpenHurricane::hdf5IO::read(List<real> &l, const string &groupName, const string &dataName) const;

template <> void OpenHurricane::hdf5IO::read(List<integer> &l, const string &dataName) const;

template <>
void OpenHurricane::hdf5IO::read(List<integer> &l, const string &groupName,
                             const string &dataName) const;

template <> void OpenHurricane::hdf5IO::read(List<std::string> &l, const string &dataName) const;

template <>
void OpenHurricane::hdf5IO::read(List<std::string> &l, const string &groupName,
                             const string &dataName) const;

template <> void OpenHurricane::hdf5IO::read(List<string> &l, const string &dataName) const;

template <>
void OpenHurricane::hdf5IO::read(List<string> &l, const string &groupName,
                             const string &dataName) const;

#include "hdf5IO.inl"