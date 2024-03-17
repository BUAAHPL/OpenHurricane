/*!
 * \file hdf5O.hpp
 * \brief Headers of output hdf5 file.
 *        The subroutines and functions are in the <i>hdf5O.cpp</i> file.
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
#include "hdf5IO.hpp"
#include <iostream>

namespace OpenHurricane {

    /*!
     * \brief The base class of output hdf5 file.
     */
    class hdf5O : public hdf5IO {
    public:
        enum openOptions : short { ONLY_MASTER, CREATE_ALL };

    private:
        openOptions openOption_;

    public:
        // Constructors

        /**\brief NUll constructor*/
        inline hdf5O();

        /**\brief Construct from file name.*/
        inline hdf5O(const fileName &fN);

        /**\brief Destructor.*/
        inline ~hdf5O() noexcept;

        inline void open(const fileName &fN);

        inline void open();
        void open(const unsigned int flg);

        inline bool onlyMaster() const noexcept;

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

        /*template<class Type>
        void writeMultipleCmpt(const List<Type>& l, const string& dataName);*/

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
    };

} // namespace OpenHurricane

template <> void OpenHurricane::hdf5O::write(const List<std::string> &l, const string &dataName);
template <>
void OpenHurricane::hdf5O::write(const List<std::string> &l, const string &groupName,
                             const string &dataName);

template <> void OpenHurricane::hdf5O::write(const List<string> &l, const string &dataName);
template <>
void OpenHurricane::hdf5O::write(const List<string> &l, const string &groupName,
                             const string &dataName);

#include "hdf5O.inl"