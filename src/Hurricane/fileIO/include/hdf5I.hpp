/*!
 * \file hdf5I.hpp
 * \brief Headers of input hdf5 file.
 *        The subroutines and functions are in the <i>hdf5I.cpp</i> file.
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
#include "hdf5IO.hpp"
#include <iostream>

namespace OpenHurricane {

    /*!
     * \brief The base class of input hdf5 file.
     */
    class hdf5I : public hdf5IO {
    public:
        enum openOptions : short { ONLY_MASTER, OPEN_ALL };

    private:
        openOptions openOption_;

    public:
        // Constructors

        /**\brief NUll constructor*/
        inline hdf5I();

        /**\brief Construct from file name.*/
        inline hdf5I(const fileName &fN);

        /**\brief Destructor.*/
        inline ~hdf5I() noexcept;

        inline void open(const fileName &fN);

        inline void open();
        void open(const unsigned int flg);

        inline bool onlyMaster() const noexcept;

        void readTypeFromAttribute(const H5::DataSet &dataset, int &dataTypr) const;

        template <class Type> void readSingleCmpt(List<Type> &l, const string &dataName) const;

        template <class Type>
        void readSingleCmpt(List<Type> &l, const string &groupName, const string &dataName) const;

        template <class Type> void readMultipleCmpt(List<Type> &l, const string &dataName) const;

        template <class Type>
        void readMultipleCmpt(List<Type> &l, const string &groupName, const string &dataName) const;

        void readString(std::string &str, const string &dataName) const;
        void readString(std::string &str, const string &groupName, const string &dataName) const;

        template <class Type> void read(List<Type> &l, const string &dataName) const;
        template <class Type>
        void read(List<Type> &l, const string &groupName, const string &dataName) const;

        template <template <typename T> class Form, typename Type>
        void readArrayArray(Form<Type> &l, const string &dataName) const;
        template <template <typename T> class Form, typename Type>
        void readArrayArray(Form<Type> &l, const string &groupName, const string &dataName) const;

        hur_nodiscard bool exist(const string &name) const;
        hur_nodiscard bool exist(const string &groupName, const string &name) const;
    };

} // namespace OpenHurricane

#include "hdf5I.inl"