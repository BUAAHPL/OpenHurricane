/*!
 * \file writeFieldVars.hpp
 * \brief Headers of class of writeFieldVars.
 *        The subroutines and functions are in the <i>writeFieldVars.cpp</i> file.
 * \author Chen Zhenyi
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

#include "writeFieldVar.hpp"

namespace OpenHurricane {
    /*!\brief The base class of writeFieldVars.*/
    class writeFieldVars {
    private:
        /**\brief The flow model.*/
        const flowModel &flows_;
        const controller &cont_;
        const iteration &iter_;

        /**\brief list of primitive field variables for time average and pulse value calculation .*/
        stringList timeSumVarList_;

        /**\brief The time of stopping calculating time average field variables.*/
        real endTimeAveCal_;

        mutable PtrList<writeFieldVar> writeFieldPtrList_;

        /**\brief Map for output variables of flow field.*/
        mutable std::map<std::string, object *> outFieldVarMap_;

    public:

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        writeFieldVars(const flowModel &flows_, const controller &cont_, const iteration &iter_);

        /**
         * \brief Destructor.
         */
        virtual ~writeFieldVars() noexcept { clearoutVarMap(); }

    public:

        inline const real &endTimeAveCal() const noexcept;
        inline const runtimeMesh &mesh() const;

        inline bool foundInMap(const string &name) const;

    protected:
        void updating() const;

        void clearoutVarMap() const;

    public:
        void setWriteFieldVarList(const controller &iterCont);

        void write() const;

        void setOutputField(const string &varName, const realArray &value) const;
        void setOutputField(const string &varName, realArray &&value) const;
        void setOutputField(const string &varName, const vectorArray &value) const;
        void setOutputField(const string &varName, vectorArray &&value) const;

    private:
        void removeDuplicate(stringList &writeList) const;
    };
} // namespace OpenHurricane

#include "writeFieldVars.inl"