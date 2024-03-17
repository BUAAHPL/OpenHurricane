/*!
 * \file profiles.hpp
 * \brief Headers of class of profiles.
 *        The subroutines and functions are in the <i>profiles.cpp</i> file.
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

#include "flowModel.hpp"
#include "runtimeMesh.hpp"

namespace OpenHurricane {
    /*!\brief The base class of profiles.*/
    class profiles {
    protected:
        fileName files_;

    public:
        declareClassNames;
        declareObjFty(profiles, controller, (const fileName &file), (file));

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        profiles(const fileName &file);

        static uniquePtr<profiles> creator(const fileName &file, const string &profileType);

        /**
         * \brief Destructor.
         */
        virtual ~profiles() noexcept {}

        void read(controller &cont) const;

        void write(const flowModel &flows, const controller &cont) const;

    protected:
        virtual void readProfiles(controller &cont) const = 0;

        virtual void writeProfiles(const controller &cont) const = 0;

        void getProfiles(const flowModel &flows, const controller &cont,
                         controller &profCont) const;

        std::string getFieldVar(const flowModel &flows, const integer fzi, const string &fieldName,
                                controller &profCont) const;

        template <class Type> void allGatherVList(const List<Type> &send, List<Type> &allL) const;
    };
} // namespace OpenHurricane

#include "profiles.inl"