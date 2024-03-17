/*!
 * \file FluentProfiles.hpp
 * \brief Headers of class of Ansys Fluent Profiles.
 *        The subroutines and functions are in the <i>FluentProfiles.cpp</i> file.
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

#include "profiles.hpp"

namespace OpenHurricane {
    /*!\brief The base class of Ansys Fluent Profiles.*/
    class FluentProfiles : public profiles {
    public:
        declareClassNames;

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        FluentProfiles(const fileName &file);

        /**
         * \brief Destructor.
         */
        virtual ~FluentProfiles() noexcept {}

    protected:
        virtual void readProfiles(controller &cont) const;

        virtual void writeProfiles(const controller &cont) const;

    private:
        std::string readFluentProfile(const std::string &str, controller &cont) const;
    };
} // namespace OpenHurricane
