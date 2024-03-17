/*!
 * \file drag.hpp
 * \brief Headers of class of drag force.
 *        The subroutines and functions are in the <i>drag.cpp</i> file.
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
#include "force.hpp"

namespace OpenHurricane {
    /*!\brief The class of drag force.*/
    class drag : public force {
    protected:
        void setReportName();

    public:
        declareClassName(drag);

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        drag(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
             const integerList &zd);

        /**
         * \brief Destructor.
         */
        virtual ~drag() noexcept {}

        inline virtual void computing() const { force::setReportArray(); }
    };
} // namespace OpenHurricane
