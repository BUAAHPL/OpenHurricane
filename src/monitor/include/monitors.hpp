/*!
 * \file monitors.hpp
 * \brief Headers of class of monitors.
 *        The subroutines and functions are in the <i>monitors.cpp</i> file.
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

#include "monitor.hpp"

namespace OpenHurricane {
    /*!\brief The base class of monitors.*/
    class monitors {
    private:
        const iteration &iter_;

        const runtimeMesh &mesh_;

        sharedPtrList<monitor> monitorList_;

        void removeDuplicate(stringList &nameList) const;

    public:

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        monitors(const iteration &iter, const runtimeMesh &mesh);

        /**
         * \brief Destructor.
         */
        virtual ~monitors() noexcept {}

        hur_nodiscard inline const iteration &iter() const noexcept { return iter_; }

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept { return mesh_; }

        void setMonitorList(const controller &monitorCont);

        void monitoring() const;
        void subMonitoring() const;
    };
} // namespace OpenHurricane
