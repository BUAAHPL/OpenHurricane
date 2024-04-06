/*!
 * \file pointMonitors.hpp
 * \brief Headers of class of monitoring points.
 *        The subroutines and functions are in the <i>pointMonitors.cpp</i> file.
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
#include "monitor.hpp"

namespace OpenHurricane {
    /*!\brief The class of monitoring points.*/
    class pointMonitors : public monitor {
    private:
        /** \brief The index of cell that contains the point. */
        integer cellIndex_;

        /** \brief The position for monitoring. */
        point position_;

        /** \brief The index of the process that holds the point. */
        integer processId_;

        stringList monitorVarName_;
        string monitorVarCmptName_;

        integer componentSize_;

        integerList componentId_;

        mutable fileOsstream fosPoint_;

        void sendToMaster(realArray &data) const;

        void writeOutScreenInMaster(const realArray &data) const;
        void writeOutFileInMaster(const realArray &data) const;

        void setMonitorVarCmpt();

        void getMonitorVarCmpt(realArray &data) const;

        real pointMa() const;

    public:
        declareClassName(pointMonitors);

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        pointMonitors(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                      const string &name);

        /**
         * \brief Destructor.
         */
        virtual ~pointMonitors() noexcept {}

        /**
         * \brief Monitoring.
         */
        virtual void monitoring() const;
        inline virtual void subMonitoring() const {}
    };
} // namespace OpenHurricane
