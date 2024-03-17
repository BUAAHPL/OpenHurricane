/*!
 * \file faceMonitors.hpp
 * \brief Headers of class of monitoring face.
 *        The subroutines and functions are in the <i>faceMonitors.cpp</i> file.
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
#include "surfaceIntegrals.hpp"

namespace OpenHurricane {
    /*!\brief The class of monitoring face.*/
    class faceMonitors : public monitor {
    private:
        /** \brief The name list of face zones monitored. */
        stringList zoneNameList_;

        /** \brief The id list of face zones monitored. */
        integerList zoneIdList_;

        sharedPtrList<surfaceIntegrals> surIntPtr_;

        bool perFaceZone_;

        bool dependOnVar_;

        string varName_;

        /**
         * \brief Read from controller.
         */
        void readFromCont(const controller &cont);

        void getZoneId();

        mutable fileOsstream fosResidual_;

        /**
         * \brief Return the mass flow rate (kg/s) of monitored faces.
         * \return The mass flow rate (kg/s)
         * \retval A real value
         */
        real getFaceMassFlowRate() const;

        void writeAndPrintReal(const cellRealArray &phi) const;
        void writeAndPrintVector(const cellVectorArray &phi) const;

        void writeAndPrint(const integer step) const;

    public:
        declareClassName(faceMonitors);

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        faceMonitors(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                     const string &name);

        /**
         * \brief Destructor.
         */
        virtual ~faceMonitors() noexcept { fosResidual_.close(); }

        /**
         * \brief Monitoring.
         */
        virtual void monitoring() const;
        virtual void subMonitoring() const;
    };
} // namespace OpenHurricane
