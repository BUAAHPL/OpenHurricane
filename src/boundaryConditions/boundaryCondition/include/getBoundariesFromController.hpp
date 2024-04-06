/*!
 * \file getBoundariesFromController.hpp
 * \brief Headers of parsing boundary condition from controller.
 *        The subroutines and functions are in the <i>getBoundariesFromController.cpp</i> file.
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
#include "parsingDirection.hpp"

namespace OpenHurricane {
    class mixture;
    /*!\brief The class for parsing boundary condition from controller.
     */
    class getBoundariesFromController {
    private:
    public:
        // Constructors

        getBoundariesFromController() {}

        ~getBoundariesFromController() noexcept {}

        static void setBoundariesController(controller &cont, const mixture &mixtures);

        /*static void setFreeStream(const controller& cont, const mixture& mixtures, basicFreeStream& fs);
        static void autoSetFreeStream(const controller& bbCont, const mixture& mixtures, basicFreeStream& fs);*/

        static realArray getSpeciesMassFractions(controller &bcCont, const mixture &mixtures,
                                                 const bool addToCont = false);

        static void getPressureFarField(const mixture &mixtures, const controller &cont,
                                        controller &bcCont, const faceZone &fz);
        static void getSupersonicInlet(const mixture &mixtures, const controller &cont,
                                       controller &bcCont, const faceZone &fz);
        static void getVelocityInlet(const mixture &mixtures, const controller &cont,
                                     controller &bcCont, const faceZone &fz);
        static void getSubsonicInlet(const mixture &mixtures, const controller &cont,
                                     controller &bcCont, const faceZone &fz);
        static void getPressureInlet(const mixture &mixtures, const controller &cont,
                                     controller &bcCont, const faceZone &fz);
        static void getMassFlowInlet(const mixture &mixtures, const controller &cont,
                                     controller &bcCont, const faceZone &fz);
        static void getDetonationInlet(const mixture &mixtures, const controller &cont,
                                       controller &bcCont, const faceZone &fz);
        static void getOutflow(const controller &cont, controller &bcCont, const faceZone &fz);
        static void getPressureOutlet(const controller &cont, controller &bcCont,
                                      const faceZone &fz);
        static void getWallCondition(const mixture &mixtures, const controller &cont,
                                     controller &bcCont, const faceZone &fz);

        static void addBcTypeToController(const string &fzName, controller &interCont,
                                          std::string &bcType);

        static void checkInitializationSetting(const controller &cont, const runtimeMesh &mesh);

        static void getSyntheticTurbInlet(const mixture &mixtures, const controller &cont,
                                          controller &bcCont, const faceZone &fz);
    };
} // namespace OpenHurricane