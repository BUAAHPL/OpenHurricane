/*!
 * \file cellBasedMethod.hpp
 * \brief Headers of class of cell-based method for computing time step.
 *        The subroutines and functions are in the <i>cellBasedMethod.cpp</i> file.
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
#include "pseudoTime.hpp"

namespace OpenHurricane {
    /*!\brief The class of cell-based method for computing time step.*/
    class cellBasedMethod : public pseudoTime {
    private:
        uniquePtr<realArray> deltaSxPtr_;
        uniquePtr<realArray> deltaSyPtr_;
        uniquePtr<realArray> deltaSzPtr_;

        realArray &deltaSx();
        realArray &deltaSy();
        realArray &deltaSz();

        bool useLowMachPrecon_;

    public:
        declareClassNames;

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        cellBasedMethod(const controller &cont, const runtimeMesh &mesh, const iteration &iter,
                        const flowModel &flows, realArray &dt, const cellRealArray &shockFactor,
                        const faceVector2DArray &rai, const faceVector2DArray &rav,
                        const integerArray &temperatureFlag, const integerArray &pressureFlag,
                        const cellIntegerArray &CFLFlag);

        /**
         * \brief Destructor.
         */
        virtual ~cellBasedMethod() noexcept;

        /**
         * \brief To compute the time-step for time marching.
         * \return The minimum timestep.
         */
        virtual void computingTimeStep();

        /**
         * \brief To compute the time-step for time marching.
         * \return The minimum timestep.
         */
        virtual void computingTimeStep(realArray &dt, const real cfl0);
    };
} // namespace OpenHurricane
