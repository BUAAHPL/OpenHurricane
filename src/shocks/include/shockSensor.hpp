/*!
 * \file shockSensor.hpp
 * \brief Headers of class of the shock sensor.
 *        The subroutines and functions are in the <i>shockSensor.cpp</i> file.
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
#include "fVArraysInclude.hpp"
#include "flowModel.hpp"

namespace OpenHurricane {
    /*!\brief The base class of shock sensor.*/
    class shockSensor {
    protected:
        /*!\brief Hold reference to the flows.*/
        const flowModel &flows_;

        hur_nodiscard inline const cellVectorArray &v() const noexcept { return flows_.v(); }

        realArray delta_;

        cellRealArray sensor_;

        /** \brief The geometrical weights. */
        vector2DArray thetaIJ_;
        void getGeometricalWeights();

        void checkDiffer(const realArray &q, const real differThr);

        real differThr_;
        real differTThr_;

        bool checkTemperature_;
        real THighLimit_;
        void checkTemperatureDiffer();

    public:
        /*!\brief Disallow null constructor.*/
        shockSensor() = delete;

        /*!\brief Disallow copy constructor.*/
        shockSensor(const shockSensor &) = delete;
        shockSensor &operator=(const shockSensor &) = delete;

        shockSensor(const flowModel &flows, const controller &cont);

        inline ~shockSensor() noexcept {}

        void calcShockSensor();

        hur_nodiscard inline const cellRealArray &sensor() const noexcept { return sensor_; }
        hur_nodiscard inline cellRealArray &sensor() noexcept { return sensor_; }

        void checkTemperatureDiffer(realArray &factor) const;
    };
} // namespace OpenHurricane
