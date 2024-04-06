/*!
 * \file calculateFaceZoneVar.hpp
 * \brief Headers of class of calculating Face-Zone-based Variables.
 *        The subroutines and functions are in the <i>calculateFaceZoneVar.cpp</i> file.
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
#include "fieldParameterFuncMap.hpp"
#include "flowModel.hpp"

namespace OpenHurricane {
    class calculateFaceZoneVar {
    private:
        void addToFuncMap(faceZoneFieldParameterFuncMap &funcMap) const;

    public:
        /**
         * \brief Calculating Mach number on the given face zone.
         * \param[in] flow - Flow
         * \param[in] fzi - The index of the face zone
         * \return The Mach number
         * \retval A real array
         */
        static hur_nodiscard realArray calcFaceZoneMachNumber(const flowModel &flow,
                                                              const integer fzi);

        /**
         * \brief Calculating stagnation pressure on the given face zone.
         * \param[in] flow - Flow
         * \param[in] fzi - The index of the face zone
         * \return The stagnation pressure
         * \retval A real array
         */
        static hur_nodiscard realArray calcFaceZoneTotalPressure(const flowModel &flow,
                                                                 const integer fzi);

        /**
         * \brief Calculating stagnation pressure on the given face zone.
         * \param[in] flow - Flow
         * \param[in] fzi - The index of the face zone
         * \return The stagnation pressure
         * \retval A real array
         */
        static hur_nodiscard realArray calcFaceZoneTotalPressureCCP(const flowModel &flow,
                                                                    const integer fzi);

        /**
         * \brief Calculating stagnation temperature on the given face zone.
         * \param[in] flow - Flow
         * \param[in] fzi - The index of the face zone
         * \return The stagnation temperature
         * \retval A real array
         */
        static hur_nodiscard realArray calcFaceZoneTotalTemperature(const flowModel &flow,
                                                                    const integer fzi);

        /**
         * \brief Calculating stagnation temperature on the given face zone.
         * \param[in] flow - Flow
         * \param[in] fzi - The index of the face zone
         * \return The stagnation temperature
         * \retval A real array
         */
        static hur_nodiscard realArray calcFaceZoneTotalTemperatureCCp(const flowModel &flow,
                                                                       const integer fzi);

        static hur_nodiscard realArray wallFaceZoneHeatFlux(const flowModel &flow,
                                                            const integer zoneId);
        static hur_nodiscard realArray wallFaceZoneFrictionCoefficient(const flowModel &flow,
                                                                       const integer zoneId);
        static hur_nodiscard realArray wallFaceZoneAbsPressCoefficient(const flowModel &flow,
                                                                       const integer zoneId);
        static hur_nodiscard realArray wallFaceZoneRelPressCoefficient(const flowModel &flow,
                                                                       const integer zoneId);
        static hur_nodiscard realArray wallFaceZoneUPlus(const flowModel &flow,
                                                         const integer zoneId);
        static hur_nodiscard realArray wallFaceZoneYPlus(const flowModel &flow,
                                                         const integer zoneId);
        static hur_nodiscard realArray wallFaceZoneHeatCoefficient(const flowModel &flow,
                                                                   const integer zoneId);

        static hur_nodiscard realArray calcFaceZoneVorticity(const flowModel &flow,
                                                             const integer zoneId);

        static hur_nodiscard realArray calcViscousRatio(const flowModel &flow,
                                                        const integer zoneId);

        /**
         * \brief Calculate x-velocity of the first layer cell near the face zone specified by zoneId.
         */
        static hur_nodiscard realArray calcU0(const flowModel &flow, const integer zoneId);

        /**
         * \brief Calculate y-velocity of the first layer cell near the face zone specified by zoneId.
         */
        static hur_nodiscard realArray calcV0(const flowModel &flow, const integer zoneId);

        /**
         * \brief Calculate z-velocity of the first layer cell near the face zone specified by zoneId.
         */
        static hur_nodiscard realArray calcW0(const flowModel &flow, const integer zoneId);

        // Constructors

        inline calculateFaceZoneVar();

        calculateFaceZoneVar(faceZoneFieldParameterFuncMap &funcMap);

        ~calculateFaceZoneVar() noexcept;

        static string setKeyName(const string &key);

        inline static string basicName();
    };
} // namespace OpenHurricane

#include "calculateFaceZoneVar.inl"