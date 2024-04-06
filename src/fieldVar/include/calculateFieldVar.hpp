/*!
 * \file calculateFieldVar.hpp
 * \brief Headers of base class of calculateFieldVar.
 *        The subroutines and functions are in the <i>calculateFieldVar.cpp</i> file.
 * \author Chen Zhenyi
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
#include "combustionModel.hpp"
#include "realArray.hpp"
#include "vectorArray.hpp"

namespace OpenHurricane {
    /*!\brief The class of functions of calculating field variables.*/
    class calculateFieldVar {
    public:
        static inline realArray calcMachNumber(const vectorArray &v, const realArray &gama,
                                               const realArray &p, const realArray &rho);

        static inline realArray calcTotalPressure(const realArray &p, const realArray &gama,
                                                  const realArray &Ma, const flowModel &flow);

        /**
         * \brief Calculating stagnation pressure on the given face zone.
         * \param[in] flow - Flow
         * \param[in] fzi - The index of the face zone
         * \return The stagnation pressure
         * \retval A real array
         */
        static hur_nodiscard realArray calcTotalPressure(const flowModel &flow);

        static inline realArray calcTotalPressure(const realArray &p, const realArray &gama,
                                                  const vectorArray &v, const realArray &rho,
                                                  const flowModel &flow);

        static inline realArray calcTotalTemperature(const realArray &T, const realArray &gama,
                                                     const realArray &Ma, const flowModel &flow);

        static inline realArray calcTotalTemperature(const realArray &T, const realArray &gama,
                                                     const vectorArray &v, const realArray &p,
                                                     const realArray &rho, const flowModel &flow);

        /**\brief calculate totalTemperature if species is not single*/
        static hur_nodiscard inline realArray
        totalTemperature(const flowModel &flow, const realArray &p, const realArray &T,
                         const vectorArray &v, const PtrList<cellRealArray> &yi);

        static hur_nodiscard realArray calcTotalEnthalpy(const flowModel &flow);

        static hur_nodiscard realArray cellVolume(const flowModel &flow);

        static hur_nodiscard inline realArray calcViscousRatio(const flowModel &flow);

        static inline realArray calcHeatReleaseRate(const flowModel &flow,
                                                    const combustionModel *chemtryPtr);
        static inline realArray calcGFO(const flowModel &flow, const combustionModel *chemtryPtr);

        static inline void calcMoleSpecies(const flowModel &flow,
                                           std::map<std::string, object *> &outFieldVarMap,
                                           const string &type);

        static inline realArray calcVorticity(const flowModel &flow);

        static inline realArray calcDamkohler(const flowModel &flow,
                                              const combustionModel *chemtryPtr);

        static hur_nodiscard realArray calcQCriterion(const flowModel &flow);
        static hur_nodiscard realArray calcOmegaCriterion(const flowModel &flow);

        static hur_nodiscard realArray calcDeltaCriterion(const flowModel &flow);
    };
} // namespace OpenHurricane

#include "calculateFieldVar.inl"