/*!
 * \file patching.hpp
 * \brief Headers of class of patching.
 *        The subroutines and functions are in the <i>patching.cpp</i> file.
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

#include "flowModel.hpp"
#include "spatialScheme.hpp"
#include "timeMarching.hpp"
#include "turbulenceModel.hpp"
#include "viscousFlux.hpp"

namespace OpenHurricane {
    class patching {
    private:
        /*!\brief Hold reference to runtime mesh.*/
        const runtimeMesh &mesh_;

        const controller &cont_;

        /*!\brief The name of the region.*/
        std::string regionName_;

        /*!\brief The type of patch method.*/
        std::string typeName_;

        /*!\brief The specific variable name.*/
        stdStringList varName_;

        /*!\brief The step to start the patching action.*/
        integer startStep_;

        /*!\brief The steps to stay the patching action.*/
        integer stayStep_;

    public:
        /*!\brief Construct from mesh and controller.*/
        patching(const runtimeMesh &mesh, const controller &cont);

        /*!\brief Destructor.*/
        inline ~patching() noexcept {}

        /*!\brief Hold reference to runtime mesh.*/
        hur_nodiscard inline const runtimeMesh &mesh() const;

        hur_nodiscard inline const controller &cont() const;

        /*!\brief The name of the region.*/
        hur_nodiscard inline const std::string &regionName() const noexcept;

        /*!\brief The type of patch method.*/
        hur_nodiscard inline const std::string &typeName() const noexcept;

        /*!\brief The specific variable name.*/
        hur_nodiscard inline const stdStringList &varName() const;

        /*!\brief The step to start the patching action.*/
        hur_nodiscard inline integer startStep() const noexcept;

        /*!\brief The steps to stay the patching action.*/
        hur_nodiscard inline integer stayStep() const noexcept;

        void patchRegion(const integer istep) const;

        void patchRegion() const;

    private:
        void uniformPatch(const integer istep) const;

        void uniformPatch() const;

        inline void distributePatch(const integer istep) const {
            if (istep >= startStep_ && istep <= (startStep_ + stayStep_)) {
                distributePatch();
            }
        }

        void distributePatch() const;
    };
} // namespace OpenHurricane

#include "patching.inl"