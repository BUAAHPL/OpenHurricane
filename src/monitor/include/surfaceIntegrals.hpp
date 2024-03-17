/*!
 * \file surfaceIntegrals.hpp
 * \brief Headers of base class of computing surface integrals.
 *        The subroutines and functions are in the <i>surfaceIntegrals.cpp</i> file.
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

#include "cellArrays.hpp"
#include "objectFactory.hpp"
#include "runtimeMesh.hpp"

namespace OpenHurricane {
    /*!\brief The base class of computing surface integrals.*/
    class surfaceIntegrals {
    private:
        const iteration &iter_;

        const runtimeMesh &mesh_;

    protected:
        /** \brief The id list of face zones for computing surface integrals. */
        integerList zoneIdList_;

    public:
        declareClassName(surfaceIntegrals);

        declareObjFty(surfaceIntegrals, controller,
                      (const iteration &_iter, const runtimeMesh &_mesh, const controller &cont,
                       const integerList &zoneIdList),
                      (_iter, _mesh, cont, zoneIdList));

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        surfaceIntegrals(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                         const integerList &zoneIdList);

        hur_nodiscard static uniquePtr<surfaceIntegrals> creator(const iteration &iter,
                                                                 const runtimeMesh &mesh,
                                                                 const controller &cont,
                                                                 const integerList &zoneIdList);
        /**
         * \brief Destructor.
         */
        virtual ~surfaceIntegrals() noexcept {}

        hur_nodiscard inline const iteration &iter() const noexcept { return iter_; }

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept { return mesh_; }

        hur_nodiscard inline virtual real surIntegral(const cellRealArray &phi) const {
            return real();
        }
        hur_nodiscard inline virtual vector surIntegral(const cellVectorArray &phi) const {
            return vector();
        }

        hur_nodiscard inline virtual bool dependOnVariables() const { return true; }

        hur_nodiscard inline virtual std::string printInfo() const {
            return std::string("undefined");
        }
    };
} // namespace OpenHurricane
