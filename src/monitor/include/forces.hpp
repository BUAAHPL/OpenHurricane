/*!
 * \file forces.hpp
 * \brief Headers of base class of computing forces.
 *        The subroutines and functions are in the <i>forces.cpp</i> file.
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

#include "cellArrays.hpp"
#include "objectFactory.hpp"
#include "runtimeMesh.hpp"

namespace OpenHurricane {
    /*!\brief The base class of computing forces.*/
    class forces {
    public:
        enum class reportType : short { force = 0, coefficient = 1 };

    private:
        const iteration &iter_;

        const runtimeMesh &mesh_;

    protected:
        /** \brief The id list of face zones for computing surface integrals. */
        integerList zoneIdList_;

        reportType repTp_;

        bool reportPressureComponent_;

        bool reportViscousComponent_;

        bool reportCoordinateComponent_;

        mutable realArray reportArray_;

        stringList reportName_;

        void setArraySize();

        /**
         * \brief The free stream dynamic pressure [Pa].
         */
        real pDynFree() const;

    public:
        declareClassName(forces);

        declareObjFty(forces, controller,
                      (const iteration &_iter, const runtimeMesh &_mesh, const controller &cont,
                       const integerList &zoneIdList),
                      (_iter, _mesh, cont, zoneIdList));

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        forces(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
               const integerList &zoneIdList);

        hur_nodiscard static uniquePtr<forces> creator(const iteration &iter,
                                                       const runtimeMesh &mesh,
                                                       const controller &cont,
                                                       const integerList &zoneIdList);
        /**
         * \brief Destructor.
         */
        virtual ~forces() noexcept {}

        hur_nodiscard inline const iteration &iter() const noexcept { return iter_; }

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept { return mesh_; }

        /**
         * \brief To get the force vector.
         * \param[out] Fp - pressure force vector
         * \param[out] Fv - viscous force vector
         */
        void getForces(vector &Fp, vector &Fv) const;

        hur_nodiscard inline const realArray &reportArray() const noexcept { return reportArray_; }

        hur_nodiscard inline const stringList &reportName() const noexcept { return reportName_; }

        virtual void computing() const = 0;

        hur_nodiscard inline integer reportSize() const noexcept { return reportArray_.size(); }
    };
} // namespace OpenHurricane
