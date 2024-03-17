/*!
 * \file pseudoTime.hpp
 * \brief Headers of class of the pseudo time.
 *        The subroutines and functions are in the <i>pseudoTime.cpp</i> file.
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
#include "CFL.hpp"
#include "fVArraysInclude.hpp"
#include "flowModel.hpp"

namespace OpenHurricane {
    /*!\brief The base class of pseudo time.*/
    class pseudoTime {
    protected:
        /*!\brief Hold reference to the runtime mesh.*/
        const runtimeMesh &mesh_;

        const iteration &iter_;

        const flowModel &flow_;

        /*!\brief Time step.*/
        realArray &dt_;

        const cellRealArray &shockFactor_;

        /*!\brief The spectral radius of the convective fluxes*/
        const faceVector2DArray &rai_;

        /*!\brief The spectral radius of the viscous fluxes.*/
        const faceVector2DArray &rav_;

        const integerArray &temperatureFlag_;
        const integerArray &pressureFlag_;
        const cellIntegerArray &CFLFlag_;

        /**\brief Timestep reduced factor.*/
        real dtReduceFactor_;

        /**\brief The pointer of Courant-Friedrichs-Lewy (CFL) numbers.*/
        uniquePtr<CFL> cflPtr_;

        /**\brief COnstant for time step computation*/
        real CForTimeStep_;

        /**\brief Note cell(n).AR is the cell aera aspect ratio. It is used to accelerate the convergence
         *  for stretched meshes.
         */
        bool isStretchAc_;

        real cflRatioMax_;

        real minStretchScale_;

        /** \brief The maximum cell orthogonality for correct the time step. */
        real minCellOrthogonality_;

    public:
        declareClassNames;
        declareObjFty(pseudoTime, controller,
                      (const controller &cont, const runtimeMesh &mesh, const iteration &iter,
                       const flowModel &flows, realArray &dt, const cellRealArray &shockFactor,
                       const faceVector2DArray &rai, const faceVector2DArray &rav,
                       const integerArray &temperatureFlag, const integerArray &pressureFlag,
                       const cellIntegerArray &CFLFlag),
                      (cont, mesh, iter, flows, dt, shockFactor, rai, rav, temperatureFlag,
                       pressureFlag, CFLFlag));

        /*!\brief Disallow null constructor.*/
        pseudoTime() = delete;

        /*!\brief Disallow copy constructor.*/
        pseudoTime(const pseudoTime &) = delete;
        pseudoTime &operator=(const pseudoTime &) = delete;

        pseudoTime(const controller &cont, const runtimeMesh &mesh, const iteration &iter,
                   const flowModel &flows, realArray &dt, const cellRealArray &shockFactor,
                   const faceVector2DArray &rai, const faceVector2DArray &rav,
                   const integerArray &temperatureFlag, const integerArray &pressureFlag,
                   const cellIntegerArray &CFLFlag);

        static uniquePtr<pseudoTime>
        creator(const controller &cont, const runtimeMesh &mesh, const iteration &iter,
                const flowModel &flows, realArray &dt, const cellRealArray &shockFactor,
                const faceVector2DArray &rai, const faceVector2DArray &rav,
                const integerArray &temperatureFlag, const integerArray &pressureFlag,
                const cellIntegerArray &CFLFlag);

        virtual ~pseudoTime() noexcept;

        /**
         * \brief To compute the time-step for time marching.
         * \return The minimum timestep.
         */
        virtual void computingTimeStep() = 0;

        /**
         * \brief To compute the time-step for time marching.
         * \return The minimum timestep.
         */
        virtual void computingTimeStep(realArray &dt, const real cfl0) = 0;

        /**\brief Note cell(n).AR is the cell aera aspect ratio. It is used to accelerate the convergence
         *  for stretched meshes.
         */
        inline bool unsetIsStretchAc();

        /**\brief The Courant-Friedrichs-Lewy (CFL) numbers.*/
        hur_nodiscard inline CFL &cfl() noexcept;
        /**\brief The Courant-Friedrichs-Lewy (CFL) numbers.*/
        hur_nodiscard inline const CFL &cfl() const noexcept;

        /**\brief COnstant for time step computation*/
        hur_nodiscard inline real CForTimeStep() const noexcept;

    protected:
        bool isShockReduce_;
        real minShockReduceFct_;

        inline void shockReduce();

        void limitFlagReduce();

        void nonorthogonalMeshReduce();

        void stretchedMeshAccelerate();

        void stretchedMeshAccelerate(realArray &dt);

        /**
         * \brief To compute the global time-step.
         * \return The minimum timestep [s].
         */
        real globalTimeStep() const;

        /**
         * \brief To compute the global time-step.
         * \return The minimum timestep [s].
         */
        real globalTimeStep(realArray &dt) const;

    public:
        real restrictTimeStep();

        real timeStep();

        /**
         * \brief To compute the global time-step.
         * \return The minimum timestep [s].
         */
        real getGlobalTimeStep(const real cfl0, const bool noScale = false);
    };
} // namespace OpenHurricane

#include "pseudoTime.inl"