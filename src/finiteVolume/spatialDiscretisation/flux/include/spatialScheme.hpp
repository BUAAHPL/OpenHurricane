/*!
 * \file spatialScheme.hpp
 * \brief Headers of spatial scheme.
 *        The subroutines and functions are in the <i>spatialScheme.cpp</i> file.
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

#include "firstOrderReconstruction.hpp"
#include "flowModel.hpp"
#include "realArray.hpp"
#include "reconstruction.hpp"
#include "rhoThermo.hpp"
#include "vectorArray.hpp"

namespace OpenHurricane {
    /*!
     * \brief The base class of spatial scheme.
     */
    class spatialScheme {
    protected:

        const runtimeMesh &mesh_;

        flowModel &flows_;

        /*!\brief Density flux field.*/
        faceRealArray rhoFlux_;

        /*!\brief Velocity field.*/
        cellVectorArray &v_;

        /*!\brief Thermo.*/
        const rhoThermo &thermo_;

        /*!\brief Reconstruction method pointer.*/
        uniquePtr<reconstruction> reconstrPtr_;

        /*!\brief The list of primitive objects.*/
        List<object *> objectList_;

        integer countParam_;

        integerVector2DList paramMap_;

        realArray flux_;

    protected:
        /**
         * \brief The function for calculating convective flux.
         * \param[in] rhol - The left state density.
         * \param[in] rhor - The right state density.
         * \param[in] VL - The left state velocity.
         * \param[in] VR - The right state velocity.
         * \param[in] pl - The left state pressure.
         * \param[in] pr - The right state pressure.
         * \param[in] gl - The left state specific heat ratio.
         * \param[in] gr - The right state specific heat ratio.
         * \param[in] el - The left state total energy.
         * \param[in] er - The right state total energy.
         * \param[in] cl - The left state acoustic speed.
         * \param[in] cr - The right state acoustic speed.
         * \param[in] MasInf - The Mach number of infinity.
         * \param[out] flux - The flux of this face.
         * \return The reaction rate conatants.
         */
        virtual void calcFlux(const real rhol, const real rhor, const vector &VL, const vector &VR,
                              const real pl, const real pr, const real gl, const real gr,
                              const real el, const real er, const real cl, const real cr,
                              const vector &faceArea, const real blend, realArray &flux) const = 0;

        cellRealArray &shockFactor_;

    public:
        declareClassNames;
        declareObjFty(spatialScheme, controller,
                      (const controller &cont, const runtimeMesh &mesh, flowModel &flow),
                      (cont, mesh, flow));

        /**\brief Construct null.*/
        inline spatialScheme() = delete;

        spatialScheme(const controller &cont, const runtimeMesh &mesh, flowModel &flow);

        static uniquePtr<spatialScheme> creator(const controller &cont, const runtimeMesh &mesh,
                                                flowModel &flow);

        static uniquePtr<spatialScheme> creator(const string &schemeType, const controller &cont,
                                                const runtimeMesh &mesh, flowModel &flow);

        /**\brief Destructor.*/
        virtual ~spatialScheme() noexcept;

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept;

        hur_nodiscard inline faceRealArray &rhoFlux() noexcept { return rhoFlux_; }

        hur_nodiscard inline const faceRealArray &rhoFlux() const noexcept { return rhoFlux_; }

        /*!\brief Calculate the convective fluxes for the continuity, momentum and energy equations.*/
        virtual void basicFlux();

        template <class Type> void invFluxTemplate(geometryArray<Type, cellMesh> &cellQ) const;

        virtual void invFlux(cellRealArray &cellQ) const;
        virtual void invFlux(cellVectorArray &cellQ) const;

        virtual void invFlux(cellRealArray &cellQ, const realBounded &bound) const;
        virtual void invFlux(cellVectorArray &cellQ, const vectorBounded &bound) const;

        /**
         * \brief Calc inviscous flux for species.
         */
        virtual void invFluxSpecies(PtrList<cellRealArray> &yi,
                                    const bool withLastSpc = false) const;

        template <class Type> void addInvFlux(geometryArray<Type, cellMesh> &cellQ);

        template <class Type> void grad(geometryArray<Type, cellMesh> &cellQ) const;

        hur_nodiscard inline const flowModel &flows() const noexcept;

        hur_nodiscard inline real MaInf() const noexcept;

    protected:
        /** \brief Is the reconstruction corrected at low mach number. Defalut is no */
        bool isLowMachCorr_;

        /**
         * \brief Dimitris's Low-Mach number correction: change velocity, keep
         *        energy constant & change static pressure in subsonic flows.
         * \param[in] rhol - density of left state
         * \param[in] rhor - density of right state
         */
        void lowMachCorrection(const real rhol, const real rhor, const real gammal,
                               const real gammar, vector &vl, vector &vr, real &pl, real &pr,
                               const vector &fA) const;
    };
} // namespace OpenHurricane

#include "spatialScheme.inl"