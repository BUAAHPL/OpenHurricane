/*!
 * \file JST.hpp
 * \brief Headers of JST central scheme.
 *        The subroutines and functions are in the <i>JST.cpp</i> file.
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
#include "shockSensor.hpp"
#include "spatialScheme.hpp"

namespace OpenHurricane {

    /*!
     * \brief The class of JST central scheme.
     */
    class JST : public spatialScheme {
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
         * \param[in] faceArea - The area vector of this face.
         * \param[in] MasInf - The Mach number of infinity.
         * \param[out] flux - The flux of this face.
         * \return The reaction rate conatants.
         */
        virtual void calcFlux(const real rhol, const real rhor, const vector &VL, const vector &VR,
                              const real pl, const real pr, const real gl, const real gr,
                              const real el, const real er, const real cl, const real cr,
                              const vector &faceArea, const real blend, realArray &flux) const;

        void spectralRadius();

        /** \brief Spectral radius at face.*/
        realArray AIJ_;

        /** \brief The geometrical weights. */
        vector2DArray thetaIJ_;

        realArray espIJ2_;
        realArray espIJ4_;

        vector2DArray rhoLR_;

        cellRealArray pgi_;

        /**
         * \brief If use a factor to modify pressure sensor to reduce numerical dissipation.
         */
        bool useModifyFactor_;

        /**
         * \brief The factor to modify pressure sensor to reduce numerical dissipation.
         */
        cellRealArray *modifyFactorPtr_;

        shockSensor shocks_;

        void calcModifyFactor();

        realArray vnf_;

        real k2_;
        real k4_;

        void getGeometricalWeights();

        void pressureSensor();

        void calcEsp();

        real artificialDissipation(const integer fi, const real QJI, const real WRL) const;

        vector artificialDissipation(const integer fi, const vector &QJI, const vector &WRL) const;

        void addArtificialDissipation(const integer fi, const real QJI, const real WRL,
                                      const integer cl, const integer cr, cellRealArray &rhs) const;

        void addArtificialDissipation(const integer fi, const vector QJI, const vector WRL,
                                      const integer cl, const integer cr,
                                      cellVectorArray &rhs) const;

    public:
        declareClassNames;

        /**\brief Construct null.*/
        JST() = delete;

        JST(const controller &cont, const runtimeMesh &mesh, flowModel &flow);

        /**\brief Destructor.*/
        virtual ~JST() noexcept { HurDelete(modifyFactorPtr_); }

        /*!\brief Calculate the convective fluxes for the continuity, momentum and energy equations.*/
        virtual void basicFlux();

        virtual void invFlux(cellRealArray &cellQ) const;
        virtual void invFlux(cellVectorArray &cellQ) const;

        virtual void invFlux(cellRealArray &cellQ, const realBounded &bound) const;
        virtual void invFlux(cellVectorArray &cellQ, const vectorBounded &bound) const;

        /**
         * \brief Calc inviscous flux for species.
         */
        virtual void invFluxSpecies(PtrList<cellRealArray> &yi,
                                    const bool withLastSpc = false) const;
    };

} // namespace OpenHurricane