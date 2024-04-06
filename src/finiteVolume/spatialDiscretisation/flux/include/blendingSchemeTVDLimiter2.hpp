/*!
 * \file blendingSchemeTVDLimiter2.hpp
 * \brief Headers of a mix of centered and upwind-biased Riemann flux scheme with a TVD limiter 2.
 *        The subroutines and functions are in the <i>blendingSchemeTVDLimiter2.cpp</i> file.
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
#include "limitersForLinear.hpp"
#include "shockSensor.hpp"
#include "spatialScheme.hpp"
#include "upwindScheme.hpp"
namespace OpenHurricane {

    /*!
     * \brief The class of the mix of centered and upwind-biased Riemann flux scheme with a TVD limiter 2.
     */
    class blendingSchemeTVDLimiter2 : public spatialScheme {
    protected:
        uniquePtr<upwindScheme> upwindMethodPtr_;

        realArray fluxC_;
        realArray fluxU_;

        mutable faceRealArray vnSf_;

        shockSensor shocks_;

        cellRealArray upwindFactor_;

        real minAlphaf_;
        real maxAlphaf_;

        mutable realArray *alphafPtr_;
        void makeAlphaf() const;

        void cntUpwdBlendFct();

        bool limitTemperature_;

        real THighLimit_;
        real TLowLimit_;

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
                              const vector &faceArea, const real blend, realArray &flux) const;

        void upwindFlux(const real rhol, const real rhor, const vector &VL, const vector &VR,
                        const real pl, const real pr, const real gl, const real gr, const real el,
                        const real er, const real cl, const real cr, const vector &faceArea,
                        const real blend, realArray &flux) const;

        void pureCentralFlux(const real rhol, const real rhor, const vector &VL, const vector &VR,
                             const real pl, const real pr, const real el, const real er,
                             const vector &faceArea, const integer fi, const real fw,
                             realArray &flux) const;

    public:
        declareClassName(blendingSchemeTVDLimiter2);

        /**\brief Construct null.*/
        blendingSchemeTVDLimiter2() = delete;

        blendingSchemeTVDLimiter2(const controller &cont, const runtimeMesh &mesh, flowModel &flow);

        /**\brief Destructor.*/
        virtual ~blendingSchemeTVDLimiter2() noexcept;

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

        hur_nodiscard inline bool switchToUpwind(const integer faceI) const {
            const auto cl = mesh().faces()[faceI].leftCell();
            const auto cr = mesh().faces()[faceI].rightCell();
            return shocks_.sensor()[cl] == 1.0 || shocks_.sensor()[cr] == 1.0;
        }

    protected:
        hur_nodiscard inline real actualUpwdFct(const integer faceI) const {
            const auto &cl = mesh().faces()[faceI].leftCell();
            const auto &cr = mesh().faces()[faceI].rightCell();
            return min(max((*alphafPtr_)[faceI], max(upwindFactor_[cl], upwindFactor_[cr])),
                       maxAlphaf_);
        }

        void extendShockSensor();
    };

    template <class Type>
    inline void
    blendingSchemeTVDLimiter2::invFluxTemplate(geometryArray<Type, cellMesh> &cellQ) const {
        LFatal("Attempt to access an empty function");
    }
} // namespace OpenHurricane
