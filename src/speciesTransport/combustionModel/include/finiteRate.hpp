/*!
 * \file finiteRate.hpp
 * \brief Header of finite-rate combustion model.
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
#include "combustionModel.hpp"

namespace OpenHurricane {
    /*!\brief The class of laminar finite rate combustion model.*/
    class finiteRate : public combustionModel {
    private:
        cellRealArray tchField_;
        cellRealArray tmpNField_;
        cellRealArray tmpYoField_;
        cellRealArray tmpYfField_;

        integer oIdx_;
        integer fIdx_;

    protected:
        real tauIntFactor_;

    public:
        declareClassNames;

        // Constructors

        /*!\brief Disallow null constructor.*/
        finiteRate() = delete;

        /*!\brief Construct from flow and controller.*/
        finiteRate(flowModel &flows, const controller &cont, turbulenceModel &turb);

        /*!\brief Disallow copy constructor.*/
        finiteRate(const finiteRate &) = delete;

        /*!\brief Destructor.*/
        virtual ~finiteRate() noexcept {}

        virtual void calcSourceTerms(const realArray &dt);

        virtual void expChemistrySource(realArray &dt, const bool isModifiedDt = false);

        virtual void impCalcSourceTerms(const realArray &dt);
        virtual void impChemistrySource(realArray &dt, const bool isModifiedDt = false);

        virtual void getChemistrySource();

        virtual void getImpChemistrySource(realArray &dt, const bool isModifiedDt = false);

        /**
         * \brief Full point-implicit for NS equations.
         *        The index of species in Jacobian matrix must be continuous.
         * \param[in,out] dt - Timestep
         * \param[out] Jac - The chemical source terms Jacobian based on conservative variables
         * \param[in] rhoId - The index of density rho in Jac
         * \param[in] rhouId - The index of rhou in Jac (rhovId=rhouid+1, rhowId=rhovId+1)
         * \param[in] rhoEId - The index of rhoE in Jac
         * \param[in] rhoYi0Id - The index of rhoYi0 in Jac (rhoYi1Id=rhoYi0Id+1,...)
         */
        virtual void fullPointImpChemistrySource(realArray &dt, cellRealSquareMatrixArray &Jac,
                                                 const integer rhoId, const integer rhouId,
                                                 const integer rhoEId, const integer rhoYi0Id);

        /**
         * \brief Return heat release rate [J/(m^3 s)].
         */
        virtual realArray heatReleaseRate();

        /*!
         * \brief Solve a constant volume reactor for all internal field.
         * The species mass fraction field and the temperature field will be changed.
         */
        virtual void constVolReactor(const realArray &dt, const real dtFactor,
                                     const bool firstCall = true);

        /*!\brief Disallow bitwise assignment.*/
        void operator=(const finiteRate &) = delete;

        /*!\brief Chemical source terms.*/
        virtual realArrayArray omegai() const;

    protected:
        /**
         * \brief Caculate chemical source terms explicitly in coupled scheme.
         */
        virtual void expChemistrySourceCoupled();

        /**
         * \brief Caculate chemical source terms implicitly with diagonal Jacobian matrices in coupled scheme.
         */
        virtual void diagImpChemistrySourceCoupled();

        virtual void fullImpChemistrySourceCoupled(cellRealSquareMatrixArray &Jac,
                                                   const integer rhoId, const integer rhouId,
                                                   const integer rhoEId, const integer rhoYi0Id);

        /**
         * \brief Caculate chemical source terms in integrated scheme.
         */
        virtual void chemistrySourceIntegrated();

        /**
         * \brief Caculate chemical source terms in strangSplitted scheme.
         */
        virtual void chemistrySourceStrangSplitted(const realArray &dt, const real dtFactor);
    };
} // namespace OpenHurricane
