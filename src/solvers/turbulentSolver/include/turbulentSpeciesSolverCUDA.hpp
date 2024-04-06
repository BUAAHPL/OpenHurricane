/*!
 * \file turbulentSpeciesSolverCUDA.hpp
 * \brief Headers of turbulent species transport Solver with CUDA accelerating.
 *        The subroutines and functions are in the <i>turbulentSpeciesSolverCUDA.cpp</i> file.
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

#include "SST.hpp"
#include "SpalartAllmaras.hpp"
#include "combustionModel.hpp"
#include "solver.hpp"
#ifdef CUDA_PARALLEL

namespace OpenHurricane {
    /*!\brief The class of turbulent species transport solver with CUDA accelerating.*/
    class turbulentSpeciesSolverCUDA : public solver {
    protected:
        uniquePtr<combustionModel> chemtryPtr_;

        integer rhoId_;
        integer rhouId_;
        integer rhoEId_;
        integer rhoYi0Id_;
        integer rhoTurb0Id_;

    public:
        /**
         * \brief The convective flux Jacobian.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[in] Vt - The contravariant velocity of the face of the control volume is set to zero for stationary grids
         * \param[out] AC - The convective flux Jacobian
         */
        virtual void Ac(const integer celli, const vector &normal, const real Vt,
                        realSquareMatrix &AC) const;

        /**
         * \brief The derivative of pressure with respect to the conservative variables.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[out] pq - The derivative of pressure with respect to the conservative variables
         */
        virtual void dpdq(const integer celli, const vector &normal, realArray &pq) const;

        /**
         * \brief The convective flux Jacobian.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[in] dq - The increment of the conservative variables
         * \param[out] adq - The convective flux Jacobian multiplied by the increment of the conservative variables
         */
        virtual void Acdq(const integer celli, const vector &normal, const realArray &dq,
                          realArray &adq) const;

    public:
        declareClassName(turbulentSpeciesSolverCUDA);

        turbulentSpeciesSolverCUDA(iteration &iter, const runtimeMesh &mesh);

        /*!\brief Destructor.*/
        virtual ~turbulentSpeciesSolverCUDA() noexcept;

        /**
         * \brief Solve the problem.
         */
        virtual void solving();

        /**
         * \brief Solve the problem with BDF method.
         */
        virtual void BDFSolve();

        /**
         * \brief Update boundary conditions.
         */
        virtual void bc();

        virtual void updateProperties();

        virtual void initialize();

        virtual void timeStep(realArray &dt);

        /**
         * \brief Evaluate the convective flux.
         */
        virtual void calculateFc();

        /**\brief Compute the viscous flux.*/
        virtual void calculateFv();

    private:
        void actualCalFv();

    public:
        /**\brief Caculate the source terms.*/
        virtual void calculateSource();

        virtual void updatePrimitives(const bool shouldUpdateTemp = false);

        virtual void calculateOtherTerms();

        virtual void updateFlowOld();

        virtual void write();

        virtual hur_deprecated void testChemical();
    };
} // namespace OpenHurricane

#include "turbulentSpeciesSolverCUDA.inl"

#endif // CUDA_PARALLEL