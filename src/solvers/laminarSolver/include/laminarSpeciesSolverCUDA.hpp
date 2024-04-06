/*!
 * \file laminarSpeciesSolverCUDA.hpp
 * \brief Headers of laminar species transport Solver using CUDA accelerating chemical source terms calculation.
 *        The subroutines and functions are in the <i>laminarSpeciesSolverCUDA.cpp</i> file.
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
#ifdef CUDA_PARALLEL
#include "combustionModel.hpp"
#include "finiteRate.hpp"
#include "solver.hpp"

namespace OpenHurricane {
    /*!\brief The class of laminar species solver using CUDA accelerating chemical source terms calculation.*/
    class laminarSpeciesSolverCUDA : public solver {
    protected:
        /*!\brief Delete the pointer.*/
        template <class Type> inline void deletePointer(geometryArray<Type, cellMesh> *_gPtr) const;

        uniquePtr<combustionModel> chemtryPtr_;

        integer rhoId_;
        integer rhouId_;
        integer rhoEId_;
        integer rhoYi0Id_;

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
        declareClassName(laminarSpeciesSolverCUDA);


        laminarSpeciesSolverCUDA(iteration &iter, const runtimeMesh &mesh);

        /*!\brief Destructor.*/
        virtual ~laminarSpeciesSolverCUDA() noexcept;

        virtual void solving();

        virtual void BDFSolve();

        

        /*!\brief Clear the solver.*/
        virtual void clear() noexcept;

        /**
         * \brief Update boundary conditions.
         */
        virtual void bc();

        void timeStep(realArray &dt);

        /**
         * \brief Update mixture physical and chemical properties.
         */
        virtual void updateProperties();

        /**
         * \brief Initialize the calculation.
         */
        virtual void initialize();

        /**
         * \brief Compute convective flux.
         */
        virtual void calculateFc();

        /**
         * \brief Compute viscous flux.
         */
        virtual void calculateFv();

        /**
         * \brief Compute viscous flux.
         */
        virtual void actualCalFv();

        /**
         * \brief Compute source terms.
         */
        virtual void calculateSource();

        virtual void updatePrimitives(const bool shouldUpdateTemp = false);

        virtual void updateFlowOld();

        // Calculate output variables

        /**
         * \brief Calculate stagnation temperature.
         */
        //virtual void calcTotalTemperature()const;

        /**
         * \brief Calculate heat release rate for chemical reactions.
         */
        //void calcHeatReleaseRate()const;

        /**
         * \brief The Takeno Flame Index.
         */
        //void calcGFO()const;

        /**
         * \brief Calculate output variables.
         */
        //virtual void calcOutput()const;

        // Write

        /**
         * \brief Write the residuals and results.
         */
        virtual void write();
    };
} // namespace OpenHurricane

#include "laminarSpeciesSolverCUDA.inl"
#endif // CUDA_PARALLEL