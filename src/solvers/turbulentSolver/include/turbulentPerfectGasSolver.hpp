/*!
 * \file turbulentPerfectGasSolver.hpp
 * \brief Headers of turbulent Perfect Gas Solver.
 *        The subroutines and functions are in the <i>turbulentPerfectGasSolver.cpp</i> file.
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

#include "SST.hpp"
#include "SpalartAllmaras.hpp"
#include "solver.hpp"

namespace OpenHurricane {
    /*!\brief The class of turbulent perfect gas solver.*/
    class turbulentPerfectGasSolver : public solver {
    private:
    protected:
        SpalartAllmaras *SAModelPtr_;

        SST *SSTModel_;
        // Protected Member Functions

        /*!\brief Delete the pointer.*/
        template <class Type> inline void deletePointer(geometryArray<Type, cellMesh> *_gPtr) const;

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
         * \brief The convective flux Jacobian.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[in] dq - The increment of the conservative variables
         * \param[out] adq - The convective flux Jacobian multiplied by the increment of the conservative variables
         */
        virtual void Acdq(const integer celli, const vector &normal, const realArray &dq,
                          realArray &adq) const;

    public:
        declareClassName(turbulentPerfectGasSolver);

        // Constructors

        turbulentPerfectGasSolver(iteration &iter, const runtimeMesh &mesh);

        /*!\brief Destructor.*/
        virtual ~turbulentPerfectGasSolver() noexcept;

        // Member Functions

        // Access

        virtual void solving();

        virtual void BDFSolve();

        

        /*!\brief Clear the solver.*/
        virtual void clear() noexcept;

        virtual void bc();

        virtual void updateProperties();

        virtual void initialize();

        virtual void timeStep(realArray &dt);

        virtual void calculateFc();

        virtual void calculateFv();

        virtual void calculateSource();

        virtual void updatePrimitives(const bool shouldUpdateTemp = false);

        virtual void calculateOtherTerms();

        virtual void updateFlowOld();

        //virtual void calcOutput()const;
        virtual void write();
    };
} // namespace OpenHurricane

#include "turbulentPerfectGasSolver.inl"