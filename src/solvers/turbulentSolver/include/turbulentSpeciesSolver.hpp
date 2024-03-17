/*!
 * \file turbulentSpeciesSolver.hpp
 * \brief Headers of turbulent species transport Solver.
 *        The subroutines and functions are in the <i>turbulentSpeciesSolver.cpp</i> file.
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
#include "combustionModel.hpp"
#include "solver.hpp"

namespace OpenHurricane {
    /*!\brief The class of turbulent species transport solver.*/
    class turbulentSpeciesSolver : public solver {
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
        declareClassName(turbulentSpeciesSolver);

        turbulentSpeciesSolver(iteration &iter, const runtimeMesh &mesh);

        /*!\brief Destructor.*/
        virtual ~turbulentSpeciesSolver() noexcept;

        virtual void solving();

        virtual void BDFSolve();

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

        virtual void write();

        virtual hur_deprecated void testChemical();

        inline virtual void previousSource();
        inline virtual void postSource();

    protected:
        // Test time

#ifdef TEST_PROCESS_TIME

        mutable fileOsstream fosTestChemTime_;
#endif
        void setCellLoadWeight(const real sorTime, const real otherTime);

        /**
         * \brief Source terms time cost.
         */
        real sorTime_;
        real otherTime_;
        bool calcLoadWeight_;

        hur_nodiscard bool isUpdatingCalcTime() const noexcept;

#ifdef CUDA_PARALLEL

        mutable uniquePtr<integerList> fzBoundStartPtr_;
        void updatePropertiesCUDA();
        void getYiTP(real *hur_restrict hostYiTPPtr_) const;
        void muKappaDiffBoundary();
        void makeFzBoundStart() const;
        void getYiTPBoundary(real *hur_restrict hostYiTPPtr_) const;
        hur_nodiscard inline const integerList &fzBoundStart() const;

        void setExtrapolate(const real *hur_restrict hostDimmMuKap,
                            const real *hur_restrict hostHaiCp);

        void transferTranP();

#endif // CUDA_PARALLEL
    };
} // namespace OpenHurricane

#include "turbulentSpeciesSolver.inl"