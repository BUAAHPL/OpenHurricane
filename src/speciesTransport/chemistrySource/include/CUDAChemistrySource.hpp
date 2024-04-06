/*!
 * \file CUDAChemistrySource.hpp
 * \brief Header of chemistry source in CUDA platform.
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
#include "EulerCUDA.hpp"
#include "chemistrySource.hpp"
#include "cudaStreams.hpp"
#include "nThreadsAndBlocks.hpp"
#include "reactionTableCUDA.hpp"

namespace OpenHurricane {
    /**
     * \brief The class of chemistry source in CUDA platform.
     */
    class CUDAChemistrySource : public chemistrySource {
    private:
        mutable cuChem::reactionTable *reactionsCUDAPtr_;

        mutable CUDAChemistrySourceODEs *CUDAChemODEsPtr_;
        void makeCUDAChemODEs() const;

        mutable EulerCUDA *EulerSolverPtr_;
        void makeEulerSolver() const;

        void getYiRhoT(real *hur_restrict hostYiRhoTPtr_) const;
        void getYiRhoTSlice(real *hur_restrict hostYiRhoTPtr_, const integer pos,
                            const integer offset) const;

        enum class ODEsSolverType : short { Euler, MTSEuler };

        ODEsSolverType ODEsSolverType_;

        nThreadsAndBlocks nTB_;
        realArray subDt_;

    public:
        declareClassNames;

        // Constructors

        /*!\brief Disallow null constructor.*/
        CUDAChemistrySource() = delete;

        /**
         * \brief Disallow bitwise copy constructor.
         */
        CUDAChemistrySource(const CUDAChemistrySource &) = delete;

        /*!\brief Construct from flow and controller.*/
        CUDAChemistrySource(flowModel &flows, const controller &cont);

        /*!\brief Destructor.*/
        virtual ~CUDAChemistrySource() noexcept;

        virtual void calculateSourceTerms(const bool withoutLastSpec = true);
        virtual void calculateSourceTermsImp(const bool withoutLastSpec = true);

        virtual void calculateSourceTermsAsync(real *hur_restrict hostYiRhoTPtr_,
                                               cu2DArray<cu_real> &dYiRhoT,
                                               real *hur_restrict RRi, cu2DArray<cu_real> &dRi,
                                               const cudaStreams &streams);

        virtual void calculateSourceTermsAsyncSlice(real *hur_restrict hostYiRhoTPtr_,
                                                    cu2DArray<cu_real> &dYiRhoT,
                                                    real *hur_restrict RRi,
                                                    cu2DArray<cu_real> &dRi,
                                                    const cudaStreams &streams, const integer pos,
                                                    const integer offset);

        virtual void calculateSourceTermsImpAsync(
            real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT,
            real *hur_restrict RRi, cu2DArray<cu_real> &dRi, real *hur_restrict dRidrhoyi,
            cu2DArray<cu_real> &ddRidrhoyi, const cudaStreams &streams);

        virtual void calculateSourceTermsImpAsyncSlice(
            real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT,
            real *hur_restrict RRi, cu2DArray<cu_real> &dRi, real *hur_restrict dRidrhoyi,
            cu2DArray<cu_real> &ddRidrhoyi, const cudaStreams &streams, const integer pos,
            const integer offset);

        virtual void calculateSourceTermsImpAsyncHybrid(
            real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT,
            real *hur_restrict RRi, cu2DArray<cu_real> &dRi, cu_float *hur_restrict dRidrhoyi,
            cu2DArray<cu_float> &ddRidrhoyi, const cudaStreams &streams);

        virtual void calculateSourceTermsImpAsyncHybridSlice(
            real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT,
            real *hur_restrict RRi, cu2DArray<cu_real> &dRi, cu_float *hur_restrict dRidrhoyi,
            cu2DArray<cu_float> &ddRidrhoyi, const cudaStreams &streams, const integer pos,
            const integer offset);

        virtual void createReactionTable() const;
        virtual void createReactionTableAsync(const cudaStreams &streams) const;

        virtual void destroyReactionTable() const;

        inline CUDAChemistrySourceODEs &CUDAChemODEs();
        inline const CUDAChemistrySourceODEs &CUDAChemODEs() const;

        inline EulerCUDA &EulerSolver();
        inline const EulerCUDA &EulerSolver() const;

        virtual real solveTest(const real t, const real dT, real &_p, real &_T, real &_rho,
                               realArray &yi, fileOsstream &fos);

        /*!\brief Solve chemical source terms.
         * \param[in] dt - timestep of every cell for reacting.
         * \param[in] dtFactor - the facetor to scale dt.
         * \return Return the minimum timestep.
         */
        virtual real solve(const realArray &dt, const real dtFactor);

        /*!\brief Solve chemical source terms and update species mass fraction and temperature.
         * \param[in] dt - timestep of every cell for reacting.
         * \param[in] dtFactor - the facetor to scale dt.
         */
        virtual void solveUpdate(const realArray &dt, const real dtFactor);

        virtual void solve(cellRealArray &rhoi, cellRealArray &pi, cellRealArray &Ti,
                           PtrList<cellRealArray> &yii, const realArray &dt, realArray &subDt,
                           const real dtFactor);

        virtual void solve(cellRealArray &rhoi, cellRealArray &pi, cellRealArray &Ti,
                           PtrList<cellRealArray> &yii, const realArray &dt, realArray &subDt,
                           const real dtFactor, const realArray &odeFactor);
    };
} // namespace OpenHurricane
#include "CUDAChemistrySource.inl"

#endif // CUDA_PARALLEL