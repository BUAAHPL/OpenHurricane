/*!
 * \file cuSparSolverTest.hpp
 * \brief Headers of class of the cuSparSolverTest.
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
#include "timeMarching.hpp"

#ifdef CUDA_PARALLEL
#include "cuSparSolver.hpp"

namespace OpenHurricane {

    /*!\brief The base class for AMG scheme.*/
    class cuSparSolverTest : public timeMarching {
    private:
        const controller &sparSolCont_;
        real omega_;

        real betaT_;

        uniquePtr<cuCRSBlockMatrix<cu_real>> APtr_;
        uniquePtr<cuBlockArray1<cu_real, cu_integer>> dqPtr_;

        uniquePtr<sparseSolver::cuRealBlockSparSolver> sparSolver_;

    protected:
        bool modifiedDiagnalTime_;

        real alphaDT_;

        real scaleTimeStep_;

        integer nIter_;

        /**
         * \brief To reserve GPU temporary memory of APtr_ and dqPtr_. Default is false.
         * If false, then sparSolver_ will also not be reserved.
         */
        bool isReserveGPUTmpMem_;

    public:
        // Static data
        declareClassName(cuSparSolverTest);

        /*!\brief Disallow null constructor.*/
        cuSparSolverTest() = delete;

        /*!\brief Disallow copy constructor.*/
        cuSparSolverTest(const cuSparSolverTest &) = delete;

        cuSparSolverTest(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
                         solver &_solver, cellVectorArray &v,
                         const bool modifiedDiagnalTime = false);

        /*!\brief Destructor.*/
        virtual ~cuSparSolverTest() noexcept;

        void update();

        virtual void initializing();

        virtual void marching();

        hur_nodiscard inline virtual bool updatedTemperature() const noexcept { return true; }

    protected:
        /**
         * \brief To construct the operator \f$A\f$ in \f[{\mathbf{A}}x = b\f].
         */
        void constructOperator();
    };

} // namespace OpenHurricane
#endif
