/*!
 * \file AMGSchemeTest.hpp
 * \brief Headers of class of the AMGSchemeTest.
 *        The subroutines and functions are in the <i>AMGSchemeTest.cpp</i> file.
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
#include "sparSolver.hpp"
#include "timeMarching.hpp"

namespace OpenHurricane {

    /*!\brief The base class for AMG scheme.*/
    class AMGSchemeTest : public timeMarching {
    private:
        real omega_;

        real betaT_;

        CRSMatrix<realSquareMatrix> A_;

        uniquePtr<sparseSolver::realBlockSparSolver> sparSolver_;

    protected:
        bool modifiedDiagnalTime_;

        real alphaDT_;

        real scaleTimeStep_;

        integer nIter_;

    public:
        declareClassName(AMGSchemeTest);

        /*!\brief Disallow null constructor.*/
        AMGSchemeTest() = delete;

        /*!\brief Disallow copy constructor.*/
        AMGSchemeTest(const AMGSchemeTest &) = delete;

        AMGSchemeTest(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
                      solver &_solver, cellVectorArray &v, const bool modifiedDiagnalTime = false);

        /*!\brief Destructor.*/
        virtual ~AMGSchemeTest() noexcept;

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
