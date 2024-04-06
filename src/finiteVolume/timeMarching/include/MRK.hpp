/*!
 * \file MRK.hpp
 * \brief Headers of class of the multi-stage Runge-Kutta.
 *        The subroutines and functions are in the <i>MRK.cpp</i> file.
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
#include "timeMarching.hpp"

namespace OpenHurricane {

    /*!\brief The base class for MRK scheme.*/
    class MRK : public timeMarching {
    private:
        /**\brief Multi-stage coefficients.*/
        realArray alpha_;

        real sigma_;

        /**\brief The index of variable of density.*/
        integer rhoOi_;

        /**\brief Smoothing coefficient. 0.3~0.8*/
        real ep_;

        /**\brief The maximum Jacobi iterations. Usually two is sufficient.*/
        integer maxJacStep_;

        /**\brief The ratio  to increase the CFl number for implicit residual smoothing.*/
        real cflRatioForImpResSmoo_;

    protected:
        /**\brief Using implicit residual smoothing.*/
        bool impResSmooth_;

    public:
        declareClassNames;

        /*!\brief Disallow null constructor.*/
        MRK() = delete;

        /*!\brief Disallow copy constructor.*/
        MRK(const MRK &) = delete;

        MRK(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
            solver &_solver, cellVectorArray &v);

        /*!\brief Destructor.*/
        virtual ~MRK() noexcept;

        void stageUpdate();

        void updatePrimitives();

        void restoreConserves();

        void setStageCoeffs(const string &w);

        virtual void timeStep();

        virtual void initializing();

        virtual void marching();

        void implicitResidualSmoothing(cellRealArray &Q) const;

        void implicitResidualSmoothing(cellVectorArray &Q) const;

        void implicitResidualSmoothing();
    };

} // namespace OpenHurricane