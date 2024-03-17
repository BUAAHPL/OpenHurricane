/*!
 * \file ESDIRK.hpp
 * \brief Headers of class of the ESDIRK dual time stepping.
 *        The subroutines and functions are in the <i>ESDIRK.cpp</i> file.
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

namespace OpenHurricane {

    /*!\brief The base class for ESDIRK scheme.
     * first-stage Explicit, Single Diagonal coefficient, Implicit Runge-Kutta.
     */
    template <class timeMethod> class ESDIRK : public timeMethod {
    public:
        enum ESDIRKType : short { ESDIRK3, ESDIRK4 };

        enum class iterationMethod : short { pseudoTimeStepping, NewtonIteration };

        enum class newTimeExterpolationMethod : short { lastTime, Lagrange };

    private:
        real gamma_;
        realArray b_;
        realSquareMatrix a_;

        /**\brief Stage of RK.*/
        integer stage_;

        /**\brief Current stage*/
        integer k_;

        ESDIRKType type_;

        integer rhoOi_;

        integer m0_;

        //realArray dtn_;

        iterationMethod iterMethod_;

        newTimeExterpolationMethod newExterpolate_;

        void setABG();

        void initialSub();

    protected:
        void setPseudoTimeStepScaleFactor(const integer subSteps);
        void setSolverWrite();
        void setESDIRKType(const controller &timeMethodCont);

    public:
        declareClassName(ESDIRK);
      

        /*!\brief Disallow null constructor.*/
        ESDIRK() = delete;

        /*!\brief Disallow copy constructor.*/
        ESDIRK(const ESDIRK &) = delete;

        ESDIRK(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
               solver &_solver, cellVectorArray &v);

        /*!\brief Destructor.*/
        virtual ~ESDIRK() noexcept;

        virtual bool explicitSource() const noexcept;

        virtual bool diagonalImpSource() const noexcept;

        virtual void initializing();

        virtual void updateOld();

        virtual void updateRhs();

        void storeRhs();

        virtual void marching();

        void stageUpdate();

        hur_nodiscard inline real getPhyTimeStep() const noexcept {
            return timeMethod::iter().pTStep().pTimeStep();
        };

        virtual hur_nodiscard bool isESDIRKScheme() const noexcept;

        virtual void timeStep();

    protected:
        void updateDyPhyTimeStep();
    };

} // namespace OpenHurricane

#include "ESDIRK.inl"