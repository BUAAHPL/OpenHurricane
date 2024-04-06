/*!
 * \file BDF123.hpp
 * \brief Headers of class of the BDF123 dual time stepping.
 *        The subroutines and functions are in the <i>BDF123.cpp</i> file.
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

/*!\brief The base class for BDF123 scheme.*/
template <class timeMethod>class BDF123 : public timeMethod {
public:
    enum BDFType : short { BDF1, BDF2, BDF3 };

    enum class newTimeExterpolationMethod : short { lastTime, Lagrange };

    enum class iterationMethod : short { pseudoTimeStepping, NewtonIteration };

private:
    realArray alpha_;

    BDFType type_;

    void setAlpha(const BDFType typ = BDF1);

    iterationMethod iterMethod_;

    newTimeExterpolationMethod newExterpolate_;

    integer m0_;

    integer storeSize_;
    integer alphaSize_;

    //realArray dtn_;

    integer rhoi_;

    void initialSub();

    /**
     * \brief If true, it is to use a linear dt or CFL for the initial stage of calculation.
     * \note Not work for restart from unsteady case.
     * \note Default is false.
     */
    bool useLinearDtCFL_;

    /**
     * \brief The steps to hold small dt or CFL for initial stage of calculation.
     * \note Not work for restart from unsteady case.
     * \note Default value is 10.
     */
    integer dtCFLConst_;

    /**
     * \brief The steps to grow dt or CFL from the small value for initial stage of calculation to the given one.
     * \note Not work for restart from unsteady case.
     * \note Default value is 10.
     */
    integer dtCFLLinear_;

    real dtCFL0_;
    real dtFinal_;

protected:
    void setPseudoTimeStepScaleFactor(const integer subSteps);

    void setSolverWrite();
    void setBDFType(const controller &timeMethodCont);

public:
    // Static data
    declareClassName(BDF123);
    // Constructors.

    /*!\brief Disallow null constructor.*/
    BDF123() = delete;

    /*!\brief Disallow copy constructor.*/
    BDF123(const BDF123 &) = delete;

    BDF123(const controller &cont, const runtimeMesh &mesh,
           const flowModel &flowM, solver &_solver, cellVectorArray &v);

    /*!\brief Destructor.*/
    virtual ~BDF123() noexcept;

    virtual bool explicitSource() const noexcept;

    virtual bool diagonalImpSource() const noexcept;

    virtual void initializing();

    virtual void updateOld();

    virtual void updateRhs();

    virtual void updateWStar();

    virtual void marching();

    hur_nodiscard inline real getPhyTimeStep() const noexcept {
        return timeMethod::iter().pTStep().pTimeStep();
    };

    virtual hur_nodiscard bool isBDFScheme() const noexcept;

protected:
    void updateDyPhyTimeStep();
};

} // namespace OpenHurricane

#include "BDF123.inl"