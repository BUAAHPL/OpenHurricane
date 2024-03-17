/*!
 * \file ODEsSolver.hpp
 * \brief Headers of ODEs solver.
 *        The subroutines and functions are in the <i>ODEsSolver.cpp</i> file.
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

#include "LUDecompose.hpp"
#include "controller.hpp"
#include "objectFactory.hpp"
#include "smartPointer.hpp"

namespace OpenHurricane {

    /**
     * \brief The base class of ODE solver.
     */
    class ODEsSolver {
    public:
        using dydtFunc_type = std::function<void(const real, const realArray &, realArray &)>;

        using JacobianFunc_type =
            std::function<void(const real, const realArray &, realArray &, realSquareMatrix &)>;

        using checkFunc_type = std::function<bool(real &, const realArray &, const realArray &)>;

    protected:
        uniquePtr<dydtFunc_type> dydtFuncPtr_;
        uniquePtr<JacobianFunc_type> JacobianFuncPtr_;
        uniquePtr<checkFunc_type> checkFuncPtr_;

    public:
        inline void setDyDtFunc(dydtFunc_type &func) {
            dydtFuncPtr_.reset(new dydtFunc_type(func));
        }

        inline void setJacobianFunc(JacobianFunc_type &func) {
            JacobianFuncPtr_.reset(new JacobianFunc_type(func));
        }

        inline void setCheckFunc(checkFunc_type &func) {
            checkFuncPtr_.reset(new checkFunc_type(func));
        }

        /**\brief Calculate the derivatives in dydt.*/
        inline void DyDt(const real t, const realArray &y, realArray &dydt) {
            (*dydtFuncPtr_)(t, y, dydt);
        }
        /**
         *\brief Calculate the Jacobian of the system y' = f(x,y)
         */
        inline void jacobian(const real t, const realArray &y, realArray &dfdt,
                             realSquareMatrix &dfdy) {
            (*JacobianFuncPtr_)(t, y, dfdt, dfdy);
        }

        /**
         * \brief Check the new y is acceptable.
         * \param[in,out] dt - The timestep. If unacceptable, then reduce dt.
         * \param[in] yOld - The old solution array
         * \param[in] yNew - The new solution array
         * \return True if acceptable.
         */
        hur_nodiscard inline bool checkNewY(real &dt, const realArray &yOld,
                                            const realArray &yNew) {
            if (checkFuncPtr_) {
                return (*checkFuncPtr_)(dt, yOld, yNew);
            }
            return true;
        }

    protected:
        /**\brief The number of equations*/
        mutable integer nEqs_;

        /**\brief The maximum number of equations*/
        const integer maxEqns_;

        /*!\brief The maximum steps for sub-iteration.*/
        integer maxSteps_;

        /*!\brief The absolute tolerance.*/
        realArray ATOL_;

        /*!\brief The relative tolerance.*/
        realArray RTOL_;

        real safeScale_;

        real minScale_;

        real maxScale_;

        real alphaInc_;
#ifdef TEST_PROCESS_TIME
        integer countIter_;

    public:
        hur_nodiscard inline integer countIter() const noexcept;

#endif
    public:
        /*!\brief To get the maximum error of the equations.
         * \param[in] y0 - The initial state solutions.
         * \param[in] y - The solutions.
         * \param[in] e - The errors of each equations.
         * \return Return the maximum error.
         */
        real maxError(const realArray &y0, const realArray &y, const realArray &e) const;

    protected:
        /*!\brief To get the maximum error of the specified equations.
         * \param[in] y0 - The initial state solutions.
         * \param[in] y - The solutions.
         * \param[in] e - The errors of each equations.
         * \param[in] Gm - The indeces of the specified equations to be solved.
         * \return Return the maximum error.
         */
        real maxError(const realArray &y0, const realArray &y, const realArray &e,
                      const integerListList &Gml, const integer im) const;

    private:
        /*!\brief Temporary array for equations.*/
        realArray yTemp_;

        realArray dydt0_;

    public:
        declareClassNames;
        declareObjFty(ODEsSolver, controller, (const integer nEqs, const controller &cont),
                      (nEqs, cont));

        /*!\brief Construct from ODEs.*/
        ODEsSolver(const integer nEqs);

        /*!\brief Construct from ODEs and by given absolute and relative torelance
         *		  as well as the maximum steps.
         */
        ODEsSolver(const integer nEqs, const real atol, const real rtol,
                   const integer maxStep = 1000);

        /*!\brief Construct from ODEs and the controller.*/
        ODEsSolver(const integer nEqs, const controller &cont);

        /*!\brief Selector for the specific ODE solver.*/
        static uniquePtr<ODEsSolver> creator(const integer nEqs, const controller &cont);

        /*!\brief Destructor.*/
        virtual ~ODEsSolver() noexcept {
            dydtFuncPtr_.clear();
            JacobianFuncPtr_.clear();
            checkFuncPtr_.clear();
        }

        /*!\brief Return the number of equations of the ODEs.*/
        hur_nodiscard inline integer nEqs() const noexcept;

        /*!\brief Return the number of equations of the ODEs.*/
        hur_nodiscard inline integer maxEqs() const noexcept;

        /*!\brief Return true if should reset the number of ODEs and reset it.*/
        inline virtual bool reset(const integer nEqsNew);

        template <typename Type> inline static void resetArray(List<Type> &ar, const integer n);

        template <typename Type> inline void resetArray(List<Type> &ar);

        /*!\brief Return the maximum steps for iteration.*/
        inline integer maxSteps() const noexcept;

        /*!\brief Return the absolute tolerance.*/
        hur_nodiscard inline const realArray &ATOL() const;

        /*!\brief Return the relative tolerance.*/
        hur_nodiscard inline const realArray &RTOL() const;

        hur_nodiscard inline real safeScale() const noexcept;

        hur_nodiscard inline real minScale() const noexcept;

        hur_nodiscard inline real maxScale() const noexcept;

        hur_nodiscard inline real alphaInc() const noexcept;

    public:
        /*!\brief Solve the ODEs.
         * \param[in] t0 - The initial time.
         * \param[in] y0 - The initial states.
         * \param[in] dydt0 - The derivatives.
         * \param[in] dt - The timestep.
         * \param[out] y - The solutions.
         * \return The maximum error by the specified tolerance.
         */
        virtual real solve(const real t0, const realArray &y0, const realArray &dydt0,
                           const real dt, realArray &y) = 0;

        /*!\brief Only solve the ODEs without computing error only for explicit solver.
         * \param[in] t0 - The initial time.
         * \param[in] y0 - The initial states.
         * \param[in] dydt0 - The derivatives.
         * \param[in] dt - The timestep.
         * \param[out] y - The solutions.
         */
        virtual void solveOnly(const real t0, const realArray &y0, const realArray &dydt0,
                               const real dt, realArray &y) = 0;

    protected:
        /*!\brief Solve the ODEs.
         * \param[in] t0 - The initial time.
         * \param[in] y0 - The initial states.
         * \param[in] dydt0 - The derivatives.
         * \param[in] dt - The timestep.
         * \param[out] y - The solutions.
         * \return The maximum error by the specified tolerance.
         */
        virtual real solve(const real t0, const realArray &y0, const realArray &dydt0,
                           const real dt, const integerListList &Gml, const integer im,
                           realArray &y) = 0;

        /*!\brief Solve MTS problem.
         * \param[in] t0 - The initial time.
         * \param[in] tn - The base physical time, and the time the solver should reach.
         * \param[out] tReached - The time reached by solver.
         * \param[in] dt - The time step to be used in im_th group.
         * \param[in] Gml - The groups of different time scale.
         * \param[in] im - The index of group of which the character time should be used.
         * \param[in] advancingEnd - Should the solver advance to the tbase.
         * \param y - The initial values and the solutions.
         */
        virtual void solve(const real t0, const real tn, real &tReached, const real &dt,
                           const integerListList &Gml, const integer im, const bool advancingEnd,
                           realArray &y);

        /*!\brief Solve MTS problem.
         * \param[in] t0 - The initial time.
         * \param[in] tn - The base physical time, and the time the solver should reach.
         * \param[out] tReached - The time reached by solver.
         * \param[in] dt - The character time scale array.
         * \param[in] Gml - The groups of different time scale.
         * \param y - The initial values and the solutions.
         */
        virtual void solveList(const real t0, const real tn, real &tReached, const realArray &dt,
                               const integerListList &Gml, realArray &y);

        /*!\brief Get error estimate.
         * \param[in] t0 - The initial time.
         * \param[in] tn - The base physical time, and the time the solver should reach.
         * \param[out] tReached - The time reached by solver.
         * \param[in] dt - The character time scale array.
         * \param[in] Gml - The groups of different time scale.
         * \param y - The initial values and the solutions.
         */
        virtual real getError(const real t0, const real tn, real &tReached, const realArray &dt,
                              const integerListList &Gml, realArray &y) {
            return Zero;
        }

        /*!\brief Solve the ODEs by giving the initial timestep
         *        and try to adjust the step according to the specified tolerance.
         * \param t - The time value.
         * \param dt - The initial timestep.
         * \param y - The solutions.
         */
        virtual void adaptiveSolve(real &t, real &dt, realArray &y);

    public:
        /*!\brief Solve the ODEs by giving the initial timestep in a period time
         *        and try to adjust the step according to the specified tolerance.
         * \param[in] t0 - The initial time value.
         * \param t0 - The end time value.
         * \param dt0 - The initial timestep.
         * \param y - The solutions.
         */
        virtual void solve(const real t0, real &tn, real &dt0, realArray &y);

        /*!\brief Get error estimate.
         * \param[in] t0 - The initial time.
         * \param[in] tn - The base physical time, and the time the solver should reach.
         * \param[out] tReached - The time reached by solver.
         * \param[in] dt - The character time scale array.
         * \param[in] Gml - The groups of different time scale.
         * \param y - The initial values and the solutions.
         */
        virtual void solve(const real t0, const real tn, real &tReached, const realArray &dt,
                           const integerListList &Gml, realArray &y);
    };

} // namespace OpenHurricane

#include "ODEsSolver.inl"
