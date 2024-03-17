/*!
 * \file ODESolver.hpp
 * \brief Headers of base class of chemistry ODEs solver.
 *        The subroutines and functions are in the <i>ODESolver.cpp</i> file.
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
#include "Euler.hpp"
#include "ODEsSolver.hpp"
#include "chemistryODESolvers.hpp"

namespace OpenHurricane {
    /*!\brief The class of ODESolver*/
    template <class ChemistrySource> class ODESolver : public chemistryODESolvers<ChemistrySource> {
    private:
        uniquePtr<ODEsSolver> odesPtr_;

        realArray cT_;

#ifdef TEST_PROCESS_TIME
    protected:
        virtual inline void setODECountIter() noexcept;
#endif
    public:
        declareClassName(ODEsSolver);

        /*!\brief Disallow null constructor.*/
        ODESolver() = delete;

        /*!\brief Construct from flow and controller.*/
        ODESolver(flowModel &flows, const controller &cont);

        /*!\brief Destructor.*/
        virtual ~ODESolver() noexcept;

        virtual void solve(realArray &yi, real &dT, real &subDT);
    };

} // namespace OpenHurricane

#include "ODESolver.inl"
