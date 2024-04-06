/*!
 * \file sparHybridLUSGS.cpp
 * \brief Main subroutines for spar Hybrid LUSGS smoother
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
#include "sparHybridLUSGS.hpp"

template <class Type>
const std::string OpenHurricane::sparseSolver::HybridLUSGS<Type>::className_("HybridLUSGS");

namespace OpenHurricane {
    createClassNameStrTmpl(sparseSolver::realHybridLUSGSSolver, "HybridLUSGS");
    createClassNameStrTmpl(sparseSolver::realBlockHybridLUSGSSolver, "HybridLUSGS");

    namespace sparseSolver {
        registerObjFty(realSparSolver, realHybridLUSGSSolver, controller);
        registerObjFty(realBlockSparSolver, realBlockHybridLUSGSSolver, controller);
        registerObjFty(realSparSolver, realHybridLUSGSSolver, Smoother);
        registerObjFty(realBlockSparSolver, realBlockHybridLUSGSSolver, Smoother);
    } // namespace sparseSolver
} // namespace OpenHurricane