/*!
 * \file sparSymGaussSeidel.cpp
 * \brief Main subroutines for sparse matrices symmetric Gauss-Seidel Smoother for \f${\bf{A}}x = b\f$
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
#include "sparSymGaussSeidel.hpp"
template <class Type>
const std::string OpenHurricane::sparseSolver::symGaussSeidel<Type>::className_ = "symGaussSeidel";

namespace OpenHurricane {
    createClassNameStrTmpl(sparseSolver::realSymGSSolver, "symGaussSeidel");
    createClassNameStrTmpl(sparseSolver::realBlockSymGSSolver, "symGaussSeidel");

    namespace sparseSolver {
        registerObjFty(realSparSolver, realSymGSSolver, controller);
        registerObjFty(realBlockSparSolver, realBlockSymGSSolver, controller);

        registerObjFty(realSparSolver, realSymGSSolver, Smoother);
        registerObjFty(realBlockSparSolver, realBlockSymGSSolver, Smoother);
    } // namespace sparseSolver
} // namespace OpenHurricane
