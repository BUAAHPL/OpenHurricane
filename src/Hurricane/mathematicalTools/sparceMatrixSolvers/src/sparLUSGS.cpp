﻿/*!
 * \file LUSGS.cpp
 * \brief Main subroutines for sparse matrices LUSGS Smoother for \f${\bf{A}}x = b\f$
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
#include "sparLUSGS.hpp"
namespace OpenHurricane {
    template <class Type> const std::string sparseSolver::LUSGS<Type>::className_("LUSGS");
    createClassNameStrTmpl(sparseSolver::realLUSGSSolver, "LUSGS");
    createClassNameStrTmpl(sparseSolver::realBlockLUSGSSolver, "LUSGS");
    namespace sparseSolver {
        registerObjFty(realSparSolver, realLUSGSSolver, controller);
        registerObjFty(realBlockSparSolver, realBlockLUSGSSolver, controller);

        registerObjFty(realSparSolver, realLUSGSSolver, Smoother);
        registerObjFty(realBlockSparSolver, realBlockLUSGSSolver, Smoother);
    } // namespace sparseSolver
} // namespace OpenHurricane
