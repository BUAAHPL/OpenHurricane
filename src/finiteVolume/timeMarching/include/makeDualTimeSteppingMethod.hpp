/*!
 * \file makeDualTimeSteppingMethod.hpp
 * \brief Headers of class of the making dual time stepping.
 *        The subroutines and functions are in the <i>makeDualTimeSteppingMethod.cpp</i> file.
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

#include "BDF123.hpp"
#include "ESDIRK.hpp"
#include "LUSGS.hpp"

namespace OpenHurricane {

    // Dual time-stepping with LU-SGS

    using BDF123LUSGS = BDF123<LUSGS>;
    using ESDIRKLUSGS = ESDIRK<LUSGS>;

    template <> void OpenHurricane::BDF123<LUSGS>::setSolverWrite();

    template <> void OpenHurricane::BDF123<LUSGS>::setBDFType(const controller &timeMethodCont);

    template <> void OpenHurricane::ESDIRK<LUSGS>::setSolverWrite();

    template <> void OpenHurricane::ESDIRK<LUSGS>::setESDIRKType(const controller &timeMethodCont);

} // namespace OpenHurricane