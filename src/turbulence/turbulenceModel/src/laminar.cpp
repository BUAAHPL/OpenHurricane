/*!
 * \file laminar.cpp
 * \brief Main subroutines for the laminar.
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

#include "laminar.hpp"
#include "ArrayInclude.hpp"

#include "fVArraysInclude.hpp"

 namespace OpenHurricane{
	createClassNameStr(laminar,"laminar");
}

namespace OpenHurricane {
 registerObjFty(turbulenceModel,laminar,controller);
}

OpenHurricane::laminar::laminar(const controller &cont, flowModel &ev)
    : turbulenceModel(cont, ev) {}

void OpenHurricane::laminar::expSource() {}

void OpenHurricane::laminar::impSource() {}

void OpenHurricane::laminar::fullImpSource(cellRealSquareMatrixArray &Jac,
                                           const integer rhoId,
                                           const integer rhoTurb0) {}

void OpenHurricane::laminar::visFlux(const faceRealArray &rhof,
                                     const faceRealArray &mulf,
                                     const faceRealArray &mutf,
                                     const cellRealArray &mul,
                                     const cellRealArray &mut,
                                     const cellRealArray &rho) {}

void OpenHurricane::laminar::update() {}

void OpenHurricane::laminar::limit() {}

OpenHurricane::realArray OpenHurricane::laminar::k() const {
    return realArray();
}

OpenHurricane::realArray OpenHurricane::laminar::epsilon() const {
    return realArray();
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::laminar::Ret() const {
    return realArray();
}

OpenHurricane::cellRealArray &OpenHurricane::laminar::var(const integer i) {
    return const_cast<cellRealArray &>(NullRefObj::nullRef<cellRealArray>());
}

const OpenHurricane::cellRealArray &
OpenHurricane::laminar::var(const integer i) const {
    return NullRefObj::nullRef<cellRealArray>();
}

OpenHurricane::faceSymmTensorArray OpenHurricane::laminar::tauEff(
    const faceRealArray &rhof, const faceRealArray &mulf,
    const faceRealArray &mutf, const faceTensorArray &deltafV) const {
    faceSymmTensorArray tau(
        object("tau", mesh(), object::NOT_WRITE, object::TEMPORARY), mesh(),
        //(muf + mutf) * (twoSymm(deltafV) - real(2.0 / 3.0) * diag(deltafV))
        (mutf) * (twoSymm(deltafV) -
                  (real(2.0 / 3.0) * (div(diagToVector(deltafV))) * I)));
    return tau;
}
