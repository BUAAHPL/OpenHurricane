/*!
 * \file viscousThermo.cpp
 * \brief Main subroutines for the viscous thermo.
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
#include "viscousThermo.hpp"

OpenHurricane::viscousThermo::viscousThermo(const runtimeMesh &mesh)
    : rhoThermo(mesh), mu_(object("mu", mesh, object::WRITE_OUTPUT), mesh),
      kappa_(object("kappa", mesh, object::WRITE_OUTPUT), mesh) {}

OpenHurricane::viscousThermo::viscousThermo(const runtimeMesh &mesh, const controller &cont)
    : rhoThermo(mesh, cont), mu_(object("mu", mesh, object::WRITE_OUTPUT), mesh),
      kappa_(object("kappa", mesh, object::WRITE_OUTPUT), mesh) {}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::viscousThermo::mu() const noexcept {
    return mu_;
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::viscousThermo::mu() noexcept {
    return mu_;
}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::viscousThermo::kappa() const noexcept {
    return kappa_;
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::viscousThermo::kappa() noexcept {
    return kappa_;
}
