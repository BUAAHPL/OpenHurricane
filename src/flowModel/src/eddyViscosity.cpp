/*!
 * \file eddyViscosity.cpp
 * \brief Main subroutines for the eddy-viscosity hypothesis flow model.
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

#include "eddyViscosity.hpp"

namespace OpenHurricane {
    createClassName(eddyViscosity);
    registerObjFty(flowModel, eddyViscosity, controller);
} // namespace OpenHurricane

OpenHurricane::eddyViscosity::eddyViscosity(const runtimeMesh &mesh)
    : flowModel(mesh), Sct_(0.5), Prt_(0.9),
      mut_(object("mut", mesh, object::WRITE_RELAY_OUTPUT), mesh), mutLow_(1.0e-15),
      mutHigh_(2.0e5) {
    thermoPtr_ = new viscousThermo(mesh);
    thermoPtr_->mixtures().transTable().setPr(Prl());
    thermoPtr_->mixtures().thermalTable().setLimitTemperature(TLow_, THigh_);
}

OpenHurricane::eddyViscosity::eddyViscosity(const runtimeMesh &mesh, const controller &cont)
    : flowModel(mesh, cont), Sct_(cont.findOrDefault("Sct", 0.5)),
      Prt_(cont.findOrDefault("Prt", 0.9)),
      mut_(object("mut", mesh, object::WRITE_RELAY_OUTPUT), mesh), mutLow_(1.0e-15),
      mutHigh_(2.0e5) {
    if (cont.found("limits")) {
        auto &limitsCont = cont.subController("limits");
        mutLow_ = limitsCont.findOrDefault<real>("mutLow", mutLow_);
        mutHigh_ = limitsCont.findOrDefault<real>("mutHigh", mutHigh_);
    } else {
        mutLow_ = cont.findOrDefault<real>("mutLow", mutLow_);
        mutHigh_ = cont.findOrDefault<real>("mutHigh", mutHigh_);
    }

    thermoPtr_ = new viscousThermo(mesh, cont);

    thermoPtr_->mixtures().transTable().setPr(Prl());

    thermoPtr_->mixtures().thermalTable().setLimitTemperature(TLow_, THigh_);
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::eddyViscosity::mul() noexcept {
    return thermoPtr_->mu();
}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::eddyViscosity::mul() const noexcept {
    return thermoPtr_->mu();
}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::eddyViscosity::mut() const noexcept {
    return mut_;
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::eddyViscosity::muEff() {
    return thermoPtr_->mu() + mut();
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::eddyViscosity::kappal() noexcept {
    return thermoPtr_->kappa();
}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::eddyViscosity::kappal() const noexcept {
    return thermoPtr_->kappa();
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::eddyViscosity::kappat() noexcept {
    return thermoPtr_->kappa();
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::eddyViscosity::kappaEff() {
    return thermoPtr_->kappa();
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::eddyViscosity::keEff() const {
    cellRealArray tmpKe(object("tmpKe", mesh(), object::NOT_WRITE), mesh(),
                        thermoPtr_->mu() / Prl_ + mut_ / Prt_);
    return tmpKe;
}

hur_nodiscard OpenHurricane::real OpenHurricane::eddyViscosity::keEff(const integer cellI) const {
    return thermoPtr_->mu()[cellI] / Prl_ + mut_[cellI] / Prt_;
}
