/*!
 * \file laminarFlow.cpp
 * \brief Main subroutines for the laminar flow model.
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

#include "laminarFlow.hpp"
namespace OpenHurricane {
    createClassName(laminarFlow);
         registerObjFty(flowModel, laminarFlow, controller);
} // namespace OpenHurricane

OpenHurricane::laminarFlow::laminarFlow(const runtimeMesh &mesh) : flowModel(mesh) {
    thermoPtr_ = new viscousThermo(mesh);

    thermoPtr_->mixtures().transTable().setPr(Prl());
    thermoPtr_->mixtures().thermalTable().setLimitTemperature(TLow_, THigh_);
}

OpenHurricane::laminarFlow::laminarFlow(const runtimeMesh &mesh, const controller &cont)
    : flowModel(mesh, cont) {
    thermoPtr_ = new viscousThermo(mesh, cont);
    thermoPtr_->mixtures().transTable().setPr(Prl());

    thermoPtr_->mixtures().thermalTable().setLimitTemperature(TLow_, THigh_);
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::laminarFlow::mul() noexcept {
    return thermoPtr_->mu();
}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::laminarFlow::mul() const noexcept {
    return thermoPtr_->mu();
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::laminarFlow::muEff() {
    return thermoPtr_->mu();
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::laminarFlow::kappal() noexcept {
    return thermoPtr_->kappa();
}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::laminarFlow::kappal() const noexcept {
    return thermoPtr_->kappa();
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::laminarFlow::kappaEff() {
    return thermoPtr_->kappa();
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::laminarFlow::keEff() const {
    cellRealArray tmpKe(object(string("tmpKe"), mesh(), object::NOT_WRITE), mesh(),
                        thermoPtr_->mu() / Prl_);
    return tmpKe;
}

hur_nodiscard OpenHurricane::real OpenHurricane::laminarFlow::keEff(const integer cellI) const {
    return thermoPtr_->mu()[cellI] / Prl_;
}
