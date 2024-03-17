/*!
 * \file flowModel.cpp
 * \brief Main subroutines for the flow model.
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

#include "flowModel.hpp"
namespace OpenHurricane {
    createClassName(flowModel);
    createObjFty(flowModel, controller);
} // namespace OpenHurricane

OpenHurricane::flowModel::flowModel(const runtimeMesh &mesh)
    : mesh_(mesh), thermoPtr_(nullptr), v_(object("v", string("\"u\",\"v\",\"w\""), mesh,
                                                  object::WRITE_RELAY_OUTPUT, object::PRIMITIVE),
                                           mesh),
      shockFactor_(object("shock_factor", mesh, object::NOT_WRITE), mesh, real(1)), Prl_(0.72),
      pLow_(0.0), pHigh_(1.0e10), TLow_(0.0), THigh_(5000.0),
      refValues_(mesh.Iteration().refValues()), temperatureFlag_(mesh.nCells(), 1),
      pressureFlag_(mesh.nCells(), 1),
      CFLFlag_(object("CFLFlag", mesh, object::NOT_WRITE), mesh, 1), freeStrPtr_(nullptr) {}

OpenHurricane::flowModel::flowModel(const runtimeMesh &mesh, const controller &cont)
    : mesh_(mesh), thermoPtr_(nullptr), v_(object("v", string("\"u\",\"v\",\"w\""), mesh,
                                                  object::WRITE_RELAY_OUTPUT, object::PRIMITIVE),
                                           mesh),
      shockFactor_(object("shock_factor", mesh, object::NOT_WRITE), mesh, real(1)),
      Prl_(cont.findOrDefault<real>("Prl", 0.72)), pLow_(0.0), pHigh_(1.0e10), TLow_(0.0),
      THigh_(5000.0), refValues_(mesh.Iteration().refValues()), temperatureFlag_(mesh.nCells(), 1),
      pressureFlag_(mesh.nCells(), 1),
      CFLFlag_(object("CFLFlag", mesh, object::NOT_WRITE), mesh, 1), freeStrPtr_(nullptr) {

    if (cont.found("limits")) {
        auto &limitsCont = cont.subController("limits");
        pLow_ = limitsCont.findOrDefault<real>("pLow", pLow_);
        pHigh_ = limitsCont.findOrDefault<real>("pHigh", pHigh_);
        TLow_ = limitsCont.findOrDefault<real>("TLow", TLow_);
        THigh_ = limitsCont.findOrDefault<real>("THigh", THigh_);
    } else {
        pLow_ = cont.findOrDefault<real>("pLow", pLow_);
        pHigh_ = cont.findOrDefault<real>("pHigh", pHigh_);
        TLow_ = cont.findOrDefault<real>("TLow", TLow_);
        THigh_ = cont.findOrDefault<real>("THigh", THigh_);
    }
}

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::flowModel>
OpenHurricane::flowModel::creator(const runtimeMesh &mesh, const controller &cont) {
    string flowType = cont.findWord(flowModel::className_);
    Pout << "    Info: setting flow model: " << flowType << std::endl;
    defineInObjCreator(flowModel, flowType, controller, (mesh, cont));
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::flowModel::mul() noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::flowModel::mul() const noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::flowModel::mul(const integer faceZoneId) {
    if (mul().size() == 0) {
        return realArray();
    }
    return fv::interpolate(mul(), faceZoneId);
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::flowModel::mul(const integer faceZoneId) const {
    if (mul().size() == 0) {
        return realArray();
    }
    return fv::interpolate(mul(), faceZoneId);
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::flowModel::nu() {
    if (mul().size() == 0) {
        return const_cast<cellRealArray &>(cellRealArray::nullObject());
    }

    cellRealArray nul(object("nul", mesh_, object::NOT_WRITE), mesh_, mul() / rho());
    return nul;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::flowModel::nu(const integer faceZoneId) {
    if (mul().size() == 0) {
        return realArray();
    }
    return fv::interpolate(mul(), faceZoneId) / fv::interpolate(rho(), faceZoneId);
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::flowModel::muEff() {
    return mul();
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::flowModel::mut() noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::flowModel::mut() const noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::flowModel::kappal() noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard const OpenHurricane::cellRealArray &OpenHurricane::flowModel::kappal() const noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::flowModel::kappaEff() {
    return kappal();
}

hur_nodiscard OpenHurricane::cellRealArray &OpenHurricane::flowModel::kappat() noexcept {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::flowModel::keEff() const {
    return const_cast<cellRealArray &>(cellRealArray::nullObject());
}

hur_nodiscard OpenHurricane::real OpenHurricane::flowModel::keEff(const integer cellI) const {
    return Zero;
}

hur_nodiscard const OpenHurricane::freeStream &OpenHurricane::flowModel::freeStr() const {
    if (freeStrPtr_ == nullptr) {
        freeStrPtr_ = new freeStream(mesh_.Iteration().cont(), *this);
    }
    return *freeStrPtr_;
}
