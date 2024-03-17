#include "combustionModel.hpp"
/*!
 * \file combustionModel.inl
 * \brief The In-Line functions of the <i>combustionModel.hpp</i> file.
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

hur_nodiscard inline OpenHurricane::realArray &OpenHurricane::combustionModel::tauInt() noexcept {
    if (tauIntPtr_ == nullptr) {
        tauIntPtr_ = new realArray(flows_.mesh().nTotalCells(), Zero);
    }
    return *tauIntPtr_;
}

hur_nodiscard inline OpenHurricane::reactionList &OpenHurricane::combustionModel::reactions() {
    return flows_.mixtures().reactions();
}

hur_nodiscard inline const OpenHurricane::reactionList &OpenHurricane::combustionModel::reactions() const {
    return flows_.mixtures().reactions();
}

inline OpenHurricane::chemistrySource &OpenHurricane::combustionModel::chemistry() noexcept {
    return *chemistryPtr_;
}

inline OpenHurricane::speciesList &OpenHurricane::combustionModel::species() {
    return flows_.mixtures().species();
}

inline const OpenHurricane::speciesList &OpenHurricane::combustionModel::species() const {
    return flows_.mixtures().species();
}

inline OpenHurricane::mixture &OpenHurricane::combustionModel::mixtures() {
    return flows_.mixtures();
}

hur_nodiscard inline const OpenHurricane::mixture &OpenHurricane::combustionModel::mixtures() const {
    return flows_.mixtures();
}

hur_nodiscard inline const OpenHurricane::stringList &
OpenHurricane::combustionModel::fuelName() const noexcept {
    return fuelName_;
}

hur_nodiscard inline const OpenHurricane::stringList &
OpenHurricane::combustionModel::oxygenName() const noexcept {
    return oxygenName_;
}

hur_nodiscard inline const OpenHurricane::stringList &
OpenHurricane::combustionModel::productionName() const noexcept {
    return productionName_;
}

inline void OpenHurricane::combustionModel::fullPointImpChemistrySource(
    realArray &dt, cellRealSquareMatrixArray &Jac, const integer rhoId, const integer rhouId,
    const integer rhoEId, const integer rhoYi0Id) {
    LFatal("This function is not implemented");
}

hur_nodiscard inline OpenHurricane::realArray OpenHurricane::combustionModel::tcFRR() {
    realArray tc(flows_.mesh().nTotalCells(), Zero);
    tcFRR(tc);
    return tc;
}

hur_nodiscard inline OpenHurricane::realArray OpenHurricane::combustionModel::tcSFR() {
    realArray tc(flows_.mesh().nTotalCells(), Zero);
    tcSFR(tc);
    return tc;
}

hur_nodiscard inline OpenHurricane::realArray OpenHurricane::combustionModel::tcGSPR() {
    realArray tc(flows_.mesh().nTotalCells(), Zero);
    tcGSPR(tc);
    return tc;
}

hur_nodiscard inline OpenHurricane::realArray OpenHurricane::combustionModel::tcJacDT() {
    realArray tc(flows_.mesh().nTotalCells(), Zero);
    tcJacDT(tc);
    return tc;
}

hur_nodiscard inline const OpenHurricane::runtimeMesh &
OpenHurricane::combustionModel::mesh() const noexcept {
    return flows_.mesh();
}

hur_nodiscard inline bool OpenHurricane::combustionModel::isStrangSplitted() const noexcept {
    return chemOptions_ == chemistryOption::strangSplitted;
}

hur_nodiscard inline bool OpenHurricane::combustionModel::isCoupled() const noexcept {
    return chemOptions_ == chemistryOption::coupled;
}

hur_nodiscard inline bool OpenHurricane::combustionModel::isIntegrated() const noexcept {
    return chemOptions_ == chemistryOption::integrated;
}

hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::combustionModel::mixZ() {
    if (mixZPtr_ == nullptr) {
        mixZPtr_ = new cellRealArray(object("mixtureFractionZ", flows_.mesh()), flows_.mesh());
    }
    return *mixZPtr_;
}

hur_nodiscard inline const OpenHurricane::cellRealArray &OpenHurricane::combustionModel::mixZ() const {
    if (mixZPtr_ == nullptr) {
        mixZPtr_ = new cellRealArray(object("mixtureFractionZ", flows_.mesh()), flows_.mesh());
    }
    return *mixZPtr_;
}
