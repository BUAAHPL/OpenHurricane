/*!
 * \file speciesList.inl
 * \brief The In-Line functions of the <i>speciesList.hpp</i> file.
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

inline OpenHurricane::speciesList::speciesList() : Base(), elementsNameList_() {}

inline OpenHurricane::speciesList::speciesList(const integer _size)
    : Base(_size), elementsNameList_() {}

inline OpenHurricane::speciesList::speciesList(const speciesList &sT)
    : Base(sT), elementsNameList_(sT.elementsNameList_) {}

inline OpenHurricane::speciesList::speciesList(speciesList &&sT) noexcept
    : Base(std::move(sT)), elementsNameList_(std::move(sT.elementsNameList_)) {}

hur_nodiscard inline const OpenHurricane::string &
OpenHurricane::speciesList::name(const integer specieI) const noexcept {
    return this->operator[](specieI).name();
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::speciesList::W(const integer specieI) const noexcept {
    return this->operator[](specieI).W();
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::speciesList::Ri(const integer specieI) const noexcept {
    return this->operator[](specieI).Ri();
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::speciesList::Rm(const realArray &yi) const noexcept {
    return constant::physicalConstant::Ru / MWbyYi(yi);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::speciesList::Z1(const real s, const real YF,
                                                                const real YO2, const real YF1,
                                                                const real YO22) const noexcept {
    return (s * YF - YO2 + YO22) / max(s * YF1 + YO22, veryTiny);
}