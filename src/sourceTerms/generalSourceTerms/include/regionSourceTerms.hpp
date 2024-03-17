﻿/*!
 * \file regionSourceTerms.hpp
 * \brief Headers of base class of region source terms.
 *        The subroutines and functions are in the <i>regionSourceTerms.cpp</i> file.
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

#include "sourceTerm.hpp"

namespace OpenHurricane {

    /*!\brief The class of regionSourceTerms.*/
   class regionSourceTerms : public sourceTerm {
    protected:
        integerArrayArray regionList_;

    public:
        declareClassName(regionSourceTerms);

        // Constructors

        regionSourceTerms() = delete;

        regionSourceTerms(const flowModel &flows, const iteration &iter, const controller &cont);

        regionSourceTerms(const regionSourceTerms &st) = delete;

        /*!\brief Destructor.*/
        virtual ~regionSourceTerms() noexcept {}

        hur_nodiscard inline const integerArrayArray &regionList() const noexcept;

        void operator=(const regionSourceTerms &st) = delete;
    };
} // namespace OpenHurricane

#include "regionSourceTerms.inl"