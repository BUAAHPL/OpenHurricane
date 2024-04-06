/*!
 * \file JANAFTransportParsing.hpp
 * \brief Header of JANAF of transport file parsing.
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
#pragma once

#include "dataStructure.hpp"
#include "kineticTheory.hpp"
#include "transportList.hpp"

namespace OpenHurricane {

    class JANAFTransportParsing {
    private:
        // Private data

        std::string transportFileString_;

        stdStringList transportFileList_;

        transportList &transport_;

        integerList speciesFlag_;

        const stringList &speciesNameList_;

        string extractSpeciesName(const std::string &speciesString) const;

        void ensureAllSpeciesFound() const;

        void ensureNoDuplicates() const;

        hur_nodiscard stdStringList getTransportSection(const stdStringList &stl) const;

    public:
        JANAFTransportParsing(const std::string &thermoFile, const stringList &speciesNameList,
                              transportList &transport);

        ~JANAFTransportParsing() noexcept {}

        void parsing(const real Prl);

    private:
        void parsing(transportList &trT, const real Prl);
    };

} // namespace OpenHurricane