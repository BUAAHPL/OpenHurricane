/*!
 * \file JANAFThermoParsing.hpp
 * \brief Header of JANAF of thermo file parsing.
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

#include "dataStructure.hpp"
#include "thermo.hpp"
#include "thermoList.hpp"
namespace OpenHurricane {

    class JANAFThermoParsing {
    private:
        // Private data

        std::string thermoFileString_;

        stdStringList therFileStringLine_;

        thermoList &thermo_;

        real globalLowT_;
        real globalCommonT_;
        real globalHighT_;

        const stringList &speciesNameList_;

        speciesList &spt_;

        integerList speciesFlag_;

        /**\brief Isotope molecular weights.*/
        std::map<std::string, real> &isotopeAtomicWts_;

    private:
        // Private functions

        void parseAllThermoData(thermoList &th);

        void ensureNoMissingSpecies() const;
        void checkDuplicate() const;

        void ensureSpeciesNamesAreValid() const;

        void getGlobalTemperatures();

        stdStringList getThermoSection(const stdStringList &stl) const;

        bool isSectionMatchedNASA(const stdStringList &stl, const integer offset) const;

        string extractSpeciesName(const std::string &speciesString) const;

        void parseNASASection(const std::string &l1, const std::string &l2, const std::string &l3,
                              const std::string &l4, thermoList &th);

        hur_nodiscard speciesElementsList parsingElements(const std::string &eleStr) const;

        bool checkInElementsList(const std::string &eleS) const;

    public:

        /**\brief Construct from components.*/
        JANAFThermoParsing(const std::string &thermoFile, const stringList &speciesNameList,
                           speciesList &spt, thermoList &thermo,
                           std::map<std::string, real> &isotopeAtomicWts);

        /**\brief Destructor.*/
        ~JANAFThermoParsing() noexcept {}

        void parsing();

    private:
        void parsing(thermoList &th);
    };

} // namespace OpenHurricane