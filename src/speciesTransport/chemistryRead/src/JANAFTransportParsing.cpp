/*!
 * \file JANAFTransportParsing.cpp
 * \brief Main subroutines for the JANAF of transport file parsing.
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

#include "JANAFTransportParsing.hpp"
#include "transport.hpp"

#include <sstream>

OpenHurricane::string
OpenHurricane::JANAFTransportParsing::extractSpeciesName(const std::string &speciesString) const {
    std::regex speciesNameRegex("\\s*(.*?)\\s+");
    string speciesName;
    std::smatch what;
   
    if (std::regex_search(speciesString, what, speciesNameRegex)) {
        speciesName = what[1];
    } else {
        errorAbortStr(("Can't find a species name. There is probably no space "
                       "after column 18. And the string is: " +
                       speciesString));
    }
    return speciesName;
}

void OpenHurricane::JANAFTransportParsing::ensureAllSpeciesFound() const {
    for (integer i = 0; i < speciesFlag_.size(); ++i) {
        if (speciesFlag_[i] == 0) {
            errorAbortStr(("Missingg species: \"" + speciesNameList_[i] +
                           "\" in transport file, please check!"));
        }
    }
}

void OpenHurricane::JANAFTransportParsing::ensureNoDuplicates() const {
    bool hasDup = false;
    for (integer i = 0; i < speciesFlag_.size(); ++i) {
        if (speciesFlag_[i] > 1) {
            hasDup = true;
            break;
        }
    }
    if (hasDup) {
        Pout << "         Please check! The following species has more than "
                "one set of transport data:"
             << std::endl;
        for (integer i = 0; i < speciesFlag_.size(); ++i) {
            if (speciesFlag_[i] > 1) {
                Pout << "           " << speciesNameList_[i] << std::endl;
            }
        }
        Pout << "         Only the first data set will be used." << std::endl;
    }
}

hur_nodiscard OpenHurricane::stdStringList
OpenHurricane::JANAFTransportParsing::getTransportSection(const stdStringList &stl) const {
    integer beginL = -1;
    integer endL = -1;
    const std::regex tranTag("\\s*TRAN(?:|SPORT).*", std::regex::icase);
    const std::regex endTag("\\s*END.*", std::regex::icase);

    for (integer i = 0; i < stl.size(); ++i) {       
        if ((beginL < 0) && std::regex_match(stringToUpperCase(stl[i]), tranTag)) {
            beginL = i;
        }

        else if ((endL < 0) && std::regex_match(stringToUpperCase(stl[i]), endTag)) {
            endL = i;
            break;
        }
    }
    if (beginL <= -1 && endL <= -1) {
        return stl;
    } else if (beginL <= -1) {
#ifdef HUR_DEBUG
        checkWarning("Transport file does not contain valid starting or ending tags");
#endif // HUR_DEBUG
    } else if (beginL >= endL) {
#ifdef HUR_DEBUG
        checkWarning("Transport file does not contain valid starting or ending tags");
#endif // HUR_DEBUG
        endL = stl.size();
    }
    stdStringList therLOrigin(endL - beginL - 1);

    integer countEffect = Zero;
    for (integer i = ++beginL; i < endL; ++i) {
        if (stl[i].find_first_of('!') == 0) {
            continue;
        }
        integer j = countEffect++;
        therLOrigin[j] = stl[i];
    }

    stdStringList therL(countEffect);
    for (integer i = 0; i < countEffect; ++i) {
        therL[i] = therLOrigin[i];
    }

    return therL;
}

OpenHurricane::JANAFTransportParsing::JANAFTransportParsing(const std::string &transportFile,
                                                        const stringList &speciesNameList,
                                                        transportList &transport)
    : transportFileString_(transportFile), transportFileList_(), transport_(transport),
      speciesFlag_(speciesNameList.size(), Zero), speciesNameList_(speciesNameList) {
    transportFileList_ = getTransportSection(stringToLineList(transportFile));
    for (integer i = 0; i < transportFileList_.size(); ++i) {
        auto exclPos = transportFileList_[i].find_first_of('!');
        if (exclPos != std::string::npos) {
            transportFileList_[i].erase(exclPos);
        }
    }
}

void OpenHurricane::JANAFTransportParsing::parsing(const real Prl) {
    parsing(transport_, Prl);
}

void OpenHurricane::JANAFTransportParsing::parsing(transportList &trT, const real Prl) {
    if (trT.tranTable().size() < speciesNameList_.size()) {
        trT.tranTable().resize(speciesNameList_.size());
    }

    for (integer i = 0; i < transportFileList_.size(); ++i) {

        if (!trimCopy(transportFileList_[i]).empty()) {
            string speciesName = extractSpeciesName(transportFileList_[i].substr(0, 18));
            for (integer j = 0; j < speciesNameList_.size(); ++j) {
                if (speciesNameList_[j] == speciesName) {
                    speciesFlag_[j]++;

                    if (speciesFlag_[j] > 1) {
                        break;
                    }
                    std::stringstream ss(transportFileList_[i].substr(18));
                    integer moleculeIndex;
                    real potentialWellDepth;
                    real collisionDiameter;
                    real dipoleMoment;
                    real polarizability;
                    real rotRelaxtionNumber;

                    auto defpre = ss.precision();

                    ss.precision(feature<real>::precision);

                    ss >> moleculeIndex >> potentialWellDepth >> collisionDiameter >>
                        dipoleMoment >> polarizability >> rotRelaxtionNumber;

                    ss.precision(defpre);

                    auto tranPtr = trT.tranTable().set(
                        j, new kineticTheory(trT.species(), j, Prl, moleculeIndex,
                                             potentialWellDepth, collisionDiameter, dipoleMoment,
                                             polarizability, rotRelaxtionNumber));

                    HurDelete(tranPtr);
                    break;
                }
            }
        }
    }

    ensureAllSpeciesFound();
    ensureNoDuplicates();
}
