/*!
 * \file JANAFThermoParsing.cpp
 * \brief Main subroutines for the JANAF of thermo file parsing.
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

#include "JANAFThermoParsing.hpp"
#include "JANAF.hpp"
#include "atomicWeights.hpp"
#include "speciesElements.hpp"

void OpenHurricane::JANAFThermoParsing::parseAllThermoData(thermoList &th) {
    stdStringList therML = getThermoSection(therFileStringLine_);

    const std::regex emptyT("\\s*");
    const std::regex globalTRegex("\\s*([0-9\\.]+)\\s+([0-9\\.]+)\\s+([0-9\\.]+)\\s*");

    for (integer i = 0; i < therML.size(); ++i) {
        if (isSectionMatchedNASA(therML, i)) {
            parseNASASection(therML[i], therML[i + 1], therML[i + 2], therML[i + 3], th);
            i += 3;
        } else if (!std::regex_match(therML[i], emptyT) &&
                   !std::regex_match(therML[i], globalTRegex)) {
            errorAbortStr(("Unmatched: " + therML[i]));
        }
    }
}

void OpenHurricane::JANAFThermoParsing::ensureNoMissingSpecies() const {
    for (integer i = 0; i < speciesFlag_.size(); ++i) {
        if (speciesFlag_[i] == 0) {
            errorAbortStr(
                ("Missingg species: \"" + spt_[i].name() + "\" in thermo file, please check!"));
        }
    }
}

void OpenHurricane::JANAFThermoParsing::checkDuplicate() const {
    bool hasDup = false;
    for (integer i = 0; i < speciesFlag_.size(); ++i) {
        if (speciesFlag_[i] > 1) {
            hasDup = true;
            break;
        }
    }

    if (hasDup) {
        Pout << "         Please check! The following species has more than "
                "one set of thermo data:"
             << std::endl;
        for (integer i = 0; i < speciesFlag_.size(); ++i) {
            if (speciesFlag_[i] > 1) {
                Pout << "           " << speciesNameList_[i] << std::endl;
            }
        }
        Pout << "         Only the first data set will be used." << std::endl;
    }
}

void OpenHurricane::JANAFThermoParsing::ensureSpeciesNamesAreValid() const {
    std::regex space("\\s");

    for (integer i = 0; i < speciesNameList_.size(); ++i) {
        if (std::regex_search(speciesNameList_[i], space)) {
            errorAbortStr(("Invalid species name \"" + speciesNameList_[i] +
                           "\". Species name must contain no space."));
        }
    }
}

void OpenHurricane::JANAFThermoParsing::getGlobalTemperatures() {
    // The regular expression for parsing global temperatures.
    std::regex globalTRegex(
        "\\s*THER(?:|MO)\\s+(?:|ALL)\\s*([0-9\\.]+)\\s+([0-9\\.]+)\\s+([0-9\\.]+)\\s*",
        std::regex::icase);

    std::smatch what;
    bool found = false;

    if (std::regex_search(thermoFileString_, what, globalTRegex)) {
        if (!what.empty()) {
            readReal(what[1].str(), globalLowT_);
            readReal(what[2].str(), globalCommonT_);
            readReal(what[3].str(), globalHighT_);
            found = true;
        }
    }
    if (!found) {
        LFatal("Can not find list of global temperatures.");
    }
}

OpenHurricane::stdStringList
OpenHurricane::JANAFThermoParsing::getThermoSection(const stdStringList &stl) const {
    integer beginL = -1;
    integer endL = -1;
    const std::regex thermoTag("\\s*THER(?:|MO).*", std::regex::icase);
    const std::regex endTag("\\s*END.*", std::regex::icase);

    for (integer i = 0; i < stl.size(); ++i) {
        if ((beginL < 0) && std::regex_match(stringToUpperCase(stl[i]), thermoTag)) {
            beginL = i;
        }

        else if ((endL < 0) && std::regex_match(stringToUpperCase(stl[i]), endTag)) {
            endL = i;
            break;
        }
    }
    if (beginL <= -1 && endL <= -1) {
        endL = stl.size();
    } else if (beginL <= -1) {
#ifdef HUR_DEBUG
        checkWarning("Thermo file does not contain valid starting or ending tags");
#endif // HUR_DEBUG
    } else if (beginL >= endL) {
#ifdef HUR_DEBUG
        checkWarning("Thermo file does not contain valid starting or ending tags");
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
    for (integer i = 0; i < therL.size(); ++i) {
        auto exclPos = therL[i].find_first_of('!');
        if (exclPos != std::string::npos) {
            therL[i].erase(exclPos);
        }
    }

    return therL;
}

bool OpenHurricane::JANAFThermoParsing::isSectionMatchedNASA(const stdStringList &stl,
                                                         const integer offset) const {
    const std::regex janafLine1("[^!]{79}1.*");
    const std::regex janafLine2("[^!]{79}2.*");
    const std::regex janafLine3("[^!]{79}3.*");
    const std::regex janafLine4("[^!]{79}4.*");
    if (stl.size() >= (offset + 4)) {
        return (std::regex_match(stl[offset], janafLine1) &&
                std::regex_match(stl[offset + 1], janafLine2) &&
                std::regex_match(stl[offset + 2], janafLine3) &&
                std::regex_match(stl[offset + 3], janafLine4));
    } else {
        return false;
    }
}

OpenHurricane::string
OpenHurricane::JANAFThermoParsing::extractSpeciesName(const std::string &speciesString) const {
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

void OpenHurricane::JANAFThermoParsing::parseNASASection(const std::string &l1, const std::string &l2,
                                                     const std::string &l3, const std::string &l4,
                                                     thermoList &th) {
    string speciesName = extractSpeciesName(l1.substr(0, 18));

    for (integer i = 0; i < speciesNameList_.size(); ++i) {
        if (speciesNameList_[i] == speciesName) {
            speciesFlag_[i]++;
            if (speciesFlag_[i] > 1) {
                break;
            }

            std::string note(trimCopy(l1.substr(18, 6)));
            std::string phase(l1.substr(44, 1));

            real speciesLowT;
            readReal(l1.substr(45, 10), speciesLowT);
            real speciesHighT;
            readReal(l1.substr(55, 10), speciesHighT);

            real speciesCommT;
            if (trimCopy(l1.substr(65, 8)) == " ") {
                speciesCommT = globalCommonT_;
            } else {
                readReal(l1.substr(65, 8), speciesCommT);
            }
            std::string elementsString = l1.substr(24, 20);

            spt_[i].elementList() = parsingElements(elementsString);
            spt_[i].getMolWeight();

            typename JANAF::coeffArray upperTempCoeff(Zero);
            typename JANAF::coeffArray lowerTempCoeff(Zero);

            readReal(l2.substr(0, 15), upperTempCoeff[0]);
            readReal(l2.substr(15, 15), upperTempCoeff[1]);
            readReal(l2.substr(30, 15), upperTempCoeff[2]);
            readReal(l2.substr(45, 15), upperTempCoeff[3]);
            readReal(l2.substr(60, 15), upperTempCoeff[4]);
            readReal(l3.substr(0, 15), upperTempCoeff[5]);
            readReal(l3.substr(15, 15), upperTempCoeff[6]);
            readReal(l3.substr(30, 15), lowerTempCoeff[0]);
            readReal(l3.substr(45, 15), lowerTempCoeff[1]);
            readReal(l3.substr(60, 15), lowerTempCoeff[2]);
            readReal(l4.substr(0, 15), lowerTempCoeff[3]);
            readReal(l4.substr(15, 15), lowerTempCoeff[4]);
            readReal(l4.substr(30, 15), lowerTempCoeff[5]);
            readReal(l4.substr(45, 15), lowerTempCoeff[6]);

            auto thPtrOld =
                th.thTable().set(i, new JANAF(th.eos(), i, speciesLowT, speciesHighT, speciesCommT,
                                              upperTempCoeff, lowerTempCoeff, phase, true));
            HurDelete(thPtrOld);

            break;
        }
    }
}
hur_nodiscard OpenHurricane::speciesElementsList
OpenHurricane::JANAFThermoParsing::parsingElements(const std::string &eleStr) const {
    speciesElementsList eleNL;

    if (eleStr.length() % 5 == 0) {
        for (integer i = 0; i < (integer)eleStr.length() / 5; ++i) {
            std::string elem = trimCopy(eleStr.substr(i * 5, 2));
            std::string countStr = trimCopy(eleStr.substr(2 + i * 5, 3));

            if (countStr.length() == 0)
                continue;

            // Check if count string is a number
            integer count = Zero;

            readInteger(countStr, count);

            if (count == 0)
                continue;

            if (!checkInElementsList(elem)) {
                errorAbortStr(("Element: \"" + elem +
                               "\" is not in elements list: " + toString(spt_.elementsNameList())));
            }
            speciesElements temElem(elem, count);

            auto isFound = isotopeAtomicWts_.find(elem);
            if (isFound != isotopeAtomicWts_.end()) {
                temElem.setAtomicWeight(isFound->second);
            } else {
                temElem.setAtomicWeight();
            }
            eleNL.append(temElem);
        }
    } else {
        errorAbortStr(("Invalid element string found for value: " + eleStr));
    }
    return eleNL;
}

bool OpenHurricane::JANAFThermoParsing::checkInElementsList(const std::string &eleS) const {
    const auto &eleList = spt_.elementsNameList();
    for (integer i = 0; i < eleList.size(); ++i) {
        if (eleS == eleList[i])
            return true;
    }
    return false;
}

OpenHurricane::JANAFThermoParsing::JANAFThermoParsing(const std::string &thermoFile,
                                                  const stringList &speciesNameList,
                                                  speciesList &spt, thermoList &thermo,
                                                  std::map<std::string, real> &isotopeAtomicWts)
    : thermoFileString_(thermoFile), therFileStringLine_(stringToLineList(thermoFile)),
      thermo_(thermo), globalLowT_(0.0), globalCommonT_(0.0), globalHighT_(0.0),
      speciesNameList_(speciesNameList), spt_(spt), speciesFlag_(speciesNameList.size(), Zero),
      isotopeAtomicWts_(isotopeAtomicWts) {
    if (speciesNameList_.size() > spt_.size()) {
        spt_.resize(speciesNameList_.size());
    }
}

void OpenHurricane::JANAFThermoParsing::parsing() {
    parsing(thermo_);
}

void OpenHurricane::JANAFThermoParsing::parsing(thermoList &th) {
    if (th.thTable().size() < spt_.size()) {
        th.thTable().resize(spt_.size());
    }
    getGlobalTemperatures();
    th.setGlobalTemperature(globalLowT_, globalCommonT_, globalHighT_);
#ifdef HUR_DEBUG
    Pout("    Info: Parsing JANAF thermo file...\n");
    if (report) {
        LInfo("Parsing NASA thermo files: %s", thermoFileString_.c_str());
    }
#endif // HUR_DEBUG
    parseAllThermoData(th);
    ensureNoMissingSpecies();
    checkDuplicate();
    ensureSpeciesNamesAreValid();
}
