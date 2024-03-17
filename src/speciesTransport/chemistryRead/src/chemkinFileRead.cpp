/*!
 * \file chemkinFileRead.cpp
 * \brief Main subroutines for the chemkin file parsing.
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

#include "chemkinFileRead.hpp"
#include "ArrheniusRR.hpp"
#include "ArrheniusThirdBodyRR.hpp"
#include "ChebyshevPolynomialsRR.hpp"
#include "JanevRR.hpp"
#include "LandauTellerRR.hpp"
#include "chemicallyActicatedBimolecularRR.hpp"
#include "collisionRR.hpp"
#include "irreversibleReactions.hpp"
#include "logInterpolationRR.hpp"
#include "powerSeriesRR.hpp"
#include "reversibleReactionWithRevPars.hpp"
#include "reversibleReactions.hpp"
#include "unimolecularFallOffRR.hpp"

hur_nodiscard std::map<std::string, OpenHurricane::chemkinFileRead::unitTypesOfA>
OpenHurricane::chemkinFileRead::unitMapOfA() const {

    std::map<std::string, unitTypesOfA> unitMap;
    unitMap.emplace("MOLE", unitTypesOfA::MOLES);
    unitMap.emplace("MOLES", unitTypesOfA::MOLES);
    unitMap.emplace("MOLEC", unitTypesOfA::MOLECULES);
    unitMap.emplace("MOLECULE", unitTypesOfA::MOLECULES);
    unitMap.emplace("MOLECULES", unitTypesOfA::MOLECULES);
    return unitMap;
}

hur_nodiscard std::map<std::string, OpenHurricane::chemkinFileRead::unitTypesOfE>
OpenHurricane::chemkinFileRead::unitMapOfE() const {
    std::map<std::string, unitTypesOfE> unitMap;
    unitMap.emplace("CAL", unitTypesOfE::CAL);
    unitMap.emplace("CAL/MOLE", unitTypesOfE::CAL);
    unitMap.emplace("EVOL", unitTypesOfE::EVOL);
    unitMap.emplace("EVOLTS", unitTypesOfE::EVOL);
    unitMap.emplace("JOUL", unitTypesOfE::JOUL);
    unitMap.emplace("JOULES/MOLE", unitTypesOfE::JOUL);
    unitMap.emplace("KCAL", unitTypesOfE::KCAL);
    unitMap.emplace("KCAL/MOLE", unitTypesOfE::KCAL);
    unitMap.emplace("KELV", unitTypesOfE::KELV);
    unitMap.emplace("KELVINS", unitTypesOfE::KELV);
    unitMap.emplace("KJOU", unitTypesOfE::KJOU);
    unitMap.emplace("KJOULES/MOLE", unitTypesOfE::KJOU);
    return unitMap;
}

void OpenHurricane::chemkinFileRead::getReacAuxKeywordMap() {
    reactionAuxiliaryKeywordMap_.emplace("NeutralThirdBody",
                                         reactionAuxiliaryKeyword::NeutralThirdBody);
    reactionAuxiliaryKeywordMap_.emplace("CHEB", reactionAuxiliaryKeyword::CHEB);
    reactionAuxiliaryKeywordMap_.emplace("COLLEFF", reactionAuxiliaryKeyword::COLLEFF);
    reactionAuxiliaryKeywordMap_.emplace("DUP", reactionAuxiliaryKeyword::DUP);
    reactionAuxiliaryKeywordMap_.emplace("DUPLICATE", reactionAuxiliaryKeyword::DUP);
    reactionAuxiliaryKeywordMap_.emplace("EXCI", reactionAuxiliaryKeyword::EXCI);
    reactionAuxiliaryKeywordMap_.emplace("FIT1", reactionAuxiliaryKeyword::FIT1);
    reactionAuxiliaryKeywordMap_.emplace("FORD", reactionAuxiliaryKeyword::FORD);
    reactionAuxiliaryKeywordMap_.emplace("HIGH", reactionAuxiliaryKeyword::HIGH);
    reactionAuxiliaryKeywordMap_.emplace("JAN", reactionAuxiliaryKeyword::JAN);
    reactionAuxiliaryKeywordMap_.emplace("LOW", reactionAuxiliaryKeyword::LOW);
    reactionAuxiliaryKeywordMap_.emplace("LT", reactionAuxiliaryKeyword::LT);
    reactionAuxiliaryKeywordMap_.emplace("MOME", reactionAuxiliaryKeyword::MOME);
    reactionAuxiliaryKeywordMap_.emplace("PCHEB", reactionAuxiliaryKeyword::PCHEB);
    reactionAuxiliaryKeywordMap_.emplace("PLOG", reactionAuxiliaryKeyword::PLOG);
    reactionAuxiliaryKeywordMap_.emplace("REV", reactionAuxiliaryKeyword::REV);
    reactionAuxiliaryKeywordMap_.emplace("RLT", reactionAuxiliaryKeyword::RLT);
    reactionAuxiliaryKeywordMap_.emplace("RORD", reactionAuxiliaryKeyword::RORD);
    reactionAuxiliaryKeywordMap_.emplace("SRI", reactionAuxiliaryKeyword::SRI);
    reactionAuxiliaryKeywordMap_.emplace("TCHEB", reactionAuxiliaryKeyword::TCHEB);
    reactionAuxiliaryKeywordMap_.emplace("TDEP", reactionAuxiliaryKeyword::TDEP);
    reactionAuxiliaryKeywordMap_.emplace("TROE", reactionAuxiliaryKeyword::TROE);
    reactionAuxiliaryKeywordMap_.emplace("UNITS", reactionAuxiliaryKeyword::UNITS);
    reactionAuxiliaryKeywordMap_.emplace("XSMI", reactionAuxiliaryKeyword::XSMI);
    reactionAuxiliaryKeywordMap_.emplace("END", reactionAuxiliaryKeyword::END);
}

hur_nodiscard bool OpenHurricane::chemkinFileRead::checkChemFile() const {

#ifdef HUR_DEBUG
    Pout << "    Info: Checking the format of chemkin file." << std::endl;
#endif // HUR_DEBUG
    if (noReaction_) {
        const std::regex fs("ELEM(?:|ENT|ENTS)([\\s\\S]*?)END([\\s\\S]*?)"
                            "SPEC(?:|IE|IES)([\\s\\S]*?)END",
                            std::regex::icase);
        if (std::regex_search(chemkinFileString_, fs)) {
            return true;
        }
    } else {
        const std::regex fs("ELEM(?:|ENT|ENTS)([\\s\\S]*?)END([\\s\\S]*?)"
                            "SPEC(?:|IE|IES)([\\s\\S]*?)END([\\s\\S]*?)"
                            "REAC(?:|TION|TIONS)([\\s\\S]*?)END",
                            std::regex::icase);
        if (std::regex_search(chemkinFileString_, fs)) {
            return true;
        }
    }

    return false;
}

hur_nodiscard std::string
OpenHurricane::chemkinFileRead::replaceComments(const std::string &fileString) const {
    std::string str = fileString;
    replaceAllMarks(str, "\r", " ");
    replaceAllMarks(str, "\t", " ");
    std::regex commentRegex("(!.*?)\\n|(!.*?)$");
    std::string formatString = " \n";
    return std::regex_replace(str, commentRegex, formatString);
}

hur_nodiscard std::string OpenHurricane::chemkinFileRead::readReactionString() const {
    std::regex reactionListRegex("REAC(?:|TION|TIONS)\\s+([\\s\\S]*?)\\s+END", std::regex::icase);
    std::smatch resultMatch;
    std::regex_search(chemkinFileString_, resultMatch, reactionListRegex);
    std::string reacString = resultMatch[0];
    eraseNullBlank(reacString);
    return reacString;
}

void OpenHurricane::chemkinFileRead::parsingReacStoich(const std::string &forwardStr,
                                                   List<reacSpcCoeffType> &forwardCoeff,
                                                   const std::string &backwardStr,
                                                   List<reacSpcCoeffType> &backwardCoeff,
                                                   const std::string &reactionName) const {
    auto readReacSpcFromMap = [&](auto &map, List<reacSpcCoeffType> &coeff) {
        integer sj = 0;
        for (auto iter = map.begin(); iter != map.end(); ++iter) {
            integer index;
            if (!checkReacSpec(iter->first, index)) {
                errorAbortStr(("Species: \"" + iter->first +
                               "\" not found in species list, in reaction:\n" + reactionName));
            }

            coeff[sj].index_ = index;
            coeff[sj].stoichCoeff_ = iter->second;
            coeff[sj].order_ = iter->second;
            sj++;
        }
    };

    // The map of reactants (forward reaction)
    auto forwardMap = parsingReacSpc(forwardStr, reactionName);
    forwardCoeff.resize(static_cast<integer>(forwardMap.size()));
    readReacSpcFromMap(forwardMap, forwardCoeff);

    // The map of products (backward reaction)
    auto backwardMap = parsingReacSpc(backwardStr, reactionName);
    backwardCoeff.resize(static_cast<integer>(backwardMap.size()));
    readReacSpcFromMap(backwardMap, backwardCoeff);
}

void OpenHurricane::chemkinFileRead::parsingReaction() {
    auto unitMapA = unitMapOfA();
    auto unitMapE = unitMapOfE();

    std::regex reactionBeginRegex("REAC(?:|TION|TIONS)", std::regex::icase);
    std::regex reactionSingleRegex("(.*?)\\s*"
                                   "(<=>|=>|=)\\s*"
                                   "(.*?)"
                                   "\\s+((?:[-+]?[0-9]+|\\.)\\.*[0-9]*(?:[eEgGdD][-+]?[0-9]+)?)"
                                   "\\s+(.*?)"
                                   "\\s+(.*?)$|\\n");
    auto reactionLineStr = stringToLineList(reactionString_);
    integer reacSize = Zero;
    boolList isDuplicate;
    boolList isPresDepFallOffList;
    boolList isThirdBodyReactionList;

    for (integer i = 0; i < reactionLineStr.size(); ++i) {
        std::string reacName = reactionLineStr[i];
        const std::string reacNameLineStr = reacName;
        std::smatch reacNameSMatch;
        if (std::regex_search(reacNameLineStr, reacNameSMatch, reactionSingleRegex)) {
            reacName = reacNameSMatch[1];
            reacName += reacNameSMatch[2];
            reacName += reacNameSMatch[3];

            integer speId = -1;
            realArray efficiences;
            bool isPresDepFallOffReac =
                checkAndReplacePressureDependency(reactionLineStr[i], speId, efficiences);
            bool isThirdBodyReaction = checkAndReplaceThirdbody(reactionLineStr[i], efficiences);
            bool hasSetPresDepFallOffReac = false;
            auto uA = globalUnitTypeOfA_;
            auto uE = globalUnitTypeOfE_;
            std::smatch what;
            reactionRateTypes *reacRatePtr = nullptr;
            reactionRateTypes *revReacRatePtr = nullptr;
            fallOffFunctions::fallOffFunction *fofFPtr = nullptr;

            isDuplicate.append(false);
            isPresDepFallOffList.append(isPresDepFallOffReac);
            isThirdBodyReactionList.append(isThirdBodyReaction);
            if (std::regex_search(reactionLineStr[i], what, reactionSingleRegex)) {
                reactionType rType = reactionType::REVERSIBLE;
                bool isReversible = true;
                if (what[2] == "=>") { // irreversible
                    isReversible = false;
                    rType = reactionType::IRREVERSIBLE;
                }
                reactionRateType rrType = reactionRateType::Arrhenius;
                fallOffFunctionType fofType = fallOffFunctionType::unknown;

                ArrheniusRR ArrCoef, *revArrCoef = nullptr;
                readReal(what[4], ArrCoef.A());
                readReal(what[5], ArrCoef.beta());
                readReal(what[6], ArrCoef.Ta());

                List<reacSpcCoeffType> forwardCoeff;
                List<reacSpcCoeffType> backwardCoeff;
                parsingReacStoich(what[1], forwardCoeff, what[3], backwardCoeff, reacName);
                if (i >= reactionLineStr.size() && isPresDepFallOffReac) {
                    errorAbortStr(
                        ("The auxiliary information line(s) must follow the reaction to identify "
                         "the fall-off formulation and parameters for reaction: \"" +
                         reacName + "\""));
                }

                while (i < reactionLineStr.size() - 1) {
                    if (std::regex_search(reactionLineStr[i + 1], reactionSingleRegex)) {
                        if (isPresDepFallOffReac && !hasSetPresDepFallOffReac) {
                            errorAbortStr(("The auxiliary information line(s) must follow the "
                                           "reaction (described below) to identify the "
                                           "fall-off formulation and parameters.\n\"" +
                                           reacName + "\""));
                        }
                        break;
                    }
                    stringList keywordList;
                    integerList posList;
                    findLineKeywordType(reactionLineStr[i + 1], keywordList, posList);
                    if (keywordList.size() == 0) {
                        errorAbortStr(("Cannot reconized line: " + reactionLineStr[i + 1]));
                    }

                    for (integer keyI = 0; keyI < keywordList.size(); ++keyI) {
                        auto keyType = reactionAuxiliaryKeywordMap_.at(keywordList[keyI]);
                        switch (keyType) {
                        case reactionAuxiliaryKeyword::NeutralThirdBody: {
                            if (!isThirdBodyReaction && !isPresDepFallOffReac) {
                                errorAbortStr(
                                    ("No third-body symbols in reaction: " + reacName +
                                     ", but has efficiences section: " + reactionLineStr[i + 1]));
                            }

                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            parsingThirdBodyCoeffs(searchStr, efficiences, isPresDepFallOffReac,
                                                   speId);
                        } break;
                        case reactionAuxiliaryKeyword::LOW: {
                            hasSetPresDepFallOffReac = true;
                            if (!isPresDepFallOffReac) {
                                errorAbortStr(("LOW keyword found: " + reactionLineStr[i + 1] +
                                               ", but has no pressure dependent species "
                                               "in: " +
                                               reacName));
                            }
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::unimolecularFallOff;
                            } else {
                                errorAbortStr(("Attempt to change the rate type to unimolecular "
                                               "fall-off in reaction: " +
                                               reacName));
                            }
                            if (fofType == fallOffFunctionType::unknown) {
                                fofType = fallOffFunctionType::Lindemann;
                            }
                            std::regex lowRegex("(LOW)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)\\s*\\/",
                                                std::regex::icase);
                            string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            ArrheniusRR lowArrCoef;
                            std::smatch lowArrCoefMatch;
                            std::regex_search(searchStr, lowArrCoefMatch, lowRegex);
                            readReal(lowArrCoefMatch[2], lowArrCoef.A());
                            readReal(lowArrCoefMatch[3], lowArrCoef.beta());
                            readReal(lowArrCoefMatch[4], lowArrCoef.Ta());
                            if (fofFPtr == nullptr) {
                                fofFPtr = new fallOffFunctions::Lindemann();
                            }
                            reacRatePtr = new unimolecularFallOffRR(lowArrCoef, ArrCoef, *fofFPtr,
                                                                    specTable_, efficiences);
                        } break;
                        case reactionAuxiliaryKeyword::HIGH: {
                            hasSetPresDepFallOffReac = true;
                            if (!isPresDepFallOffReac) {
                                errorAbortStr(("HIGH keyword found: " + reactionLineStr[i + 1] +
                                               ", but has no pressure dependent "
                                               "species in: " +
                                               reacName));
                            }
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::chemicallyActivatedBimolecular;
                            } else {
                                errorAbortStr(("Attempt to change the rate type to chemically "
                                               "activated bimolecular "
                                               "fall-off in reaction: " +
                                               reacName));
                            }
                            if (fofType == fallOffFunctionType::unknown) {
                                fofType = fallOffFunctionType::Lindemann;
                            }
                            std::regex highRegex("(HIGH)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)\\s*\\/",
                                                 std::regex::icase);
                            string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            ArrheniusRR highArrCoef;
                            std::smatch highArrCoefMatch;
                            std::regex_search(searchStr, highArrCoefMatch, highRegex);
                            readReal(highArrCoefMatch[2], highArrCoef.A());
                            readReal(highArrCoefMatch[3], highArrCoef.beta());
                            readReal(highArrCoefMatch[4], highArrCoef.Ta());
                            if (fofFPtr == nullptr) {
                                fofFPtr = new fallOffFunctions::Lindemann();
                            }
                            reacRatePtr = new chemicallyActicatedBimolecularRR(
                                ArrCoef, highArrCoef, *fofFPtr, specTable_, efficiences);
                        } break;
                        case reactionAuxiliaryKeyword::COLLEFF: {
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::collision;
                            } else {
                                errorAbortStr(
                                    ("Attempt to change the rate type to collision in reaction: " +
                                     reacName));
                            }
                            if (forwardCoeff.size() != 2) {
                                errorAbortStr(("The collision reaction should be "
                                               "bimolecular, but the reaction is: " +
                                               reacName));
                            }

                            //Using the Lennard-Jones diameter as an approximation for the spherical diameter of a species.
                            realArray collisionCoeffs(3);

                            //  The average diameter of the two spherical particles
                            real dd = Zero;
                            for (integer collI = 0; collI < forwardCoeff.size(); ++collI) {
                                dd += transport_.tranTable()[forwardCoeff[collI].index_].sigma();
                            }
                            dd /= 2.0;

                            reacRatePtr = new collisionRR(forwardCoeff.first().index_,
                                                          forwardCoeff.last().index_, specTable_,
                                                          dd, ArrCoef);
                        } break;
                        case reactionAuxiliaryKeyword::TROE: {
                            hasSetPresDepFallOffReac = true;
                            if (!isPresDepFallOffReac) {
                                errorAbortStr(("TROE keyword found: " + reactionLineStr[i + 1] +
                                               ", but has no pressure dependent "
                                               "species in: " +
                                               reacName));
                            }
                            if (fofType == fallOffFunctionType::unknown ||
                                fofType == fallOffFunctionType::Lindemann) {
                                fofType = fallOffFunctionType::Troe;
                                HurDelete(fofFPtr);
                            } else {
                                errorAbortStr(
                                    ("Attempt to change the fall-off type to Troe in reaction: " +
                                     reacName));
                            }
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            parsingTroe((void **)&fofFPtr, searchStr);
                        } break;
                        case reactionAuxiliaryKeyword::SRI: {
                            hasSetPresDepFallOffReac = true;
                            if (!isPresDepFallOffReac) {
                                errorAbortStr(("SRI keyword found: " + reactionLineStr[i + 1] +
                                               ", but has no pressure dependent species "
                                               "in: " +
                                               reacName));
                            }
                            if (fofType == fallOffFunctionType::unknown ||
                                fofType == fallOffFunctionType::Lindemann) {
                                fofType = fallOffFunctionType::SRI;
                                HurDelete(fofFPtr);
                            } else {
                                errorAbortStr(
                                    ("Attempt to change the fall-off type to SRI in reaction: " +
                                     reacName));
                            }
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            parsingSRI((void **)&fofFPtr, searchStr);
                        } break;
                        case reactionAuxiliaryKeyword::LT: {
                            if (isPresDepFallOffReac) {
                                errorAbortStr(("LT keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot be used in "
                                               "pressure-dependent reaction in: " +
                                               reacName));
                            }
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::LandauTeller;
                            } else if (rrType != reactionRateType::LandauTeller) {
                                errorAbortStr(("Attempt to change the rate type to Landau-Teller "
                                               "in reaction: " +
                                               reacName));
                            }
                            realArray LandauTellerCoeffs(2);
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex LTRegex("(LT)\\s*\\/\\s*(.*?)\\s+(.*?)\\s*\\/",
                                               std::regex::icase);
                            if (!std::regex_search(searchStr, what, LTRegex)) {
                                errorAbortStr(("LT keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot find paremeters in: " + reacName));
                            }
                            readReal(what[2], LandauTellerCoeffs[0]);
                            readReal(what[3], LandauTellerCoeffs[1]);
                            reacRatePtr =
                                new LandauTellerRR(ArrCoef.A(), ArrCoef.beta(), ArrCoef.Ta(),
                                                   LandauTellerCoeffs[0], LandauTellerCoeffs[1]);
                        } break;
                        case reactionAuxiliaryKeyword::RLT: {
                            if (isPresDepFallOffReac) {
                                errorAbortStr(("RLT keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot be used in "
                                               "pressure-dependent reaction in: " +
                                               reacName));
                            }
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::LandauTeller;
                            } else if (rrType != reactionRateType::LandauTeller) {
                                errorAbortStr(("Attempt to change the rate type to reverse "
                                               "Landau-Teller in reaction: " +
                                               reacName));
                            }
                            if (rType != reactionType::NONEQUILIBRIUMREVERSIBLE) {
                                LFatal(
                                    "The keyword RLT must follow the explicit REV parameters");
                            }
                            realArray revLandauTellerCoeffs(2);
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex RLTRegex("(RLT)\\s*\\/\\s*(.*?)\\s+(.*?)\\s*\\/",
                                                std::regex::icase);
                            if (!std::regex_search(searchStr, what, RLTRegex)) {
                                errorAbortStr(("RLT keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot find paremeters in: " + reacName));
                            }
                            readReal(what[2], revLandauTellerCoeffs[0]);
                            readReal(what[3], revLandauTellerCoeffs[1]);

                            if (revArrCoef == nullptr) {
                                revReacRatePtr = new LandauTellerRR(
                                    ArrCoef.A(), ArrCoef.beta(), ArrCoef.Ta(),
                                    revLandauTellerCoeffs[0], revLandauTellerCoeffs[1]);
                            } else {
                                revReacRatePtr = new LandauTellerRR(
                                    revArrCoef->A(), revArrCoef->beta(), revArrCoef->Ta(),
                                    revLandauTellerCoeffs[0], revLandauTellerCoeffs[1]);
                            }
                        } break;
                        case reactionAuxiliaryKeyword::JAN: {
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::Janev;
                            } else {
                                errorAbortStr(
                                    ("Attempt to change the rate type to Janev in reaction: " +
                                     reacName));
                            }
                            FixedList<real, 9> janevList;
                            janevList.setZero();
                            std::smatch JanWhat;
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex JANRegex(
                                "(JAN)\\s*\\/"
                                "\\s*(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(."
                                "*?)\\s+(.*?)\\s+(.*?)\\s*\\/",
                                std::regex::icase);
                            if (!std::regex_search(searchStr, JanWhat, JANRegex)) {
                                errorAbortStr(("JAN keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot find paremeters in: " + reacName));
                            }
                            for (integer jli = 0; jli < 9; jli++) {
                                readReal(JanWhat[jli + 2], janevList[jli]);
                            }
                            reacRatePtr =
                                new JanevRR(ArrCoef.A(), ArrCoef.beta(), ArrCoef.Ta(), janevList);
                        } break;
                        case reactionAuxiliaryKeyword::FIT1: {
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::powerSeries;
                            } else {
                                errorAbortStr(("Attempt to change the rate type to power-series in "
                                               "reaction: " +
                                               reacName));
                            }
                            FixedList<real, 4> powerSeriesCoeffs;
                            powerSeriesCoeffs.setZero();
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex FIT1Regex(
                                "(FIT1)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s*\\/",
                                std::regex::icase);
                            if (!std::regex_search(searchStr, what, FIT1Regex)) {
                                errorAbortStr(("FIT1 keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot find paremeters in: " + reacName));
                            }
                            for (integer jli = 0; jli < 4; jli++) {
                                readReal(what[jli + 2], powerSeriesCoeffs[jli]);
                            }
                            reacRatePtr =
                                new powerSeriesRR(ArrCoef.A(), ArrCoef.beta(), powerSeriesCoeffs);
                        } break;
                        case reactionAuxiliaryKeyword::TDEP: {
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex TDEPRegex("(TDEP)\\s*\\/\\s*(.*?)\\s*\\/",
                                                 std::regex::icase);
                            std::regex_search(searchStr, what, TDEPRegex);
                            std::string speN = what[2];
                            integer TDepId;
                            if (!checkReacSpec(speN, TDepId)) {
                                errorAbortStr(("Unknown species: " + speN + " in: " + reacName));
                            }
                            errorAbortStr(("Reaction type: Species Temperature "
                                           "Dependence not yet supported in:" +
                                           reacName));
                        } break;
                        case reactionAuxiliaryKeyword::XSMI: {
                            errorAbortStr(("Reaction type: Collision Cross Section "
                                           "not yet supported in:" +
                                           reacName));
                        } break;
                        case reactionAuxiliaryKeyword::MOME: {
                            errorAbortStr(("Reaction type: Plasma Momentum-Transfer Collision "
                                           "Frequency not yet supported in:" +
                                           reacName));
                        } break;
                        case reactionAuxiliaryKeyword::EXCI: {
                            errorAbortStr(("Reaction type: Energy Loss Parameter "
                                           "not yet supported in:" +
                                           reacName));
                        } break;
                        case reactionAuxiliaryKeyword::REV: {
                            rType = reactionType::NONEQUILIBRIUMREVERSIBLE;
                            realArray reverseArrheniusCoeffs(3, Zero);
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex RevRegex("(REV)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)\\s*\\/",
                                                std::regex::icase);
                            std::smatch resultMatch;
                            if (std::regex_search(searchStr, resultMatch, RevRegex)) {
                                revArrCoef = new ArrheniusRR();
                                readReal(resultMatch[2], revArrCoef->A());
                                readReal(resultMatch[3], revArrCoef->beta());
                                readReal(resultMatch[4], revArrCoef->Ta());
                            }
                        } break;
                        case reactionAuxiliaryKeyword::DUP: {
                            isDuplicate[reacSize] = true;
                        } break;
                        case reactionAuxiliaryKeyword::FORD: {
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex fordRegex("(FORD)\\s*\\/\\s*(.*?)\\s+(.*?)\\s*\\/",
                                                 std::regex::icase);
                            parsingReacSpcOrder(reacName, searchStr, fordRegex, forwardCoeff);
                        } break;
                        case reactionAuxiliaryKeyword::RORD: {
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex rordRegex("(RORD)\\s*\\/\\s*(.*?)\\s+(.*?)\\s*\\/",
                                                 std::regex::icase);
                            parsingReacSpcOrder(reacName, searchStr, rordRegex, backwardCoeff);
                        } break;
                        case reactionAuxiliaryKeyword::UNITS: {
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex unitRegex("(UNITS)\\s*\\/\\s*(.*?)\\s*\\/",
                                                 std::regex::icase);
                            if (std::regex_search(searchStr, what, unitRegex)) {
                                std::string newUnit = what[2];
                                auto iterA = unitMapA.find(newUnit);
                                auto iterE = unitMapE.find(newUnit);
                                if (iterA != unitMapA.end()) {
                                    uA = iterA->second;
                                } else if (iterE != unitMapE.end()) {
                                    uE = iterE->second;
                                } else {
                                    errorAbortStr(
                                        ("Wrong units: " + newUnit + " in reaction: " + reacName));
                                }
                            }
                        } break;
                        case reactionAuxiliaryKeyword::PLOG: {
                            if (isPresDepFallOffReac) {
                                errorAbortStr(("Plog keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot be used in "
                                               "third-body reaction in: " +
                                               reacName));
                            }
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::
                                    pressureDependenceUsingLogarithmicInterpolation;
                            } else if (rrType !=
                                       reactionRateType::
                                           pressureDependenceUsingLogarithmicInterpolation) {
                                errorAbortStr(("Attempt to change the rate type to pressure "
                                               "dependence using logarithmic "
                                               "interpolation in reaction: " +
                                               reacName));
                            }
                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex plogRegex(
                                "(PLOG)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s*\\/",
                                std::regex::icase);
                            std::smatch pArrCoefMatch;
                            if (!std::regex_search(searchStr, pArrCoefMatch, plogRegex)) {
                                errorAbortStr(("Plog keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot find paremeters in: " + reacName));
                            }
                            ArrheniusRR pArrCoef;
                            if (reacRatePtr == nullptr) {
                                reacRatePtr = new logInterpolationRR();
                            }
                            real p;
                            readReal(pArrCoefMatch[2], p);
                            readReal(pArrCoefMatch[3], pArrCoef.A());
                            readReal(pArrCoefMatch[4], pArrCoef.beta());
                            readReal(pArrCoefMatch[5], pArrCoef.Ta());

                            dynamic_cast<logInterpolationRR *>(reacRatePtr)->append(p, pArrCoef);
                        } break;
                        case reactionAuxiliaryKeyword::CHEB: {
                            hasSetPresDepFallOffReac = true;
                            if (!isPresDepFallOffReac) {
                                errorAbortStr(("CHEB keyword found: " + reactionLineStr[i + 1] +
                                               ", but has no pressure dependent "
                                               "species in: " +
                                               reacName));
                            }
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::ChebyshevReactionRate;
                            } else if (rrType != reactionRateType::ChebyshevReactionRate) {
                                errorAbortStr(("Attempt to change the rate type to Chebyshev "
                                               "polynomial rate expression in reaction: " +
                                               reacName));
                            }

                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex chebRegex("(CHEB)\\s*\\/\\s*(.*?)\\s*\\/",
                                                 std::regex::icase);
                            if (!std::regex_search(searchStr, what, chebRegex)) {
                                errorAbortStr(("CHEB keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot find paremeters in: " + reacName));
                            }
                            std::string dataStr = what[2];
                            std::regex dataRegex("\\S+");
                            std::string::const_iterator start = dataStr.begin();
                            std::string::const_iterator end = dataStr.end();
                            realArray dataList;
                            while (std::regex_search(start, end, what, dataRegex)) {
                                real datai;
                                readReal(what[0], datai);
                                dataList.append(datai);
                                start = what[0].second;
                            }
                            if (reacRatePtr == nullptr) {
                                reacRatePtr = new ChebyshevPolynomialsRR();
                            }
                            dynamic_cast<ChebyshevPolynomialsRR *>(reacRatePtr)
                                ->tmpAnm()
                                .append(dataList);
                        } break;
                        case reactionAuxiliaryKeyword::PCHEB: {
                            hasSetPresDepFallOffReac = true;
                            if (!isPresDepFallOffReac) {
                                errorAbortStr(("PCHEB keyword found: " + reactionLineStr[i + 1] +
                                               ", but has no pressure dependent "
                                               "species in: " +
                                               reacName));
                            }
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::ChebyshevReactionRate;
                            } else if (rrType != reactionRateType::ChebyshevReactionRate) {
                                errorAbortStr(("Attempt to change the rate type to Chebyshev "
                                               "polynomial rate expression in reaction: " +
                                               reacName));
                            }

                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex pchebRegex("(PCHEB)\\s*\\/\\s*(.*?)\\s+(.*?)\\s*\\/",
                                                  std::regex::icase);
                            if (!std::regex_search(searchStr, what, pchebRegex)) {
                                errorAbortStr(("PCHEB keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot find paremeters in: " + reacName));
                            }
                            realArray pMinMax(2);
                            readReal(what[2], pMinMax[0]);
                            readReal(what[3], pMinMax[1]);
                            if (reacRatePtr == nullptr) {
                                reacRatePtr = new ChebyshevPolynomialsRR();
                            }
                            dynamic_cast<ChebyshevPolynomialsRR *>(reacRatePtr)->Pmin() =
                                pMinMax[0];
                            dynamic_cast<ChebyshevPolynomialsRR *>(reacRatePtr)->Pmax() =
                                pMinMax[1];
                        } break;
                        case reactionAuxiliaryKeyword::TCHEB: {
                            hasSetPresDepFallOffReac = true;
                            if (!isPresDepFallOffReac) {
                                errorAbortStr(("TCHEB keyword found: " + reactionLineStr[i + 1] +
                                               ", but has no pressure dependent "
                                               "species in: " +
                                               reacName));
                            }
                            if (rrType == reactionRateType::Arrhenius) {
                                rrType = reactionRateType::ChebyshevReactionRate;
                            } else if (rrType != reactionRateType::ChebyshevReactionRate) {
                                errorAbortStr(("Attempt to change the rate type to Chebyshev "
                                               "polynomial rate expression in reaction: " +
                                               reacName));
                            }

                            std::string searchStr = reactionLineStr[i + 1].substr(posList[keyI]);
                            std::regex tchebRegex("(TCHEB)\\s*\\/\\s*(.*?)\\s+(.*?)\\s*\\/",
                                                  std::regex::icase);
                            if (!std::regex_search(searchStr, what, tchebRegex)) {
                                errorAbortStr(("TCHEB keyword found: " + reactionLineStr[i + 1] +
                                               ", which cannot find paremeters in: " + reacName));
                            }
                            realArray tMinMax(2);
                            readReal(what[2], tMinMax[0]);
                            readReal(what[3], tMinMax[1]);
                            if (reacRatePtr == nullptr) {
                                reacRatePtr = new ChebyshevPolynomialsRR();
                            }
                            dynamic_cast<ChebyshevPolynomialsRR *>(reacRatePtr)->Tmin() =
                                tMinMax[0];
                            dynamic_cast<ChebyshevPolynomialsRR *>(reacRatePtr)->Tmax() =
                                tMinMax[1];
                        } break;
                        case reactionAuxiliaryKeyword::END:
                            break;
                        default: {
                            errorAbortStr(("Unknown chemkin keyword: " + keywordList[keyI] +
                                           " in: " + reacName));
                        } break;
                        }
                    }
                    i++;
                }

                if (rrType == reactionRateType::Arrhenius && efficiences.size() != 0) {
                    rrType = reactionRateType::thirdBodyArrhenius;
                    reacRatePtr = new ArrheniusThirdBodyRR(ArrCoef, specTable_, efficiences);
                } else if (reacRatePtr == nullptr) {
                    reacRatePtr = new ArrheniusRR(ArrCoef);
                }
                if (checkReactionDuplicate(forwardCoeff, backwardCoeff, isDuplicate,
                                           isPresDepFallOffList, isThirdBodyReactionList,
                                           reacSize)) {
                    errorAbortStr(("Reaction: " + reacName +
                                   " is duplicate, but with no DUPLICATE keyword."));
                }

                addReactionToList(reacName, forwardCoeff, backwardCoeff, isReversible, rType,
                                  rrType, fofType, reacRatePtr, revReacRatePtr, fofFPtr, uA, uE,
                                  efficiences, &ArrCoef, revArrCoef, reacSize);
                forwardCoeff.clear();
                backwardCoeff.clear();
                HurDelete(revArrCoef);
                HurDelete(reacRatePtr);
                HurDelete(revReacRatePtr);
                HurDelete(fofFPtr);
                reacSize++;
            }
        } else {
            if (!std::regex_search(reactionLineStr[i], reactionBeginRegex)) {
                errorAbortStr(("Unknown Reaction line: " + reactionLineStr[i]));
            }
        }
    }
}

void OpenHurricane::chemkinFileRead::parsingThirdBodyCoeffs(const std::string &thirdBodyLineStr,
                                                        realArray &efficiens,
                                                        bool isPressureDependent,
                                                        const integer enhancedSpeId) const {
    const std::regex thirdBodyRegex("(\\S+?)\\s*\\/\\s*(.*?)\\s*\\/");
    std::string::const_iterator start = thirdBodyLineStr.begin();
    auto end = thirdBodyLineStr.end();
    std::smatch what;
    while (std::regex_search(start, end, what, thirdBodyRegex)) {
        integer id;
        if (checkReacSpec(what[1], id)) {
            if (isPressureDependent) {
                if (enhancedSpeId != -1 && enhancedSpeId != id) {
                    errorAbortStr(("The third body in reaction should be: " +
                                   specTable_[enhancedSpeId].name() + ", but " + thirdBodyLineStr +
                                   " were found."));
                }
            }
            readReal(what[2], efficiens[id]);
        } else {
            errorAbortStr(("Unknown species: " + std::string(what[1])));
        }
        start = what[0].second;
    }
}

hur_nodiscard OpenHurricane::stdStringList
OpenHurricane::chemkinFileRead::splitReacSpecToStrList(const std::string &reacSpecies,
                                                   const std::string &reactionName) const {
    stdStringList splitRS;
    auto trimSpecStr = trimCopy(reacSpecies);
    if (trimSpecStr.empty()) {
        errorAbortStr(("Empty species in reaction: \"" + reactionName + "\""));
    }

    auto pos = trimSpecStr.find_first_of("+", 0);
    if (pos == 0) {
        errorAbortStr(
            ("Wrong species form: \"" + reacSpecies + "\" in reaction: \"" + reactionName + "\""));
    }
    if (pos == std::string::npos || pos == trimSpecStr.size() - 1) {
        splitRS.push_back(trimSpecStr);
        return splitRS;
    }

    auto lastPos = trimSpecStr.find_first_not_of("+", 0);
    integer count = 0;
    while (pos != std::string::npos) {
        splitRS.push_back(trimCopy(trimSpecStr.substr(lastPos, pos - lastPos)));
        auto nextPos = trimSpecStr.find_first_of("+", pos + 1);
        lastPos = trimSpecStr.find_first_not_of("+", pos + 1);
        if (nextPos == std::string::npos && lastPos == std::string::npos) {
            splitRS[count] += "+";
            count++;
            pos = nextPos;
        } else if (nextPos == std::string::npos) {
            count++;
            pos = nextPos;
            splitRS.push_back(trimCopy(trimSpecStr.substr(lastPos, pos - lastPos)));
            count++;
        } else if (lastPos == std::string::npos) {
            splitRS[count] += "++";
            count++;
            pos = nextPos;
        } else if (nextPos < lastPos) {
            splitRS[count] += "+";
            count++;
            nextPos = trimSpecStr.find_first_of("+", nextPos + 1);
            if (nextPos == std::string::npos) {
                splitRS.push_back(trimCopy(trimSpecStr.substr(lastPos, nextPos - lastPos)));
                count++;
            }
            pos = nextPos;
        } else if (trimCopy(trimSpecStr.substr(lastPos, nextPos - lastPos)).empty()) {
            splitRS[count] += "+";
            lastPos = trimSpecStr.find_first_not_of("+", nextPos + 1);
            nextPos = trimSpecStr.find_first_of("+", nextPos + 1);

            if (nextPos == std::string::npos && lastPos == std::string::npos) {
                splitRS[count] += "+";
            } else if (nextPos == std::string::npos) {
                splitRS.push_back(trimCopy(trimSpecStr.substr(lastPos, nextPos - lastPos)));
            }
            count++;
            pos = nextPos;
        } else {
            count++;
            pos = nextPos;
        }
    }
    return splitRS;
}

bool OpenHurricane::chemkinFileRead::checkReacSpec(const std::string &rsp, integer &id) const {
    for (integer i = 0; i < specTable_.size(); ++i) {
        if (rsp == specTable_[i].name()) {
            id = i;
            return true;
        }
    }
    id = -1;
    return false;
}

hur_nodiscard std::map<std::string, OpenHurricane::real>
OpenHurricane::chemkinFileRead::parsingReacSpc(const std::string &reacSpc,
                                           const std::string &reactionName) const {
    auto strl = splitReacSpecToStrList(reacSpc, reactionName);
    std::map<std::string, real> reacSpcM;
    std::regex stoichiometryRegex("((?:[0-9]*)\\.*[0-9]*)\\s*(\\S+)");

    for (auto &e : strl) {
        std::smatch resultMatch;
        auto found = std::regex_search(e, resultMatch, stoichiometryRegex);
        auto spcName = trimCopy(resultMatch[2]);
        real stc = 0;
        if (resultMatch[1] == "") {
            stc = 1;
        } else {
            readReal(resultMatch[1], stc);
        }

        if (reacSpcM.find(spcName) != reacSpcM.end()) {
            reacSpcM[spcName] += stc;
        } else {
            reacSpcM.emplace(spcName, stc);
        }
    }

    return reacSpcM;
}

void OpenHurricane::chemkinFileRead::parsingTroe(void **fofPtr, const string &TroeStr) const {
    std::regex TroeRegex("(TROE)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)(?:|\\s+(.*?))\\s*\\/",
                         std::regex::icase);
    std::smatch resultMatch;
    if (std::regex_search(TroeStr, resultMatch, TroeRegex)) {
        realList Coeffs;
        for (size_t jj = 2; jj < resultMatch.size(); ++jj) {
            if (resultMatch[jj] != "") {
                real val = 0;
                readReal(resultMatch[jj], val);
                Coeffs.push_back(val);
            }
        }

        if (Coeffs.size() == 3) {
            *fofPtr = new fallOffFunctions::Troe(Coeffs[0], Coeffs[1], Coeffs[2]);
        } else if (Coeffs.size() == 4) {
            *fofPtr = new fallOffFunctions::Troe(Coeffs[0], Coeffs[1], Coeffs[2], Coeffs[3]);
        }
    }
}

void OpenHurricane::chemkinFileRead::parsingSRI(void **fofPtr, const string &SRIStr) const {
    std::regex SRIRegex("(SRI)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)(?:|\\s+(.*?)\\s+(.*?))\\s*\\/",
                        std::regex::icase);
    std::smatch resultMatch;
    if (std::regex_search(SRIStr, resultMatch, SRIRegex)) {
        realList Coeffs;
        for (size_t jj = 2; jj < resultMatch.size(); ++jj) {
            if (resultMatch[jj] != "") {
                real val = 0;
                readReal(resultMatch[jj], val);
                Coeffs.push_back(val);
            }
        }

        if (Coeffs.size() == 3) {
            *fofPtr = new fallOffFunctions::SRI(Coeffs[0], Coeffs[1], Coeffs[2]);
        } else if (Coeffs.size() == 5) {
            *fofPtr =
                new fallOffFunctions::SRI(Coeffs[0], Coeffs[1], Coeffs[2], Coeffs[3], Coeffs[4]);
        }
    }
}

void OpenHurricane::chemkinFileRead::parsingReacSpcOrder(const std::string &reacName,
                                                     const std::string &str, std::regex &reg,
                                                     List<reacSpcCoeffType> &coef) const {
    std::smatch resultMatch;
    std::regex_search(str, resultMatch, reg);
    std::string speN = resultMatch[2];
    real newOrder;
    readReal(resultMatch[3], newOrder);
    integer newId;
    if (!checkReacSpec(speN, newId)) {
        errorAbortStr(("Unknown species name: " + speN + " for setting new order in " + reacName));
    }
    bool foundInLRHS = false;
    for (integer i = 0; i < coef.size(); ++i) {
        integer index = coef[i].index_;
        if (newId == index) {
            foundInLRHS = true;
            coef[i].order_ = newOrder;
            break;
        }
    }
    if (!foundInLRHS) {
        errorAbortStr(
            ("The species : " + speN +
             " to be set to new order is not found in reactants or products list in: " + reacName));
    }
}

bool OpenHurricane::chemkinFileRead::checkAndReplacePressureDependency(std::string &reacStr,
                                                                   integer &id,
                                                                   realArray &efficiences) const {
    std::regex pressureDependentRegex("\\(\\s*?\\+\\s*?(.*?)\\s*?\\)");
    std::smatch what;
    std::string::const_iterator start = reacStr.begin();
    std::string::const_iterator end = reacStr.end();
    integer count = Zero;
    stringList whatStr(3);
    while (std::regex_search(start, end, what, pressureDependentRegex)) {
        whatStr[count] = what[1];
        count++;
        if (count > 2) {
            id = -2;
            break;
        }
        start = what[0].second;
    }

    if (count == 0) {
        id = -1;
        return false;
    } else if (count == 1) {
        id = -2;
        errorAbortStr(
            ("Reaction: \"" + trimCopy(reacStr) + "\" is wrong for pressure-dependent reaction."));
    } else if (count == 2) {
        trim(whatStr[0]);
        trim(whatStr[1]);
        if (whatStr[0] != whatStr[1]) {
            id = -2;
            checkWarning("The species acting as the third body are not the same.");
        } else {

            if (whatStr[0] == "M") {
                reacStr = std::regex_replace(reacStr, pressureDependentRegex, "");
                id = -1;
            } else {
                if (!checkReacSpec(whatStr[0], id)) {
                    id = -2;
                    checkWarningStr(
                        ("The species: \"" + whatStr[0] + "\" is not in the species list."));
                } else {
                    reacStr = std::regex_replace(reacStr, pressureDependentRegex, "");
                }
            }
        }
    }

    if (checkAndReplaceThirdbody(reacStr, efficiences)) {
        errorAbortStr(("A species or all species are a third-body in a fall-off reaction, "
                       "and +M also appears in the reaction.\nReaction: \"" +
                       trimCopy(reacStr) + "\""));
    }
    efficiences.resize(specTable_.size(), 1.0);
    if (id != -1) {
        efficiences = Zero;
        efficiences[id] = 1.0;
    }
    return true;
}

bool OpenHurricane::chemkinFileRead::checkAndReplaceThirdbody(std::string &reacStr,
                                                          realArray &efficiences) const {
    const std::regex thirdBodyRegex("\\+\\s*?M");
    std::string::const_iterator start = reacStr.begin();
    std::string::const_iterator end = reacStr.end();
    integer count = Zero;
    std::smatch what;
    while (std::regex_search(start, end, what, thirdBodyRegex)) {
        count++;
        start = what[0].second;
    }

    if (count == 0) {
        return false;
    } else if (count != 2) {
        errorAbortStr(
            ("Reaction: \"" + trimCopy(reacStr) + "\" is wrong for third body reaction."));
    } else {
        reacStr = std::regex_replace(reacStr, thirdBodyRegex, "");
        efficiences.resize(specTable_.size(), 1.0);
    }

    return true;
}

hur_nodiscard bool OpenHurricane::chemkinFileRead::checkReactionDuplicate(
    const List<reacSpcCoeffType> &forwardCoeff, const List<reacSpcCoeffType> &backwardCoeff,
    const boolList &isDuplicate, const boolList &isPressureDependent, const boolList &isThirdBody,
    const integer currentId) const {
    for (integer irc = 0; irc < reactions_.size(); ++irc) {
        if (checkReacSpecList(forwardCoeff, reactions_[irc].forwardCoeffs()) &&
            checkReacSpecList(backwardCoeff, reactions_[irc].backwardCoeffs())) {
            integer id1 = reactions_[irc].index();
            if (isDuplicate[id1] && isDuplicate[currentId]) {
                return false;
            } else {
                if (!(isPressureDependent[id1] && isPressureDependent[currentId])) {
                    return false;
                } else if (!(isThirdBody[id1] && isThirdBody[currentId])) {
                    return false;
                } else if (!(isPressureDependent[id1] && isThirdBody[currentId])) {
                    return false;
                } else if (!(isPressureDependent[currentId] && isThirdBody[id1])) {
                    return false;
                }
                return true;
            }
        }
    }
    return false;
}

hur_nodiscard bool
OpenHurricane::chemkinFileRead::checkReacSpecList(const List<reacSpcCoeffType> &coeff1,
                                              const List<reacSpcCoeffType> &coeff2) const {
    if (coeff1.size() != coeff2.size()) {
        return false;
    } else {
        for (integer i = 0; i < coeff1.size(); ++i) {
            if (coeff1[i].index_ != coeff2[i].index_) {
                return false;
            } else if (coeff1[i].stoichCoeff_ != coeff2[i].stoichCoeff_) {
                return false;
            }
        }
    }
    return true;
}

OpenHurricane::stringList OpenHurricane::chemkinFileRead::parsingElements() {
    //	An isotope name (that is, a name not on the periodic chart) must be followed by its atomic weight (in grams per mole) delimited by slashes.
    std::regex isotopeAtomicWeightRegex("\\/\\s*(.*?)\\s*\\/");

    std::regex elementsListRegex("ELEM(?:|ENT|ENTS)\\s+([\\s\\S]*?)\\s+END", std::regex::icase);
    std::smatch resultMatch;
    std::regex_search(chemkinFileString_, resultMatch, elementsListRegex);
    std::string elementString = resultMatch[1];

    std::string::const_iterator start = elementString.begin();
    std::string::const_iterator end = elementString.end();
    std::match_results<std::string::const_iterator> matchs;
    stringList tempEle;
    std::regex elementRegex("(\\w+)");

    while (std::regex_search(start, end, matchs, elementRegex)) {
        string tempEn = string(matchs[1]);
        tempEle.push_back(tempEn);
        start = matchs[0].second;
        if (std::regex_search(start, end, matchs, isotopeAtomicWeightRegex)) {
            real weights;
            readReal(std::string(matchs[1]), weights);
            if (weights <= real(0.0)) {
                errorAbortStr(("The molecular weight: " + toString(weights) +
                               " for isotope atomic: " + tempEn + " is negative"));
            }
            isotopeAtomicWts_.emplace(tempEn, weights);
            start = matchs[0].second;
        }
    }

    specTable_.elementsNameList() = tempEle;

    for (integer i = 0; i < tempEle.size(); ++i) {
        elementsIndex_.emplace(tempEle[i], i);
    }
    return tempEle;
}

OpenHurricane::stringList OpenHurricane::chemkinFileRead::parsingSpecies() {
    std::regex speciesListRegex("SPEC(?:|IE|IES)\\s+([\\s\\S]*?)\\s+END", std::regex::icase);
    std::smatch resultMatch;
    std::regex_search(chemkinFileString_, resultMatch, speciesListRegex);
    std::string speciesString = resultMatch[1];

    std::string::const_iterator start = speciesString.begin();
    std::string::const_iterator end = speciesString.end();
    std::match_results<std::string::const_iterator> matchs;
    stringList tempSp;
    std::regex speciesRegex("\\S+");
    while (std::regex_search(start, end, matchs, speciesRegex)) {
        tempSp.append(string(matchs[0]));
        start = matchs[0].second;
    }

    specTable_.resize(tempSp.size());
    for (integer i = 0; i < tempSp.size(); ++i) {
        specTable_[i].name() = tempSp[i];
    }

    return tempSp;
}

void OpenHurricane::chemkinFileRead::parsingGlobalUnits() {
    std::regex globalUnitsRegex("REAC(?:|TION|TIONS)\\s*"
                                "\\b(\\w+)?\\b\\s*?"
                                "\\b(\\w+)?\\b\\s*?\\n",
                                std::regex::icase);

    std::smatch units;
    std::string::const_iterator start = chemkinFileString_.begin();
    std::string::const_iterator end = chemkinFileString_.end();
    bool matched = false;
    while (std::regex_search(start, end, units, globalUnitsRegex)) {
        if (matched) {
            LFatal("Reaction section error: multiple.");
        }
        if (units[1].matched) {
            std::string unit1;
            unit1 = units[1];
            std::string unit2 = nullString;
            stringToUpperCase(unit1);

            auto unitMapA = unitMapOfA();
            auto unitMapE = unitMapOfE();
            auto iterUA = unitMapA.find(unit1);
            auto iterUE = unitMapE.find(unit1);
            integer flag12 = 0;
            if (iterUA == unitMapA.end() && iterUE == unitMapE.end()) {
                errorAbortStr(("Unknown unite type: \"" + unit1 + "\" in chemkin file."));
            } else if (iterUA == unitMapA.end()) {
                flag12 = 1;
                globalUnitTypeOfE_ = iterUE->second;
            } else {
                flag12 = 2;
                globalUnitTypeOfA_ = iterUA->second;
            }

            if (units[2].matched) {
                unit2 = units[2];
                stringToUpperCase(unit2);
                if (unit1 == unit2) {
                    checkWarningStr(("Units are the same in chemkin file for "
                                     "reaction. Unit1 = \"" +
                                     unit1 + "\"; Unit2 = \"" + unit2 + "\"."));
                    start = units[0].second;
                    continue;
                }

                auto iterA2 = unitMapA.find(unit2);
                auto iterE2 = unitMapE.find(unit2);
                if (iterA2 == unitMapA.end() && iterE2 == unitMapE.end()) {
                    errorAbortStr(("Unknown unite type: \"" + unit2 + "\" in chemkin file."));
                } else if (iterA2 == unitMapA.end()) {
                    if (flag12 == 1) {
                        checkWarningStr(("Unit for E has been set twice. unit1 = \"" + unit1 +
                                         "\"; unit2 = \"" + unit2 + "\""));
                    }
                    globalUnitTypeOfE_ = iterE2->second;
                } else {
                    if (flag12 == 2) {
                        checkWarningStr(("Unit for A has been set twice. unit1 = \"" + unit1 +
                                         "\"; unit2 = \"" + unit2 + "\""));
                    }
                    globalUnitTypeOfA_ = iterA2->second;
                }
            }

        } else {
            if (units[2].matched) {
                std::string unit2;
                unit2 = units[2];
                stringToUpperCase(unit2);

                auto unitMapA = unitMapOfA();
                auto unitMapE = unitMapOfE();
                auto iterUA = unitMapA.find(unit2);
                auto iterUE = unitMapE.find(unit2);
                if (iterUA == unitMapA.end() && iterUE == unitMapE.end()) {
                    errorAbortStr(("Unknown unite type: \"" + unit2 + "\" in chemkin file."));
                } else if (iterUA == unitMapA.end()) {
                    globalUnitTypeOfE_ = iterUE->second;
                } else {
                    globalUnitTypeOfA_ = iterUA->second;
                }
            }
        }
        start = units[0].second;
        matched = true;
    }
}

void OpenHurricane::chemkinFileRead::findLineKeywordType(const std::string &lineStr,
                                                     stringList &lineKeyword,
                                                     integerList &pos) const {
    std::regex keywordRegex("(?:(\\S+?)\\s*\\/[\\s\\S]*?\\/|(\\S+?)(?:\\b|$))");
    std::smatch resultMatch;
    auto start = lineStr.begin();
    const std::string::const_iterator start0 = start;
    auto end = lineStr.end();

    while (std::regex_search(start, end, resultMatch, keywordRegex)) {
        pos.append(integer(std::distance(start0, start)));
        string keyw = string(resultMatch[1]);
        if (!resultMatch[1].matched) {
            keyw = string(resultMatch[2]);
        }
        stringToUpperCase(keyw);
        integer id;
        auto iter = reactionAuxiliaryKeywordMap_.find(keyw);
        if (iter != reactionAuxiliaryKeywordMap_.end()) {
            lineKeyword.append(keyw);
        } else {
            if (!checkReacSpec(resultMatch[1], id)) {
                errorAbortStr(("Unknown keyword: " + keyw));
            }
            lineKeyword.append("NeutralThirdBody");
            return;
        }
        start = resultMatch[0].second;
    }
}

hur_nodiscard bool
OpenHurricane::chemkinFileRead::checkElementBalance(const List<reacSpcCoeffType> &fwdCoeff,
                                                const List<reacSpcCoeffType> &bwdCoeff) const {
    realList nAtoms(specTable_.elementsNameList().size());
    for (auto &iA : nAtoms) {
        iA = 0;
    }
    for (auto &e : fwdCoeff) {
        const auto &specEle = specTable_[e.index_].elementList();
        for (auto &jSE : specEle) {
            const auto labelI = elementsIndex_.at(jSE.name());
            nAtoms[labelI] += e.stoichCoeff_ * jSE.nAtoms();
        }
    }

    for (auto &e : fwdCoeff) {
        const auto &specEle = specTable_[e.index_].elementList();
        for (auto &jSE : specEle) {
            const auto labelI = elementsIndex_.at(jSE.name());
            nAtoms[labelI] -= e.stoichCoeff_ * jSE.nAtoms();
        }
    }

    for (const auto &iA : nAtoms) {
        if (fabs(iA) > tiny) {
            return false;
        }
    }

    return true;
}

void OpenHurricane::chemkinFileRead::addReactionToList(
    const string &reactionName, const List<reacSpcCoeffType> &forwardCoeff,
    const List<reacSpcCoeffType> &backwardCoeff, const bool isReversible, const reactionType rType,
    const reactionRateType rrType, const fallOffFunctionType fType, void *reacRatePtr,
    void *revReacRatePtr, void *fallOffPtr, const unitTypesOfA iuA, const unitTypesOfE iuE,
    const realArray &effi, void *ArrCoeffPtr, void *revArrCoeffPtr, const integer irc) {

    if (!checkElementBalance(forwardCoeff, backwardCoeff)) {
        errorAbortStr(("The reaction does not balance in elements, in " + reactionName));
    }
    if (reacRatePtr == nullptr) {
        errorAbortStr(("The pointer of reaction rate is null for reaction: " + reactionName));
    }
    if (ArrCoeffPtr == nullptr) {
        errorAbortStr(
            ("The pointer of forward Arrhenius parameter is null for reaction: " + reactionName));
    }
    auto reactionOrder = [](const List<reacSpcCoeffType> &coeff) -> real {
        real sumOrd = 0;
        for (auto &e : coeff) {
            sumOrd += e.order_;
        }
        return sumOrd;
    };

    real Afactor;
    real Efactor;
    changeUnit(reactionOrder(forwardCoeff), iuA, iuE, !(effi.empty()), Afactor, Efactor);
    real AfactorRev = Afactor;
    real EfactorRev = Efactor;

    if (revArrCoeffPtr != nullptr) {
        changeUnit(reactionOrder(backwardCoeff), iuA, iuE, !(effi.empty()), AfactorRev, EfactorRev);
    }
    const auto &ArrheniusCoeff = *(ArrheniusRR *)ArrCoeffPtr;
    auto &reacRate = *(reactionRateTypes *)reacRatePtr;
    constexpr real baseFactor = 1e-3;
    if (rrType == reactionRateType::Arrhenius) {
        if (rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
            const auto &revArrheniusCoeff = *(ArrheniusRR *)revArrCoeffPtr;
            reactions_.append(new reversibleReactionWithRevPars(
                reactionTypes(reactionName, specTable_, reactions_, irc, forwardCoeff,
                              backwardCoeff),
                *(ArrheniusRR(ArrheniusCoeff.A() * Afactor, ArrheniusCoeff.beta(),
                              ArrheniusCoeff.Ta() / Efactor)
                      .clone()),
                *(ArrheniusRR(revArrheniusCoeff.A() * AfactorRev, revArrheniusCoeff.beta(),
                              revArrheniusCoeff.Ta() / EfactorRev))
                     .clone()));
        } else if (rType == reactionType::IRREVERSIBLE) {
            reactions_.append(new irreversibleReactions(
                reactionTypes(reactionName, specTable_, reactions_, irc, forwardCoeff,
                              backwardCoeff),
                *(ArrheniusRR(ArrheniusCoeff.A() * Afactor, ArrheniusCoeff.beta(),
                              ArrheniusCoeff.Ta() / Efactor)
                      .clone())));
        } else if (rType == reactionType::REVERSIBLE) {
            reactions_.append(new reversibleReactions(
                reactionTypes(reactionName, specTable_, reactions_, irc, forwardCoeff,
                              backwardCoeff),
                *(ArrheniusRR(ArrheniusCoeff.A() * Afactor, ArrheniusCoeff.beta(),
                              ArrheniusCoeff.Ta() / Efactor)
                      .clone())));
        } else {
            errorAbortStr(("Unkonwn reaction type for:" + reactionName));
        }
    } else if (rrType == reactionRateType::thirdBodyArrhenius) {
        if (rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
            const auto &revArrheniusCoeff = *(ArrheniusRR *)revArrCoeffPtr;
            reactions_.append(new reversibleReactionWithRevPars(
                reactionTypes(reactionName, specTable_, reactions_, irc, forwardCoeff,
                              backwardCoeff),
                *(ArrheniusThirdBodyRR(baseFactor * ArrheniusCoeff.A() * Afactor,
                                       ArrheniusCoeff.beta(), ArrheniusCoeff.Ta() / Efactor,
                                       specTable_, effi)
                      .clone()),
                *(ArrheniusThirdBodyRR(baseFactor * revArrheniusCoeff.A() * AfactorRev,
                                       revArrheniusCoeff.beta(),
                                       revArrheniusCoeff.Ta() / EfactorRev, specTable_, effi))
                     .clone()));
        } else if (rType == reactionType::IRREVERSIBLE) {
            reactions_.append(new irreversibleReactions(
                reactionTypes(reactionName, specTable_, reactions_, irc, forwardCoeff,
                              backwardCoeff),
                *(ArrheniusThirdBodyRR(baseFactor * ArrheniusCoeff.A() * Afactor,
                                       ArrheniusCoeff.beta(), ArrheniusCoeff.Ta() / Efactor,
                                       specTable_, effi)
                      .clone())));
        } else if (rType == reactionType::REVERSIBLE) {
            reactions_.append(new reversibleReactions(
                reactionTypes(reactionName, specTable_, reactions_, irc, forwardCoeff,
                              backwardCoeff),
                *(ArrheniusThirdBodyRR(baseFactor * ArrheniusCoeff.A() * Afactor,
                                       ArrheniusCoeff.beta(), ArrheniusCoeff.Ta() / Efactor,
                                       specTable_, effi)
                      .clone())));
        } else {
            errorAbortStr(("Unkonwn reaction type for:" + reactionName));
        }
    } else if (rrType == reactionRateType::unimolecularFallOff) {
        auto &uofReacRate = dynamic_cast<unimolecularFallOffRR &>(reacRate);
        auto &k0 = uofReacRate.k0();
        k0.A() *= (baseFactor * Afactor);
        k0.Ta() /= Efactor;
        auto &kinf = uofReacRate.kinf();
        kinf.A() *= Afactor;
        kinf.Ta() /= Efactor;
        uofReacRate.resetFallOff(*(fallOffFunctions::fallOffFunction *)fallOffPtr);
        uofReacRate.resetThirdBodyEff(thirdBodyEfficiency(specTable_, effi));
        if (isReversible) {
            if (revArrCoeffPtr != nullptr || rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
                errorAbortStr(("\"reversibleReactionWithRevParameter\" cannot use with "
                               "\"unimolecularFallOff\" rate in: " +
                               reactionName));
            }

            reactions_.append(
                new reversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                      forwardCoeff, backwardCoeff),
                                        reacRate));
        } else {
            reactions_.append(
                new irreversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                        forwardCoeff, backwardCoeff),
                                          reacRate));
        }
    } else if (rrType == reactionRateType::chemicallyActivatedBimolecular) {
        auto &cheBiReacRate = dynamic_cast<chemicallyActicatedBimolecularRR &>(reacRate);
        auto &k0CB = cheBiReacRate.k0();
        k0CB.A() *= Afactor;
        k0CB.Ta() /= Efactor;
        auto &kinfCB = cheBiReacRate.kinf();
        kinfCB.A() *= (Afactor / baseFactor);
        kinfCB.Ta() /= Efactor;
        cheBiReacRate.resetFallOff(*(fallOffFunctions::fallOffFunction *)fallOffPtr);
        cheBiReacRate.resetThirdBodyEff(thirdBodyEfficiency(specTable_, effi));
        if (isReversible) {
            if (revArrCoeffPtr != nullptr || rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
                errorAbortStr(("\"reversibleReactionWithRevParameter\" cannot use with "
                               "\"chemicallyActivatedBimolecular\" rate in: " +
                               reactionName));
            }

            reactions_.append(
                new reversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                      forwardCoeff, backwardCoeff),
                                        reacRate));
        } else {
            reactions_.append(
                new irreversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                        forwardCoeff, backwardCoeff),
                                          reacRate));
        }
    } else if (rrType == reactionRateType::LandauTeller) {
        auto &LandauTellerCoeffs = dynamic_cast<LandauTellerRR &>(reacRate);
        LandauTellerCoeffs.A() = ArrheniusCoeff.A() * Afactor;
        LandauTellerCoeffs.beta() = ArrheniusCoeff.beta();
        LandauTellerCoeffs.Ta() = ArrheniusCoeff.Ta() / Efactor;
        if (rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
            const auto &revArrheniusCoeff = *(ArrheniusRR *)revArrCoeffPtr;
            auto &revReacRate = *(reactionRateTypes *)revReacRatePtr;
            auto &reverseLandauTellerCoeffs = dynamic_cast<LandauTellerRR &>(revReacRate);
            reverseLandauTellerCoeffs.A() = revArrheniusCoeff.A() * AfactorRev;
            reverseLandauTellerCoeffs.Ta() = revArrheniusCoeff.Ta() / EfactorRev;
            reactions_.append(new reversibleReactionWithRevPars(
                reactionTypes(reactionName, specTable_, reactions_, irc, forwardCoeff,
                              backwardCoeff),
                reacRate, revReacRate));
        } else if (rType == reactionType::IRREVERSIBLE) {
            reactions_.append(
                new irreversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                        forwardCoeff, backwardCoeff),
                                          reacRate));
        } else if (rType == reactionType::REVERSIBLE) {
            reactions_.append(
                new reversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                      forwardCoeff, backwardCoeff),
                                        reacRate));
        } else {
            errorAbortStr(("Unkonwn reaction type for:" + reactionName));
        }
    } else if (rrType == reactionRateType::Janev) {
        auto &JanevRRCoeffs = dynamic_cast<JanevRR &>(reacRate);
        JanevRRCoeffs.A() = ArrheniusCoeff.A() * Afactor;
        JanevRRCoeffs.beta() = ArrheniusCoeff.beta();
        JanevRRCoeffs.Ta() = ArrheniusCoeff.Ta() / Efactor;
        if (isReversible) {
            if (revArrCoeffPtr != nullptr || rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
                errorAbortStr(("\"reversibleReactionWithRevParameter\" cannot use with "
                               "\"Janev\" rate in: " +
                               reactionName));
            }

            reactions_.append(
                new reversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                      forwardCoeff, backwardCoeff),
                                        reacRate));
        } else {
            reactions_.append(
                new irreversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                        forwardCoeff, backwardCoeff),
                                          reacRate));
        }
    } else if (rrType == reactionRateType::collision) {
        auto &collisionRRCoeffs = dynamic_cast<collisionRR &>(reacRate);
        collisionRRCoeffs.kArr().A() = ArrheniusCoeff.A() * Afactor;
        collisionRRCoeffs.kArr().beta() = ArrheniusCoeff.beta();
        collisionRRCoeffs.kArr().Ta() = ArrheniusCoeff.Ta() / Efactor;
        if (isReversible) {
            if (revArrCoeffPtr != nullptr || rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
                errorAbortStr(("\"reversibleReactionWithRevParameter\" cannot use with "
                               "\"collision\" rate in: " +
                               reactionName));
            }

            reactions_.append(
                new reversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                      forwardCoeff, backwardCoeff),
                                        reacRate));
        } else {
            reactions_.append(
                new irreversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                        forwardCoeff, backwardCoeff),
                                          reacRate));
        }
    } else if (rrType == reactionRateType::powerSeries) {
        auto &powerSeriesRRCoeffs = dynamic_cast<powerSeriesRR &>(reacRate);
        powerSeriesRRCoeffs.A() = ArrheniusCoeff.A() * Afactor;
        powerSeriesRRCoeffs.beta() = ArrheniusCoeff.beta();
        if (isReversible) {
            if (revArrCoeffPtr != nullptr || rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
                errorAbortStr(("\"reversibleReactionWithRevParameter\" cannot use with "
                               "\"powerSeries\" rate in: " +
                               reactionName));
            }

            reactions_.append(
                new reversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                      forwardCoeff, backwardCoeff),
                                        reacRate));
        } else {
            reactions_.append(
                new irreversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                        forwardCoeff, backwardCoeff),
                                          reacRate));
        }
    } else if (rrType == reactionRateType::pressureDependenceUsingLogarithmicInterpolation) {
        auto &logInterpolationRRCoeffs = dynamic_cast<logInterpolationRR &>(reacRate);
        auto &kList = logInterpolationRRCoeffs.ks();
        logInterpolationRRCoeffs.ps() *= constant::physicalConstant::Patm;
        for (integer i = 0; i < kList.size(); ++i) {
            kList[i].A() *= Afactor;
            kList[i].Ta() /= Efactor;
        }
        if (!effi.empty()) {
            logInterpolationRRCoeffs.resetThirdBodyEff(thirdBodyEfficiency(specTable_, effi));
        }
        if (isReversible) {
            if (revArrCoeffPtr != nullptr || rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
                errorAbortStr(("\"reversibleReactionWithRevParameter\" cannot use with "
                               "\"pressureDependenceUsingLogarithmicInterpolation\" rate in: " +
                               reactionName));
            }

            reactions_.append(
                new reversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                      forwardCoeff, backwardCoeff),
                                        reacRate));
        } else {
            reactions_.append(
                new irreversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                        forwardCoeff, backwardCoeff),
                                          reacRate));
        }
    } else if (rrType == reactionRateType::ChebyshevReactionRate) {
        auto &ChebyshevPolynomialsRRCoeffs = dynamic_cast<ChebyshevPolynomialsRR &>(reacRate);
        // Change the pressure unit from atm to Pa.
        ChebyshevPolynomialsRRCoeffs.Pmax() *= constant::physicalConstant::Patm;
        ChebyshevPolynomialsRRCoeffs.Pmin() *= constant::physicalConstant::Patm;
        ChebyshevPolynomialsRRCoeffs.restructFromTmpAnm();
        if (!effi.empty()) {
            ChebyshevPolynomialsRRCoeffs.resetThirdBodyEff(thirdBodyEfficiency(specTable_, effi));
        }
        if (isReversible) {
            if (revArrCoeffPtr != nullptr || rType == reactionType::NONEQUILIBRIUMREVERSIBLE) {
                errorAbortStr(("\"reversibleReactionWithRevParameter\" cannot use with "
                               "\"ChebyshevReactionRate\" rate in: " +
                               reactionName));
            }

            reactions_.append(
                new reversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                      forwardCoeff, backwardCoeff),
                                        reacRate));
        } else {
            reactions_.append(
                new irreversibleReactions(reactionTypes(reactionName, specTable_, reactions_, irc,
                                                        forwardCoeff, backwardCoeff),
                                          reacRate));
        }
    } else if (rrType == reactionRateType::unknown) {
        LFatal("Reaction rate type has not been set.");
    }
}

void OpenHurricane::chemkinFileRead::changeUnit(const real reactionOrder, const unitTypesOfA uA,
                                            const unitTypesOfE uE, const bool isThirdBody,
                                            real &Afactor, real &Efactor) const {
    const real RRjoule = constant::physicalConstant::R;             // J/mol-K
    const real RRcal = RRjoule / constant::physicalConstant::calTh; // cal/mol-K

    constexpr real concFactor = 0.001;

    // Change unit for A.
    if (uA == unitTypesOfA::MOLECULES) {
        // First change from MOLECULES to MOLE
        if (isThirdBody) {
            Afactor = pow(constant::physicalConstant::NA, reactionOrder);
        } else {
            Afactor = pow(constant::physicalConstant::NA, reactionOrder - 1);
        }
        // Second, change from mole/cm^3 to kmole/m^3 concentraction units
        Afactor *= pow(concFactor, reactionOrder - 1);
    } else {
        Afactor = pow(concFactor, reactionOrder - 1);
    }

    // Change unit for E.
    switch (uE) {
    case unitTypesOfE::CAL: // Cal/mole
        Efactor = RRcal;
        break;
    case unitTypesOfE::KCAL: // KCAL/MOLE
        Efactor = RRcal / 1000.0;
        break;
    case unitTypesOfE::EVOL: // EVOLTS
        Efactor = RRjoule / constant::physicalConstant::eV;
        break;
    case unitTypesOfE::JOUL: // JOULES/MOLE
        Efactor = RRjoule;
        break;
    case unitTypesOfE::KJOU: // KJOULES/MOLE
        Efactor = RRjoule / 1000.0;
        break;
    case unitTypesOfE::KELV: // KELVINS
        Efactor = 1.0;
        break;
    default:
        break;
    }
}

hur_nodiscard std::string
OpenHurricane::chemkinFileRead::getThermoString(const std::string &chemFileStr) const {
    auto chemkinFileStr = replaceComments(chemFileStr);
    std::regex thermoListRegex("THER(?:|MO)\\s+([\\s\\S]*?)\\s+END", std::regex::icase);
    std::smatch resultMatch;
    std::regex_search(chemkinFileStr, resultMatch, thermoListRegex);
    std::string thermoString = resultMatch[0];
    return thermoString;
}

hur_nodiscard std::string
OpenHurricane::chemkinFileRead::getTransportString(const std::string &chemFileStr) const {
    auto chemkinFileStr = replaceComments(chemFileStr);
    std::regex thransportListRegex("TRAN(?:|SPORT)\\s+([\\s\\S]*?)\\s+END", std::regex::icase);
    std::smatch resultMatch;
    std::regex_search(chemkinFileStr, resultMatch, thransportListRegex);
    std::string transportString = resultMatch[0];
    return transportString;
}

OpenHurricane::chemkinFileRead::chemkinFileRead(const controller &cont, reactionList &rsT,
                                            speciesList &spT, thermoList &ther,
                                            transportList &tran, const std::string &chemFileString,
                                            const std::string &thermoFileString,
                                            const std::string &transportFileString,
                                            const bool noReaction, const bool inviscous)
    : reactions_(rsT), thermo_(ther), transport_(tran), cont_(cont), noReaction_(noReaction),
      inviscous_(inviscous), chemkinFileString_(replaceComments(chemFileString)), reactionString_(),
      thermoFileString_(), transportFileString_(), specTable_(spT), imbalanceTol_(rootTiny),
      globalUnitTypeOfA_(unitTypesOfA::MOLES), globalUnitTypeOfE_(unitTypesOfE::CAL),
      reactionAuxiliaryKeywordMap_(), isotopeAtomicWts_(), elementsIndex_() {
    if (thermoFileString.empty()) {
        thermoFileString_ = getThermoString(chemFileString);
    } else {
        thermoFileString_ = thermoFileString;
    }

    if (transportFileString.empty()) {
        transportFileString_ = getTransportString(chemFileString);
    } else {
        transportFileString_ = transportFileString;
    }

    if (!checkChemFile()) {
        LFatal("The format of chemkin file is wrong");
    }

    getReacAuxKeywordMap();

    if (!noReaction_) {
        reactionString_ = readReactionString();
    }
}

void OpenHurricane::chemkinFileRead::parsing() {
    Pout << "    Info: Parsing elements: " << std::endl;
    parsingElements();
    Pout << "         " << specTable_.elementsNameList().size()
         << " elements including: " << std::endl;
    Pout << "         "
         << "------------------------------------------------------------------"
            "--------------------------------------";
    for (integer i = 0; i < specTable_.elementsNameList().size(); i++) {
        if (i % 8 == 0 || i == 0) {
            Pout << std::endl << "         ";
        }
        Pout << std::left << std::setfill(' ') << std::setw(7) << specTable_.elementsNameList()[i];
    }
    Pout << std::endl
         << "         "
         << "------------------------------------------------------------------"
            "--------------------------------------"
         << std::endl;
    Pout << std::right;
    Pout << "    Info: Parsing species: " << std::endl;
    auto speNList = parsingSpecies();
    Pout << "         " << speNList.size() << " species including: " << std::endl;
    Pout << "         "
         << "------------------------------------------------------------------"
            "--------------------------------------";
    for (integer i = 0; i < speNList.size(); i++) {
        if (i % 8 == 0 || i == 0) {
            Pout << std::endl << "         ";
        }
        Pout << std::left << std::setfill(' ') << std::setw(13) << speNList[i];
    }
    Pout << std::endl
         << "         "
         << "------------------------------------------------------------------"
            "--------------------------------------"
         << std::endl;
    Pout << std::right;
    parsingGlobalUnits();
    Pout << "    Info: Parsing thermo and transport file: " << std::endl;
    JANAFThermoParsing myTher(thermoFileString_, speNList, specTable_, thermo_, isotopeAtomicWts_);
    myTher.parsing();
    if (!inviscous_) {
        JANAFTransportParsing myTran(transportFileString_, speNList, transport_);
        myTran.parsing(cont_.findOrDefault<real>("Prl", 0.72, true));
    }

    if (!noReaction_) {
        Pout << "    Info: Parsing reactions: " << std::endl;
        parsingReaction();
        Pout << "    Reactions: " << reactions_.size() << std::endl;
    }
}