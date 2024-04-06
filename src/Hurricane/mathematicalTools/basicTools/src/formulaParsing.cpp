/*!
 * \file formulaParsing.cpp
 * \brief Main subroutines for formula parsing.
 * \author Peng Jian
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
#include "formulaParsing.hpp"
#include <forward_list>

std::map<std::string, OpenHurricane::formulaParsing::validOperators>
    OpenHurricane::formulaParsing::validOperatorMaps_;

std::map<OpenHurricane::formulaParsing::validOperators,
         OpenHurricane::formulaParsing::operatorTypePriority>
    OpenHurricane::formulaParsing::operatorPriority_;

OpenHurricane::formulaParsing::initValidOperatorMaps::initValidOperatorMaps() {
    validOperatorMaps_.emplace("(", validOperators::leftPare);
    validOperatorMaps_.emplace(")", validOperators::rightPare);
    validOperatorMaps_.emplace("~", validOperators::nagetive);
    validOperatorMaps_.emplace("abs", validOperators::abs);
    validOperatorMaps_.emplace("sqrt", validOperators::sqrt);
    validOperatorMaps_.emplace("exp", validOperators::exp);
    validOperatorMaps_.emplace("ln", validOperators::ln);
    validOperatorMaps_.emplace("log10", validOperators::log10);
    validOperatorMaps_.emplace("sin", validOperators::sin);
    validOperatorMaps_.emplace("cos", validOperators::cos);
    validOperatorMaps_.emplace("tanh", validOperators::tanh);
    validOperatorMaps_.emplace("+", validOperators::plus);
    validOperatorMaps_.emplace("-", validOperators::minus);
    validOperatorMaps_.emplace("*", validOperators::times);
    validOperatorMaps_.emplace("/", validOperators::divide);
    validOperatorMaps_.emplace("%", validOperators::mod);
    validOperatorMaps_.emplace("^", validOperators::power);

    operatorPriority_.emplace(validOperators::leftPare,
                              operatorTypePriority(operatorType::leftParen, 0));
    operatorPriority_.emplace(validOperators::rightPare,
                              operatorTypePriority(operatorType::rightParen, 0));
    operatorPriority_.emplace(validOperators::nagetive,
                              operatorTypePriority(operatorType::unary, 3));
    operatorPriority_.emplace(validOperators::abs, operatorTypePriority(operatorType::unary, 3));
    operatorPriority_.emplace(validOperators::sqrt, operatorTypePriority(operatorType::unary, 3));
    operatorPriority_.emplace(validOperators::exp, operatorTypePriority(operatorType::unary, 3));
    operatorPriority_.emplace(validOperators::ln, operatorTypePriority(operatorType::unary, 3));
    operatorPriority_.emplace(validOperators::log10, operatorTypePriority(operatorType::unary, 3));
    operatorPriority_.emplace(validOperators::sin, operatorTypePriority(operatorType::unary, 3));
    operatorPriority_.emplace(validOperators::cos, operatorTypePriority(operatorType::unary, 3));
    operatorPriority_.emplace(validOperators::tanh, operatorTypePriority(operatorType::unary, 3));
    operatorPriority_.emplace(validOperators::plus, operatorTypePriority(operatorType::binary, 1));
    operatorPriority_.emplace(validOperators::minus, operatorTypePriority(operatorType::binary, 1));
    operatorPriority_.emplace(validOperators::times, operatorTypePriority(operatorType::binary, 2));
    operatorPriority_.emplace(validOperators::divide,
                              operatorTypePriority(operatorType::binary, 2));
    operatorPriority_.emplace(validOperators::mod, operatorTypePriority(operatorType::binary, 2));
    operatorPriority_.emplace(validOperators::power, operatorTypePriority(operatorType::unary, 3));
}

OpenHurricane::formulaParsing::initValidOperatorMaps
    OpenHurricane::formulaParsing::dummyInitValidOperatorMaps_;

void OpenHurricane::formulaParsing::clear() noexcept {
    formulaStr_.clear();
    priorityIndex_.clear();
    value_ = 0;
    isValid_ = invalidType::valid;
}

void OpenHurricane::formulaParsing::checkValid() {

    replaceAllMarks(formulaStr_, " ", "");
    replaceAllMarks(formulaStr_, "\n", "");
    replaceAllMarks(formulaStr_, "\r", "");

    if (formulaStr_.empty()) {
        isValid_ = invalidType::empty;
        return;
    }
    std::forward_list<char> eleL;
    for (size_t m = 0; m < formulaStr_.size(); ++m) {
        if (!isValidChar(formulaStr_[m])) {
            isValid_ = invalidType::invalidChar;
            return;
        }
        if (formulaStr_[m] == '(') {
            eleL.emplace_front('(');
        } else if (formulaStr_[m] == ')') {
            if (!eleL.empty()) {
                eleL.pop_front();
            } else {
                isValid_ = invalidType::mismatchedParen;
                return;
            }
        }
    }
    if (!eleL.empty()) {
        isValid_ = invalidType::mismatchedParen;
    }
}

void OpenHurricane::formulaParsing::preParsing() {
    errorAndExit();
    for (std::string::iterator iter = formulaStr_.begin(); iter != formulaStr_.end();) {
        if (isdigit(*iter)) {
            parsingDigit(formulaStr_, iter);
            if (isInvalid()) {
                return;
            }
        } else if (isalpha(*iter)) {
            parsingChar(formulaStr_, iter);
            if (isInvalid()) {
                return;
            }
        } else if (*iter == '.') {
            isValid_ = invalidType::invalidDecimal;
            return;
        } else {
            std::string str(iter, iter + 1);
            if (str == "-") {
                if (currentIdentifiers_.empty()) {
                    str = "~";
                } else if (currentIdentifiers_.last().operator_ != validOperators::Operand) {
                    if (currentIdentifiers_.last().operator_ != validOperators::rightPare) {
                        str = "~";
                    }
                }
            } else if (str == "+") {
                if (currentIdentifiers_.empty() ||
                    (!currentIdentifiers_.empty() &&
                     currentIdentifiers_.last().operator_ == validOperators::leftPare)) {
                    ++iter;
                    continue;
                }
            }
            if (!checkOperator(str)) {
                return;
            }
            identifierType idet;
            idet.name_ = str;
            auto citer = validOperatorMaps_.find(str);
            idet.val_ = 0;
            if (citer != validOperatorMaps_.end()) {
                auto &op = operatorPriority_.at(citer->second);
                idet.operator_ = citer->second;
                idet.operatorType_ = op.operatorType_;
                idet.priority_ = op.priority_;
            } else {
                idet.operator_ = validOperators::unknown;
                idet.priority_ = 0;
                idet.operatorType_ = operatorType::end;
                isValid_ = invalidType::unknown;
            }
            currentIdentifiers_.append(idet);
            ++iter;
        }
    }
}

void OpenHurricane::formulaParsing::parsingDigit(const std::string &str,
                                                 std::string::iterator &iter) {
    std::string::iterator endIter = iter;
    while (endIter != str.end() && ((*endIter) == '.' || isdigit(*endIter))) {
        ++endIter;
    }
    if (*(endIter - 1) == '.') {
        isValid_ = invalidType::invalidDecimal;
        return;
    }
    std::string name(iter, endIter);
    identifierType idet;
    idet.name_ = name;
    idet.operatorType_ = operatorType::operand;
    idet.val_ = atof(name.c_str());
    idet.priority_ = 0;
    idet.operator_ = validOperators::Operand;
    currentIdentifiers_.append(idet);
    iter = endIter;
}

void OpenHurricane::formulaParsing::parsingChar(const std::string &str,
                                                std::string::iterator &iter) {
    std::string::iterator endIter = iter;
    while (endIter != str.end() && isalpha(*endIter)) {
        ++endIter;
    }

    std::string name(iter, endIter);
    stringToLowerCase(name);

    identifierType idet;
    idet.name_ = name;
    auto citer = validOperatorMaps_.find(name);
    idet.val_ = 0;
    if (citer != validOperatorMaps_.end()) {
        auto &op = operatorPriority_.at(citer->second);
        idet.operator_ = citer->second;
        idet.operatorType_ = op.operatorType_;
        idet.priority_ = op.priority_;
    } else {
        idet.operator_ = validOperators::unknown;
        idet.priority_ = 0;
        idet.operatorType_ = operatorType::end;
        isValid_ = invalidType::unknown;
    }
    currentIdentifiers_.append(idet);
    iter = endIter;
}

bool OpenHurricane::formulaParsing::checkOperator(const std::string &str) {
    const auto &iter = validOperatorMaps_.find(str);
    if (iter != validOperatorMaps_.end()) {
        const auto &op = operatorPriority_.at(iter->second);
        if (currentIdentifiers_.empty() && op.operatorType_ == operatorType::binary) {
            isValid_ = invalidType::mismatchedNumChar;
            return false;
        } else if (op.operatorType_ == operatorType::unary) {
            if (!currentIdentifiers_.empty() &&
                (currentIdentifiers_.last().operatorType_ != operatorType::binary)) {
                isValid_ = invalidType::mismatchedNumChar;
                return false;
            }
        } else if (op.operatorType_ == operatorType::binary ||
                   op.operatorType_ == operatorType::rightParen) {
            const auto &lci = currentIdentifiers_.last();
            if (lci.operatorType_ != operatorType::operand &&
                lci.operatorType_ != operatorType::rightParen) {
                isValid_ = invalidType::mismatchedNumChar;
                return false;
            }
        } else if (op.operatorType_ == operatorType::leftParen) {
            const auto &lci = currentIdentifiers_.last();
            if (lci.operatorType_ == operatorType::operand ||
                lci.operatorType_ == operatorType::rightParen) {
                isValid_ = invalidType::mismatchedNumChar;
                return false;
            }
        }
        return true;
    }
    isValid_ = invalidType::unknown;
    return false;
}

void OpenHurricane::formulaParsing::adjust() {
    errorAndExit();
    integer ef = 0;
    std::forward_list<integer> opp;
    for (integer i = 0; i < currentIdentifiers_.size(); ++i) {
        const auto &e = currentIdentifiers_[i];
        if (e.operatorType_ == operatorType::operand) {
            priorityIndex_.push_back(i);
        } else {
            if (opp.empty()) {
                opp.push_front(i);
            } else {
                ef = opp.front();

                if (e.operatorType_ == operatorType::leftParen ||
                    e.operatorType_ == operatorType::rightParen) {
                    priorityIndex_.push_front(i);
                    continue;
                }
                if (e.priority_ <= currentIdentifiers_[ef].priority_) {
                    ef = opp.front();
                    opp.pop_front();
                    priorityIndex_.push_back(ef);
                    --i;
                    continue;
                } else {
                    opp.push_front(i);
                }
            }
        }
    }
    while (!opp.empty()) {
        ef = opp.front();
        opp.pop_front();
        priorityIndex_.push_back(ef);
    }
}

void OpenHurricane::formulaParsing::errorAndExit() const {
    std::string errMsg;
    if (isValid_ == invalidType::valid) {
        return;
    } else if (isValid_ == invalidType::unknown) {
        errMsg = "unknown operators";
    } else if (isValid_ == invalidType::invalidChar) {
        errMsg = "invalid character ";
    } else if (isValid_ == invalidType::mismatchedParen) {
        errMsg = "mismatched parenthesis";
    } else if (isValid_ == invalidType::empty) {
        errMsg = "empty formular";
    } else if (isValid_ == invalidType::invalidDecimal) {
        errMsg = "invalid decimal point";
    } else if (isValid_ == invalidType::mismatchedNumChar) {
        errMsg = "mismatched number and character";
    }
    LFatal("Invalid formular: %s for %s", formulaStr_.c_str(), errMsg.c_str());
}

void OpenHurricane::formulaParsing::runFormular() {
    real re;
    real re1;
    real re2;
    std::forward_list<real> results;
    for (integer i = 0; i < priorityIndex_.size(); ++i) {
        const auto m = priorityIndex_[i];
        const auto op = currentIdentifiers_[m].operator_;
        switch (op) {
        case validOperators::leftPare:
        case validOperators::rightPare:
            break;
        case validOperators::nagetive:
            re = results.front();
            results.pop_front();
            results.push_front(-re);
            break;
        case validOperators::abs:
            re = results.front();
            results.pop_front();
            results.push_front(abs(re));
            break;
        case validOperators::sqrt:
            re = results.front();
            results.pop_front();
            if (re < 0) {
                LFatal("Negative real number");
            }
            results.push_front(sqrt(re));
            break;
        case validOperators::exp:
            re = results.front();
            results.pop_front();
            results.push_front(exp(re));
            break;
        case validOperators::ln:
            re = results.front();
            results.pop_front();
            if (re < 0) {
                LFatal("Negative real number");
            }
            results.push_front(log(re));
            break;
        case validOperators::log10:
            re = results.front();
            results.pop_front();
            if (re < 0) {
                LFatal("Negative real number");
            }
            results.push_front(log10(re));
            break;
        case validOperators::sin:
            re = results.front();
            results.pop_front();
            results.push_front(sin(re));
            break;
        case validOperators::cos:
            re = results.front();
            results.pop_front();
            results.push_front(cos(re));
            break;
        case validOperators::tanh:
            re = results.front();
            results.pop_front();
            results.push_front(tanh(re));
            break;
        case validOperators::plus:
            re = results.front();
            results.pop_front();
            re2 = re;
            re = results.front();
            results.pop_front();
            re1 = re;
            results.push_front(re1 + re2);
            break;
        case validOperators::minus:
            re = results.front();
            results.pop_front();
            re2 = re;
            re = results.front();
            results.pop_front();
            re1 = re;
            results.push_front(re1 - re2);
            break;
        case validOperators::times:
            re = results.front();
            results.pop_front();
            re2 = re;
            re = results.front();
            results.pop_front();
            re1 = re;
            results.push_front(re1 * re2);
            break;
        case validOperators::divide:
            re = results.front();
            results.pop_front();
            re2 = re;
            if (re2 == 0) {
                LFatal("The denominator can not be zero");
            }
            re = results.front();
            results.pop_front();
            re1 = re;
            results.push_front(re1 / re2);
            break;
        case validOperators::mod:
            re = results.front();
            results.pop_front();
            re2 = re;
            if (re2 != int(re2)) {
                LFatal("The denominator can not be zero");
            }
            re = results.front();
            results.pop_front();
            re1 = re;
            results.push_front(integer(re1) % integer(re2));
            break;
        case validOperators::power:
            re = results.front();
            results.pop_front();
            re2 = re;

            re = results.front();
            results.pop_front();
            re1 = re;
            results.push_front(pow(re1, re2));
            break;
        default:
            results.push_front(currentIdentifiers_[m].val_);
            break;
        }
    }
    value_ = results.front();
}

hur_nodiscard OpenHurricane::real OpenHurricane::formulaParsing::parsing() {
    checkValid();
    preParsing();
    adjust();
    runFormular();
    return value_;
}
