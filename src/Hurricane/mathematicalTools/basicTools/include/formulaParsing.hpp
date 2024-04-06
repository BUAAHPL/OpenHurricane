/*!
 * \file formulaParsing.hpp
 * \brief Header of data formula parsing.
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
#pragma once
#include "OpenHurricane.hpp"

namespace OpenHurricane {

    namespace formulaParsingSection {}

    /*!
     * \breif The class of formula Parsing
     */
    class formulaParsing {
    public:
        /**
         * \brief The enum of arithmetic operator type.
         */
        enum class operatorType : short {
            unary = 0,  //!< Unary operator
            binary,     //!< Binary operator
            operand,    //!< Operand
            leftParen,  //!< Left paren
            rightParen, //!< Right paren
            end         //!< Ending of formula parsing
        };

        enum class invalidType : short {
            valid = 0,             /**< Valid. */
            invalidChar = 1,       /**< Invalid character. */
            mismatchedParen = 2,   /**< Mismatched paren. */
            invalidDecimal = 3,    /**< Invalid decimal point. */
            mismatchedNumChar = 4, /**< Mismatched number and character. */
            empty = 5,             /**< Empty formular. */
            unknown
        };

        enum class validOperators : short {
            leftPare,  //!< Left paren
            rightPare, //!< right paren
            nagetive,  //!< Unary nagetive operator
            abs,
            sqrt,
            exp,
            ln,
            log10,
            sin,
            cos,
            tanh,
            plus,
            minus,
            times,
            divide,
            mod,
            power,
            Operand,
            unknown
        };

        static std::map<std::string, validOperators> validOperatorMaps_;

        struct operatorTypePriority {
            operatorType operatorType_;
            short priority_;

            inline operatorTypePriority(operatorType operatorType, short priority)
                : operatorType_(operatorType), priority_(priority) {}
        };

        static std::map<validOperators, operatorTypePriority> operatorPriority_;
        /**
         * \brief The class of arithmetic identifier.
         */
        class identifierType {
        public:
            /**\brief The name of arithmetic identifier.*/
            std::string name_;

            /**\brief The type of arithmetic identifier.*/
            operatorType operatorType_;

            validOperators operator_;

            /**\brief The value of arithmetic identifier.*/
            real val_;

            short priority_;
        };

        //! \cond internalClass
        class initValidOperatorMaps {
        public:
            initValidOperatorMaps();
        };

        static initValidOperatorMaps dummyInitValidOperatorMaps_;

    private:
        /**\brief Formula expression.*/
        std::string formulaStr_;

        /**\brief Priority_ index.*/
        integerList priorityIndex_;

        /**\brief The value of formula.*/
        real value_;

        /**\brief Is valid?*/
        invalidType isValid_;

        List<identifierType> currentIdentifiers_;

    public:
        inline formulaParsing()
            : formulaStr_(), priorityIndex_(), value_(), isValid_(invalidType::valid),
              currentIdentifiers_() {
            priorityIndex_.reserve(100);
            currentIdentifiers_.reserve(100);
        }

        inline formulaParsing(const std::string &formularStr)
            : formulaStr_(formularStr), priorityIndex_(), value_(), isValid_(invalidType::valid),
              currentIdentifiers_() {
            checkValid();
            priorityIndex_.reserve(100);
            currentIdentifiers_.reserve(100);
        }

        inline ~formulaParsing() noexcept {}

        void clear() noexcept;

        hur_nodiscard inline bool isValid() const noexcept {
            return isValid_ == invalidType::valid;
        }

        hur_nodiscard inline bool isInvalid() const noexcept {
            return isValid_ != invalidType::valid;
        }

    protected:
        hur_nodiscard inline bool isValidChar(const char ch) const {
            return isdigit(ch) || isalpha(ch) || isOperator(ch) || (ch == '.');
        }

        hur_nodiscard inline bool isOperator(const char ch) const {
            return (ch == '+') || (ch == '-') || (ch == '/') || (ch == '*') || (ch == '=') ||
                   (ch == '%') || (ch == '!') || (ch == '^') || (ch == '(') || (ch == ')');
        }

        void checkValid();

        void preParsing();
        void parsingDigit(const std::string &str, std::string::iterator &iter);
        void parsingChar(const std::string &str, std::string::iterator &iter);

        bool checkOperator(const std::string &str);

        void adjust();

        void errorAndExit() const;

        void runFormular();
    public:
        hur_nodiscard real parsing();
    };

} //  namespace OpenHurricane
