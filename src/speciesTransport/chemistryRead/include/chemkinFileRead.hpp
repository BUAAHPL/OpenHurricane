/*!
 * \file chemkinFileRead.hpp
 * \brief Header of chemkin reaction file parsing.
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

#include "JANAFThermoParsing.hpp"
#include "JANAFTransportParsing.hpp"
#include "controller.hpp"
#include "dataStructure.hpp"
#include "reactionList.hpp"
#include "speciesList.hpp"

#include "reactionTypes.hpp"
#include "thermoList.hpp"
#include "transportList.hpp"

namespace OpenHurricane {

    /*!\breif Class of chemkin file parsing.*/
    class chemkinFileRead {
        using reacSpcCoeffType = reactionTypes::reactionSpeciesCoeffs;

    private:
        // Private data

        /**\brief Reaction table.*/
        reactionList &reactions_;

        /**\brief Thermo table.*/
        thermoList &thermo_;

        /**\brief The transport table.*/
        transportList &transport_;

        const controller &cont_;

        bool noReaction_;

        bool inviscous_;

        /**\brief Chemkin file string.*/
        std::string chemkinFileString_;

        /**\brief Reaction section string.*/
        std::string reactionString_;

        /**\brief Thermo file string.*/
        std::string thermoFileString_;

        /**\brief Transport file string.*/
        std::string transportFileString_;

        /**\brief Hold reference to the species list.*/
        speciesList &specTable_;

        /**\brief Imbalance tolerance.*/
        real imbalanceTol_;

        /**\brief Reaction bcType.*/
        enum class reactionType : short {
            IRREVERSIBLE,
            REVERSIBLE,
            NONEQUILIBRIUMREVERSIBLE,
            UNKNOWN
        };

        static const char *reactionTypeNames[4];

        /**\brief The bcType of reaction rate.*/
        enum class reactionRateType : short {
            Arrhenius,
            thirdBodyArrhenius,
            unimolecularFallOff,
            chemicallyActivatedBimolecular,
            collision,
            Janev,
            LandauTeller,
            powerSeries,
            pressureDependenceUsingLogarithmicInterpolation,
            ChebyshevReactionRate,
            unknown
        };

        static const char *reactionRateTypeNames[11];

        /**\brief The bcType of fall off function.*/
        enum class fallOffFunctionType : short { Lindemann, Troe, SRI, unknown };

        static const char *fallOffFunctionNames[4];

        /**
         * \brief The unit types of pre-exponential factor A.
         * \note If units are not specified, A must be in cm, mole, sec, K.
         */
        enum class unitTypesOfA : short {
            MOLES = 0, /// A is in cm, mole, sec, and K.
            MOLECULES  /// A is in cm, molecules, sec, and K.
        } globalUnitTypeOfA_;

        hur_nodiscard std::map<std::string, unitTypesOfA> unitMapOfA() const;

        /**
         * \brief The unit types of energy E.
         * \note If units are not specified, E must be in cal[/mole]. Chemkin uses thermal calorie, 4.184 Joules.
         */
        enum class unitTypesOfE : short {
            CAL = 0, ///CAL[/MOLE], in calorie per mole
            EVOL,    ///EVOL[TS]
            JOUL,    ///JOUL[ES/MOLE]
            KCAL,    ///KCAL[/MOLE]
            KELV,    ///KELV[INS]
            KJOU     ///KJOU[LES/MOLE]
        } globalUnitTypeOfE_;

        hur_nodiscard std::map<std::string, unitTypesOfE> unitMapOfE() const;

        enum class reactionAuxiliaryKeyword : short {
            NeutralThirdBody = 0, /**< Neutral Third Body Efficiency.*/
            CHEB,                 /**<Chebyshev Polynomial Rate Expressions.*/
            COLLEFF,              /**<Efficiency of Collision Frequency Expression.*/
            DUP,                  /**<Duplicate Reactions.*/
            EXCI,                 /**<Energy Loss Parameter.*/
            FIT1,
            FORD, /**<Forward Reaction Order Parameter.*/
            HIGH, /**<Defines the high-pressure limit for pressure-dependent chemically activated bimolecular reactions.*/
            JAN,  /**<Optional Rate Fit Expressions .*/
            LOW, /**<Defines the low-pressure limit for pressure-dependent unimolecular fall-off reactions .*/
            LT,   /**<Landau-Teller Reactions .*/
            MOME, /**<Plasma Momentum-Transfer Collision Frequency Options .*/
            PCHEB, /**<Supersedes the default pressure limits for a Chebyshev polynomial rate expression .*/
            PLOG, /**<Pressure Dependence Through Logarithmic Interpolation .*/
            REV,  /**<Reverse Rate Parameters .*/
            RLT, /**<Supersedes the default reverse reaction rate expression by the Landau-Teller reaction rate .*/
            RORD, /**<Reverse Reaction Order Parameter.*/
            SRI,  /**<Defines the SRI pressure-dependent reaction rate.*/
            TCHEB, /**<Supersedes the default temperature limits for a Chebyshev polynomial rate expression.*/
            TDEP,  /**<Species Temperature Dependence.*/
            TROE,  /**<Defines the Troe pressure-dependent reaction rate.*/
            UNITS, /**<Reaction Units.*/
            XSMI, /**<Flags a reaction as representing collision cross-section information for the determination of ion momentum-transfer collision frequencies in a plasma simulation.*/
            END
        };

        std::map<std::string, reactionAuxiliaryKeyword> reactionAuxiliaryKeywordMap_;

        /**\brief Isotope molecular weights.*/
        std::map<std::string, real> isotopeAtomicWts_;

        /**\brief The map of elements index.*/
        std::map<std::string, integer> elementsIndex_;

        void getReacAuxKeywordMap();

        /**\brief Check the chemkin file.*/
        hur_nodiscard bool checkChemFile() const;

        /**
         *\brief Replace the comments in file string.
         *\param[in] fileString - The string to replace comments.
         */
        hur_nodiscard std::string replaceComments(const std::string &fileString) const;

        /**\brief Read reaction section string from chemkin file.*/
        hur_nodiscard std::string readReactionString() const;

        void parsingReacStoich(const std::string &forwardStr, List<reacSpcCoeffType> &forwardCoeff,
                               const std::string &backwardStr,
                               List<reacSpcCoeffType> &backwardCoeff,
                               const std::string &reactionName) const;

        /**\brief Parsing reaction.*/
        void parsingReaction();

        void parsingThirdBodyCoeffs(const std::string &thirdBodyLineStr, realArray &efficiens,
                                    bool isPressureDependent, const integer enhancedSpeId) const;

        /**\brief Split reaction species string to species string list.*/
        hur_nodiscard stdStringList splitReacSpecToStrList(const std::string &reacSpecies,
                                                        const std::string &reactionName) const;

        /**
         *\brief Check if the reaction species is within the species list.
         *       And find the index of species.
         */
        bool checkReacSpec(const std::string &rsp, integer &id) const;

        hur_nodiscard std::map<std::string, real>
        parsingReacSpc(const std::string &reacSpc, const std::string &reactionName) const;

        void parsingTroe(void **fofPtr, const string &TroeStr) const;
        void parsingSRI(void **fofPtr, const string &SRIStr) const;

        void parsingReacSpcOrder(const std::string &reacName, const std::string &str,
                                 std::regex &reg, List<reacSpcCoeffType> &coef) const;

        /**
         * \brief Check if the reaction is a pressure-dependent reaction.
         *       The integer id indicates that the idth species is acting as the third body in the fall-off region, not the total concentration M.
         *       If id = -1, then all species act as third bodies.
         *       If id = -2, then error.
         *       An M as a reactant and product surrounded by parentheses indicates that the reaction is a pressure-dependent reaction,
         *       in which case auxiliary information line(s) (described below) must follow the reaction to identify the fall-off formulation and parameters.
         *       A species may also be enclosed in parenthesis. Here, for example, (+H2O) indicates that water is acting as
         *       the third body in the fall-off region, not the total concentration M.
         * \param[in] reacStr - The reaction line string.
         * \param[out] id - The index of the species which is acting as the third body in the fall-off region.
         * \param[out] efficiences - The enhanced third body efficiencies.
         * \return Return true if the reaction is a pressure-dependent reaction.
         */
        bool checkAndReplacePressureDependency(std::string &reacStr, integer &id,
                                               realArray &efficiences) const;

        /**
         * \brief Check if the reaction contains an M as a reactant and product stands for an arbitrary third body.
         *        An M in the reaction description indicates that a third body is participating in the reaction.
         *        In a reaction containing an M, species can be specified to have enhanced third body efficiencies,
         *        in which case auxiliary information (described below) must follow the reaction line.
         *        If no enhanced third body efficiencies are specified, then all species act equally as third bodies
         *        and the effective concentration of the third body is the total concentration of the mixture.
         * \param[in] reacStr - The reaction line string.
         * \param[out] efficiences - The enhanced third body efficiencies.
         * \return Return true if the reaction contains an M as a reactant and product stands for an arbitrary third body.
         */
        bool checkAndReplaceThirdbody(std::string &reacStr, realArray &efficiences) const;

        hur_nodiscard bool checkReactionDuplicate(const List<reacSpcCoeffType> &forwardCoeff,
                                                  const List<reacSpcCoeffType> &backwardCoeff,
                                                  const boolList &isDuplicate,
                                                  const boolList &isPressureDependent,
                                                  const boolList &isThirdBody,
                                                  const integer currentId) const;

        hur_nodiscard bool checkReacSpecList(const List<reacSpcCoeffType> &coeff1,
                                             const List<reacSpcCoeffType> &coeff2) const;

        /**
         * \brief Read elements name list.
         * \return The elements name list.
         */
        stringList parsingElements();

        /**
         * \brief Read species name list.
         * \return The species name list.
         */
        stringList parsingSpecies();

        /**\brief Parsing global units.*/
        void parsingGlobalUnits();

        /**\brief Find the keyword bcType of the line string from the reaction section.*/
        void findLineKeywordType(const std::string &lineStr, stringList &lineKeyword,
                                 integerList &pos) const;

        hur_nodiscard bool checkElementBalance(const List<reacSpcCoeffType> &fwdCoeff,
                                               const List<reacSpcCoeffType> &bwdCoeff) const;
        void addReactionToList(const string &reactionName, const List<reacSpcCoeffType> &forwardCoeff,
                               const List<reacSpcCoeffType> &backwardCoeff, const bool isReversible,
                               const reactionType rType, const reactionRateType rrType,
                               const fallOffFunctionType fType, void *reacRatePtr,
                               void *revReacRatePtr, void *fallOffPtr, const unitTypesOfA iuA,
                               const unitTypesOfE iuE, const realArray &effi, void *ArrCoeffPtr,
                               void *revArrCoeffPtr, const integer irc);

        void changeUnit(const real reactionOrder, const unitTypesOfA uA, const unitTypesOfE uE,
                        const bool isThirdBody, real &Afactor, real &Efactor) const;

        hur_nodiscard std::string getThermoString(const std::string &chemFileStr) const;
        hur_nodiscard std::string getTransportString(const std::string &chemFileStr) const;

    public:
        // Constructors

        /**\brief Disallow null constructor.*/
        chemkinFileRead() = delete;

        /**\brief Construct from components.*/
        chemkinFileRead(const controller &cont, reactionList &rsT, speciesList &spT,
                        thermoList &ther, transportList &tran, const std::string &chemFileString,
                        const std::string &thermoFileString = nullString,
                        const std::string &transportFileString = nullString,
                        const bool noReaction = false, const bool inviscous = false);

        /**\brief Disallow copy constructor.*/
        chemkinFileRead(const chemkinFileRead &) = delete;

        inline ~chemkinFileRead() noexcept {}

        /**\brief Parsing the chemkin file.*/
        void parsing();

        /*!\brief Reaction table.*/
        hur_nodiscard inline reactionList &reactions() { return reactions_; }

        // Operators

        /**\brief Disallow bitwise assignment.*/
        void operator=(const chemkinFileRead &) = delete;
    };

} // namespace OpenHurricane