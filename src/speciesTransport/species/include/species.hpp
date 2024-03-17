/*!
 * \file species.hpp
 * \brief All the information about the definition of the species.
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
#include "OpenHurricane.hpp"
#include "speciesElements.hpp"

namespace OpenHurricane {
    class species;

    /*!\brief The single species: air (mixture)*/
    extern species air;

    /**
     * \brief The class of species.
     */
    class species {
    private:
        /** \brief The name of specie*/
        string name_;

        /** \brief Molecular weight of specie [g/mol or kg/kmol]*/
        mutable real molecularWeight_;

        /** \brief Gas constant [J/(kg K)]*/
        mutable real Ri_;

        /** \brief The list of elements of the specie*/
        speciesElementsList elementList_;

    public:
        species();

        /**
         * \brief Construct from name and molecular weight.
         * For example, for air.
         * \param[in] names - The molecular name.
         * \param[in] molWeight - The molecular weight [g/mol or kg/kmol].
         * \param[in] Rr - The gas constant [J/(kg K)].
         */
        species(const string &names, const real molWeight, const real Rr);

        /**
         * \brief Construct from name and molecular weight.
         * For example, for air.
         * \param[in] names - The molecular name.
         * \param[in] molWeight - The molecular weight [g/mol or kg/kmol].
         */
        species(const string &names, const real molWeight);

        /**
         * \brief Construct from components.
         * \param[in] names - The molecular name.
         * \param[in] molWeight - The molecular weight [g/mol or kg/kmol].
         * \param[in] eleList - The list of elements of the specie.
         */
        species(const string &name, const real molWeight, const speciesElementsList &eleList);

        /**
         * \brief Construct from components without molecular weight.
         * \param[in] names - The molecular name.
         * \param[in] eleList - The list of elements of the specie.
         */
        species(const string &name, const speciesElementsList &eleList);

        /** \brief Construct as copy*/
        species(const species &sp);

        /**
         * \brief Construct as copy by given new specie name.
         */
        species(const species &sp, const string &name);

        /**\brief Assignment to the given species.*/
        species &operator=(const species &);

        /*!\brief Destructor.*/
        inline ~species() noexcept {}

        /** \brief Return the name of the specie*/
        hur_nodiscard inline const string &name() const noexcept { return name_; }
        hur_nodiscard inline string &name() noexcept { return name_; }

        /** \brief Return the list of the elements of the specie*/
        hur_nodiscard inline const speciesElementsList &elementList() const noexcept {
            return elementList_;
        }
        hur_nodiscard inline speciesElementsList &elementList() noexcept { return elementList_; }

        /** \brief Molecular weight [g/mol or kg/kmol]*/
        hur_nodiscard inline real W() const noexcept {
            if (molecularWeight_ == 0.0) {
                getMolWeight();
            }
            return molecularWeight_;
        }

        /*!\brief Change the molecular weight [g/mol or kg/kmol].*/
        inline void changeW(const real newMW) noexcept {
            molecularWeight_ = newMW;
            Ri_ = constant::physicalConstant::Ru / molecularWeight_;
        }

        /** \brief Gas constant [J/(kg K)]*/
        hur_nodiscard inline real Ri() const noexcept { return Ri_; }

        /** \brief get molecular weight [g/mol or kg/kmol]*/
        void getMolWeight() const;

        /**
         * \brief If found the element by the given name.
         * \param[in] elementName - The name of the element.
         * \param[out] nAtoms - The number of atoms of the element.
         * \param[OUT] elementWeight - tHE Atomic weight of the given element [g/mol or kg/kmol].
         * \return True if found.
         */
        bool foundElement(const string &elementName, integer &nAtoms, real &elementWeight) const;
    };

} // namespace OpenHurricane
