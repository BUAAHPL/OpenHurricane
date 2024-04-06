/*!
 * \file speciesList.hpp
 * \brief All the information about the definition of the specie list.
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

#include "cellArrays.hpp"
#include "realArray.hpp"
#include "species.hpp"

namespace OpenHurricane {
    /**
     * \brief The class of species list.
     */
    class speciesList : public List<species> {
    public:
        using Base = List<species>;

    private:
        /*!\brief Elements name list*/
        stringList elementsNameList_;

    public:
        hur_nodiscard static inline const speciesList &nullObject() {
            return NullRefObj::nullRef<speciesList>();
        }

        /**\brief Construct null.*/
        inline speciesList();

        /**\brief Construct by given size.*/
        inline speciesList(const integer _size);

        /**\brief Construct as copy.*/
        inline speciesList(const speciesList &sT);

        /**\brief Construct as copy.*/
        inline speciesList(speciesList &&sT) noexcept;

        /**\brief Destructor.*/
        inline ~speciesList() noexcept {}

        /**\brief Return the name list of species.*/
        hur_nodiscard string nameList() const;

        /**\brief Return the elements name list*/
        hur_nodiscard const stringList &elementsNameList() const noexcept;

        hur_nodiscard stringList &elementsNameList() noexcept;

        /**\brief Return the name of species i.*/
        hur_nodiscard inline const string &name(const integer specieI) const noexcept;

        /**\brief Return the molecular weight of species i.*/
        hur_nodiscard inline real W(const integer specieI = 0) const noexcept;

        /** \brief Gas constant [J/(kg K)] of species i.*/
        hur_nodiscard inline real Ri(const integer specieI = 0) const noexcept;

        /**\brief Does the species list contain this species?*/
        hur_nodiscard bool contains(const string &_spName) const;

        /**
         * \brief Does the species list contain this species?
         * If true, return the index of this species in the table.
         * If false, _id = -1.
         */
        bool contains(const string &_spName, integer &_id) const;

        /**
         * \brief Return the index of this species in the table by given name.
         * If not exit, throw a fatal error.
         */
        hur_nodiscard integer index(const string &_spName) const;

        /*!
         * \brief Return the gas constant of the mixture specified by this species list [J/kg/K].
         * \param[in] yi - the mass fraction.
         */
        hur_nodiscard inline real Rm(const realArray &yi) const noexcept;

        /*!
         * \brief Return the gas constant of the mixture specified by this species list.
         * \param[in] yi - the mass fraction.
         */
        hur_nodiscard real Rm(const PtrList<cellRealArray> &yi, const integer cellI) const noexcept;

        /*!\brief Return the mean molecular weight by given molar fraction*/
        hur_nodiscard real MWbyXi(const realArray &xi) const noexcept;

        /*!\brief Return the mean molecular weight by given mass fraction*/
        hur_nodiscard real MWbyYi(const realArray &yi) const noexcept;

        /*!\brief Return the mean molecular weight by given mass fraction*/
        hur_nodiscard real MWbyYi(const PtrList<cellRealArray> &yi,
                                  const integer cellI) const noexcept;

        /*!\brief Change mass fraction yi to molar fraction xi.*/
        hur_nodiscard realArray Yi2Xi(const realArray &yi) const;

        /*!\brief Change mass fraction yi to molar fraction xi.*/
        hur_nodiscard realArray Yi2Xi(const PtrList<cellRealArray> &yi, const integer cellI) const;

        /*!\brief Change mass fraction yi to molar fraction xi.*/
        void Yi2Xi(const PtrList<cellRealArray> &yi, const integer cellI, realArray &xi) const;

        /*!\brief Change mass fraction yi to molar fraction xi.*/
        void Yi2Xi(const PtrList<cellRealArray> &yi, PtrList<cellRealArray> &xi,
                   const bool isOnlyInternal = true) const;

        /*!\brief Change mass fraction yi to molar fraction xi.*/
        void Yi2Xi(const realArray &yi, realArray &xi) const;

        /*!\brief Change molar fraction xi to mass fraction yi.*/
        hur_nodiscard realArray Xi2Yi(const realArray &xi) const;

        /*!\brief Change molar fraction xi to mass fraction yi.*/
        void Xi2Yi(const realArray &xi, realArray &yi) const;

        /**
         * \brief To get the element mass fraction.
         * \param[in] elementName - The name of the element.
         * \param[in] yi - The species mass fraction array.
         * \return The element mass fraction.
         * \retval A real value.
         */
        hur_nodiscard real Zj(const string &elementName, const realArray &yi) const;

        /**
         * \brief To get the element mass fraction.
         * \param[in] elementName - The name of the element.
         * \param[in] yi - The species mass fraction array.
         * \param[in] cellI - The index of cell.
         * \return The element mass fraction.
         * \retval A real value.
         */
        hur_nodiscard real Zj(const string &elementName, const PtrList<cellRealArray> &yi,
                              const integer cellI) const;

        /**
         * \brief The mixture fraction given by:
         * \f[Z = \frac{{s{Y_F} - {Y_{{O_2}}} + {Y_{{O_2},2}}}}{{s{Y_{F,1}} + {Y_{{O_2},2}}}}\f]
         * \param[in] s - The stoichiometric ratio.
         * \param[in] YF - The mass fraction of fuel.
         * \param[in] YO2 - The mass fraction of oxidizer.
         * \param[in] YF1 - The mass fraction of fule in the fuel stream.
         * \param[in] YO22 - The mass fraction of oxidizer in the oxidizer stream.
         */
        hur_nodiscard inline real Z1(const real s, const real YF, const real YO2, const real YF1,
                                     const real YO22) const noexcept;
    };

} // namespace OpenHurricane

#include "speciesList.inl"