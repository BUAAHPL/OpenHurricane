/*!
 * \file rhoThermo.hpp
 * \brief Headers of base class of rho thermo.
 *        The subroutines and functions are in the <i>rhoThermo.cpp</i> file.
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
#include "cellArrays.hpp"
#include "getBoundariesFromController.hpp"
#include "mixture.hpp"
#include "runtimeMesh.hpp"

namespace OpenHurricane {
    /*!\brief The class of rho thermo*/
    class rhoThermo {
    private:
        /*!\brief The pressure field.*/
        mutable cellRealArray p_;

        /*!\brief The temperature field.*/
        mutable cellRealArray T_;

        /*!\brief The density field.*/
        mutable cellRealArray rho_;

        /*!\brief The energy field.*/
        mutable cellRealArray E_;

        /*!\brief Ratio of specific heat capacities.*/
        mutable cellRealArray gamma_;

        /**\brief The heat capacity at constant pressure.*/
        mutable cellRealArray cp_;

        uniquePtr<mixture> mixPtr_;

    public:
        /*!\brief Construct from mesh*/
        rhoThermo(const runtimeMesh &mesh);

        /*!\brief Construct from mesh*/
        rhoThermo(const runtimeMesh &mesh, const controller &cont, const bool inviscous = false);

        /*!\brief Destructor.*/
        virtual ~rhoThermo() noexcept { mixPtr_.clear(); }

        inline const runtimeMesh &mesh() const noexcept;

        /*!\brief Return const-reference to the density field.*/
        virtual hur_nodiscard const cellRealArray &rho() const noexcept;

        /*!\brief Return access to the density field.*/
        virtual hur_nodiscard cellRealArray &rho() noexcept;

        /*!\brief Return const-reference to the pressure field.*/
        virtual hur_nodiscard const cellRealArray &p() const noexcept;

        /*!\brief Return access to the pressure field.*/
        virtual hur_nodiscard cellRealArray &p() noexcept;

        /*!\brief Return const-reference to the temperature field.*/
        virtual hur_nodiscard const cellRealArray &T() const noexcept;

        /*!\brief Return access to the temperature field.*/
        virtual hur_nodiscard cellRealArray &T() noexcept;

        /*!\brief Return const-reference to the energy field.*/
        virtual hur_nodiscard const cellRealArray &E() const noexcept;

        /*!\brief Return access to the energy field.*/
        virtual hur_nodiscard cellRealArray &E() noexcept;

        /*!\brief Ratio of specific heat capacities.*/
        virtual hur_nodiscard const cellRealArray &gamma() const noexcept;

        /*!\brief Ratio of specific heat capacities.*/
        virtual hur_nodiscard cellRealArray &gamma() noexcept;

        /**\brief The heat capacity at constant pressure.*/
        hur_nodiscard cellRealArray &cp() noexcept;

        /**\brief The heat capacity at constant pressure.*/
        hur_nodiscard const cellRealArray &cp() const noexcept;

        /*!\brief Ratio of specific heat capacities.*/
        virtual hur_nodiscard const cellRealArray &mu() const noexcept;

        /*!\brief Ratio of specific heat capacities.*/
        virtual hur_nodiscard cellRealArray &mu() noexcept;

        /*!\brief Ratio of thermo conductivity.*/
        virtual hur_nodiscard const cellRealArray &kappa() const noexcept;

        /*!\brief Ratio of thermo conductivity.*/
        virtual hur_nodiscard cellRealArray &kappa() noexcept;

        /*!\brief Return access to mixture pointer.*/
        hur_nodiscard mixture &mixtures() noexcept;

        /*!\brief Return true if mixtures_ is a null pointer.*/
        hur_nodiscard inline bool isMixturesPtrNUll() const noexcept;

        /*!\brief Return True if mixture is a single specie.*/
        hur_nodiscard inline bool isSingleSpecie() const noexcept;

        /*
         * \brief Compute temperature by given specific enthalpy.
         * \param[in] Ha - Specific absolute enthalpy (dimensionless).
         * \param[in] v - Velocity (dimensionless).
         * \param[in] rho - Density (dimensionless).
         * \param[out] T - Temperature (dimensionless).
         */
        void THa(const cellRealArray &ha, const cellRealArray &rho, cellRealArray &T) const;

        /*
         * \brief Compute temperature by given specific enthalpy with flags.
         * \param[in] Ha - Specific absolute enthalpy (dimensionless).
         * \param[in] v - Velocity (dimensionless).
         * \param[in] rho - Density (dimensionless).
         * \param[out] T - Temperature (dimensionless).
         * \param[out] _flag - The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        void THa(const cellRealArray &ha, const cellRealArray &rho, cellRealArray &T,
                 integerArray &_flag) const;

        /*
         * \brief Compute temperature by given absolute internal energy.
         * \param[in] ea - Absolute internal energy (dimensionless).
         * \param[in] v - Velocity (dimensionless).
         * \param[in] rho - Density (dimensionless).
         * \param[out] T - Temperature (dimensionless).
         */
        void TEa(const cellRealArray &ea, const cellRealArray &rho, cellRealArray &T) const;

        /*
         * \brief Compute temperature by given absolute internal energy with flags.
         * \param[in] ea - Absolute internal energy (dimensionless).
         * \param[in] v - Velocity (dimensionless).
         * \param[in] rho - Density (dimensionless).
         * \param[out] T - Temperature (dimensionless).
         * \param[out] _flag - The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        void TEa(const cellRealArray &ea, const cellRealArray &rho, cellRealArray &T,
                 integerArray &_flag) const;

        /*
         * \brief Compute temperature by given total enthalpy.
         * \param[in] Ha - Total enthalpy, not specific absolute enthalpy (dimensionless).
         * \param[in] v - Velocity (dimensionless).
         * \param[in] rho - Density (dimensionless).
         * \param[out] T - Temperature (dimensionless).
         */
        void THa(const cellRealArray &Ha, const cellVectorArray &v, const cellRealArray &rho,
                 cellRealArray &T) const;

        /*
         * \brief Compute temperature by given total enthalpy with flags.
         * \param[in] Ha - Total enthalpy, not specific absolute enthalpy (dimensionless).
         * \param[in] v - Velocity (dimensionless).
         * \param[in] rho - Density (dimensionless).
         * \param[out] T - Temperature (dimensionless).
         * \param[out] _flag - The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        void THa(const cellRealArray &Ha, const cellVectorArray &v, const cellRealArray &rho,
                 cellRealArray &T, integerArray &_flag) const;

        /*
         * \brief Compute temperature by given total energy.
         * \param[in] ea - Total energy, not absolute internal energy (dimensionless).
         * \param[in] v - Velocity (dimensionless).
         * \param[in] rho - Density (dimensionless).
         * \param[out] T - Temperature (dimensionless).
         */
        void TEa(const cellRealArray &ea, const cellVectorArray &v, const cellRealArray &rho,
                 cellRealArray &T) const;

        /*
         * \brief Compute temperature by given total energy with flags.
         * \param[in] ea - Total energy, not absolute internal energy (dimensionless).
         * \param[in] v - Velocity (dimensionless).
         * \param[in] rho - Density (dimensionless).
         * \param[out] T - Temperature (dimensionless).
         * \param[out] _flag - The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        void TEa(const cellRealArray &ea, const cellVectorArray &v, const cellRealArray &rho,
                 cellRealArray &T, integerArray &_flag) const;

        /*!\brief Correct the array field that exceed the limits by using average method.
         * \param[in] lowLmt - The low limit.
         * \param[in] highLmt - The high limit.
         * \param[in] flag - The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        void correctLimits(cellRealArray &T, const integer cellI, const real lowLmt,
                           const real highLmt, const integer flag) const;

        void p(const cellRealArray &rho, const cellRealArray &T, cellRealArray &p) const;

        void mu(const cellRealArray &p, const cellRealArray &T, cellRealArray &mu) const;
    };
} // namespace OpenHurricane

#include "rhoThermo.inl"