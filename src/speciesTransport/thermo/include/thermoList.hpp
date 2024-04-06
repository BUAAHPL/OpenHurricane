/*!
 * \file thermoList.hpp
 * \brief Header of thermo properties list.
 *       The subroutines and functions are in the <i>thermoList.cpp</i> file.
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

#include "smartPointerList.hpp"
#include "thermo.hpp"

namespace OpenHurricane {
    /*!\brief The class of thermo properties list.*/
    class thermoList : public object {
    private:
        /*!\brief Hold const-reference to the list of species.*/
        const speciesList &species_;

        /*!\brief The thermo list.*/
        PtrList<thermo> thList_;

        uniquePtr<equationOfState> eosPtr_;

        // Temperature limits of applicability of functions
        real TLow_, THigh_;

        real globalLowT_, globalCommonT_, globaHighT_;

        /*!\brief Convergence tolerance of energy to temperature inversion functions.*/
        static const real tol_;

        /*!\brief Max number of iterations in energy to temperature inversion functions.*/
        static const int maxIter_;

    public:
        thermoList();

        thermoList(const speciesList &species);

        thermoList(const speciesList &species, const runtimeMesh &mesh);

        thermoList(const speciesList &species, const controller &cont);

        thermoList(const speciesList &species, const controller &cont, const runtimeMesh &mesh);

        thermoList(const thermoList &);

        thermoList(thermoList &&) noexcept;

        ~thermoList() noexcept;

        /*!\brief Return the table of species.*/
        hur_nodiscard inline const speciesList &species() const noexcept { return species_; }

        /*!\brief The transport table.*/
        const hur_nodiscard PtrList<thermo> &thTable() const noexcept;

        hur_nodiscard inline PtrList<thermo> &thTable() noexcept { return thList_; }

        hur_nodiscard const equationOfState &eos() const noexcept;

        hur_nodiscard inline equationOfState *eosPtr() noexcept { return eosPtr_.get(); }

        hur_nodiscard inline bool isNUllEos() const noexcept { return eosPtr_.isNull(); }

        void setEos(const controller &cont);

        /*!\brief Return temperature lower limit.*/
        hur_nodiscard inline real TLow() const noexcept { return TLow_; }

        /*!\brief Return temperature higher limit.*/
        hur_nodiscard inline real THigh() const noexcept { return THigh_; }

        inline void setLimitTemperature(const real LT, const real HT) noexcept {
            TLow_ = LT;
            THigh_ = HT;
        }

        inline void setGlobalTemperature(const real gLT, const real gCT, const real gHT) noexcept {
            globalLowT_ = gLT;
            globalCommonT_ = gCT;
            globaHighT_ = gHT;
        }

        hur_nodiscard inline real globalLowT() const noexcept { return globalLowT_; }

        hur_nodiscard inline real globaHighT() noexcept { return globaHighT_; }

        hur_nodiscard inline real globalCommonT() noexcept { return globalCommonT_; }

        /*!\brief Limit the temperature to be in the range TLow_ to THigh_.*/
        hur_nodiscard real limit(const real T) const;

        /*!
         * \brief Limit the temperature to be in the range TLow_ to THigh_.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real limit(const real T, integer &iFlag) const;

        /**\brief Return the pressure [Pa]. Note: all dimensional.
         * \param[in] rho - Density [kg/m^3].
         * \param[in] T - temperature [K].
         * \param[in] yi - species mass fractions
         * \return pressure [Pa].
         */
        hur_nodiscard inline real p(const real rho, const real T, const realArray &yi) const;

        /**\brief Return the pressure [Pa]. Note: all dimensional.
         * \param[in] rho - Density [kg/m^3].
         * \param[in] T - temperature [K].
         * \param[in] yi - species mass fractions
         * \param[in] cellI - The index of cell
         * \return pressure [Pa].
         */
        hur_nodiscard inline real p(const real rho, const real T, const PtrList<cellRealArray> &yi,
                                    const integer cellI) const;

        /**\brief Return the density [kg/m^3]. Note: all dimensional.
         * \param[in] p - Pressure [Pa].
         * \param[in] T - temperature [K].
         * \param[in] yi - species mass fractions
         * \param[in] cellI - The index of cell
         * \return density [kg/m^3].
         */
        hur_nodiscard inline real rho(const real p, const real T, const realArray &yi) const;

        /**\brief Return the density [kg/m^3]. Note: all dimensional.
         * \param[in] p - Pressure [Pa].
         * \param[in] T - temperature [K].
         * \param[in] yi - species mass fractions
         * \return density [kg/m^3].
         */
        hur_nodiscard inline real rho(const real p, const real T, const PtrList<cellRealArray> &yi,
                                      const integer cellI) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] at standard pressure (1 atm).
         * \param[in] T - temperature [K].
         * \param[in] yi - species mass fractions
         * \return Heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard real cp0(const real T, const realArray &yi) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] at standard pressure (1 atm).
         * \param[in] T - temperature [K].
         * \param[in] yi - species mass fractions
         * \param[in] cellI - The index of cell
         * \return Heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard real cp0(const real T, const PtrList<cellRealArray> &yi,
                               const integer cellI) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] at standard pressure (1 atm).*/
        hur_nodiscard inline real cp0i(const real T, const integer i) const;

        /*!\brief Heat capacity at constant pressure [J/(kmol K)] at standard pressure (1 atm).*/
        hur_nodiscard real Cp0(const real T, const realArray &xi) const;

        /*!\brief Heat capacity at constant pressure [J/(kmol K)] at standard pressure (1 atm).*/
        hur_nodiscard inline real Cp0i(const real T, const integer i) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given density rho.*/
        hur_nodiscard inline real cp_rho(const real rho, const real T, const realArray &yi) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given density rho.*/
        hur_nodiscard inline real cp_rho(const real rho, const real T,
                                         const PtrList<cellRealArray> &yi,
                                         const integer cellI) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given pressure p.*/
        hur_nodiscard inline real cp_p(const real p, const real T, const realArray &yi) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given pressure p.*/
        hur_nodiscard inline real cp_p(const real p, const real T, const PtrList<cellRealArray> &yi,
                                       const integer cellI) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given density rho.*/
        hur_nodiscard inline real cpi_rho(const real rho, const real T, const integer i) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given pressure p.*/
        hur_nodiscard inline real cpi_p(const real p, const real T, const integer i) const;

        /*!\brief Heat capacity at constant pressure [J/(kmol K)] by given density rho.*/
        hur_nodiscard inline real Cp_rho(const real rho, const real T, const realArray &xi) const;

        /*!\brief Heat capacity at constant pressure [J/(kmol K)] by given pressure p.*/
        hur_nodiscard inline real Cp_p(const real p, const real T, const realArray &xi) const;

        /*!\brief Heat capacity at constant pressure [J/(kmol K)] by given density rho.*/
        hur_nodiscard inline real Cpi_rho(const real rho, const real T, const integer i) const;

        /*!\brief Heat capacity at constant pressure [J/(kmol K)] by given pressure p.*/
        hur_nodiscard inline real Cpi_p(const real p, const real T, const integer i) const;

        /*!\brief Heat capacity at constant volume [J/(kg K)] at standard pressure (1 atm).*/
        hur_nodiscard real cv0(const real T, const realArray &yi) const;

        /*!\brief Heat capacity at constant volume [J/(kg K)] at standard pressure (1 atm).*/
        hur_nodiscard real cv0(const real T, const PtrList<cellRealArray> &yi,
                               const integer cellI) const;

        /*!\brief Heat capacity at constant volume [J/(kmol K)] at standard pressure (1 atm).*/
        hur_nodiscard inline real Cv0i(const real T, const integer i) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given pressure p.*/
        hur_nodiscard inline real cv_p(const real p, const real T, const realArray &yi) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given pressure p.*/
        hur_nodiscard inline real cv_p(const real p, const real T, const PtrList<cellRealArray> &yi,
                                       const integer cellI) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given density rho.*/
        hur_nodiscard inline real cv_rho(const real rho, const real T, const realArray &yi) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] by given density rho.*/
        hur_nodiscard inline real cv_rho(const real rho, const real T,
                                         const PtrList<cellRealArray> &yi,
                                         const integer cellI) const;

        /*!\brief Ratio of specific heat capacities.*/
        hur_nodiscard inline real gamma(const real p, const real T, const realArray &yi) const;

        /*!\brief Ratio of specific heat capacities.*/
        hur_nodiscard inline real gamma(const real p, const real T,
                                        const PtrList<cellRealArray> &yi,
                                        const integer cellI) const;

        /*!\brief Ratio of specific heat capacities.*/
        hur_nodiscard inline real gammai(const real p, const real T, const integer i) const;

        /*!\brief Ratio of specific heat capacities for ascoutic speed calculation.*/
        hur_nodiscard inline real gammaTh(const real p, const real rho, const real T,
                                          const realArray &yi) const;

        /*!\brief Ratio of specific heat capacities for ascoutic speed calculation.*/
        hur_nodiscard inline real gammaTh(const real p, const real rho, const real T,
                                          const PtrList<cellRealArray> &yi,
                                          const integer cellI) const;

        /*!\brief Ratio of specific heat capacities for ascoutic speed calculation.*/
        hur_nodiscard inline real gammaThi(const real p, const real rho, const real T,
                                           const integer i) const;

        /*!\brief Absolute Enthalpy [J/kg] at standard pressure.*/
        hur_nodiscard real ha0(const real T, const realArray &yi) const;

        /*!\brief Absolute Enthalpy [J/kg] at standard pressure.*/
        hur_nodiscard real ha0(const real T, const PtrList<cellRealArray> &yi,
                               const integer cellI) const;

        /*!\brief Absolute Enthalpy field [J/kg] at standard-state.*/
        hur_nodiscard realArray ha0(const real T) const;

        /*!\brief Absolute Enthalpy [J/kg] by given pressure p.*/
        hur_nodiscard inline real ha_p(const real p, const real T, const realArray &yi) const;

        /*!\brief Absolute Enthalpy [J/kg] by given pressure p.*/
        hur_nodiscard inline real ha_p(const real p, const real T, const PtrList<cellRealArray> &yi,
                                       const integer cellI) const;

        /*!\brief Absolute Enthalpy [J/kg] by given pressure p.*/
        hur_nodiscard inline real ha_p(const real p, const real T, const integer i) const;

        /*!\brief Absolute Enthalpy [J/kg] by given density rho.*/
        hur_nodiscard inline real ha_rho(const real rho, const real T, const realArray &yi) const;

        /*!\brief Absolute Enthalpy [J/kg] by given density rho.*/
        hur_nodiscard inline real ha_rho(const real rho, const real T,
                                         const PtrList<cellRealArray> &yi,
                                         const integer cellI) const;

        /*!\brief Chemical enthalpy [J/kg].*/
        hur_nodiscard realArray hc() const;

        /*!\brief Chemical enthalpy [J/kg].*/
        hur_nodiscard real hc(const realArray &yi) const;

        /*!\brief Chemical enthalpy [J/kg].*/
        hur_nodiscard real hc(const PtrList<cellRealArray> &yi, const integer cellI) const;

        /*!\brief Sensible Enthalpy [J/kg] by given pressure p.*/
        hur_nodiscard inline real hs_p(const real p, const real T, const realArray &yi) const;

        /*!\brief Sensible Enthalpy [J/kg] by given pressure p.*/
        hur_nodiscard inline real hs_p(const real p, const real T, const PtrList<cellRealArray> &yi,
                                       const integer cellI) const;

        /*!\brief Sensible Enthalpy [J/kg] by given density rho.*/
        hur_nodiscard inline real hs_rho(const real rho, const real T, const realArray &yi) const;

        /*!\brief Sensible Enthalpy [J/kg] by given density rho.*/
        hur_nodiscard inline real hs_rho(const real rho, const real T,
                                         const PtrList<cellRealArray> &yi,
                                         const integer cellI) const;

        /*!\brief Absolute internal energy [J/kg] by given pressure p.*/
        hur_nodiscard inline real ea_p(const real p, const real T, const realArray &yi) const;

        /*!\brief Absolute internal energy [J/kg] by given pressure p.*/
        hur_nodiscard inline real ea_p(const real p, const real T, const PtrList<cellRealArray> &yi,
                                       const integer cellI) const;

        /*!\brief Absolute internal energy [J/kg] by given density rho.*/
        hur_nodiscard inline real ea_rho(const real rho, const real T, const realArray &yi) const;

        /*!\brief Absolute internal energy [J/kg] by given density rho.*/
        hur_nodiscard inline real ea_rho(const real rho, const real T,
                                         const PtrList<cellRealArray> &yi,
                                         const integer cellI) const;

        /*!\brief Entropy [J/(kg K)].*/
        hur_nodiscard real s_p(const real pm, const real T, const realArray &yi) const;

        /*!\brief Entropy [J/(kg K)].*/
        hur_nodiscard real s_p(const real pm, const real T, const PtrList<cellRealArray> &yi,
                               const integer cellI) const;

        /*!\brief Entropy [J/(kg K)].*/
        hur_nodiscard inline real s_rho(const real rhom, const real T, const realArray &yi) const;

        /*!\brief Entropy [J/(kg K)].*/
        hur_nodiscard inline real s_rho(const real rhom, const real T,
                                        const PtrList<cellRealArray> &yi,
                                        const integer cellI) const;

        virtual real inteCp0dT(const real T1, const real T2, const realArray &yi) const;
        virtual real inteCp0dT(const real T1, const real T2, const PtrList<cellRealArray> &yi,
                               const integer cellI) const;

        /*!\brief The standard-state Gibbs free energy field [J/kg].*/
        hur_nodiscard realArray g0(const real T) const;

        /*!\brief The standard-state Gibbs free energy field [J/kg].*/
        void g0(const real T, realArray &gf) const;

        /*!\brief The standard-state molar Gibbs free energy field [J/kmol].*/
        hur_nodiscard realArray G0(const real T) const;

        /*!\brief The standard-state molar Gibbs free energy field [J/kmol].*/
        void G0(const real T, realArray &Gf) const;

        /*!\brief The standard-state molar Gibbs free energy field [J/kmol].*/
        void DG0DT(const real T, realArray &DG0DTf) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real T(real f, real p, real T0, integer &iFlag, const PtrList<cellRealArray> &yi,
               const integer cellI,
               real (thermoList::*F)(const real, const real, const PtrList<cellRealArray> &,
                                      const integer) const,
               real (thermoList::*dFdT)(const real, const real, const PtrList<cellRealArray> &,
                                         const integer) const,
               real (thermoList::*limit)(const real, integer &) const) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real T(real f, real p, real T0, integer &iFlag, const realArray &yi,
               real (thermoList::*F)(const real, const real, const realArray &) const,
               real (thermoList::*dFdT)(const real, const real, const realArray &) const,
               real (thermoList::*limit)(const real, integer &) const) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real THa_rho(real ha, real rho, real T0, integer &iFlag, const realArray &yi) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real THa_rho(real ha, real rho, real T0, integer &iFlag, const PtrList<cellRealArray> &yi,
                     const integer cellI) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real THa_p(real ha, real p, real T0, integer &iFlag, const realArray &yi) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real THa_p(real ha, real p, real T0, integer &iFlag, const PtrList<cellRealArray> &yi,
                   const integer cellI) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real TEa_rho(real ea, real rho, real T0, integer &iFlag, const realArray &yi) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real TEa_rho(real ea, real rho, real T0, integer &iFlag, const PtrList<cellRealArray> &yi,
                     const integer cellI) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real TEa_p(real ea, real p, real T0, integer &iFlag, const realArray &yi) const;

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p/rho, T)
         *  and dF(p/rho, T)/dT.
         * \param[in]  f - enthalpy or internal energy.
         * \param[in]  p - static pressure/density.
         * \param[in]  T0 - initial temperature.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        real TEa_p(real ea, real p, real T0, integer &iFlag, const PtrList<cellRealArray> &yi,
                   const integer cellI) const;

        hur_nodiscard inline thermo &operator[](const integer i) { return thList_[i]; }

        hur_nodiscard inline const thermo &operator[](const integer i) const { return thList_[i]; }
    };
} // namespace OpenHurricane

#include "thermoList.inl"