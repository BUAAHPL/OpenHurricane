/*!
 * \file transportList.hpp
 * \brief Header of transport properties list.
 *       The subroutines and functions are in the <i>transportList.cpp</i> file.
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

#include "smartPointerList.hpp"
#include "transport.hpp"

namespace OpenHurricane {
    /*!\brief The class of transport properties list.*/
    class transportList {
    private:

        /*!\brief Hold const-reference to the list of species.*/
        const speciesList &species_;

        /*!\brief The transport list.*/
        PtrList<transport> tranList_;

        /** \brief Temporary array. */
        mutable realArray txi_, tmuf_, tkppaf_, tphi_;

        // Wilke
        mutable realSquareMatrix *phiPrePtr_;
        mutable realSquareMatrix *WiWjPtr_;

        hur_nodiscard const realSquareMatrix &phiPre() const;
        hur_nodiscard const realSquareMatrix &WiWj() const;

    public:

        transportList();

        /*!\brief Construct from components.*/
        transportList(const speciesList &species);

        /*!\brief Construct from components.*/
        transportList(const speciesList &species, const controller &cont);

        /*!\brief Copy constructor.*/
        transportList(const transportList &);

        /*!\brief Copy constructor.*/
        transportList(transportList &&) noexcept;

        /*!\brief Destructor.*/
        ~transportList() noexcept;


        /*!\brief Return the table of species.*/
        hur_nodiscard inline const speciesList &species() const noexcept { return species_; }

        /*!\brief The transport table.*/
        hur_nodiscard inline const PtrList<transport> &tranTable() const noexcept;

        /*!\brief The transport table.*/
        hur_nodiscard inline PtrList<transport> &tranTable() noexcept { return tranList_; }

        /*!\brief Return laminar Prandtl number.*/
        hur_nodiscard inline real Pr() const { return tranList_[0].Pr(); }

        /*!\brief Return laminar Prandtl number.*/
        inline void setPr(const real newPr) {
            for (integer i = 0; i < tranList_.size(); ++i) {
                tranList_[i].setPr(newPr);
            }
        }

        /*!
         * \brief Return the dynamic viscosity [kg/(m s)] of species i.
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        hur_nodiscard real mu(const real p, const real T, const integer i) const;

        /*!
         * \brief Return the dynamic viscosity [kg/(m s)] of mixtures.
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        hur_nodiscard real mu(const real p, const real T, const realArray &yi) const;

        /*!
         * \brief Return the dynamic viscosity [kg/(m s)] of mixtures.
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        hur_nodiscard real mu(const real p, const real T, const PtrList<cellRealArray> &yi,
                              const integer cellI) const;

        /*!
         * \brief Return the dynamic viscosity [kg/(m s)] of every species.
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        hur_nodiscard realArray muf(const real p, const real T) const;

        /*!
         * \brief Return the dynamic viscosity [kg/(m s)] of every species.
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        void muf(const real p, const real T, realArray &mui) const;

        /*!
         * \brief Return the linematic viscosity [kg/(m s)] of species i.
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        hur_nodiscard real nu(const integer i, const real p, const real T, const real rho) const;

        /*!
         * \brief Return the thermal conductivity [W/mK] of species i.
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard real kappa(const real p, const real T, const real cpi, const integer i) const;

        /*!
         * \brief Return the thermal conductivity [W/mK] of species i.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] mui - Dynamic viscosity [kg/(m s)].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard real kappa(const real p, const real T, const real mui, const real cpi,
                                 const integer i) const;

        /*!
         * \brief Return the thermal conductivity [W/mK] for every species.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] mui - Dynamic viscosity [kg/(m s)].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard realArray kappaf(const real p, const real T, const realArray &mui,
                                       const realArray &cpi) const;

        /*!
         * \brief Return the thermal conductivity [W/mK] for every species.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] mui - Dynamic viscosity [kg/(m s)].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        void kappaf(const real p, const real T, const realArray &mui, const realArray &cpi,
                    realArray &kappaf) const;

        /*!
         * \brief Return the thermal conductivity [W/mK] for mixtures.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] mui - Dynamic viscosity [kg/(m s)].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard real kappa(const real p, const real T, const realArray &mui,
                                 const realArray &cpi, const realArray &yi) const;

        /*!
         * \brief The dynamic viscosity and thermal conductivity [W/mK] for mixtures.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        void muKappa(const real p, const real T, const realArray &cpi, const realArray &yi,
                     real &mum, real &kappam) const;

        /*!
         * \brief The dynamic viscosity and thermal conductivity [W/mK] for mixtures.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        void muKappa(const real p, const real T, const realArray &cpi,
                     const PtrList<cellRealArray> &yi, const integer cellI, real &mum,
                     real &kappam) const;

        /*!
         * \brief Return the Species mass diffusivity of specie i.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] yi - Mass fraction of species.
         */
        hur_nodiscard real D(const real p, const real T, const realArray &yi,
                             const integer i) const;

        /**
         * \brief The mixture diffusivity.
         */
        hur_nodiscard real DiffMix(const real p, const real T, const PtrList<cellRealArray> &yi,
                                   const PtrList<cellRealArray> &Dim, const integer cellI) const;

        /*!
         * \brief Return the Species mass diffusivity of specie i.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] xi - Mole fraction of species.
         */
        hur_nodiscard real D_Xi(const real p, const real T, const realArray &xi,
                                const integer i) const;

        /*!
         * \brief Return the Species mass diffusivity of specie i.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] yi - Mass fraction of species.
         */
        hur_nodiscard real D(const real p, const real T, const realArray &yi, const integer i,
                             realArray &Dij) const;

        /*!
         * \brief Return the Species mass diffusivity of specie i.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] xi - Mole fraction of species.
         */
        hur_nodiscard inline real D_Xi(const real p, const real T, const realArray &xi,
                                       const integer i, realArray &Dij) const;

        /*!
         * \brief Return the Species mass diffusivity of specie i.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] xi - Mole fraction of species.
         */
        hur_nodiscard inline real D_Xi(const real p, const real T, const real Tv1p5,
                                       const realArray &xi, const integer i, realArray &Dij) const;

        /*!
         * \brief Return the Species mass diffusivity of specie i.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] xi - Mole fraction of species.
         */
        hur_nodiscard inline real D_Xi(const real p, const real T, const real Tv1p5,
                                       const realArray &xi, const integer i) const;

        /*!
         * \brief Return the Species mass diffusivity of specie i.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] yi - Mass fraction of species.
         */
        hur_nodiscard real D(const real p, const real T, const PtrList<cellRealArray> &yi,
                             const integer cellI, const integer i) const;

        /*!
         * \brief Return the Species mass diffusivity of specie i.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] yi - Mass fraction of species.
         */
        hur_nodiscard real D(const real p, const real T, const PtrList<cellRealArray> &yi,
                             const integer cellI, const integer i, realArray &Dij) const;

        /*!
         * \brief Return the Species mass diffusivity for every species.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] yi - Mass fraction of species.
         */
        hur_nodiscard realArray Df(const real p, const real T, const realArray &yi) const;

        /*!
         * \brief Return the Species mass diffusivity for every species.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] yi - Mass fraction of species.
         */
        hur_nodiscard realArray Df(const real p, const real T, const PtrList<cellRealArray> &yi,
                                   const integer cellI) const;

        /*!
         * \brief Return the Species mass diffusivity for every species.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] yi - Mass fraction of species.
         */
        realArray Df(const real p, const real T, const PtrList<cellRealArray> &yi,
                     const integer cellI, realArray &Dij) const;

        /*!
         * \brief Return the Species mass diffusivity for every species.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] yi - Mass fraction of species.
         */
        void Df(realArray &Dff, const real p, const real T, const PtrList<cellRealArray> &yi,
                const integer cellI) const;

        /*!
         * \brief Return the Species mass diffusivity for every species.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] yi - Mass fraction of species.
         */
        void Df(realArray &Dff, const real p, const real T, const PtrList<cellRealArray> &yi,
                const integer cellI, realArray &Dij) const;
    };
} // namespace OpenHurricane

#include "transportList.inl"