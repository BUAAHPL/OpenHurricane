/*!
 * \file kineticTheory.hpp
 * \brief Header of transport properties by kinetic theory.
 *       The subroutines and functions are in the <i>kineticTheory.cpp</i> file.
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
    /*!\brief The class of kinetic theory transport properties.*/
    class kineticTheory : public transport {
    private:
        static constexpr real muPre_ = 2.6693e-6;

        static constexpr real omegaPre_ = 1.147;

        static constexpr real omegaExponent_ = -0.145;

        /*!
         *\brief An index indicating whether the molecule has a monatomic, linear or nonlinear geometrical configuration.
         *       If the index is 0, the molecule is a single atom. If the index is 1 the molecule is linear, and if it is 2, the molecule
         *       is nonlinear.
         */
        integer geometryIndex_;

        /*!\brief The Lennard-Jones potential well depth in Kelvins.*/
        real ekb_;

        /*!\brief The reverse of Lennard-Jones potential well depth in Kelvins.*/
        real rekb_;

        /*!\brief The Lennard-Jones collision diameter in angstroms.*/
        real sigma_;

        /*!\brief The square of Lennard-Jones collision diameter in angstroms.*/
        real sigma2_;

        /*!\brief The dipole moment in Debye.*/
        real mmu_;

        /*!\brief The polarizability in cubic angstroms.*/
        real alpha_;

        /*!\brief The rotational relaxation collision number at 298 K.*/
        real Zrot_;

        mutable realArray *DijConstPartPtr_;
        mutable realArray *rekbij05Ptr_;
        // Private member functions

        /*!\brief The collision integral omega(2,2).*/
        hur_nodiscard inline real omega22(const real Tstar) const;

    public:
        declareClassNames;

        inline kineticTheory(const speciesList &sp, const integer index, const real Prl,
                             const integer gi, const real ekb, const real sig, const real mmu,
                             const real alp, const real Zr);

        inline kineticTheory(const speciesList &sp, const integer index, const controller &cont);

        kineticTheory(const kineticTheory &tra);

        /*!\brief Construct as copy and given new species list.*/
        kineticTheory(const kineticTheory &tra, const speciesList &sp);

        /*!\brief Return a clone.*/
        virtual hur_nodiscard uniquePtr<transport> clone() const {
            return uniquePtr<transport>(new kineticTheory(*this));
        }

        /*!\brief Return a clone by given new species list.*/
        virtual hur_nodiscard uniquePtr<transport> clone(const speciesList &sp) const {
            return uniquePtr<transport>(new kineticTheory(*this, sp));
        }

        /*!\brief Destructor.*/
        virtual ~kineticTheory() noexcept;

        /*!
         *\brief An index indicating whether the molecule has a monatomic, linear or nonlinear geometrical configuration.
         *       If the index is 0, the molecule is a single atom. If the index is 1 the molecule is linear, and if it is 2, the molecule
         *       is nonlinear.
         */
        hur_nodiscard inline integer geometryIndex() const noexcept;

        /*!\brief The Lennard-Jones potential well depth in Kelvins.*/
        hur_nodiscard inline real ekb() const noexcept;

        /*!\brief The Lennard-Jones collision diameter in angstroms.*/
        hur_nodiscard inline virtual real sigma() const noexcept;

        /*!\brief The dipole moment in Debye.*/
        hur_nodiscard inline real mmu() const noexcept;

        /*!\brief The polarizability in cubic angstroms.*/
        hur_nodiscard inline real alpha() const noexcept;

        /*!\brief The rotational relaxation collision number at 298 K.*/
        hur_nodiscard inline real Zrot() const noexcept;

        /*!
         * \brief Return the dynamic viscosity [kg/(m s)].
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        hur_nodiscard virtual real mu(const real p, const real T) const;

        /*!
         * \brief Return the thermal conductivity [W/mK].
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard virtual real kappa(const real p, const real T, const real cpi) const;

        /*!
         * \brief Return the thermal conductivity [W/mK].
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] mui - Dynamic viscosity [kg/(m s)].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard virtual real kappa(const real p, const real T, const real mui,
                                         const real cpi) const;

        hur_nodiscard virtual real Di(const real p, const real T, const realArray &xi,
                                      const PtrList<transport> &tranPtr) const;

        hur_nodiscard inline virtual real Di(const real p, const real T, const realArray &xi,
                                             const PtrList<transport> &tranPtr,
                                             realArray &Dij) const;

        hur_nodiscard virtual real Di(const real p, const real T, const real Tv1p5,
                                      const realArray &xi, const PtrList<transport> &tranPtr,
                                      realArray &Dij) const;

        hur_nodiscard virtual real Di(const real p, const real T, const real Tv1p5,
                                      const realArray &xi, const PtrList<transport> &tranPtr) const;

        kineticTheory &operator=(const kineticTheory &);
    };
} // namespace OpenHurricane

#include "kineticTheory.inl"