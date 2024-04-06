/*!
 * \file SpalartAllmaras.hpp
 * \brief Header of Spalart-Allmaras turbulence model
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

#include "RANSModel.hpp"

namespace OpenHurricane {

    class SpalartAllmaras : public RANSModel {
    protected:
        inline real fw(const real g) const {
            return g * pow((1.0 + pow6(cw3_)) / (pow6(g) + pow6(cw3_)), 1.0 / 6.0);
        }

        inline real g(const real r) const { return r + cw2_ * (pow6(r) - r); }

        inline real ft2(const real chi) const { return ct3_ * exp(-ct4_ * sqr(chi)); }

    protected:
        enum SAVersion : short {
            SA,   /**<"Standard" Spalart-Allmaras One-Equation Model.*/
            SANEG /**<"Negative" Spalart-Allmaras One-Equation Model.*/
        };

        cellRealArray nut_;

        realArray nutLastTime_;

        // Model coefficients

        real cb1_;
        real cb2_;
        real sigma_;
        real cv1_;
        real kappa_;
        real cw1_;
        real cw2_;
        real cw3_;

        real ct3_;
        real ct4_;

        real cv2_;
        real cv3_;

        SAVersion SAVer_;

        /**\brief Relaxation factor of nut equation for splitting solver.*/
        real nutRelax_;

    public:
        declareClassNames;

        /*!\brief Construct from components.*/
        SpalartAllmaras(const controller &cont, flowModel &ev);

        SpalartAllmaras(const SpalartAllmaras &) = delete;
        SpalartAllmaras &operator=(const SpalartAllmaras &) = delete;
        /*!\brief Destructor.*/
        virtual ~SpalartAllmaras() noexcept {}

        // Member Functions

        /**\brief Turbulent parameters initialization function*/
        virtual void turbParaInitialize();

        /**\brief mut at boudaryface calculation function*/
        virtual void mutBoundaryUpdate();

        /*!\brief Set turbulent complete boundary condition.*/
        virtual void bndValueSetting(const controller &cont);

        /*!\brief explicit source.*/
        virtual void expSource();

        /*!\brief implicit source.*/
        virtual void impSource();

        /*!\brief implicit source.*/
        virtual void fullImpSource(cellRealSquareMatrixArray &Jac, const integer rhoId,
                                   const integer rhoTurb0);

        /*!\brief viscous flux.*/
        virtual void visFlux(const faceRealArray &rhof, const faceRealArray &mulf,
                             const faceRealArray &mutf, const cellRealArray &mul,
                             const cellRealArray &mut, const cellRealArray &rho);

        virtual void update();

        virtual void limit();

        inline cellRealArray &nut() { return nut_; }

        virtual realArray k() const;

        virtual realArray epsilon() const;

        /**
         * \brief Turbulent Reynolds number [dimensionless].
         */
        virtual hur_nodiscard realArray Ret() const;

        virtual cellRealArray &var(const integer i);

        virtual const cellRealArray &var(const integer i) const;

        virtual void solving(const realArray &dt);
        virtual void initialize() { turbParaInitialize(); }

        virtual void initializeRestart() {
            this->updateBoundary();
            this->mut().updateBoundary();
        }

        virtual void updateBoundary();
        virtual void limitAndUpdateBoundary();

        /**
         * \brief Only available for RANS models.
         */
        virtual symmTensorArray ReynoldsStressTensor() const;

        virtual faceSymmTensorArray tauEff(const faceRealArray &rhof, const faceRealArray &muf,
                                           const faceRealArray &mutf,
                                           const faceTensorArray &deltafV) const;

        virtual void calcGrad(const spatialScheme &sps);
    };

} // namespace OpenHurricane
