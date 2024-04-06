/*!
 * \file SST.hpp
 * \brief Header of SST turbulence model
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

    class SST : public RANSModel {
    protected:
        enum SSTVersion : short {
            standardSST, /**<"Standard" Menter SST Two-Equation Model.*/
            SST2003,     /**<Menter SST Two-Equation Model from 2003*/
            SSTSUST,     /**<Menter SST Two-Equation Model with Controlled Decay.*/
            SSTV         /**<Menter SST Two-Equation Model with Vorticity Source Term.*/
        };

        /*!\brief Turbulent kinetic energy.*/
        cellRealArray k_;

        /*!\brief The specific turbulent dissipation rate.*/
        cellRealArray w_;

        realArray kLastTime_;
        realArray wLastTime_;

        cellRealArray F1_;

        cellRealArray kSource;
        cellRealArray wSource;

        cellRealArray kInvFlux;
        cellRealArray wInvFlux;

        cellRealArray kVisFluxTmp;
        cellRealArray wVisFluxTmp;

        cellRealArray F2_;
        cellRealArray SS_;

        // Model coefficients

        real a1_, a2_, a3_;
        real betas_;

        /*!\brief Von Karman constant*/
        real kappa_;

        real sigmak1_, sigmaw1_, beta1_;
        real sigmak2_, sigmaw2_, beta2_;
        real gam1_;
        real gam2_;

        real Cmu_;
        real E_;

        real kamb_;
        real wamb_;

        real F1(const integer celli, const real cdsst) const;

        real F2(const integer celli) const;

        SSTVersion SSTVer_;

        /**\brief Relaxation factor of k equation for splitting solver.*/
        real kRelax_;

        /**\brief Relaxation factor of w equation for splitting solver.*/
        real wRelax_;

        vector2DArray rak_;
        vector2DArray raw_;

        /*!\brief Viscous flux of k equation.*/
        realArray kVisFlux_;
        realArray Pk_;
        realArray dqk_;
        realArray dqw_;

        real yplusLam_;
        static real yPlusLam(const real kappa, const real E);

        bool omegaWallFunction_;
        integerList omegaWallFunctionFaceZoneList_;

        real minK_;
        real minW_;

        real CDkOmega(const real cdsst) const;

        /*!\brief In order to avoid the buildup of turbulent kinetic energy in the stagnation regions,
         * the production term in the turbulence equations should be limited.
         */
        real minP(const real P, const real brwk) const;

        real dpkdrk(const real p, const real brwk, const real rhok) const;
        real dpkdrk(const real p, const real brwk, const real rhok, const real dpkdrkt) const;

        real dpkdrw(const real P, const real brwk, const real rhow) const;

    public:
        declareClassNames;

        /*!\brief Construct from components.*/
        SST(const controller &cont, flowModel &ev);

        SST(const SST &) = delete;
        SST &operator=(const SST &) = delete;

        /*!\brief Destructor.*/
        virtual ~SST() noexcept {}

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

        void updateF1();

        virtual void limit();

        inline cellRealArray &k() noexcept { return k_; }

        inline cellRealArray &w() noexcept { return w_; }

        /*!\brief Turbulent kinetic energy.*/
        virtual realArray k() const;

        /*!\brief The turbulent dissipation rate.*/
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

    private:
        void limitKAndW();

    public:
        /**
         * \brief Only available for RANS models.
         */
        virtual symmTensorArray ReynoldsStressTensor() const;

        /**
         * \brief The form of SST model is given on the linear eddy viscosity models.
         * Linear models use the Boussinesq assumption.
         */
        virtual faceSymmTensorArray tauEff(const faceRealArray &rhof, const faceRealArray &mulf,
                                           const faceRealArray &mutf,
                                           const faceTensorArray &deltafV) const;

        virtual void calcGrad(const spatialScheme &sps);


        /*!\brief Correct energy equation (molecular diffusion and turbulent transport in the energy equation)*/
        virtual void correctEnergyEqVisFlux(cellRealArray &E) const;

    protected:
        /**
         * \brief Get the production of k.
         * \param[out] Pk - The production of k
         * \param[out] dPkdrhokPtr - The pointer of Jacobian of the production of k.
         */
        void getPk(realArray &Pk, realArray *dPkdrhokPtr = nullptr);

        virtual void correctOmegaSource() {}
        virtual void correctOmegaSourceImp() {}

        void correctOmegaRHS();
        void correctOmegaRHSDiagImp();
        void correctOmegaRHSImp(cellRealSquareMatrixArray &Jac, const integer rhoId,
                                const integer rhoTurb0);
    };

} // namespace OpenHurricane
