/*!
 * \file turbulenceModel.hpp
 * \brief Header of turbulence model
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

#include "fVArraysInclude.hpp"
#include "faceGrad.hpp"
#include "faceInterpolation.hpp"
#include "gradients.hpp"
#include "runtimeMesh.hpp"
#include "spatialScheme.hpp"
#include "wallDistance.hpp"

namespace OpenHurricane {
    /*!\brief The basic class of turbulence model*/
    class turbulenceModel {
    public:
        enum solverType : short { coupled, splitting };

    protected:
        /*!\brief Wall distance.*/
        uniquePtr<wallDistance> yPtr_;

        /**\brief Reference to flow model.*/
        flowModel &flowM_;

        /**\brief The number of turbulence model equations.*/
        integer nEq_;

        solverType solveType_;

        integer cStep_;

        integer maxStep_;

        cellRealArray &mul_;
        cellRealArray &mut_;

    public:
        declareClassNames;
        declareObjFty(turbulenceModel, controller, (const controller &cont, flowModel &flowMod),
                      (cont, flowMod));

        /*!\brief Construct from components.*/
        turbulenceModel(const controller &cont, flowModel &flowMod);

        /*!\brief Disallow default copy constructor.*/
        turbulenceModel(const turbulenceModel &) = delete;

        /*!\brief Disallow default bitwise assignment.*/
        void operator=(const turbulenceModel &) = delete;

        static uniquePtr<turbulenceModel> creator(const controller &cont, flowModel &flowMod);

        /*!\brief Destructor.*/
        virtual ~turbulenceModel() noexcept {}

        /*!\brief Hold const reference to runtime mesh.*/
        hur_nodiscard inline const runtimeMesh &mesh() const noexcept { return flowM_.mesh(); }

        hur_nodiscard inline const controller &cont() const noexcept {
            return mesh().Iteration().cont();
        }

        /*!\brief Wall distance.*/
        hur_nodiscard inline const wallDistance &y() const noexcept {
            if (!yPtr_) {
                LFatal("Attempt to access to the null pointer of wall distance method.");
            }
            return *yPtr_;
        }

        /*!\brief Wall distance.*/
        hur_nodiscard inline wallDistance &y() noexcept {
            if (!yPtr_) {
                LFatal("Attempt to access to the null pointer of wall distance method.");
            }
            return *yPtr_;
        }

        /*!\brief The distance to wall.*/
        hur_nodiscard inline const cellRealArray &wallDist() const {
            if (!yPtr_) {
                LFatal("Attempt to access to the null pointer of wall distance method.");
                return cellRealArray::nullObject();
            } else {
                return yPtr_->wallDist();
            }
        }

        hur_nodiscard inline cellRealArray &mul() noexcept { return mul_; }

        hur_nodiscard inline const cellRealArray &mul() const noexcept { return mul_; }

        hur_nodiscard inline realArray mul(const integer zoneid) { return flowM_.mul(zoneid); }

        hur_nodiscard inline cellRealArray nu() { return flowM_.nu(); }

        hur_nodiscard inline realArray nu(const integer zoneid) { return flowM_.nu(zoneid); }

        /*!\brief Return turbulent viscosity field.*/
        hur_nodiscard inline cellRealArray &mut() noexcept { return mut_; }

        hur_nodiscard inline const cellRealArray &mut() const noexcept { return mut_; }

        /*!\brief Return access to the velocity field.*/
        hur_nodiscard inline cellVectorArray &v() noexcept { return flowM_.v(); }

        /*!\brief Return access to the velocity field.*/
        hur_nodiscard inline const cellVectorArray &v() const noexcept { return flowM_.v(); }

        /*!\brief Return access to the density field.*/
        hur_nodiscard inline cellRealArray &rho() noexcept { return flowM_.rho(); }

        /*!\brief Return const access to the density field.*/
        hur_nodiscard inline const cellRealArray &rho() const noexcept { return flowM_.rho(); }

        /*!\brief Lower limit of mut*/
        hur_nodiscard inline real mutLow() const noexcept { return flowM_.mutLow(); }

        /*!\brief Higher limit of mut*/
        hur_nodiscard inline real mutHigh() const noexcept { return flowM_.mutHigh(); }

        hur_nodiscard inline solverType solveType() const noexcept { return solveType_; }

        hur_nodiscard inline bool isCoupled() const noexcept { return solveType_ == coupled; }

        hur_nodiscard inline bool isSplitting() const noexcept { return solveType_ == splitting; }

        /**\brief The dimensionless velocity: u plus*/
        hur_nodiscard realArray up(const realArray &muw, const integer zoneId);

        /**\brief The dimensionless velocity: u plus*/
        hur_nodiscard realArray up(const integer zoneId);

        /**\brief The dimensionless wall distance: y plus*/
        hur_nodiscard realArray yPlus(const integer zoneId);

        /*!\brief explicit source.*/
        virtual void expSource() = 0;

        /*!\brief implicit source.*/
        virtual void impSource() = 0;

        /*!\brief implicit source.*/
        virtual void fullImpSource(cellRealSquareMatrixArray &Jac, const integer rhoId,
                                   const integer rhoTurb0) = 0;

        /*!\brief viscous flux.*/
        virtual void visFlux(const faceRealArray &rhof, const faceRealArray &mulf,
                             const faceRealArray &mutf, const cellRealArray &mul,
                             const cellRealArray &mut, const cellRealArray &rho) = 0;

        virtual void update() = 0;

        virtual void limit() = 0;

        /**
         * \brief Turbulence kinetic energy [m^2/s^2].
         */
        hur_nodiscard virtual realArray k() const = 0;

        /**
         * \brief Dissipation rate [m^2/s^3].
         */
        hur_nodiscard virtual realArray epsilon() const = 0;

        /**
         * \brief Turbulent Reynolds number [dimensionless].
         */
        hur_nodiscard virtual hur_nodiscard realArray Ret() const = 0;

        /**
         * \brief Kolmogorov length scale.
         */
        virtual hur_nodiscard realArray KolmogorovLengthScale() const;

        /**
         * \brief Kolmogorov time scale.
         */
        virtual hur_nodiscard realArray KolmogorovTimeScale() const;

        /**
         * \brief Kolmogorov velocity scale.
         */
        virtual hur_nodiscard realArray KolmogorovVelocityScale() const;

        /**
         * \brief Integral length scale.
         */
        virtual hur_nodiscard realArray integralLengthScale() const;

        /**
         * \brief Integral time scale.
         */
        virtual hur_nodiscard realArray integralTimeScale() const;

        /**
         * \brief Integral velocity scale.
         */
        virtual hur_nodiscard realArray integralVelocityScale() const;

        /**\brief The number of turbulence model equations.*/
        hur_nodiscard inline integer nEq() const noexcept { return nEq_; }

        virtual cellRealArray &var(const integer i) = 0;

        virtual const cellRealArray &var(const integer i) const = 0;

        virtual void solving(const realArray &dt);

        bool loop();

        virtual void initialize() {}

        virtual void initializeRestart() {}

        virtual void updateBoundary() {}
        virtual void limitAndUpdateBoundary() {}

        virtual void calcGrad(const spatialScheme &sps) {}

        virtual faceSymmTensorArray tauEff(const faceRealArray &rhof, const faceRealArray &mulf,
                                           const faceRealArray &mutf,
                                           const faceTensorArray &deltafV) const;

        /**
         * \brief Only available for RANS models.
         */
        virtual symmTensorArray ReynoldsStressTensor() const;

        /*!\brief Correct energy equation (molecular diffusion and turbulent transport in the energy equation)*/
        virtual void correctEnergyEqVisFlux(cellRealArray &E) const {}

    };
} // namespace OpenHurricane