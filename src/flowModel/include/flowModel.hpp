/*!
 * \file flowModel.hpp
 * \brief Header of flow model
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
#include "faceInterpolation.hpp"
#include "freeStream.hpp"
#include "objectFactory.hpp"
#include "rhoThermo.hpp"

namespace OpenHurricane {
    /*!\brief The class of flow model.*/
    class flowModel {
    private:
        // Private data

        /*!\brief Hold const reference to mesh*/
        const runtimeMesh &mesh_;

    protected:
        // Protected data

        /*!\brief Thermo pointer.*/
        rhoThermo *thermoPtr_;

        /*!\brief Velocity field.*/
        cellVectorArray v_;

        cellRealArray shockFactor_;

        /*!\brief Laminar Prandtl number.*/
        real Prl_;

        /**\brief Pressure lower limit.*/
        real pLow_;

        /**\brief Pressure lower limit.*/
        real pHigh_;

        /**\brief Temperature lower limit.*/
        real TLow_;

        /**\brief Temperature lower limit.*/
        real THigh_;

        /**\brief Reference values.*/
        const referenceValues &refValues_;

        /**\brief The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        integerArray temperatureFlag_;

        /**\brief The flag of state to evaluate pressure.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        integerArray pressureFlag_;

        /*!\brief The flag to adjust CFL number.
         * -# 1 - Do not adjust CFL number.
         * -# 0 - Reduce the CFL number.
         */
        cellIntegerArray CFLFlag_;

        mutable freeStream *freeStrPtr_;

    public:
        declareClassName(flowModel);

        declareObjFty(flowModel, controller, (const runtimeMesh &mesh, const controller &cont),
                      (mesh, cont));

        /*!\brief Null constructor.*/
        flowModel(const runtimeMesh &mesh);

        flowModel(const runtimeMesh &mesh, const controller &cont);

        hur_nodiscard static uniquePtr<flowModel> creator(const runtimeMesh &mesh,
                                                       const controller &cont);

        /*!\brief Destructor.*/
        inline virtual ~flowModel() noexcept {
            HurDelete(thermoPtr_);
            HurDelete(freeStrPtr_);
        }

        // Member functions

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept { return mesh_; }

        hur_nodiscard inline rhoThermo &thermo() noexcept { return *thermoPtr_; }

        hur_nodiscard inline const referenceValues &refValues() const noexcept {
            return refValues_;
        }

        /*!\brief Return const-reference to the density field.*/
        hur_nodiscard inline const cellRealArray &rho() const noexcept { return thermoPtr_->rho(); }

        /*!\brief Return access to the density field.*/
        hur_nodiscard inline cellRealArray &rho() noexcept { return thermoPtr_->rho(); }

        /*!\brief Return const-reference to the pressure field.*/
        hur_nodiscard inline const cellRealArray &p() const noexcept { return thermoPtr_->p(); }

        /*!\brief Return access to the pressure field.*/
        hur_nodiscard inline cellRealArray &p() noexcept { return thermoPtr_->p(); }

        /*!\brief Return const-reference to the temperature field.*/
        hur_nodiscard inline const cellRealArray &T() const noexcept { return thermoPtr_->T(); }

        /*!\brief Return access to the temperature field.*/
        hur_nodiscard inline cellRealArray &T() noexcept { return thermoPtr_->T(); }

        /*!\brief Return const-reference to the energy field.*/
        hur_nodiscard inline const cellRealArray &E() const noexcept { return thermoPtr_->E(); }

        /*!\brief Return access to the energy field.*/
        hur_nodiscard inline cellRealArray &E() noexcept { return thermoPtr_->E(); }

        /*!\brief Return const-reference to the velocity field.*/
        hur_nodiscard inline const cellVectorArray &v() const noexcept { return v_; }

        /*!\brief Return access to the velocity field.*/
        hur_nodiscard inline cellVectorArray &v() noexcept { return v_; }

        /*!\brief Return const-reference to the shock factor field.*/
        hur_nodiscard inline const cellRealArray &shockFactor() const noexcept {
            return shockFactor_;
        }

        /*!\brief Return access to the shock factor field.*/
        hur_nodiscard inline cellRealArray &shockFactor() noexcept { return shockFactor_; }

        hur_nodiscard inline cellRealArray &gama() noexcept { return thermoPtr_->gamma(); }

        hur_nodiscard inline const cellRealArray &gama() const noexcept {
            return thermoPtr_->gamma();
        }

        /**\brief The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        hur_nodiscard inline integerArray &temperatureFlag() noexcept { return temperatureFlag_; }

        /**\brief The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        hur_nodiscard inline const integerArray &temperatureFlag() const noexcept {
            return temperatureFlag_;
        }

        /**\brief The flag of state to evaluate pressure.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        hur_nodiscard inline integerArray &pressureFlag() noexcept { return pressureFlag_; }

        /**\brief The flag of state to evaluate pressure.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        hur_nodiscard inline const integerArray &pressureFlag() const noexcept {
            return pressureFlag_;
        }

        /*!\brief The flag to adjust CFL number.
         * -# 1 - Do not adjust CFL number.
         * -# 0 - Reduce the CFL number.
         */
        hur_nodiscard inline cellIntegerArray &CFLFlag() noexcept { return CFLFlag_; }

        /*!\brief The flag to adjust CFL number.
         * -# 1 - Do not adjust CFL number.
         * -# 0 - Reduce the CFL number.
         */
        hur_nodiscard inline const cellIntegerArray &CFLFlag() const noexcept { return CFLFlag_; }

        /*!\brief Return laminar dynamic viscosity field.*/
        virtual hur_nodiscard cellRealArray &mul() noexcept;

        /*!\brief Return laminar dynamic viscosity field.*/
        virtual hur_nodiscard const cellRealArray &mul() const noexcept;

        /*!\brief Return laminar dynamic viscosity field.*/
        virtual hur_nodiscard realArray mul(const integer faceZoneId);

        /*!\brief Return laminar dynamic viscosity field.*/
        virtual hur_nodiscard realArray mul(const integer faceZoneId) const;

        /*!\brief Return laminar kinematic viscosity field.*/
        virtual hur_nodiscard cellRealArray nu();

        /*!\brief Return laminar kinematic viscosity field.*/
        virtual hur_nodiscard realArray nu(const integer faceZoneId);

        /*!\briefReturn the effective dynamic viscosity.*/
        virtual hur_nodiscard cellRealArray muEff();

        /*!\brief Return turbulent dynamic viscosity field.*/
        virtual hur_nodiscard cellRealArray &mut() noexcept;

        /*!\brief Return turbulent dynamic viscosity field.*/
        virtual hur_nodiscard const cellRealArray &mut() const noexcept;

        /*!\brief Return laminar thermo conductivity field.*/
        virtual hur_nodiscard cellRealArray &kappal() noexcept;

        /*!\brief Return laminar thermo conductivity field.*/
        virtual hur_nodiscard const cellRealArray &kappal() const noexcept;

        /*!\briefReturn the effective thermo conductivity field.*/
        virtual hur_nodiscard cellRealArray kappaEff();

        /*!\brief Return turbulent thermo conductivity field.*/
        virtual hur_nodiscard cellRealArray &kappat() noexcept;

        /*!\briefReturn the effective mu/Pr.*/
        virtual hur_nodiscard cellRealArray keEff() const;

        /*!\briefReturn the effective mu/Pr.*/
        virtual hur_nodiscard real keEff(const integer cellI) const;

        /*!\brief Return access to mixture pointer.*/
        hur_nodiscard inline mixture &mixtures() noexcept { return thermoPtr_->mixtures(); }

        /*!\brief Return access to mixture pointer.*/
        hur_nodiscard inline const mixture &mixtures() const noexcept {
            return thermoPtr_->mixtures();
        }

        /*!
         * \brief The Takeno Flame Index.
         * \param[in] fuelName - The name list of fuel
         * \param[in] oxygenName - The name list of oxygen
         * \return The Takeno Flame Index
         * \retval Return a real array
         */
        hur_nodiscard inline realArray GFO(const stringList &fuelName,
                                           const stringList &oxygenName) const {
            return mixtures().GFO(fuelName, oxygenName);
        }

        /*!
         * \brief The normalizing, Takeno Flame Index.
         * \param[in] fuelName - The name list of fuel
         * \param[in] oxygenName - The name list of oxygen
         * \return The Takeno Flame Index: 1.0 - indicates premixed combustion; -1.0 indicates non-premixed combustion; 0 - indicates no flame
         * \retval Return a real array
         */
        hur_nodiscard inline realArray nGFO(const stringList &fuelName,
                                            const stringList &oxygenName) const {
            return mixtures().nGFO(fuelName, oxygenName);
        }

        /*!\brief Lower limit of mut.*/
        hur_nodiscard inline virtual real mutLow() const noexcept { return real(0.0); }

        /*!\brief Higher limit of mut.*/
        hur_nodiscard inline virtual real mutHigh() const noexcept { return real(0.0); }

        /*!\brief Laminar Prandtl number.*/
        hur_nodiscard inline real Prl() const noexcept { return Prl_; }

        /*!\brief Schmidt number.*/
        virtual hur_nodiscard real Sct() const noexcept { return real(0.0); }

        /*!\brief Turbulent Prandtl number*/
        virtual hur_nodiscard real Prt() const noexcept { return real(0.0); }

        /**\brief Pressure lower limit (dimensionless).*/
        hur_nodiscard inline real pLow() const noexcept { return pLow_; }

        /**\brief Pressure lower limit (dimensionless).*/
        hur_nodiscard inline real pHigh() const noexcept { return pHigh_; }

        /**\brief Temperature lower limit (dimensionless).*/
        hur_nodiscard inline real TLow() const noexcept { return TLow_; }

        /**\brief Temperature lower limit (dimensionless).*/
        hur_nodiscard inline real THigh() const noexcept { return THigh_; }

        /**
         * \brief Free stream.
         */
        hur_nodiscard const freeStream &freeStr() const;

        hur_nodiscard inline virtual string typeName() const { return className_; }
    };
} // namespace OpenHurricane