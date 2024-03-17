/*!
 * \file RANSModel.hpp
 * \brief Header of RANS turbulence model
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

#include "LUSGS.hpp"
#include "eddyViscosity.hpp"
#include "turbulenceModel.hpp"

namespace OpenHurricane {
    /*!\brief The class of RANS turbulence model.*/
    class RANSModel : public turbulenceModel {
    public:
        // Boundary specification method
        enum specificationMethod : short {
            viscosityRatio,
            intensityAndLength,
            intensityAndHYdraulicD,
            origianalTurbEquation,
            givenDirectly
        };

    private:
        /*!\brief Turbulent intensity.*/
        real intensity_;

        /*!\brief Turbulent viscosity ratio.*/
        real viscosityRatio_;

        /*!\brief Turbulent length scale.*/
        real lengthScale_;

        /*!\brief Turbulent hydraulic diameter.*/
        real hydraulicD_;

        /*!\brief The coefficient for limiting turbulence parameter increase.*/
        real coefTurb_;

    protected:
        void limitTurbVarIncrease(cellRealArray &cellQ, realArray &oldQ) const;

        /*!\brief The definition mehtod of turbulent boundary condition.*/
        specificationMethod bndSpMethod_;

    public:
        /*!\brief Construct from components.*/
        RANSModel(const controller &cont, flowModel &ev);

        /*!\brief Disallow default copy constructor.*/
        RANSModel(const RANSModel &) = delete;

        /*!\brief Disallow default bitwise assignment.*/
        void operator=(const RANSModel &) = delete;

        /*!\brief Destructor.*/
        virtual ~RANSModel() noexcept {}

        /*!\brief Return turbulent intensity.*/
        inline real intensity() { return intensity_; }

        /*!\brief Return viscosity ratio.*/
        inline real getViscosityRatio() { return viscosityRatio_; }

        /*!\brief Return length scale.*/
        inline real getLengthScale() { return lengthScale_; }

        /*!\brief Return hydraulic diameter.*/
        inline real getHydraulicD() { return hydraulicD_; }

        /*!\brief Return the specification method of turbulent boundary condition.*/
        inline specificationMethod bndSpMethod() { return bndSpMethod_; }

        /*!\brief Specific turbulence-relevant parameters specofication method.*/
        virtual void eqBndSpMehtod(){};

        /*!\brief Set turbulent inlet boundary condition.*/
        void inletBndSetting(controller &turbCont);

        /*!\brief Set turbulent complete boundary condition.*/
        virtual void bndValueSetting(controller &cont){};

    protected:
        bool averageLimit_;

        /*!\brief Limit cellQ to positive.
         * If the value is negative in the grid, then the average value of its adjacent grid is taken
         */
        void limitNegative(cellRealArray &cellQ, const real lowestV = tiny,
                           const bool average = true, const bool reportToScreen = false);
    };
} // namespace OpenHurricane
