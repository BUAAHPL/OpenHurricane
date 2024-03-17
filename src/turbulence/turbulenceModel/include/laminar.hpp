/*!
 * \file laminar.hpp
 * \brief Header of laminar
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

#include "turbulenceModel.hpp"

namespace OpenHurricane {

class laminar : public turbulenceModel {
private:
    // Private Member Functions

    // Disallow default bitwise copy construct and assignment
    laminar(const laminar &) = delete;
    void operator=(const laminar &) = delete;

public:
    declareClassNames;
    // Constructors

    /*!\brief Construct from components.*/
    laminar(const controller &cont, flowModel &ev);

    /*!\brief Destructor.*/
    virtual ~laminar() noexcept {}

    // Member Functions

    /*!\brief explicit source.*/
    virtual void expSource();

    /*!\brief implicit source.*/
    virtual void impSource();

    /*!\brief implicit source.*/
    virtual void fullImpSource(cellRealSquareMatrixArray &Jac,
                               const integer rhoId, const integer rhoTurb0);

    /*!\brief viscous flux.*/
    virtual void visFlux(const faceRealArray &rhof, const faceRealArray &mulf,
                         const faceRealArray &mutf, const cellRealArray &mul,
                         const cellRealArray &mut, const cellRealArray &rho);

    virtual void update();

    virtual void limit();

    virtual realArray k() const;

    virtual realArray epsilon() const;

    /**
     * \brief Turbulent Reynolds number [dimensionless].
     */
    virtual hur_nodiscard realArray Ret() const;

    virtual cellRealArray &var(const integer i);

    virtual const cellRealArray &var(const integer i) const;

    virtual faceSymmTensorArray tauEff(const faceRealArray &rhof,
                                       const faceRealArray &mulf,
                                       const faceRealArray &mutf,
                                       const faceTensorArray &deltafV) const;
};

} // namespace OpenHurricane
