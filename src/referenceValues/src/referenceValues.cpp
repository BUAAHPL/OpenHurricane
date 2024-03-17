/*!
 * \file referenceValues.cpp
 * \brief The subroutines and functions of reference values
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

#include "referenceValues.hpp"
#include "controller.hpp"

OpenHurricane::referenceValues::referenceValues(const controller &cont)
    : Ma_(1.0), u_(1.0), p_(0.0), T_(288.16), rho_(1.225), mu_(1.7894e-5), gamma_(1.4), area_(1),
      length_(1), Tt_(300.0), R_(287.06) {
    reset(cont);
}

void OpenHurricane::referenceValues::reset(const controller &cont) {
    if (cont.found("givenBy")) {
        string gw = cont.findWord("givenBy");

        const auto &pcont = cont.parent();
        if (pcont.found("boundaryCondition")) {
            const auto &bcont = pcont.subController("boundaryCondition");
            if (bcont.found(gw)) {
                const auto &gbcont = bcont.subController(gw);
                gamma_ = gbcont.findOrDefault<real>("gamma", gamma_);
                p_ = gbcont.findOrDefault<real>("p", p_);
                rho_ = gbcont.findOrDefault<real>("rho", rho_);
                u_ = gbcont.findOrDefault<real>("v", u_);
                T_ = gbcont.findOrDefault<real>("T", T_);
                mu_ = gbcont.findOrDefault<real>("mu", mu_);

                R_ = gbcont.findOrDefault<real>("Rgas", R_);

                Ma_ = u_ / sqrt(gamma_ * R_ * T_);
            } else {
                LFatal("Cannot find boundary conditions: \"%s\"in ", gw.c_str(),
                       bcont.name().c_str());
            }
        } else {
            LFatal("Cannot find boundary conditions in %s", pcont.name().c_str());
        }
    } else {
        rho_ = cont.findOrDefault<real>("density", rho_);
        p_ = cont.findOrDefault<real>("pressure", p_);
        T_ = cont.findOrDefault<real>("temperature", T_);
        u_ = cont.findOrDefault<real>("velocity", u_);
        mu_ = cont.findOrDefault<real>("viscosity", mu_);
        gamma_ = cont.findOrDefault<real>("specificHeatRatio", gamma_);
        R_ = cont.findOrDefault<real>("gasConstantNumber", R_);
        Ma_ = u_ / sqrt(gamma_ * R_ * T_);
    }

    area_ = cont.findOrDefault<real>("area", area_);
    length_ = cont.findOrDefault<real>("length", length_);

    Tt_ = T_ * (1 + (gamma_ - 1) / 2.0 * sqr(Ma_));
}