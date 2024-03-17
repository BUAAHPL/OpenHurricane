/*!
 * \file perfectGas.cpp
 * \brief The subroutines and functions of class of equation of state of perfect gas.
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
#include "perfectGas.hpp"
#include "commonInclude.hpp"
namespace OpenHurricane {
    createClassNameStr(perfectGas, "perfectGas");
    registerObjFty(equationOfState, perfectGas, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::rhom(const real Pm, const real Tm,
                                                          const realArray &_yi) const noexcept {
    real Rm = species().Rm(_yi);
    return Pm / (Rm * Tm);
}

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::rhom(const real pm, const real Tm,
                                                          const PtrList<cellRealArray> &_yi,
                                                          const integer celli) const noexcept {
    real Rm = species().Rm(_yi, celli);
    return pm / (Rm * Tm);
}

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::pm(const real Rhom, const real Tm,
                                                        const realArray &_yi) const noexcept {
    real Rm = species().Rm(_yi);
    return Rhom * Rm * Tm;
}

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::pm(const real rhom, const real Tm,
                                                        const PtrList<cellRealArray> &_yi,
                                                        const integer celli) const noexcept {
    real Rm = species().Rm(_yi, celli);
    return rhom * Rm * Tm;
}

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::sm_p(const real p, const real T,
                                                          const realArray &yi) const {
    return Sm_p(p, T, species().Yi2Xi(yi)) / species().MWbyYi(yi);
}

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::sm_p(const real p, const real T,
                                                          const PtrList<cellRealArray> &yi,
                                                          const integer celli) const {
    return Sm_p(p, T, species().Yi2Xi(yi, celli)) / species().MWbyYi(yi, celli);
}

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::Sm_p(const real p, const real T,
                                                          const realArray &xi) const {
#ifdef HUR_DEBUG

    if (xi.size() != species().size()) {
        LFatal("The size of species list: %d is not equal to Xi(%d)", species().size(), xi.size());
    }
#endif // HUR_DEBUG
    if (species().size() == 1) {
        return Si_p(p, T, 0);
    }
    real Sm = 0.0;
    for (integer i = 0; i < xi.size(); ++i) {
        Sm += Si_p(p * xi[i], T, i) * xi[i];
    }
    return Sm;
}

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::sm_rho(const real rho, const real T,
                                                            const realArray &yi) const {
    return Sm_rho(rho, T, species().Yi2Xi(yi)) / species().MWbyYi(yi);
}

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::sm_rho(const real rho, const real T,
                                                            const PtrList<cellRealArray> &yi,
                                                            const integer celli) const {
    return Sm_rho(rho, T, species().Yi2Xi(yi, celli)) / species().MWbyYi(yi, celli);
}

hur_nodiscard OpenHurricane::real OpenHurricane::perfectGas::Sm_rho(const real rho, const real T,
                                                            const realArray &xi) const {
#ifdef HUR_DEBUG

    if (xi.size() != species().size()) {
        LFatal("The size of species list: %d is not equal to Xi(%d)", species().size(), xi.size());
    }

#endif // HUR_DEBUG

    if (species().size() == 1) {
        return Si_rho(rho, T, 0);
    }
    real Sm = 0.0;
    const auto yi = species().Xi2Yi(xi);
    for (integer i = 0; i < xi.size(); ++i) {
        Sm += Si_rho(rho * yi[i], T, i) * xi[i];
    }
    return Sm;
}