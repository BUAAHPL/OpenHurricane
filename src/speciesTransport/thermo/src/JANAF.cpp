/*!
 * \file JANAF.cpp
 * \brief Main subroutines for JANAF.
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
#include "JANAF.hpp"
#include "commonInclude.hpp"
using namespace OpenHurricane::constant::physicalConstant;
namespace OpenHurricane {
    createClassNameStr(JANAF, "JANAF");
    registerObjFty(thermo, JANAF, controller);
} // namespace OpenHurricane

template <>
const std::streamsize
    OpenHurricane::JANAF::coeffArray::precision(OpenHurricane::feature<OpenHurricane::real>::precision);

OpenHurricane::JANAF::JANAF(const equationOfState &st, const integer id, const real Tlow,
                        const real Thigh, const real Tcommon,
                        const typename JANAF::coeffArray &highCpCoeffs,
                        const typename JANAF::coeffArray &lowCpCoeffs, const bool convertCoeffs)
    : thermo(st, id), TLow_(Tlow), THigh_(Thigh), TCommon_(Tcommon) {
    if (convertCoeffs) {
        for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel] * this->Ri();
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel] * this->Ri();
        }
    } else {
        for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
        }
    }
}

OpenHurricane::JANAF::JANAF(const equationOfState &st, const integer id, const real Tlow,
                        const real Thigh, const real Tcommon, const coeffArray &highCpCoeffs,
                        const coeffArray &lowCpCoeffs, const phaseType pt, const bool convertCoeffs)
    : thermo(st, id, pt), TLow_(Tlow), THigh_(Thigh), TCommon_(Tcommon) {
    if (convertCoeffs) {
        for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel] * this->Ri();
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel] * this->Ri();
        }
    } else {
        for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
        }
    }
}

OpenHurricane::JANAF::JANAF(const equationOfState &st, const integer id, const real Tlow,
                        const real Thigh, const real Tcommon, const coeffArray &highCpCoeffs,
                        const coeffArray &lowCpCoeffs, const std::string &pt,
                        const bool convertCoeffs)
    : thermo(st, id, pt), TLow_(Tlow), THigh_(Thigh), TCommon_(Tcommon) {
    if (convertCoeffs) {
        for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel] * this->Ri();
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel] * this->Ri();
        }
    } else {
        for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
        }
    }
}

OpenHurricane::JANAF::JANAF(const controller &cont, const equationOfState &st, const integer id)
    : thermo(cont, st, id),
      TLow_(cont.subController("thermodynamics").findType<real>("TLow", real(0.0))),
      THigh_(cont.subController("thermodynamics").findType<real>("THigh", real(0.0))),
      TCommon_(cont.subController("thermodynamics").findType<real>("TCommon", real(0.0))),
      highCpCoeffs_(
          cont.subController("thermodynamics").findType<coeffArray>("highCpCoeffs", coeffArray())),
      lowCpCoeffs_(
          cont.subController("thermodynamics").findType<coeffArray>("lowCpCoeffs", coeffArray())) {
    // Convert coefficients to mass-basis
    for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
        highCpCoeffs_[coefLabel] *= this->Ri();
        lowCpCoeffs_[coefLabel] *= this->Ri();
    }
}

OpenHurricane::JANAF::JANAF(const JANAF &jt)
    : thermo(jt), TLow_(jt.TLow_), THigh_(jt.THigh_), TCommon_(jt.TCommon_) {
    for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
        highCpCoeffs_[coefLabel] = jt.highCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = jt.lowCpCoeffs_[coefLabel];
    }
}

OpenHurricane::JANAF &OpenHurricane::JANAF::operator=(const JANAF &jt) {
    if (this != std::addressof(jt)) {
        thermo::operator=(jt);
        TLow_ = jt.TLow_;
        THigh_ = jt.THigh_;
        TCommon_ = jt.TCommon_;

        for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
            highCpCoeffs_[coefLabel] = jt.highCpCoeffs_[coefLabel];
            lowCpCoeffs_[coefLabel] = jt.lowCpCoeffs_[coefLabel];
        }
    }
    return *this;
}


OpenHurricane::JANAF::JANAF(const JANAF &jt, const string &name)
    : thermo(jt), TLow_(jt.TLow_), THigh_(jt.THigh_), TCommon_(jt.TCommon_) {
    for (integer coefLabel = 0; coefLabel < coeffArray::nElements_; coefLabel++) {
        highCpCoeffs_[coefLabel] = jt.highCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = jt.lowCpCoeffs_[coefLabel];
    }
}

OpenHurricane::real OpenHurricane::JANAF::limit(const real T, integer &iFlag) const noexcept {
    if (T < TLow_) {
        iFlag = 0;
        return TLow_;
    } else if (T > THigh_) {
        iFlag = 2;
        return THigh_;
    }

    iFlag = 1;
    return T;
}

hur_nodiscard OpenHurricane::real OpenHurricane::JANAF::cp0(const real T) const noexcept {
    const coeffArray &a = coeffs(T);
    const real T0 = limit(T);
    return ((((a[4] * T0 + a[3]) * T0 + a[2]) * T0 + a[1]) * T0 + a[0]);
}

hur_nodiscard OpenHurricane::real OpenHurricane::JANAF::ha0(const real T) const noexcept {
    const coeffArray &a = coeffs(T);

    if (T < TLow_) {
        real T0 = TLow_;

        real h0 =
            ((((a[4] / 5.0 * T0 + a[3] / 4.0) * T0 + a[2] / 3.0) * T0 + a[1] / 2.0) * T0 + a[0]) *
                T0 +
            a[5];
        real cp0 = ((((a[4] * T0 + a[3]) * T0 + a[2]) * T0 + a[1]) * T0 + a[0]);
        h0 += cp0 * (T - T0);

        return h0;
    } else if (T > THigh_) {
        real T0 = THigh_;

        real h0 =
            ((((a[4] / 5.0 * T0 + a[3] / 4.0) * T0 + a[2] / 3.0) * T0 + a[1] / 2.0) * T0 + a[0]) *
                T0 +
            a[5];
        real cp0 = ((((a[4] * T0 + a[3]) * T0 + a[2]) * T0 + a[1]) * T0 + a[0]);
        h0 += cp0 * (T - T0);

        return h0;
    } else {
        return (((((a[4] / 5.0 * T + a[3] / 4.0) * T + a[2] / 3.0) * T + a[1] / 2.0) * T + a[0]) *
                    T +
                a[5]);
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::JANAF::hc() const noexcept {
    const coeffArray &a = lowCpCoeffs_;
    return (((((a[4] / 5.0 * Tstd + a[3] / 4.0) * Tstd + a[2] / 3.0) * Tstd + a[1] / 2.0) * Tstd +
             a[0]) *
                Tstd +
            a[5]);
}

hur_nodiscard OpenHurricane::real OpenHurricane::JANAF::s0(const real T) const {
    const coeffArray &a = coeffs(T);
    return ((((a[4] / 4.0 * T + a[3] / 3.0) * T + a[2] / 2.0) * T + a[1]) * T +
            a[0] * log(max(T, tiny)) + a[6]);
}

hur_nodiscard OpenHurricane::real OpenHurricane::JANAF::g0(const real T) const {
    const coeffArray &a = coeffs(T);

    return (((((-a[4] / 20.0 * T - a[3] / 12.0) * T - a[2] / 6.0) * T - a[1] / 2.0) * T +
             a[0] * (1 - log(max(T, tiny))) - a[6]) *
                T +
            a[5]);
}

hur_nodiscard OpenHurricane::real OpenHurricane::JANAF::Dcp_pDT(const real pi, const real T)const noexcept {
    const coeffArray &a = coeffs(T);
    return (a[1] + T * (2.0 * a[2] + T * (3.0 * a[3] + 4.0 * a[4] * T)));
}
