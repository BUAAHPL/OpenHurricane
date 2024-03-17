#include "reactionTableCUDA.hpp"
/*!
 * \file reactionTableCUDA.inl
 * \brief The In-Line functions of the <i>reactionTableCUDA.hpp</i> file.
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

#ifdef CUDA_PARALLEL

inline cu_dual const OpenHurricane::cuChem::speciesTable &
OpenHurricane::cuChem::reactionTable::species() const {
    return species_;
}

inline OpenHurricane::cuChem::reactionTable::reactionTable()
    : nrc_(0), nsp_(0), fc_(), rc_(), sumStoi_(), reactionType_(), thirdBodyType_(), kf_(), kr_(),
      tbc_(), pdc_(), species_() {}

inline cu_host
OpenHurricane::cuChem::reactionTable::reactionTable(const cuChemInterface &reacTableInt,
                                                const speciesTable spt)
    : nrc_(reacTableInt.nrc_), nsp_(reacTableInt.nsp_),
      fc_(reacTableInt.nrc_, reacTableInt.reacInt_),
      rc_(reacTableInt.nrc_, reacTableInt.revReacInt_),
      sumStoi_(reacTableInt.nrc_, reacTableInt.sumStoiIntPtr_),
      reactionType_(reacTableInt.nrc_, reacTableInt.reactionTypeIntPtr_),
      thirdBodyType_(reacTableInt.nrc_, reacTableInt.thirdBodyTypeIntPtr_),
      kf_(reacTableInt.nrc_, reacTableInt.kfInt_), kr_(reacTableInt.nrc_, reacTableInt.rfInt_),
      tbc_(reacTableInt.nsp_, reacTableInt.nrc_, reacTableInt.thirdBInt_),
      pdc_(reacTableInt.nrc_, reacTableInt.pressDepInt_), species_(spt) {
    bool hasNonEqReac = false;
    for (cu_integer i = 0; i < nrc_; ++i) {
        if (reacTableInt.reactionTypeIntPtr_[i] == 2) {
            hasNonEqReac = true;
            break;
        }
    }
    if (!hasNonEqReac) {
        kr_.destroy();
    }
}

inline OpenHurricane::cuChem::reactionTable::reactionTable(const cuChemInterface &reacTableInt,
                                                       const speciesTable spt,
                                                       const cudaStreams &streams)
    : nrc_(reacTableInt.nrc_), nsp_(reacTableInt.nsp_),
      fc_(reacTableInt.nrc_, reacTableInt.reacInt_, streams),
      rc_(reacTableInt.nrc_, reacTableInt.revReacInt_, streams),
      sumStoi_(reacTableInt.nrc_, reacTableInt.sumStoiIntPtr_, streams()),
      reactionType_(reacTableInt.nrc_, reacTableInt.reactionTypeIntPtr_, streams()),
      thirdBodyType_(reacTableInt.nrc_, reacTableInt.thirdBodyTypeIntPtr_, streams()),
      kf_(reacTableInt.nrc_, reacTableInt.kfInt_, streams),
      kr_(reacTableInt.nrc_, reacTableInt.rfInt_, streams),
      tbc_(reacTableInt.nsp_, reacTableInt.nrc_, reacTableInt.thirdBInt_, streams),
      pdc_(reacTableInt.nrc_, reacTableInt.pressDepInt_, streams), species_(spt) {}

inline cu_dual OpenHurricane::cuChem::reactionTable::reactionTable(const reactionTable &rt)
    : nrc_(rt.nrc_), nsp_(rt.nsp_), fc_(rt.fc_), rc_(rt.rc_), sumStoi_(rt.sumStoi_),
      reactionType_(rt.reactionType_), thirdBodyType_(rt.thirdBodyType_), kf_(rt.kf_), kr_(rt.kr_),
      tbc_(rt.tbc_), pdc_(rt.pdc_), species_(rt.species_) {}

inline cu_dual OpenHurricane::cuChem::reactionTable::~reactionTable() noexcept {}

inline cu_host void OpenHurricane::cuChem::reactionTable::destroy() {
    fc_.destroy();
    rc_.destroy();
    destroyCuArray(sumStoi_);
    destroyCuArray(reactionType_);
    destroyCuArray(thirdBodyType_);
    kf_.destroy();
    kr_.destroy();
    tbc_.destroy();
    pdc_.destroy();
}

inline cu_device cu_real OpenHurricane::cuChem::reactionTable::Kc(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ GdRT) const {
    return Kp(irc, T, GdRT) * pow(CUDA_Patm / (CUDA_Ru * T), sumStoi_(irc));
}

inline cu_device cu_real OpenHurricane::cuChem::reactionTable::kf(const cu_ushort irc,
                                                                       const cu_real T) const {
    return kf_.k(irc, T);
}

inline cu_device cu_real OpenHurricane::cuChem::reactionTable::kr(const cu_ushort irc,
                                                                       const cu_real T) const {
    return kr_.k(irc, T);
}

inline cu_device cu_real OpenHurricane::cuChem::reactionTable::kr(const cu_real kf,
                                                                       const cu_real kc) const {
    return kf / cu_max(kc, cu_veryTiny);
}

inline cu_device cu_real OpenHurricane::cuChem::reactionTable::qfr(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real &kff, cu_real &krr) const {
    cu_real pf, pr;
    return qfr(irc, T, ci, GdRT, kff, krr, pf, pr);
}

inline cu_device cu_real OpenHurricane::cuChem::reactionTable::qfr(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT) const {
    cu_real kff, krr;
    return qfr(irc, T, ci, GdRT, kff, krr);
}

inline cu_device cu_ushort OpenHurricane::cuChem::reactionTable::nsp() const noexcept {
    return nsp_;
}

inline cu_device cu_ushort OpenHurricane::cuChem::reactionTable::nrc() const noexcept {
    return nrc_;
}

template <class dOmegaType>
inline cu_device cu_real OpenHurricane::cuChem::reactionTable::qfr(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real &kff, cu_real &krr, cu_real &qf,
    cu_real &qr, dOmegaType *__restrict__ dOmegaidCi) const {
    kfAndKr(irc, T, ci, GdRT, kff, krr);

    qf = kff;
    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto id = fc_.index_(i, irc);
        qf *= pow(ci[id], fc_.order_(i, irc));
    }
    qr = krr;
    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto id = rc_.index_(i, irc);
        qr *= pow(ci[id], rc_.order_(i, irc));
    }

    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto spi = fc_.index_(i, irc);
        const dOmegaType stoichi = static_cast<dOmegaType>(fc_.stoichCoeff_(i, irc));
        dOmegaType pf = static_cast<dOmegaType>(kff);
        for (cu_ushort k = 0; k < fc_.size(irc); ++k) {
            const auto sk = fc_.index_(k, irc);
            const dOmegaType orderk = static_cast<dOmegaType>(fc_.order_(k, irc));

            if (i == k) {
                if (orderk < dOmegaType(1)) {
                    if (static_cast<dOmegaType>(ci[sk]) > cu_tiny) {
                        pf *= orderk * pow(static_cast<dOmegaType>(ci[sk]), orderk - dOmegaType(1));
                    } else {
                        pf = dOmegaType(0);
                    }
                } else if (orderk == dOmegaType(1)) {
                    continue;
                } else {
                    pf *= orderk * pow(static_cast<dOmegaType>(ci[sk]), orderk - dOmegaType(1));
                }
            } else {
                pf *= pow(static_cast<dOmegaType>(ci[sk]), orderk);
            }
        }
        if (spi != nsp() - 1) {
            atomicAdd(&dOmegaidCi[spi], -stoichi * pf);
            for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
                const auto id1 = rc_.index_(j, irc);
                if (id1 == spi) {
                    const dOmegaType stoic1 = static_cast<dOmegaType>(rc_.stoichCoeff_(j, irc));
                    atomicAdd(&dOmegaidCi[spi], stoic1 * pf);
                }
            }
        } else //spi=Ns
        {
            const dOmegaType Wns = static_cast<dOmegaType>(species_.Wi(spi));
            for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
                const auto id1 = fc_.index_(j, irc);
                const dOmegaType stoic1 = static_cast<dOmegaType>(fc_.stoichCoeff_(j, irc));
                if (id1 != spi) {
                    atomicAdd(&dOmegaidCi[id1],
                              static_cast<dOmegaType>(species_.Wi(id1)) / Wns * stoic1 * pf);
                }
            }
            for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
                const auto id1 = rc_.index_(j, irc);
                const dOmegaType stoic1 = static_cast<dOmegaType>(rc_.stoichCoeff_(j, irc));
                if (id1 != spi) {
                    atomicAdd(&dOmegaidCi[id1],
                              -static_cast<dOmegaType>(species_.Wi(id1)) / Wns * stoic1 * pf);
                }
            }
        }
    }

    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto spi = rc_.index_(i, irc);
        const dOmegaType stoichi = static_cast<dOmegaType>(rc_.stoichCoeff_(i, irc));
        dOmegaType pr = static_cast<dOmegaType>(krr);
        for (cu_ushort k = 0; k < rc_.size(irc); ++k) {
            const auto sk = rc_.index_(k, irc);
            const dOmegaType orderk = static_cast<dOmegaType>(rc_.order_(k, irc));

            if (i == k) {
                if (orderk < dOmegaType(1)) {
                    if (static_cast<dOmegaType>(ci[sk]) > cu_tiny) {
                        pr *= orderk * pow(static_cast<dOmegaType>(ci[sk]), orderk - dOmegaType(1));
                    } else {
                        pr = dOmegaType(0);
                    }
                } else if (orderk == dOmegaType(1)) {
                    continue;
                } else {
                    pr *= orderk * pow(static_cast<dOmegaType>(ci[sk]), orderk - dOmegaType(1));
                }
            } else {
                pr *= pow(static_cast<dOmegaType>(ci[sk]), orderk);
            }
        }

        if (spi != nsp() - 1) {
            atomicAdd(&dOmegaidCi[spi], -stoichi * pr);
            for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
                const auto id1 = fc_.index_(j, irc);
                if (id1 == spi) {
                    const dOmegaType stoic1 = static_cast<dOmegaType>(fc_.stoichCoeff_(j, irc));
                    atomicAdd(&dOmegaidCi[spi], stoic1 * pr);
                }
            }
        } else //spi=Ns
        {
            const dOmegaType Wns = static_cast<dOmegaType>(species_.Wi(spi));
            for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
                const auto id1 = fc_.index_(j, irc);
                const dOmegaType stoic1 = static_cast<dOmegaType>(fc_.stoichCoeff_(j, irc));
                if (id1 != spi) {
                    atomicAdd(&dOmegaidCi[id1],
                              -static_cast<dOmegaType>(species_.Wi(id1)) / Wns * stoic1 * pr);
                }
            }
            for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
                const auto id1 = rc_.index_(j, irc);
                const dOmegaType stoic1 = static_cast<dOmegaType>(rc_.stoichCoeff_(j, irc));
                if (id1 != spi) {
                    atomicAdd(&dOmegaidCi[id1],
                              static_cast<dOmegaType>(species_.Wi(id1)) / Wns * stoic1 * pr);
                }
            }
        }
    }

    return qf - qr;
}

inline cu_device void
OpenHurricane::cuChem::reactionTable::Ri(const cu_ushort isp, const cu_real *__restrict__ omegai,
                                     cu_real *__restrict__ Rii) const {
    Rii[isp] = omegai[isp] * species_.Wi(isp);
}

template <class dOmegaType>
inline cu_device void OpenHurricane::cuChem::reactionTable::omegaCoupled3(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai,
    dOmegaType *__restrict__ dOmegaidCi) const {
    cu_real kf, kr, qf, qr;
    const auto qi = qfr<dOmegaType>(irc, T, ci, GdRT, kf, kr, qf, qr, dOmegaidCi);

    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto id = fc_.index_(i, irc);
        const auto stoic = fc_.stoichCoeff_(i, irc);
        atomicAdd(&omegai[id], -stoic * qi);
    }

    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto id = rc_.index_(i, irc);
        const auto stoic = rc_.stoichCoeff_(i, irc);
        const auto ord = rc_.order_(i, irc);
        atomicAdd(&omegai[id], stoic * qi);
    }
}
#endif // CUDA_PARALLEL