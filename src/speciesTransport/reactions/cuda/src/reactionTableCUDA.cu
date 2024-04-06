/*!
 * \file reactionTableCUDA.cu
 * \brief The subroutines and functions of reaction Table in CUDA platform.
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

#include "reactionTableCUDA.hpp"
#include <cmath>
#ifdef CUDA_PARALLEL

cu_device cu_real OpenHurricane::cuChem::reactionTable::Kp(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ GdRT) const {
    cu_real Kptemp = 0;
    for (cu_ushort i = 0; i < fc_.size_(irc); ++i) {
        const cu_ushort index = fc_.index_(i, irc);
        Kptemp -= GdRT[index] * fc_.stoichCoeff_(i, irc);
    }
    for (cu_ushort i = 0; i < rc_.size_(irc); ++i) {
        const cu_ushort index = rc_.index_(i, irc);
        Kptemp += GdRT[index] * rc_.stoichCoeff_(i, irc);
    }

    if (Kptemp > 600) {
        return cu_rootLarge;
    } else {
        return exp(Kptemp);
    }
}

cu_device void OpenHurricane::cuChem::reactionTable::kfAndKr(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real &kff, cu_real &krr) const {
    kff = kf(irc, T);
    krr = 0;
    if (reactionType_(irc) == irreversibleReaction) {
        krr = 0;
    } else if (reactionType_(irc) == nonEquilibriumReversibleReaction) {
        krr = kr(irc, T);
    } else if (reactionType_(irc) == reversibleReaction) {
        krr = kr(kff, Kc(irc, T, GdRT));
    }

    if (thirdBodyType_(irc) == noThirdBody) {
        // Do nothing
    } else if (thirdBodyType_(irc) == thirdBody) {
        const auto M = tbc_.thirdBodyEfficiencies(ci, irc, nsp());
        kff *= M;
        krr *= M;
    } else if (thirdBodyType_(irc) == UnimolecularFallOff) {
        const auto k0 = pdc_.k(irc, T);
        const auto kinf = kff;
        const auto M = tbc_.thirdBodyEfficiencies(ci, irc, nsp());
        const auto pr = pdc_.Pr(k0, kinf, M);

        const auto Mf = pdc_.UnimolecularFallOffFactor(pr, pdc_.F(irc, T, pr));
        kff *= Mf;
        krr *= Mf;
    } else if (thirdBodyType_(irc) == bimolecularFallOff) {
        const auto kinf = pdc_.k(irc, T);
        const auto k0 = kff;
        const auto M = tbc_.thirdBodyEfficiencies(ci, irc, nsp());
        const auto pr = pdc_.Pr(k0, kinf, M);

        const auto Mf = pdc_.BimolecularFallOffFactor(pr, pdc_.F(irc, T, pr));
        kff *= Mf;
        krr *= Mf;
    }
}

cu_device cu_real OpenHurricane::cuChem::reactionTable::qfr(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real &kff, cu_real &krr, cu_real &qf,
    cu_real &qr) const {
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

    return qf - qr;
}

//cu_device cu_real OpenHurricane::cuChem::reactionTable::qfr(
//    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
//    const cu_real *__restrict__ GdRT, cu_real &kff, cu_real &krr, cu_real &qf,
//    cu_real &qr, cu_real *__restrict__ dOmegaidCi) const {
//    kfAndKr(irc, T, ci, GdRT, kff, krr);
//
//    qf = kff;
//    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
//        const auto id = fc_.index_(i, irc);
//        qf *= pow(ci[id], fc_.order_(i, irc));
//    }
//
//    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
//        const auto spi = fc_.index_(i, irc);
//        const auto stoichi = fc_.stoichCoeff_(i, irc);
//        auto pf = kff;
//        for (cu_ushort k = 0; k < fc_.size(irc); ++k) {
//            const auto sk = fc_.index_(k, irc);
//            const auto orderk = fc_.order_(k, irc);
//
//            if (i == k) {
//                if (orderk < 1.0) {
//                    if (ci[sk] > cu_tiny) {
//                        pf *= orderk * pow(ci[sk] + cu_veryTiny, orderk - cu_real(1.0));
//                    } else {
//                        pf = cu_real(0.0);
//                    }
//                } else if (orderk == 1.0) {
//                    continue;
//                } else {
//                    pf *= orderk * pow(ci[sk], orderk - cu_real(1.0));
//                }
//            } else {
//                pf *= pow(ci[sk], orderk);
//            }
//        }
//        if (spi != nsp() - 1) {
//            atomicAdd(&dOmegaidCi[spi], -stoichi * pf);
//            for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
//                const auto id1 = rc_.index_(j, irc);
//                if (id1 == spi) {
//                    const auto stoic1 = rc_.stoichCoeff_(j, irc);
//                    atomicAdd(&dOmegaidCi[spi], stoic1 * pf);
//                }
//            }
//        } else //spi=Ns
//        {
//            const cu_real Wns = species_.Wi(spi);
//            for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
//                const auto id1 = fc_.index_(j, irc);
//                const auto stoic1 = fc_.stoichCoeff_(j, irc);
//                if (id1 != spi) {
//                    atomicAdd(&dOmegaidCi[id1], species_.Wi(id1) / Wns * stoic1 * pf);
//                }
//            }
//            for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
//                const auto id1 = rc_.index_(j, irc);
//                const auto stoic1 = rc_.stoichCoeff_(j, irc);
//                if (id1 != spi) {
//                    atomicAdd(&dOmegaidCi[id1], -species_.Wi(id1) / Wns * stoic1 * pf);
//                }
//            }
//        }
//    }
//
//    qr = krr;
//    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
//        const auto id = rc_.index_(i, irc);
//        qr *= pow(ci[id], rc_.order_(i, irc));
//    }
//
//    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
//        const auto spi = rc_.index_(i, irc);
//        const auto stoichi = rc_.stoichCoeff_(i, irc);
//        auto pr = krr;
//        for (cu_ushort k = 0; k < rc_.size(irc); ++k) {
//            const auto sk = rc_.index_(k, irc);
//            const auto orderk = rc_.order_(k, irc);
//
//            if (i == k) {
//                if (orderk < 1.0) {
//                    if (ci[sk] > cu_tiny) {
//                        pr *= orderk * pow(ci[sk] + cu_veryTiny, orderk - cu_real(1.0));
//                    } else {
//                        pr = cu_real(0.0);
//                    }
//                } else if (orderk == 1.0) {
//                    continue;
//                } else {
//                    pr *= orderk * pow(ci[sk], orderk - cu_real(1.0));
//                }
//            } else {
//                pr *= pow(ci[sk], orderk);
//            }
//        }
//
//        if (spi != nsp() - 1) {
//            atomicAdd(&dOmegaidCi[spi], -stoichi * pr);
//            for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
//                const auto id1 = fc_.index_(j, irc);
//                if (id1 == spi) {
//                    const auto stoic1 = fc_.stoichCoeff_(j, irc);
//                    atomicAdd(&dOmegaidCi[spi], stoic1 * pr);
//                }
//            }
//        } else //spi=Ns
//        {
//            const cu_real Wns = species_.Wi(spi);
//            for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
//                const auto id1 = fc_.index_(j, irc);
//                const auto stoic1 = fc_.stoichCoeff_(j, irc);
//                if (id1 != spi) {
//                    atomicAdd(&dOmegaidCi[id1], -species_.Wi(id1) / Wns * stoic1 * pr);
//                }
//            }
//            for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
//                const auto id1 = rc_.index_(j, irc);
//                const auto stoic1 = rc_.stoichCoeff_(j, irc);
//                if (id1 != spi) {
//                    atomicAdd(&dOmegaidCi[id1], species_.Wi(id1) / Wns * stoic1 * pr);
//                }
//            }
//        }
//    }
//
//    return qf - qr;
//}

cu_device cu_real OpenHurricane::cuChem::reactionTable::qfrAllSpecies(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real &kff, cu_real &krr, cu_real &qf,
    cu_real &qr, cu_real *__restrict__ dOmegaidCi) const {
    kfAndKr(irc, T, ci, GdRT, kff, krr);

    qf = kff;
    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto id = fc_.index_(i, irc);
        qf *= pow(ci[id], fc_.order_(i, irc));
    }

    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto spi = fc_.index_(i, irc);
        const auto stoichi = fc_.stoichCoeff_(i, irc);
        auto pf = kff;
        for (cu_ushort k = 0; k < fc_.size(irc); ++k) {
            const auto sk = fc_.index_(k, irc);
            const auto orderk = fc_.order_(k, irc);

            if (i == k) {
                if (orderk < 1.0) {
                    if (ci[sk] > cu_tiny) {
                        pf *= orderk * pow(ci[sk] + cu_veryTiny, orderk - cu_real(1.0));
                    } else {
                        pf = cu_real(0.0);
                    }
                } else if (orderk == 1.0) {
                    continue;
                } else {
                    pf *= orderk * pow(ci[sk], orderk - cu_real(1.0));
                }
            } else {
                pf *= pow(ci[sk], orderk);
            }
        }
        atomicAdd(&dOmegaidCi[spi], -stoichi * pf);
        for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
            const auto id1 = rc_.index_(j, irc);
            if (id1 == spi) {
                const auto stoic1 = rc_.stoichCoeff_(j, irc);
                atomicAdd(&dOmegaidCi[spi], stoic1 * pf);
            }
        }
    }

    qr = krr;
    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto id = rc_.index_(i, irc);
        qr *= pow(ci[id], rc_.order_(i, irc));
    }

    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto spi = rc_.index_(i, irc);
        const auto stoichi = rc_.stoichCoeff_(i, irc);
        auto pr = krr;
        for (cu_ushort k = 0; k < rc_.size(irc); ++k) {
            const auto sk = rc_.index_(k, irc);
            const auto orderk = rc_.order_(k, irc);

            if (i == k) {
                if (orderk < 1.0) {
                    if (ci[sk] > cu_tiny) {
                        pr *= orderk * pow(ci[sk] + cu_veryTiny, orderk - cu_real(1.0));
                    } else {
                        pr = cu_real(0.0);
                    }
                } else if (orderk == 1.0) {
                    continue;
                } else {
                    pr *= orderk * pow(ci[sk], orderk - cu_real(1.0));
                }
            } else {
                pr *= pow(ci[sk], orderk);
            }
        }
        atomicAdd(&dOmegaidCi[spi], -stoichi * pr);
        for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
            const auto id1 = fc_.index_(j, irc);
            if (id1 == spi) {
                const auto stoic1 = fc_.stoichCoeff_(j, irc);
                atomicAdd(&dOmegaidCi[spi], stoic1 * pr);
            }
        }
    }

    return qf - qr;
}

cu_device void OpenHurricane::cuChem::reactionTable::dWdciAllSpecies(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real &kff, cu_real &krr,
    cu_real *__restrict__ dOmegaidCi) const {
    kfAndKr(irc, T, ci, GdRT, kff, krr);

    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto spi = fc_.index_(i, irc);
        const auto stoichi = fc_.stoichCoeff_(i, irc);
        auto pf = kff;
        for (cu_ushort k = 0; k < fc_.size(irc); ++k) {
            const auto sk = fc_.index_(k, irc);
            const auto orderk = fc_.order_(k, irc);

            if (i == k) {
                if (orderk < 1.0) {
                    if (ci[sk] > cu_tiny) {
                        pf *= orderk * pow(ci[sk] + cu_veryTiny, orderk - cu_real(1.0));
                    } else {
                        pf = cu_real(0.0);
                    }
                } else if (orderk == 1.0) {
                    continue;
                } else {
                    pf *= orderk * pow(ci[sk], orderk - cu_real(1.0));
                }
            } else {
                pf *= pow(ci[sk], orderk);
            }
        }
        atomicAdd(&dOmegaidCi[spi], -stoichi * pf);
        for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
            const auto id1 = rc_.index_(j, irc);
            if (id1 == spi) {
                const auto stoic1 = rc_.stoichCoeff_(j, irc);
                atomicAdd(&dOmegaidCi[spi], stoic1 * pf);
            }
        }
    }

    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto spi = rc_.index_(i, irc);
        const auto stoichi = rc_.stoichCoeff_(i, irc);
        auto pr = krr;
        for (cu_ushort k = 0; k < rc_.size(irc); ++k) {
            const auto sk = rc_.index_(k, irc);
            const auto orderk = rc_.order_(k, irc);

            if (i == k) {
                if (orderk < 1.0) {
                    if (ci[sk] > cu_tiny) {
                        pr *= orderk * pow(ci[sk] + cu_veryTiny, orderk - cu_real(1.0));
                    } else {
                        pr = cu_real(0.0);
                    }
                } else if (orderk == 1.0) {
                    continue;
                } else {
                    pr *= orderk * pow(ci[sk], orderk - cu_real(1.0));
                }
            } else {
                pr *= pow(ci[sk], orderk);
            }
        }
        atomicAdd(&dOmegaidCi[spi], -stoichi * pr);
        for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
            const auto id1 = fc_.index_(j, irc);
            if (id1 == spi) {
                const auto stoic1 = fc_.stoichCoeff_(j, irc);
                atomicAdd(&dOmegaidCi[spi], stoic1 * pr);
            }
        }
    }
}

cu_device void OpenHurricane::cuChem::reactionTable::omegaCoupled(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai) const {
    const auto qi = qfr(irc, T, ci, GdRT);
    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto id = fc_.index_(i, irc);
        const auto stoic = fc_.stoichCoeff_(i, irc);
        if (id < nsp() - 1) {
            atomicAdd(&omegai[id], -stoic * qi);
        }
    }

    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto id = rc_.index_(i, irc);
        const auto stoic = rc_.stoichCoeff_(i, irc);
        if (id < nsp() - 1) {
            atomicAdd(&omegai[id], stoic * qi);
        }
    }
}

cu_device void OpenHurricane::cuChem::reactionTable::omega(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai) const {
    const auto qi = qfr(irc, T, ci, GdRT);

    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto id = fc_.index_(i, irc);
        const auto stoic = fc_.stoichCoeff_(i, irc);
        atomicAdd(&omegai[id], -stoic * qi);
    }

    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto id = rc_.index_(i, irc);
        const auto stoic = rc_.stoichCoeff_(i, irc);
        atomicAdd(&omegai[id], stoic * qi);
    }
}

cu_device void OpenHurricane::cuChem::reactionTable::omegaCoupled(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai,
    cu_real *__restrict__ dOmegaidCi) const {
    cu_real kf, kr, qf, qr;
    const auto qi = qfr(irc, T, ci, GdRT, kf, kr, qf, qr);

    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto id = fc_.index_(i, irc);
        const auto stoic = fc_.stoichCoeff_(i, irc);
        const auto ord = fc_.order_(i, irc);
        atomicAdd(&omegai[id], -stoic * qi);
        auto cii = ci[id];
        if (cii != 0) {
            if (id < nsp_ - 1) {
                atomicAdd(&dOmegaidCi[id], -stoic * ord * qf / cii);
            }
        }
    }

    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto id = rc_.index_(i, irc);
        const auto stoic = rc_.stoichCoeff_(i, irc);
        const auto ord = rc_.order_(i, irc);
        atomicAdd(&omegai[id], stoic * qi);
        auto cii = ci[id];
        if (cii != 0) {
            if (id < nsp_ - 1) {
                atomicAdd(&dOmegaidCi[id], -stoic * ord * qr / cii);
            }
        }
    }
}

cu_device void OpenHurricane::cuChem::reactionTable::omegaCoupled2(
    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
    const cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai,
    cu_real *__restrict__ dOmegaidCi) const {
    cu_real kf, kr, qf, qr;
    const auto qi = qfr(irc, T, ci, GdRT, kf, kr, qf, qr);

    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
        const auto id = fc_.index_(i, irc);
        const auto stoic = fc_.stoichCoeff_(i, irc);
        const auto ord = fc_.order_(i, irc);

        atomicAdd(&omegai[id], -stoic * qi);
        auto cii = ci[id];
        if (cii != 0) {
            if (id < nsp_ - 1) {
                atomicAdd(&dOmegaidCi[id], -stoic * ord * qf / cii);
                for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
                    const auto id1 = rc_.index_(j, irc);
                    if (id1 == id) {
                        const auto stoic1 = rc_.stoichCoeff_(j, irc);
                        atomicAdd(&dOmegaidCi[id], stoic1 * ord * qf / cii);
                    }
                }
            } else {
                const cu_real dwidcns = ord * qf / cii;
                const cu_real Wns = species_.Wi(id);
                for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
                    const auto id1 = fc_.index_(j, irc);
                    const auto stoic1 = fc_.stoichCoeff_(j, irc);
                    if (id1 != id) {
                        atomicAdd(&dOmegaidCi[id1], species_.Wi(id1) / Wns * stoic1 * dwidcns);
                    }
                }
                for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
                    const auto id1 = rc_.index_(j, irc);
                    const auto stoic1 = rc_.stoichCoeff_(j, irc);
                    if (id1 != id) {
                        atomicAdd(&dOmegaidCi[id1], -species_.Wi(id1) / Wns * stoic1 * dwidcns);
                    }
                }
            }
        }
    }

    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
        const auto id = rc_.index_(i, irc);
        const auto stoic = rc_.stoichCoeff_(i, irc);
        const auto ord = rc_.order_(i, irc);

        atomicAdd(&omegai[id], stoic * qi);
        auto cii = ci[id];
        if (cii != 0) {
            if (id < nsp_ - 1) {
                atomicAdd(&dOmegaidCi[id], -stoic * ord * qr / cii);
                for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
                    const auto id1 = fc_.index_(j, irc);
                    if (id1 == id) {
                        const auto stoic1 = fc_.stoichCoeff_(j, irc);
                        atomicAdd(&dOmegaidCi[id], stoic1 * ord * qr / cii);
                    }
                }
            } else {
                const cu_real dwidcns = ord * qr / cii;
                const cu_real Wns = species_.Wi(id);
                for (cu_ushort j = 0; j < fc_.size(irc); ++j) {
                    const auto id1 = fc_.index_(j, irc);
                    const auto stoic1 = fc_.stoichCoeff_(j, irc);
                    if (id1 != id) {
                        atomicAdd(&dOmegaidCi[id1], -species_.Wi(id1) / Wns * stoic1 * dwidcns);
                    }
                }
                for (cu_ushort j = 0; j < rc_.size(irc); ++j) {
                    const auto id1 = rc_.index_(j, irc);
                    const auto stoic1 = rc_.stoichCoeff_(j, irc);
                    if (id1 != id) {
                        atomicAdd(&dOmegaidCi[id1], species_.Wi(id1) / Wns * stoic1 * dwidcns);
                    }
                }
            }
        }
    }
}

//cu_device void OpenHurricane::cuChem::reactionTable::omegaCoupled3(
//    const cu_ushort irc, const cu_real T, const cu_real *__restrict__ ci,
//    const cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai,
//    cu_real *__restrict__ dOmegaidCi) const {
//    cu_real kf, kr, qf, qr;
//    const auto qi = qfr(irc, T, ci, GdRT, kf, kr, qf, qr, dOmegaidCi);
//
//    for (cu_ushort i = 0; i < fc_.size(irc); ++i) {
//        const auto id = fc_.index_(i, irc);
//        const auto stoic = fc_.stoichCoeff_(i, irc);
//        atomicAdd(&omegai[id], -stoic * qi);
//    }
//
//    for (cu_ushort i = 0; i < rc_.size(irc); ++i) {
//        const auto id = rc_.index_(i, irc);
//        const auto stoic = rc_.stoichCoeff_(i, irc);
//        const auto ord = rc_.order_(i, irc);
//        atomicAdd(&omegai[id], stoic * qi);
//    }
//}

#endif // CUDA_PARALLEL