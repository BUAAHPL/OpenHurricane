/*!
 * \file reactionTypes.cpp
 * \brief The subroutines and functions of base class of reaction.
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
#include "reactionTypes.hpp"
#include "reactionList.hpp"
namespace OpenHurricane {
    createClassName(reactionTypes);
         createObjFty(reactionTypes, controller);
} // namespace OpenHurricane

OpenHurricane::reactionTypes::reactionTypes(const string &names, const speciesList &specie,
                                        const reactionList &rt, const integer id,
                                        const List<reactionSpeciesCoeffs> &forwardCoe,
                                        const List<reactionSpeciesCoeffs> &backwardCoe)
    : name_(names), species_(specie), G0_(rt.G0()), DG0DT_(rt.DG0DT()), index_(id),
      forwardCoeffs_(forwardCoe), backwardCoeffs_(backwardCoe), bfCoeffs_(), sumStoi_(Zero),
      reacType_(types::irreversible) {
    setBFCoeffs();
}

OpenHurricane::reactionTypes::reactionTypes(const string &names, const speciesList &specie,
                                        const reactionList &rt, const integer id)
    : name_(names), species_(specie), G0_(rt.G0()), DG0DT_(rt.DG0DT()), index_(id), sumStoi_(Zero),
      reacType_(types::irreversible) {}

OpenHurricane::reactionTypes::reactionTypes(const speciesList &sp, const reactionList &rt,
                                        const controller &cont)
    : name_(cont.findWord("name")), species_(sp), G0_(rt.G0()), DG0DT_(rt.DG0DT()),
      index_(cont.findType<integer>("index", -1)), sumStoi_(Zero), reacType_(types::irreversible) {
    realList forCoefTemp = cont.findType<realList>("forwardCoeffs", realList());
    realList bacCoefTemp = cont.findType<realList>("backwardCoeffs", realList());

    if (forCoefTemp.size() != sp.size() && bacCoefTemp.size() != sp.size()) {
        errorAbortStr(("The size of stoichiometry coeffs list in reaction: " + toString(index_) +
                       " is not equal to the species list!"));
    }
    integer countFC = Zero;
    integer countBC = Zero;
    for (integer i = 0; i < sp.size(); ++i) {
        if (forCoefTemp[i] != real(0.0)) {
            countFC++;
        }

        if (bacCoefTemp[i] != real(0.0)) {
            countBC++;
        }
    }
    forwardCoeffs_.resize(countFC);
    backwardCoeffs_.resize(countBC);

    countFC = Zero;
    countBC = Zero;
    for (integer i = 0; i < sp.size(); ++i) {
        if (forCoefTemp[i] != real(0.0)) {
            forwardCoeffs_[countFC].index_ = i;
            forwardCoeffs_[countFC].stoichCoeff_ = forCoefTemp[i];
            countFC++;
        }

        if (bacCoefTemp[i] != real(0.0)) {
            backwardCoeffs_[countBC].index_ = i;
            backwardCoeffs_[countBC].stoichCoeff_ = bacCoefTemp[i];
            countBC++;
        }
    }

    if (cont.found("forwardOrder")) {
        realList forwardOrderTemp = cont.findType<realList>("forwardOrder", realList());
        if (forwardOrderTemp.size() != sp.size()) {
            errorAbortStr(("The size of forward stoichiometry order list in reaction: " +
                           toString(index_) + " is not equal to the species list!"));
        }

        for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
            integer id = forwardCoeffs_[i].index_;
            forwardCoeffs_[i].order_ = forwardOrderTemp[i];
        }
    }

    if (cont.found("backwardOrder")) {
        realList backwardOrderTemp = cont.findType<realList>("backwardOrder", realList());
        if (backwardOrderTemp.size() != sp.size()) {
            errorAbortStr(("The size of backward stoichiometry order list in reaction: " +
                           toString(index_) + " is not equal to the species list!"));
        }

        for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
            integer id = backwardCoeffs_[i].index_;
            backwardCoeffs_[i].order_ = backwardOrderTemp[i];
        }
    }

    setBFCoeffs();
}

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::reactionTypes>
OpenHurricane::reactionTypes::creator(const speciesList &sp, const reactionList &rt,
                                  const controller &cont) {
    string rType = cont.findWord(reactionTypes::className_);
    defineInObjCreator(reactionTypes, rType, controller, (sp, rt, cont));
}

void OpenHurricane::reactionTypes::setBFCoeffs() {
    List<reactionSpeciesCoeffs> tempBFL(species_.size());

    for (integer i = 0; i < tempBFL.size(); ++i) {
        tempBFL[i].index_ = i;
        tempBFL[i].stoichCoeff_ = Zero;
        tempBFL[i].order_ = Zero;
    }

    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        integer index = backwardCoeffs_[i].index_;

        tempBFL[index].stoichCoeff_ = backwardCoeffs_[i].stoichCoeff_;
        tempBFL[index].order_ = backwardCoeffs_[i].order_;
    }

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        integer index = forwardCoeffs_[i].index_;

        tempBFL[index].stoichCoeff_ -= forwardCoeffs_[i].stoichCoeff_;
        tempBFL[index].order_ -= forwardCoeffs_[i].order_;
    }

    integer countBF = Zero;

    for (integer i = 0; i < tempBFL.size(); ++i) {
        if (tempBFL[i].stoichCoeff_ != real(0.0)) {
            countBF++;
        }
    }

    sumStoi_ = Zero;
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        sumStoi_ -= forwardCoeffs_[i].stoichCoeff_;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        sumStoi_ += backwardCoeffs_[i].stoichCoeff_;
    }
    if (countBF > 0) {
        bfCoeffs_.resize(countBF);
        countBF = Zero;
        for (integer i = 0; i < tempBFL.size(); ++i) {
            if (tempBFL[i].stoichCoeff_ != real(0.0)) {
                bfCoeffs_[countBF] = tempBFL[i];
                countBF++;
            }
        }
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::reactionTypes::Kp(const real T) const {
    real Kptemp = Zero;
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        integer index = forwardCoeffs_[i].index_;
        Kptemp -= G0_[index] * forwardCoeffs_[i].stoichCoeff_;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        integer index = backwardCoeffs_[i].index_;
        Kptemp += G0_[index] * backwardCoeffs_[i].stoichCoeff_;
    }
    Kptemp /= (max(T, tiny) * constant::physicalConstant::Ru);

    const real expKp = exp(-Kptemp);

    if (std::isnan(expKp) || isinf(expKp)) {
        return rootVeryLarge;
    }
    return expKp;
}

hur_nodiscard OpenHurricane::real OpenHurricane::reactionTypes::DKpDT(const real T) const {
    const real kp = Kp(T);
    real dKptempDT = Zero;
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        integer index = forwardCoeffs_[i].index_;
        dKptempDT -= (DG0DT_[index] - G0_[index] / max(T, tiny)) * forwardCoeffs_[i].stoichCoeff_;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        integer index = backwardCoeffs_[i].index_;
        dKptempDT += (DG0DT_[index] - G0_[index] / max(T, tiny)) * backwardCoeffs_[i].stoichCoeff_;
    }
    dKptempDT /= (max(T, tiny) * constant::physicalConstant::Ru);
    const real dkpdt = kp * dKptempDT;
    if (isnan(dkpdt) || isinf(dkpdt)) {
        return veryLarge;
    }
    return dkpdt;
}

hur_nodiscard OpenHurricane::real OpenHurricane::reactionTypes::Kc(const real T) const {
    if (fabs(sumStoi_) < veryTiny) {
        return Kp(T);
    } else {
        return min(Kp(T) *
                       pow(constant::physicalConstant::Patm / (constant::physicalConstant::Ru * T),
                           sumStoi_),
                   veryLarge);
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::reactionTypes::DKcDT(const real T) const {
    if (fabs(sumStoi_) < veryTiny) {
        return DKpDT(T);
    } else {
        const real PatmdR0TStoi =
            pow(constant::physicalConstant::Patm / (constant::physicalConstant::Ru * T), sumStoi_);
        const real kp = Kp(T);
        const real dkcdt = (DKpDT(T) + kp * (-sumStoi_) / T) * PatmdR0TStoi;
        if (isnan(dkcdt) || isinf(dkcdt)) {
            return veryLarge;
        }
        return dkcdt;
    }
}

OpenHurricane::real OpenHurricane::reactionTypes::q(const real p, const real T, const realArray &c,
                                            real &qf, real &qr) const {
    qf = this->kf(p, T, c);
    qr = this->kr(qf, p, T, c);
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;
        qf *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;
        qr *= pow(max(c[si], real(0)), ord);
    }
    return qf - qr;
}

OpenHurricane::real OpenHurricane::reactionTypes::q(const real p, const real T, const realArray &c,
                                            real &qf, real &qr, real &kf, real &kr) const {
    qf = this->kf(p, T, c);
    kf = qf;
    qr = this->kr(qf, p, T, c);
    kr = qr;
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;
        qf *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;
        qr *= pow(max(c[si], real(0)), ord);
    }
    return qf - qr;
}

OpenHurricane::real OpenHurricane::reactionTypes::q(const real p, const real T, const realArray &c,
                                            realArray &dcidt) const {
    real qf, qr;
    const auto qfr = q(p, T, c, qf, qr);
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        dcidt[si] -= stoi * qfr;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        dcidt[si] += stoi * qfr;
    }
    return qfr;
}

OpenHurricane::real OpenHurricane::reactionTypes::q(const real p, const real T, const realArray &c,
                                            realArray &dcidt, real &qf, real &qr, real &kf,
                                            real &kr) const {
    const auto qfr = q(p, T, c, qf, qr, kf, kr);
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        dcidt[si] -= stoi * qfr;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        dcidt[si] += stoi * qfr;
    }
    return qfr;
}

OpenHurricane::real OpenHurricane::reactionTypes::q(const real p, const real T, const realArray &c,
                                            const bool reduced,
                                            const integerArray &fullToSimplified,
                                            realArray &dcidt) const {
    real qf, qr;
    const auto qfr = q(p, T, c, qf, qr);
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        dcidt[id] -= stoi * qfr;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        dcidt[id] += stoi * qfr;
    }
    return qfr;
}

OpenHurricane::real OpenHurricane::reactionTypes::q(const real p, const real T, const realArray &c,
                                            const bool reduced,
                                            const integerArray &fullToSimplified, realArray &dcidt,
                                            real &qf, real &qr, real &kf, real &kr) const {
    const auto qfr = q(p, T, c, qf, qr, kf, kr);
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        dcidt[id] -= stoi * qfr;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        dcidt[id] += stoi * qfr;
    }
    return qfr;
}

void OpenHurricane::reactionTypes::DqDci(const real kf, const real kr, const real qfr, const real p,
                                     const real T, const realArray &c,
                                     realSquareMatrix &Jac) const {
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        real dqfdci = kf;
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
            const integer sj = forwardCoeffs_[j].index_;
            const real ord = forwardCoeffs_[j].order_;
            if (i == j) {
                if (ord == real(1.0)) {
                    continue;
                } else if (ord < real(1)) {
                    if (c[sj] > tiny) {
                        dqfdci *= ord * pow(c[sj], ord - real(1));
                    } else {
                        dqfdci *= real(0);
                    }
                } else {
                    dqfdci *= ord * pow(max(c[sj], real(0)), ord - real(1));
                }
            } else {
                dqfdci *= pow(max(c[sj], real(0)), ord);
            }
        }
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
            const integer sj = forwardCoeffs_[j].index_;
            const real stoi = forwardCoeffs_[j].stoichCoeff_;
            Jac(sj, si) -= stoi * dqfdci;
        }
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
            const integer sj = backwardCoeffs_[j].index_;
            const real stoi = backwardCoeffs_[j].stoichCoeff_;
            Jac(sj, si) += stoi * dqfdci;
        }
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        real dqrdci = kr;
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
            const integer sj = backwardCoeffs_[j].index_;
            const real ord = backwardCoeffs_[j].order_;
            if (i == j) {
                if (ord == real(1.0)) {
                    continue;
                } else if (ord < real(1)) {
                    if (c[sj] > tiny) {
                        dqrdci *= ord * pow(c[sj], ord - real(1));
                    } else {
                        dqrdci *= real(0);
                    }
                } else {
                    dqrdci *= ord * pow(max(c[sj], real(0)), ord - real(1));
                }
            } else {
                dqrdci *= pow(max(c[sj], real(0)), ord);
            }
        }
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
            const integer sj = forwardCoeffs_[j].index_;
            const real stoi = forwardCoeffs_[j].stoichCoeff_;
            Jac(sj, si) += stoi * dqrdci;
        }
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
            const integer sj = backwardCoeffs_[j].index_;
            const real stoi = backwardCoeffs_[j].stoichCoeff_;
            Jac(sj, si) -= stoi * dqrdci;
        }
    }

    if (hasThirdBodyTerms()) {
        realArray dggdci(species_.size());
        this->DGGDc(p, T, c, dggdci);
        for (integer i = 0; i < dggdci.size(); ++i) {
            const integer si = i;
            for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
                const integer sj = forwardCoeffs_[j].index_;
                const real stoi = forwardCoeffs_[j].stoichCoeff_;
                if (!isnan(dggdci[i]) && !isinf(dggdci[i])) {
                    Jac(sj, si) -= stoi * dggdci[i] * qfr;
                }
            }
            for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
                const integer sj = backwardCoeffs_[j].index_;
                const real stoi = backwardCoeffs_[j].stoichCoeff_;
                if (!isnan(dggdci[i]) && !isinf(dggdci[i])) {
                    Jac(sj, si) += stoi * dggdci[i] * qfr;
                }
            }
        }
    }
}

void OpenHurricane::reactionTypes::DqDci(const real kf, const real kr, const real qfr, const real p,
                                     const real T, const realArray &c, realArray &diagJac,
                                     const bool calcLastSpc) const {
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        real dqfdci = kf;
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
            const integer sj = forwardCoeffs_[j].index_;
            const real ord = forwardCoeffs_[j].order_;
            if (i == j) {
                if (ord == real(1.0)) {
                    continue;
                } else if (ord < real(1)) {
                    if (c[sj] > tiny) {
                        dqfdci *= ord * pow(c[sj], ord - real(1));
                    } else {
                        dqfdci *= real(0);
                    }
                } else {
                    dqfdci *= ord * pow(max(c[sj], real(0)), ord - real(1));
                }
            } else {
                dqfdci *= pow(max(c[sj], real(0)), ord);
            }
        }
        if (calcLastSpc || (!calcLastSpc && si != species_.size() - 1)) {
            const real stoif = forwardCoeffs_[i].stoichCoeff_;
            diagJac[si] -= stoif * dqfdci;

            for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
                const integer sj = backwardCoeffs_[j].index_;
                if (si == sj) {
                    const real stoi = backwardCoeffs_[j].stoichCoeff_;
                    diagJac[si] += stoi * dqfdci;
                }
            }
        } else {
            for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
                const integer sj = forwardCoeffs_[j].index_;
                if (sj != si) {
                    const real stoi = forwardCoeffs_[j].stoichCoeff_;
                    diagJac[sj] += stoi * dqfdci * species().W(sj) / species().W(si);
                }
            }
            for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
                const integer sj = backwardCoeffs_[j].index_;
                if (sj != si) {
                    const real stoi = backwardCoeffs_[j].stoichCoeff_;

                    diagJac[sj] -= stoi * dqfdci * species().W(sj) / species().W(si);
                }
            }
        }
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        real dqrdci = kr;
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
            const integer sj = backwardCoeffs_[j].index_;
            const real ord = backwardCoeffs_[j].order_;
            if (i == j) {
                if (ord == real(1.0)) {
                    continue;
                } else if (ord < real(1)) {
                    if (c[sj] > tiny) {
                        dqrdci *= ord * pow(c[sj], ord - real(1));
                    } else {
                        dqrdci *= real(0);
                    }
                } else {
                    dqrdci *= ord * pow(max(c[sj], real(0)), ord - real(1));
                }
            } else {
                dqrdci *= pow(max(c[sj], real(0)), ord);
            }
        }
        if (calcLastSpc || (!calcLastSpc && si != species_.size() - 1)) {
            for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
                const integer sj = forwardCoeffs_[j].index_;
                if (si == sj) {
                    const real stoi = forwardCoeffs_[j].stoichCoeff_;
                    diagJac[si] += stoi * dqrdci;
                }
            }
            const real stoif = backwardCoeffs_[i].stoichCoeff_;
            diagJac[si] -= stoif * dqrdci;
        } else {
            for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
                const integer sj = forwardCoeffs_[j].index_;
                if (sj != si) {
                    const real stoi = forwardCoeffs_[j].stoichCoeff_;
                    diagJac[sj] -= stoi * dqrdci * species().W(sj) / species().W(si);
                }
            }
            for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
                const integer sj = backwardCoeffs_[j].index_;
                if (sj != si) {
                    const real stoi = backwardCoeffs_[j].stoichCoeff_;

                    diagJac[sj] += stoi * dqrdci * species().W(sj) / species().W(si);
                }
            }
        }
    }

    if (hasThirdBodyTerms()) {
        realArray dggdci(species_.size());
        this->DGGDc(p, T, c, dggdci);
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
            const integer sj = forwardCoeffs_[j].index_;
            const real stoi = forwardCoeffs_[j].stoichCoeff_;
            if (calcLastSpc) {
                diagJac[sj] -= stoi * dggdci[sj] * qfr;
            } else {
                if (sj != species_.size() - 1) {
                    const integer nsm1 = species_.size() - 1;
                    diagJac[sj] -= stoi *
                                   (dggdci[sj] - dggdci[nsm1] * species_.W(sj) / species_.W(nsm1)) *
                                   qfr;
                }
            }
        }
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
            const integer sj = backwardCoeffs_[j].index_;
            const real stoi = backwardCoeffs_[j].stoichCoeff_;
            if (calcLastSpc) {
                diagJac[sj] += stoi * dggdci[sj] * qfr;
            } else {
                if (sj != species_.size() - 1) {
                    const integer nsm1 = species_.size() - 1;
                    diagJac[sj] += stoi *
                                   (dggdci[sj] - dggdci[nsm1] * species_.W(sj) / species_.W(nsm1)) *
                                   qfr;
                }
            }
        }
    }
}

OpenHurricane::real OpenHurricane::reactionTypes::q(const real p, const real T, const realArray &c,
                                            realArray &dcidt, realArray &diagJac,
                                            const bool calcLastSpc, real &kf, real &kr) const {
    real pf, pr;
    real qfr = q(p, T, c, dcidt, pf, pr, kf, kr);
    DqDci(kf, kr, qfr, p, T, c, diagJac, calcLastSpc);
    return qfr;
}
void OpenHurricane::reactionTypes::DqDci(const real kf, const real kr, const real qfr, const real p,
                                     const real T, const realArray &c, const bool reduced,
                                     const integerArray &fullToSimplified,
                                     realSquareMatrix &Jac) const {
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si =
            reduced ? fullToSimplified[forwardCoeffs_[i].index_] : forwardCoeffs_[i].index_;
        real dqfdci = kf;
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
            const integer sj = forwardCoeffs_[j].index_;
            const real ord = forwardCoeffs_[j].order_;
            if (i == j) {
                dqfdci *= ord * pow(max(c[sj], real(0)), ord - real(1));

            } else {
                dqfdci *= pow(max(c[sj], real(0)), ord);
            }
        }
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {

            const integer sj =
                reduced ? fullToSimplified[forwardCoeffs_[j].index_] : forwardCoeffs_[j].index_;
            const real stoi = forwardCoeffs_[j].stoichCoeff_;
            Jac(sj, si) -= stoi * dqfdci;
        }
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
            const integer sj =
                reduced ? fullToSimplified[backwardCoeffs_[j].index_] : backwardCoeffs_[j].index_;
            const real stoi = backwardCoeffs_[j].stoichCoeff_;
            Jac(sj, si) += stoi * dqfdci;
        }
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {

        const integer si =
            reduced ? fullToSimplified[backwardCoeffs_[i].index_] : backwardCoeffs_[i].index_;
        real dqrdci = kr;
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
            const integer sj = backwardCoeffs_[j].index_;
            const real ord = backwardCoeffs_[j].order_;
            if (i == j) {

                dqrdci *= ord * pow(max(c[sj], real(0)), ord - real(1));

            } else {
                dqrdci *= pow(max(c[sj], real(0)), ord);
            }
        }
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
            const integer sj =
                reduced ? fullToSimplified[forwardCoeffs_[j].index_] : forwardCoeffs_[j].index_;
            const real stoi = forwardCoeffs_[j].stoichCoeff_;
            Jac(sj, si) += stoi * dqrdci;
        }
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {

            const integer sj =
                reduced ? fullToSimplified[backwardCoeffs_[j].index_] : backwardCoeffs_[j].index_;
            const real stoi = backwardCoeffs_[j].stoichCoeff_;
            Jac(sj, si) -= stoi * dqrdci;
        }
    }

    if (hasThirdBodyTerms()) {
        realArray dggdci(species_.size());
        this->DGGDc(p, T, c, dggdci);
        for (integer i = 0; i < dggdci.size(); ++i) {
            const integer si = reduced ? fullToSimplified[i] : i;
            if (si != -1) {
                for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
                    const integer sj = reduced ? fullToSimplified[forwardCoeffs_[j].index_]
                                               : forwardCoeffs_[j].index_;
                    const real stoi = forwardCoeffs_[j].stoichCoeff_;
                    Jac(sj, si) -= stoi * dggdci[i] * qfr;
                }
                for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
                    const integer sj = reduced ? fullToSimplified[backwardCoeffs_[j].index_]
                                               : backwardCoeffs_[j].index_;
                    const real stoi = backwardCoeffs_[j].stoichCoeff_;
                    Jac(sj, si) += stoi * dggdci[i] * qfr;
                }
            }
        }
    }
}

void OpenHurricane::reactionTypes::DqDci(const real kf, const real kr, const real qfr, const real p,
                                     const real T, const realArray &c, const bool reduced,
                                     const integerArray &fullToSimplified, realArray &diagJac,
                                     const bool calcLastSpc) const {
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        real dqfdci = kf;
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
            const integer sj = forwardCoeffs_[j].index_;
            const real ord = forwardCoeffs_[j].order_;
            if (i == j) {
                if (ord == real(1.0)) {
                    continue;
                } else if (ord < real(1)) {
                    if (c[sj] > tiny) {
                        dqfdci *= ord * pow(c[sj], ord - real(1));
                    } else {
                        dqfdci *= real(0);
                    }
                } else {
                    dqfdci *= ord * pow(max(c[sj], real(0)), ord - real(1));
                }
            } else {
                dqfdci *= pow(max(c[sj], real(0)), ord);
            }
        }
        if (calcLastSpc || (!calcLastSpc && si != species_.size() - 1)) {
            const real stoif = forwardCoeffs_[i].stoichCoeff_;
            const integer sii = reduced ? fullToSimplified[si] : si;
            diagJac[sii] -= stoif * dqfdci;

            for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
                const integer sj = backwardCoeffs_[j].index_;
                if (si == sj) {
                    const real stoi = backwardCoeffs_[j].stoichCoeff_;
                    const integer sii = reduced ? fullToSimplified[si] : si;
                    diagJac[sii] += stoi * dqfdci;
                }
            }
        } else {
            for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
                const integer sj = forwardCoeffs_[j].index_;
                if (sj != si) {
                    const real stoi = forwardCoeffs_[j].stoichCoeff_;
                    const integer sjj = reduced ? fullToSimplified[sj] : sj;
                    diagJac[sjj] += stoi * dqfdci * species().W(sj) / species().W(si);
                }
            }
            for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
                const integer sj = backwardCoeffs_[j].index_;
                if (sj != si) {
                    const real stoi = backwardCoeffs_[j].stoichCoeff_;
                    const integer sjj = reduced ? fullToSimplified[sj] : sj;

                    diagJac[sjj] -= stoi * dqfdci * species().W(sj) / species().W(si);
                }
            }
        }
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        real dqrdci = kr;
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
            const integer sj = backwardCoeffs_[j].index_;
            const real ord = backwardCoeffs_[j].order_;
            if (i == j) {
                if (ord == real(1.0)) {
                    continue;
                } else if (ord < real(1)) {
                    if (c[sj] > tiny) {
                        dqrdci *= ord * pow(c[sj], ord - real(1));
                    } else {
                        dqrdci *= real(0);
                    }
                } else {
                    dqrdci *= ord * pow(max(c[sj], real(0)), ord - real(1));
                }
            } else {
                dqrdci *= pow(max(c[sj], real(0)), ord);
            }
        }
        if (calcLastSpc || (!calcLastSpc && si != species_.size() - 1)) {
            for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
                const integer sj = forwardCoeffs_[j].index_;
                if (si == sj) {
                    const real stoi = forwardCoeffs_[j].stoichCoeff_;
                    const integer sii = reduced ? fullToSimplified[si] : si;
                    diagJac[sii] += stoi * dqrdci;
                }
            }
            const real stoif = backwardCoeffs_[i].stoichCoeff_;
            const integer sii = reduced ? fullToSimplified[si] : si;
            diagJac[sii] -= stoif * dqrdci;
        } else {
            for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
                const integer sj = forwardCoeffs_[j].index_;
                if (sj != si) {
                    const real stoi = forwardCoeffs_[j].stoichCoeff_;
                    const integer sjj = reduced ? fullToSimplified[sj] : sj;
                    diagJac[sjj] -= stoi * dqrdci * species().W(sj) / species().W(si);
                }
            }
            for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
                const integer sj = backwardCoeffs_[j].index_;
                if (sj != si) {
                    const real stoi = backwardCoeffs_[j].stoichCoeff_;
                    const integer sjj = reduced ? fullToSimplified[sj] : sj;

                    diagJac[sjj] += stoi * dqrdci * species().W(sj) / species().W(si);
                }
            }
        }
    }

    if (hasThirdBodyTerms()) {
        realArray dggdci(species_.size());
        this->DGGDc(p, T, c, dggdci);
        for (integer j = 0; j < forwardCoeffs_.size(); ++j) {
            const integer sj = forwardCoeffs_[j].index_;
            const real stoj = forwardCoeffs_[j].stoichCoeff_;
            if (calcLastSpc) {
                const integer sjj = reduced ? fullToSimplified[sj] : sj;
                diagJac[sjj] -= stoj * dggdci[sj] * qfr;
            } else {
                if (sj != species_.size() - 1) {
                    const integer sjj = reduced ? fullToSimplified[sj] : sj;
                    const integer nsm1 = species_.size() - 1;
                    diagJac[sjj] -=
                        stoj * (dggdci[sj] - dggdci[nsm1] * species_.W(sj) / species_.W(nsm1)) *
                        qfr;
                }
            }
        }
        for (integer j = 0; j < backwardCoeffs_.size(); ++j) {
            const integer sj = backwardCoeffs_[j].index_;
            const real stoj = backwardCoeffs_[j].stoichCoeff_;
            if (calcLastSpc) {
                const integer sjj = reduced ? fullToSimplified[sj] : sj;
                diagJac[sjj] += stoj * dggdci[sj] * qfr;
            } else {
                if (sj != species_.size() - 1) {
                    const integer sjj = reduced ? fullToSimplified[sj] : sj;
                    const integer nsm1 = species_.size() - 1;
                    diagJac[sjj] +=
                        stoj * (dggdci[sj] - dggdci[nsm1] * species_.W(sj) / species_.W(nsm1)) *
                        qfr;
                }
            }
        }
    }
}

OpenHurricane::real OpenHurricane::reactionTypes::q(const real p, const real T, const realArray &c,
                                            const bool reduced,
                                            const integerArray &fullToSimplified, realArray &dcidt,
                                            realArray &diagJac, const bool calcLastSpc, real &kf,
                                            real &kr) const {
    real pf, pr;
    real qfr = q(p, T, c, reduced, fullToSimplified, dcidt, pf, pr, kf, kr);
    DqDci(kf, kr, qfr, p, T, c, reduced, fullToSimplified, diagJac, calcLastSpc);
    return qfr;
}

void OpenHurricane::reactionTypes::DqDT(const real kf, const real kr, const real qfr, const real p,
                                    const real T, const realArray &c, realSquareMatrix &Jac) const {
    real dkfdT = this->DkfDT(kf, p, T, c);
    real dkrdT = this->DkrDT(kr, p, T, c, dkfdT);

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;
        dkfdT *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;
        dkrdT *= pow(max(c[si], real(0)), ord);
    }

    real dqjdT = dkfdT - dkrdT;
    real dgdT = qfr * this->DGGDT(p, T, c);
    if (isnan(dqjdT) || isinf(dqjdT)) {
        dqjdT = 0;
    }
    if (isnan(dgdT) || isinf(dgdT)) {
        dgdT = 0;
    }

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        Jac(si, species_.size()) -= stoi * (dqjdT + dgdT);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        Jac(si, species_.size()) += stoi * (dqjdT + dgdT);
    }
}

void OpenHurricane::reactionTypes::DqDT(const real kf, const real kr, const real qfr, const real p,
                                    const real T, const realArray &c, realArray &dwdT) const {
    real dkfdT = this->DkfDT(kf, p, T, c);
    real dkrdT = this->DkrDT(kr, p, T, c, dkfdT);

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;
        dkfdT *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;
        dkrdT *= pow(max(c[si], real(0)), ord);
    }

    const real dqjdT = dkfdT - dkrdT;
    const real dgdT = qfr * this->DGGDT(p, T, c);

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        dwdT[si] -= stoi * (dqjdT + dgdT);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        dwdT[si] += stoi * (dqjdT + dgdT);
    }
}

void OpenHurricane::reactionTypes::DqDT(const real kf, const real kr, const real qfr, const real p,
                                    const real T, const realArray &c, const bool reduced,
                                    const integerArray &fullToSimplified, const integer Tid,
                                    realSquareMatrix &Jac) const {
    real dkfdT = this->DkfDT(kf, p, T, c);
    real dkrdT = this->DkrDT(kr, p, T, c, dkfdT);
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;

        dkfdT *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;
        dkrdT *= pow(max(c[si], real(0)), ord);
    }

    const real dqjdT = dkfdT - dkrdT;
    const real dgdT = qfr * this->DGGDT(p, T, c);

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        Jac(id, Tid) -= stoi * (dqjdT + dgdT);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        Jac(id, Tid) += stoi * (dqjdT + dgdT);
    }
}

void OpenHurricane::reactionTypes::DqDT(const real kf, const real kr, const real qfr, const real p,
                                    const real T, const realArray &c, const bool reduced,
                                    const integerArray &fullToSimplified, realArray &dwdT) const {
    real dkfdT = this->DkfDT(kf, p, T, c);
    real dkrdT = this->DkrDT(kr, p, T, c, dkfdT);

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;
        dkfdT *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;

        dkrdT *= pow(max(c[si], real(0)), ord);
    }

    const real dqjdT = dkfdT - dkrdT;
    const real dgdT = qfr * this->DGGDT(p, T, c);

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        dwdT[id] -= stoi * (dqjdT + dgdT);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        dwdT[id] += stoi * (dqjdT + dgdT);
    }
}

void OpenHurricane::reactionTypes::DqDp(const real kf, const real kr, const real qfr, const real p,
                                    const real T, const realArray &c, realSquareMatrix &Jac) const {
    real dkfdp = this->DkfDp(kf, p, T, c);
    real dkrdp = this->DkrDp(kr, p, T, c, dkfdp);

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;
        dkfdp *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;
        dkrdp *= pow(max(c[si], real(0)), ord);
    }

    real dqjdT = dkfdp - dkrdp;

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        Jac(si, species_.size() + 1) -= stoi * dqjdT;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        Jac(si, species_.size() + 1) += stoi * dqjdT;
    }
}

void OpenHurricane::reactionTypes::DqDp(const real kf, const real kr, const real qfr, const real p,
                                    const real T, const realArray &c, realArray &dwdp) const {
    if (!isPressureDenpendent()) {
        return;
    }
    real dkfdp = this->DkfDp(kf, p, T, c);
    real dkrdp = this->DkrDp(kr, p, T, c, dkfdp);

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;
        dkfdp *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;
        dkrdp *= pow(max(c[si], real(0)), ord);
    }

    const real dqjdp = dkfdp - dkrdp;

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        dwdp[si] -= stoi * dqjdp;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        dwdp[si] += stoi * dqjdp;
    }
}

void OpenHurricane::reactionTypes::DqDp(const real kf, const real kr, const real qfr, const real p,
                                    const real T, const realArray &c, const bool reduced,
                                    const integerArray &fullToSimplified, const integer pid,
                                    realSquareMatrix &Jac) const {
    if (!isPressureDenpendent()) {
        return;
    }
    real dkfdp = this->DkfDp(kf, p, T, c);
    real dkrdp = this->DkrDp(kr, p, T, c, dkfdp);
    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;

        dkfdp *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;

        dkrdp *= pow(max(c[si], real(0)), ord);
    }

    const real dqjdp = dkfdp - dkrdp;

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        Jac(id, pid) -= stoi * dqjdp;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        Jac(id, pid) += stoi * dqjdp;
    }
}

void OpenHurricane::reactionTypes::DqDp(const real kf, const real kr, const real qfr, const real p,
                                    const real T, const realArray &c, const bool reduced,
                                    const integerArray &fullToSimplified, realArray &dwdp) const {
    if (!isPressureDenpendent()) {
        return;
    }
    real dkfdp = this->DkfDp(kf, p, T, c);
    real dkrdp = this->DkrDp(kr, p, T, c, dkfdp);

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real ord = forwardCoeffs_[i].order_;

        dkfdp *= pow(max(c[si], real(0)), ord);
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real ord = backwardCoeffs_[i].order_;

        dkrdp *= pow(max(c[si], real(0)), ord);
    }

    const real dqjdp = dkfdp - dkrdp;

    for (integer i = 0; i < forwardCoeffs_.size(); ++i) {
        const integer si = forwardCoeffs_[i].index_;
        const real stoi = forwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        dwdp[id] -= stoi * dqjdp;
    }
    for (integer i = 0; i < backwardCoeffs_.size(); ++i) {
        const integer si = backwardCoeffs_[i].index_;
        const real stoi = backwardCoeffs_[i].stoichCoeff_;
        const integer id = reduced ? fullToSimplified[si] : si;
        dwdp[id] += stoi * dqjdp;
    }
}

#ifdef CUDA_PARALLEL
void OpenHurricane::reactionTypes::cuSetReactionCoeff(
    const integer reacId, const integer nrc,
    cuChem::cuChemInterface::reactionCoeffsInterface &reacInt,
    const List<reactionSpeciesCoeffs> &coef) const {
    reacInt.sizePtr_[reacId] = coef.size();
    if (coef.size() > reacInt.maxSpcSizeInRec_) {
        errorAbortStr(
            ("The size: " + toString(coef.size()) +
             " of reaction coefficients of reaction: " + toString(reacId) +
             " is larger than maximum size: " + toString(integer(reacInt.maxSpcSizeInRec_))));
    }

    for (integer i = 0; i < coef.size(); ++i) {
        reacInt.reactionCoeffsIndexPtr_[i * nrc + reacId] = coef[i].index_;
        reacInt.stoichCoeffPtr_[i * nrc + reacId] = coef[i].stoichCoeff_;
        reacInt.orderPtr_[i * nrc + reacId] = coef[i].order_;
    }
}

void OpenHurricane::reactionTypes::cuSetReactionType(const integer reacId, const integer nrc,
                                                 cuChem::cuChemInterface &reacInt) const {
    cuSetReactionCoeff(reacId, nrc, reacInt.reacInt_, forwardCoeffs_);
    cuSetReactionCoeff(reacId, nrc, reacInt.revReacInt_, backwardCoeffs_);
}

#endif // CUDA_PARALLEL