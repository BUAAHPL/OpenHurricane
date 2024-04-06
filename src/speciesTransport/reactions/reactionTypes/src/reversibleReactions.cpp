/*!
 * \file reversibleReactions.cpp
 * \brief The subroutines and functions of reversible reaction.
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
#include "reversibleReactions.hpp"
#include "reactionList.hpp"

#ifdef CUDA_PARALLEL
#include "ArrheniusRR.hpp"
#include "ArrheniusThirdBodyRR.hpp"
#include "chemicallyActicatedBimolecularRR.hpp"
#include "fallOffFunctions.hpp"
#include "unimolecularFallOffRR.hpp"
#endif // CUDA_PARALLEL

namespace OpenHurricane {
    createClassName(reversibleReactions);
         registerObjFty(reactionTypes, reversibleReactions, controller);
} // namespace OpenHurricane

OpenHurricane::reversibleReactions::reversibleReactions(const reactionTypes &reac,
                                                    const reactionRateTypes &rrf)
    : reactionTypes(reac), kfPtr_(rrf.clone()) {
    reactionTypes::reacType_ = reactionTypes::types::reversible;
}

OpenHurricane::reversibleReactions::reversibleReactions(const speciesList &sp, const reactionList &rt,
                                                    const controller &cont)
    : reactionTypes(sp, rt, cont),
      kfPtr_(reactionRateTypes::creator(sp, cont.subController("forwardRate"))) {
    reactionTypes::reacType_ = reactionTypes::types::reversible;
}

OpenHurricane::reversibleReactions::reversibleReactions(const reversibleReactions &other)
    : reactionTypes(other), kfPtr_(other.kfPtr_->clone()) {
    reactionTypes::reacType_ = reactionTypes::types::reversible;
}

#ifdef CUDA_PARALLEL
void OpenHurricane::reversibleReactions::cuSetReactionType(const integer reacId, const integer nrc,
                                                       cuChem::cuChemInterface &reacInt) const {
    if (!kfPtr_) {
        LFatal("The pointer of reaction rate is null");
    }
    reactionTypes::cuSetReactionType(reacId, nrc, reacInt);
    if (kfPtr_->isArrhenius()) {
        const auto &ArrK = static_cast<const ArrheniusRR &>(*kfPtr_);
        reacInt.kfInt_.aPtr_[reacId] = ArrK.A();
        reacInt.kfInt_.bPtr_[reacId] = ArrK.beta();
        reacInt.kfInt_.TaPtr_[reacId] = ArrK.Ta();
        reacInt.thirdBodyTypeIntPtr_[reacId] = cuChem::cuChemInterface::noThirdBody;
    } else if (kfPtr_->isArrheniusThirdBody()) {
        const auto &ArrK = static_cast<const ArrheniusThirdBodyRR &>(*kfPtr_);
        reacInt.kfInt_.aPtr_[reacId] = ArrK.A();
        reacInt.kfInt_.bPtr_[reacId] = ArrK.beta();
        reacInt.kfInt_.TaPtr_[reacId] = ArrK.Ta();
        reacInt.thirdBodyTypeIntPtr_[reacId] = cuChem::cuChemInterface::thirdBody;

        reacInt.thirdBInt_.thidrBodyIndexPtr_[reacId] = reacId;

        for (integer i = 0; i < species_.size(); ++i) {
            reacInt.thirdBInt_.coefThirdPtr_[i * nrc + reacId] = ArrK.thirdBodyEff()[i];
        }
    } else if (kfPtr_->isUnimolecularFallOff()) {
        const auto &ArrK = static_cast<const unimolecularFallOffRR &>(*kfPtr_);

        reacInt.thirdBodyTypeIntPtr_[reacId] = cuChem::cuChemInterface::UnimolecularFallOff;
        reacInt.thirdBInt_.thidrBodyIndexPtr_[reacId] = reacId;
        for (integer i = 0; i < species_.size(); ++i) {
            reacInt.thirdBInt_.coefThirdPtr_[i * nrc + reacId] = ArrK.thirdBodyEff()[i];
        }

        // High-pressure limit Arrhenius parameters
        reacInt.kfInt_.aPtr_[reacId] = ArrK.kinf().A();
        reacInt.kfInt_.bPtr_[reacId] = ArrK.kinf().beta();
        reacInt.kfInt_.TaPtr_[reacId] = ArrK.kinf().Ta();
        ;

        // Low-pressure limit Arrhenius parameters
        reacInt.pressDepInt_.indexPtr_[reacId] = reacId;
        reacInt.pressDepInt_.aPtr_[reacId] = ArrK.k0().A();
        reacInt.pressDepInt_.bPtr_[reacId] = ArrK.k0().beta();
        reacInt.pressDepInt_.TaPtr_[reacId] = ArrK.k0().Ta();
        if (ArrK.fallOff().isLindemann()) {
            reacInt.pressDepInt_.fallOffTypePtr_[reacId] =
                cuChem::cuChemInterface::pressureDepCoeffInterface::Lindmann;
        } else if (ArrK.fallOff().isTroe()) {
            const auto &troeCoeff = static_cast<const fallOffFunctions::Troe &>(ArrK.fallOff());
            if (troeCoeff.isTssOmitted()) {
                reacInt.pressDepInt_.fallOffCoeffPtr_[reacId] = troeCoeff.alpha();
                reacInt.pressDepInt_.fallOffCoeffPtr_[nrc + reacId] = troeCoeff.Tsss();
                reacInt.pressDepInt_.fallOffCoeffPtr_[2 * nrc + reacId] = troeCoeff.Ts();
                reacInt.pressDepInt_.fallOffTypePtr_[reacId] =
                    cuChem::cuChemInterface::pressureDepCoeffInterface::TroeWithoutLastTerm;
            } else {
                reacInt.pressDepInt_.fallOffCoeffPtr_[reacId] = troeCoeff.alpha();
                reacInt.pressDepInt_.fallOffCoeffPtr_[nrc + reacId] = troeCoeff.Tsss();
                reacInt.pressDepInt_.fallOffCoeffPtr_[2 * nrc + reacId] = troeCoeff.Ts();
                reacInt.pressDepInt_.fallOffCoeffPtr_[3 * nrc + reacId] = troeCoeff.Tss();
                reacInt.pressDepInt_.fallOffTypePtr_[reacId] =
                    cuChem::cuChemInterface::pressureDepCoeffInterface::TroeWithLastTerm;
            }
        } else if (ArrK.fallOff().isSRI()) {
            reacInt.pressDepInt_.fallOffTypePtr_[reacId] =
                cuChem::cuChemInterface::pressureDepCoeffInterface::SRI;
            const auto &SRICoeff = static_cast<const fallOffFunctions::SRI &>(ArrK.fallOff());
            reacInt.pressDepInt_.fallOffCoeffPtr_[reacId] = SRICoeff.a();
            reacInt.pressDepInt_.fallOffCoeffPtr_[nrc + reacId] = SRICoeff.b();
            reacInt.pressDepInt_.fallOffCoeffPtr_[2 * nrc + reacId] = SRICoeff.c();
            reacInt.pressDepInt_.fallOffCoeffPtr_[3 * nrc + reacId] = SRICoeff.d();
            reacInt.pressDepInt_.fallOffCoeffPtr_[4 * nrc + reacId] = SRICoeff.e();
        }
    } else if (kfPtr_->isChemicallyActicatedBimolecular()) {
        const auto &ArrK = static_cast<const chemicallyActicatedBimolecularRR &>(*kfPtr_);
        reacInt.thirdBodyTypeIntPtr_[reacId] = cuChem::cuChemInterface::bimolecularFallOff;
        reacInt.thirdBInt_.thidrBodyIndexPtr_[reacId] = reacId;
        for (integer i = 0; i < species_.size(); ++i) {
            reacInt.thirdBInt_.coefThirdPtr_[i * nrc + reacId] = ArrK.thirdBodyEff()[i];
        }

        // Low-pressure limit Arrhenius parameters
        reacInt.kfInt_.aPtr_[reacId] = ArrK.k0().A();
        reacInt.kfInt_.bPtr_[reacId] = ArrK.k0().beta();
        reacInt.kfInt_.TaPtr_[reacId] = ArrK.k0().Ta();

        // High-pressure limit Arrhenius parameters
        reacInt.pressDepInt_.indexPtr_[reacId] = reacId;
        reacInt.pressDepInt_.aPtr_[reacId] = ArrK.kinf().A();
        reacInt.pressDepInt_.bPtr_[reacId] = ArrK.kinf().beta();
        reacInt.pressDepInt_.TaPtr_[reacId] = ArrK.kinf().Ta();

        if (ArrK.fallOff().isLindemann()) {
            reacInt.pressDepInt_.fallOffTypePtr_[reacId] =
                cuChem::cuChemInterface::pressureDepCoeffInterface::Lindmann;
        } else if (ArrK.fallOff().isTroe()) {
            const auto &troeCoeff = static_cast<const fallOffFunctions::Troe &>(ArrK.fallOff());
            if (troeCoeff.isTssOmitted()) {
                reacInt.pressDepInt_.fallOffCoeffPtr_[reacId] = troeCoeff.alpha();
                reacInt.pressDepInt_.fallOffCoeffPtr_[nrc + reacId] = troeCoeff.Tsss();
                reacInt.pressDepInt_.fallOffCoeffPtr_[2 * nrc + reacId] = troeCoeff.Ts();
                reacInt.pressDepInt_.fallOffTypePtr_[reacId] =
                    cuChem::cuChemInterface::pressureDepCoeffInterface::TroeWithoutLastTerm;
            } else {
                reacInt.pressDepInt_.fallOffCoeffPtr_[reacId] = troeCoeff.alpha();
                reacInt.pressDepInt_.fallOffCoeffPtr_[nrc + reacId] = troeCoeff.Tsss();
                reacInt.pressDepInt_.fallOffCoeffPtr_[2 * nrc + reacId] = troeCoeff.Ts();
                reacInt.pressDepInt_.fallOffCoeffPtr_[3 * nrc + reacId] = troeCoeff.Tss();
                reacInt.pressDepInt_.fallOffTypePtr_[reacId] =
                    cuChem::cuChemInterface::pressureDepCoeffInterface::TroeWithLastTerm;
            }
        } else if (ArrK.fallOff().isSRI()) {
            reacInt.pressDepInt_.fallOffTypePtr_[reacId] =
                cuChem::cuChemInterface::pressureDepCoeffInterface::SRI;
            const auto &SRICoeff = static_cast<const fallOffFunctions::SRI &>(ArrK.fallOff());
            reacInt.pressDepInt_.fallOffCoeffPtr_[reacId] = SRICoeff.a();
            reacInt.pressDepInt_.fallOffCoeffPtr_[nrc + reacId] = SRICoeff.b();
            reacInt.pressDepInt_.fallOffCoeffPtr_[2 * nrc + reacId] = SRICoeff.c();
            reacInt.pressDepInt_.fallOffCoeffPtr_[3 * nrc + reacId] = SRICoeff.d();
            reacInt.pressDepInt_.fallOffCoeffPtr_[4 * nrc + reacId] = SRICoeff.e();
        }
    } else {
        LFatal("Unsupported reaction rate type");
    }
    reacInt.reactionTypeIntPtr_[reacId] = cuChem::cuChemInterface::reversibleReaction;
}
#endif // CUDA_PARALLEL