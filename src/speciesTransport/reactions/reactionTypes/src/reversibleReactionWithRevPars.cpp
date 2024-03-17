/*!
 * \file reversibleReactionWithRevPars.cpp
 * \brief The subroutines and functions of non-equilibrium reversible reaction.
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
#include "reversibleReactionWithRevPars.hpp"
#include "reactionList.hpp"

#ifdef CUDA_PARALLEL
#include "ArrheniusRR.hpp"
#include "ArrheniusThirdBodyRR.hpp"
#include "chemicallyActicatedBimolecularRR.hpp"
#include "fallOffFunctions.hpp"
#include "unimolecularFallOffRR.hpp"
#endif // CUDA_PARALLEL

namespace OpenHurricane {
    createClassName(reversibleReactionWithRevPars);
         registerObjFty(reactionTypes, reversibleReactionWithRevPars, controller);
} // namespace OpenHurricane

OpenHurricane::reversibleReactionWithRevPars::reversibleReactionWithRevPars(
    const reactionTypes &reac, const reactionRateTypes &rrf, const reactionRateTypes &rrr)
    : reactionTypes(reac), kfPtr_(rrf.clone()), krPtr_(rrr.clone()) {
    reactionTypes::reacType_ = reactionTypes::types::nonEqReversible;
}

OpenHurricane::reversibleReactionWithRevPars::reversibleReactionWithRevPars(
    const speciesList &sp, const reactionList &rt, const controller &cont)
    : reactionTypes(sp, rt, cont),
      kfPtr_(reactionRateTypes::creator(sp, cont.subController("forwardRate"))),
      krPtr_(reactionRateTypes::creator(sp, cont.subController("backwardRate"))) {
    reactionTypes::reacType_ = reactionTypes::types::nonEqReversible;
}

OpenHurricane::reversibleReactionWithRevPars::reversibleReactionWithRevPars(
    const reversibleReactionWithRevPars &other)
    : reactionTypes(other), kfPtr_(other.kfPtr_->clone()), krPtr_(other.krPtr_->clone()) {
    reactionTypes::reacType_ = reactionTypes::types::nonEqReversible;
}

#ifdef CUDA_PARALLEL
void OpenHurricane::reversibleReactionWithRevPars::cuSetReactionType(
    const integer reacId, const integer nrc, cuChem::cuChemInterface &reacInt) const {
    if (!kfPtr_) {
        LFatal("The pointer of forward reaction rate is null");
    }
    if (!krPtr_) {
        LFatal("The pointer of backward reaction rate is null");
    }
    reactionTypes::cuSetReactionType(reacId, nrc, reacInt);
    if (kfPtr_->isArrhenius()) {
        const auto &ArrK = static_cast<const ArrheniusRR &>(*kfPtr_);
        reacInt.kfInt_.aPtr_[reacId] = ArrK.A();
        reacInt.kfInt_.bPtr_[reacId] = ArrK.beta();
        reacInt.kfInt_.TaPtr_[reacId] = ArrK.Ta();
        reacInt.thirdBodyTypeIntPtr_[reacId] = cuChem::cuChemInterface::noThirdBody;
        const auto &ArrKr = static_cast<const ArrheniusRR &>(*krPtr_);
        reacInt.rfInt_.aPtr_[reacId] = ArrKr.A();
        reacInt.rfInt_.bPtr_[reacId] = ArrKr.beta();
        reacInt.rfInt_.TaPtr_[reacId] = ArrKr.Ta();
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

        const auto &ArrKr = static_cast<const ArrheniusRR &>(*krPtr_);
        reacInt.rfInt_.aPtr_[reacId] = ArrKr.A();
        reacInt.rfInt_.bPtr_[reacId] = ArrKr.beta();
        reacInt.rfInt_.TaPtr_[reacId] = ArrKr.Ta();
    } else if (kfPtr_->isUnimolecularFallOff()) {
        LFatal("The unimolecular fall-off reaction is not consistent with reversible reaction "
                   "with rev parameter.");
    } else if (kfPtr_->isChemicallyActicatedBimolecular()) {
        LFatal("The chemically activated bimolecular reaction is not consistent with "
                   "reversible reaction with rev parameter.");
    } else {
        LFatal("Unsupported reaction rate type");
    }
    reacInt.reactionTypeIntPtr_[reacId] = cuChem::cuChemInterface::reversibleReactionWithRevPars;
}
#endif // CUDA_PARALLEL