#include "turbulentSpeciesSolver.hpp"
/*!
 * \file turbulentSpeciesSolver.inl
 * \brief The In-Line functions of the <i>turbulentSpeciesSolver.hpp</i> file.
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

inline OpenHurricane::turbulentSpeciesSolver::~turbulentSpeciesSolver() noexcept {
    chemtryPtr_.clear();
#ifdef TEST_PROCESS_TIME
    fosTestChemTime_.close();
#endif // TEST_PROCESS_TIME
#ifdef CUDA_PARALLEL
    fzBoundStartPtr_.clear();
#endif //CUDA_PARALLEL
}

inline void OpenHurricane::turbulentSpeciesSolver::previousSource() {
    if (!mixtures().noReaction()) {
        hrClock myClocks;
        if (chemtryPtr_->isStrangSplitted()) {
            chemtryPtr_->strangSplittedChemistrySource(*timeMarcingPtr_);
            mixtures().E(p(), T(), v(), E(), true, true);
            updatePrimitives();
        } else if (chemtryPtr_->isIntegrated()) {
            chemtryPtr_->integratedChemistrySource();
        }
        if (isUpdatingCalcTime()) {
            HurMPI::barrier();
            sorTime_ = myClocks.elapsedClockTime();
        }
    }
}

inline void OpenHurricane::turbulentSpeciesSolver::postSource() {
    if (!mixtures().noReaction()) {
        hrClock myClocks;
        if (chemtryPtr_->isStrangSplitted()) {
            chemtryPtr_->strangSplittedChemistrySource(*timeMarcingPtr_);
            mixtures().E(p(), T(), v(), E(), true, true);
            updatePrimitives();
        }
        if (isUpdatingCalcTime()) {
            HurMPI::barrier();
            sorTime_ += myClocks.elapsedClockTime();
        }
    }
}

#ifdef CUDA_PARALLEL
hur_nodiscard inline const OpenHurricane::integerList &
OpenHurricane::turbulentSpeciesSolver::fzBoundStart() const {
    if (!fzBoundStartPtr_) {
        makeFzBoundStart();
    }
    return *fzBoundStartPtr_;
}
#endif // CUDA_PARALLEL