#include "ODESolver.hpp"
/*!
 * \file ODESolver.inl
 * \brief In-Line subroutines of the <i>ODESolver.hpp</i> file.
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
#pragma once

#ifdef TEST_PROCESS_TIME
template <class ChemistrySource>
inline void OpenHurricane::ODESolver<ChemistrySource>::setODECountIter() noexcept {
    if (odesPtr_) {
        this->odeCountIter_ += odesPtr_->countIter();
    }
}
#endif

template <class ChemistrySource>
inline OpenHurricane::ODESolver<ChemistrySource>::ODESolver(flowModel &flows,
                                                            const controller &cont)
    : chemistryODESolvers<ChemistrySource>(flows, cont), odesPtr_(nullptr),
      cT_(this->nEqns(), Zero) {
    if (cont.found("solver")) {
        const auto &solverCont = cont.subController("solver");
        odesPtr_ = ODEsSolver::creator(this->nEqns(), solverCont);
    } else {
        odesPtr_.reset(new Euler(this->nEqns(), cont));
    }

    ODEsSolver::dydtFunc_type func = std::bind(&ChemistrySource::DyDt, this, std::placeholders::_1,
                                               std::placeholders::_2, std::placeholders::_3);
    ODEsSolver::JacobianFunc_type JacFunc =
        std::bind(&ChemistrySource::jacobian, this, std::placeholders::_1, std::placeholders::_2,
                  std::placeholders::_3, std::placeholders::_4);

    ODEsSolver::checkFunc_type checkFunc =
        std::bind(&ChemistrySource::checkNewY, this, std::placeholders::_1, std::placeholders::_2,
                  std::placeholders::_3);

    odesPtr_->setDyDtFunc(func);
    odesPtr_->setJacobianFunc(JacFunc);
    odesPtr_->setCheckFunc(checkFunc);
}

template <class ChemistrySource>
inline OpenHurricane::ODESolver<ChemistrySource>::~ODESolver() noexcept {
    odesPtr_.clear();
}

template <class ChemistrySource>
inline void OpenHurricane::ODESolver<ChemistrySource>::solve(realArray &c, real &dT, real &subDT) {
    if (odesPtr_->reset(this->nEqns())) {
        odesPtr_->resetArray(cT_);
    }

    for (integer i = 0; i < this->nSpc(); ++i) {
        cT_[i] = c[i];
    }
    cT_[this->nSpc()] = this->T_;

    odesPtr_->solve(0, dT, subDT, cT_);
    this->T_ = cT_[this->nSpc()];

    for (integer i = 0; i < this->nSpc(); ++i) {
        c[i] = cT_[i];
    }
}
