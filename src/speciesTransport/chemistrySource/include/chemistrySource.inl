#include "chemistrySource.hpp"
/*!
 * \file chemistrySource.inl
 * \brief The In-Line functions of the <i>chemistrySource.hpp</i> file.
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

//hur_nodiscard inline bool OpenHurricane::chemistrySource::hasODEFactor() const noexcept
//{
//	return odeFactorPtr_ != nullptr;
//}
//
//hur_nodiscard inline OpenHurricane::real OpenHurricane::chemistrySource::odeFactor() const
//{
//	if (!hasODEFactor()) {
//		errorAbortStr("The odeFactor is not set", HUR_FUNCTION);
//	}
//	return *odeFactorPtr_;
//}
//
//hur_nodiscard inline void OpenHurricane::chemistrySource::setODEFactor(const real odef)
//{
//	HurDelete(odeFactorPtr_);
//
//	odeFactorPtr_ = new real(odef);
//}
//
//hur_nodiscard inline void OpenHurricane::chemistrySource::unsetODEFactor()
//{
//	HurDelete(odeFactorPtr_);
//}
#ifdef TEST_PROCESS_TIME
inline void OpenHurricane::chemistrySource::setODECountIter() noexcept {
    odeCountIter_ = 0;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::chemistrySource::odeCountIter() const noexcept {
    return odeCountIter_;
}
#endif

inline OpenHurricane::realArrayArray &OpenHurricane::chemistrySource::Ri() {
    if (RiPtr_ == nullptr) {
        RiPtr_ = new realArrayArray(nsp_);
        for (integer isp = 0; isp < RiPtr_->size(); ++isp) {
            (*RiPtr_)[isp].resize(mesh_.nCells(), Zero);
        }
    }
    return *RiPtr_;
}

inline const OpenHurricane::realArrayArray &OpenHurricane::chemistrySource::Ri() const {
    if (RiPtr_ == nullptr) {
        RiPtr_ = new realArrayArray(nsp_);
        for (integer isp = 0; isp < RiPtr_->size(); ++isp) {
            (*RiPtr_)[isp].resize(mesh_.nCells(), Zero);
        }
    }
    return *RiPtr_;
}

inline OpenHurricane::realArray &OpenHurricane::chemistrySource::dtInit(const realArray &dt,
                                                                const real dtFactor) {
    if (dtInitPtr_ == nullptr) {
        dtInitPtr_ = new realArray(mesh_.nCells());
        for (integer i = 0; i < mesh_.nCells(); ++i) {
            (*dtInitPtr_)[i] = dt[i] * dtFactor;
        }
    }
    return *dtInitPtr_;
}

#ifdef TEST_PROCESS_TIME
hur_nodiscard inline OpenHurricane::real
OpenHurricane::chemistrySource::totalSolveIngTime() const noexcept {
    return calcTime_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::chemistrySource::reducionTime() const noexcept {
    return reducionTime_;
}

inline void OpenHurricane::chemistrySource::writeTimeTitle(fileOsstream &fos) const {
    if (HurMPI::master()) {
        fos.os() << "variables = "
                    "\"step\",\t\"aveSolveODETime[s]\",\t\"minSolveODETime[s]"
                    "\",\t\"maxSolveODETime[s]\""
                 << std::endl;
    }
}

inline void OpenHurricane::chemistrySource::writeSolveODETime(fileOsstream &fos) const {
    real maxCalcTime = calcTime_;
    real minCalcTime = calcTime_;
    real sumCalcTime = calcTime_;
    HurMPI::reduce(maxCalcTime, MPI_MAX);
    HurMPI::reduce(minCalcTime, MPI_MIN);
    HurMPI::reduce(sumCalcTime, MPI_SUM);
    if (HurMPI::master()) {
        fos.setRealPrecision();
        fos.os() << flows_.mesh().Iteration().cStep() << "\t" << sumCalcTime / HurMPI::getProcSize()
                 << "\t" << minCalcTime << "\t" << maxCalcTime << std::endl;
        fos.unsetRealPrecision();
    }
}

#endif // TEST_PROCESS_TIME

inline const OpenHurricane::runtimeMesh &OpenHurricane::chemistrySource::mesh() const noexcept {
    return mesh_;
}

inline OpenHurricane::PtrList<OpenHurricane::cellRealArray> &OpenHurricane::chemistrySource::yi() noexcept {
    return yi_;
}

inline const OpenHurricane::PtrList<OpenHurricane::cellRealArray> &
OpenHurricane::chemistrySource::yi() const noexcept {
    return yi_;
}

inline OpenHurricane::reactionList &OpenHurricane::chemistrySource::reactions() noexcept {
    return reactions_;
}

inline const OpenHurricane::reactionList &OpenHurricane::chemistrySource::reactions() const noexcept {
    return reactions_;
}

inline const OpenHurricane::speciesList &OpenHurricane::chemistrySource::species() const noexcept {
    return species_;
}

inline OpenHurricane::integer OpenHurricane::chemistrySource::nsp() const noexcept {
    return nsp_;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::chemistrySource::nSpc() const noexcept {
    return nSpc_;
}

inline OpenHurricane::integer OpenHurricane::chemistrySource::nrc() const noexcept {
    return nrc_;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::chemistrySource::nRact() const noexcept {
    return nRact_;
}

inline bool OpenHurricane::chemistrySource::isReacting() const noexcept {
    return isReacting_;
}

inline OpenHurricane::integer OpenHurricane::chemistrySource::nEqns() const {
    return nsp() + 1;
}

inline void OpenHurricane::chemistrySource::setReactionFlowType(const ReactionFlowType t) noexcept {
    reactionFlowType_ = t;
}

inline void OpenHurricane::chemistrySource::setConstantPressure() noexcept {
    reactionFlowType_ = ReactionFlowType::ConstantPressure;
}

inline void OpenHurricane::chemistrySource::setConstantVolume() noexcept {
    reactionFlowType_ = ReactionFlowType::ConstantVolume;
}

inline bool OpenHurricane::chemistrySource::isConstPressure() const noexcept {
    return this->reactionFlowType_ == ReactionFlowType::ConstantPressure;
}

inline bool OpenHurricane::chemistrySource::isConstVolume() const noexcept {
    return this->reactionFlowType_ == ReactionFlowType::ConstantVolume;
}

inline void OpenHurricane::chemistrySource::init(const real hOre, const real T0, const real p0,
                                             const real rho0) noexcept {
    hea0_ = hOre;
    T_ = T0;
    p0_ = p0;
    p_ = p0;

    rho0_ = rho0;
    rho_ = rho0;
}

inline OpenHurricane::real OpenHurricane::chemistrySource::hea0() const noexcept {
    return hea0_;
}

inline OpenHurricane::real &OpenHurricane::chemistrySource::hea0() noexcept {
    return hea0_;
}

inline OpenHurricane::real OpenHurricane::chemistrySource::p0() const noexcept {
    return p0_;
}

inline OpenHurricane::real &OpenHurricane::chemistrySource::p0() noexcept {
    return p0_;
}

inline OpenHurricane::real OpenHurricane::chemistrySource::T() const noexcept {
    return T_;
}

inline OpenHurricane::real &OpenHurricane::chemistrySource::T() noexcept {
    return T_;
}

inline OpenHurricane::real OpenHurricane::chemistrySource::p() const noexcept {
    return p_;
}

inline OpenHurricane::real &OpenHurricane::chemistrySource::p() noexcept {
    return p_;
}

inline OpenHurricane::real OpenHurricane::chemistrySource::rho0() const noexcept {
    return rho0_;
}

inline OpenHurricane::real &OpenHurricane::chemistrySource::rho0() noexcept {
    return rho0_;
}

inline OpenHurricane::real OpenHurricane::chemistrySource::rho() const noexcept {
    return rho_;
}

inline OpenHurricane::real &OpenHurricane::chemistrySource::rho() noexcept {
    return rho_;
}

//inline OpenHurricane::integer OpenHurricane::chemistrySource::speciesSize() const noexcept
//{
//	return nsp_;
//}

inline const OpenHurricane::mixture &OpenHurricane::chemistrySource::mixtures() const {
    return flows_.mixtures();
}

inline OpenHurricane::mixture &OpenHurricane::chemistrySource::mixtures() {
    return flows_.mixtures();
}

inline const OpenHurricane::thermoList &OpenHurricane::chemistrySource::therm() const {
    return flows_.mixtures().thermalTable();
}

inline OpenHurricane::thermoList &OpenHurricane::chemistrySource::therm() {
    return flows_.mixtures().thermalTable();
}

inline OpenHurricane::realArray &OpenHurricane::chemistrySource::diagJacPerReac() {
    if (diagJacPerReacPtr_ == nullptr) {
        diagJacPerReacPtr_ = new realArray();
    }
    return *diagJacPerReacPtr_;
}

inline OpenHurricane::realArray &OpenHurricane::chemistrySource::Ri(const integer i) {
    return Ri()[i];
}

inline const OpenHurricane::realArray &OpenHurricane::chemistrySource::Ri(const integer i) const {
    return Ri()[i];
}

inline bool OpenHurricane::chemistrySource::isReduced() const noexcept {
    return false;
}

#ifdef CUDA_PARALLEL

inline void OpenHurricane::chemistrySource::calculateSourceTermsAsync(real *hur_restrict hostYiRhoTPtr_,
                                                                  cu2DArray<cu_real> &dYiRhoT,
                                                                  real *hur_restrict RRi,
                                                                  cu2DArray<cu_real> &dRi,
                                                                  const cudaStreams &streams) {
    LFatal("This function is not implemented");
}

inline void OpenHurricane::chemistrySource::calculateSourceTermsImpAsync(
    real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT, real *hur_restrict RRi,
    cu2DArray<cu_real> &dRi, real *hur_restrict dRidrhoyi, cu2DArray<cu_real> &ddRidrhoyi,
    const cudaStreams &streams) {
    LFatal("This function is not implemented");
}

inline void OpenHurricane::chemistrySource::calculateSourceTermsImpAsyncHybrid(
    real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT, real *hur_restrict RRi,
    cu2DArray<cu_real> &dRi, cu_float *hur_restrict dRidrhoyi,
    cu2DArray<cu_float> &ddRidrhoyi, const cudaStreams &streams) {
    LFatal("This function is not implemented");
}

inline void OpenHurricane::chemistrySource::createReactionTable() const {
    LFatal("This function is not implemented");
}

inline void OpenHurricane::chemistrySource::createReactionTableAsync(const cudaStreams &streams) const {
    LFatal("This function is not implemented");
}

inline void OpenHurricane::chemistrySource::destroyReactionTable() const {
    LFatal("This function is not implemented");
}

#endif // CUDA_PARALLEL

inline const OpenHurricane::cellRealArray &OpenHurricane::chemistrySource::nSPInEveryCell() const {
    if (nSPInEveryCellPtr_ == nullptr) {
        nSPInEveryCellPtr_ = new cellRealArray(
            object("nSPInEveryCell", mesh(), object::WRITE_OUTPUT), mesh(), nsp());
    }
    return *nSPInEveryCellPtr_;
}

inline OpenHurricane::cellRealArray &OpenHurricane::chemistrySource::nSPInEveryCell() {
    if (nSPInEveryCellPtr_ == nullptr) {
        nSPInEveryCellPtr_ = new cellRealArray(
            object("nSPInEveryCell", mesh(), object::WRITE_OUTPUT), mesh(), nsp());
    }
    return *nSPInEveryCellPtr_;
}

inline const OpenHurricane::cellRealArray &OpenHurricane::chemistrySource::nRCInEveryCell() const {
    if (nRCInEveryCellPtr_ == nullptr) {
        nRCInEveryCellPtr_ = new cellRealArray(
            object("nRCInEveryCell", mesh(), object::WRITE_OUTPUT), mesh(), nrc());
    }
    return *nRCInEveryCellPtr_;
}

inline OpenHurricane::cellRealArray &OpenHurricane::chemistrySource::nRCInEveryCell() {
    if (nRCInEveryCellPtr_ == nullptr) {
        nRCInEveryCellPtr_ = new cellRealArray(
            object("nRCInEveryCell", mesh(), object::WRITE_OUTPUT), mesh(), nrc());
    }
    return *nRCInEveryCellPtr_;
}

hur_nodiscard inline const OpenHurricane::cellRealArray &
OpenHurricane::chemistrySource::cellSourceCalTime() const {
    if (cellSourceCalTimePtr_ == nullptr) {
        cellSourceCalTimePtr_ = new cellRealArray(
            object("cellSourceCalTime", mesh(), object::WRITE_OUTPUT), mesh(), Zero);
    }
    return *cellSourceCalTimePtr_;
}
hur_nodiscard inline OpenHurricane::cellRealArray &OpenHurricane::chemistrySource::cellSourceCalTime() {
    if (cellSourceCalTimePtr_ == nullptr) {
        cellSourceCalTimePtr_ = new cellRealArray(
            object("cellSourceCalTime", mesh(), object::WRITE_OUTPUT), mesh(), Zero);
    }
    return *cellSourceCalTimePtr_;
}

#ifdef TEST_PROCESS_TIME

hur_nodiscard inline OpenHurricane::integer OpenHurricane::chemistrySource::nSpcGroup() const noexcept {
    return 1;
}

hur_nodiscard inline OpenHurricane::integerList
OpenHurricane::chemistrySource::nEachSpcGroup() const noexcept {
    integerList nespg(1);
    nespg[0] = nsp();
    return nespg;
}

#endif // TEST_PROCESS_TIME