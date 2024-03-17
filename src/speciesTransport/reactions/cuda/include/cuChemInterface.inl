#include "cuChemInterface.hpp"
/*!
 * \file cuChemInterface.inl
 * \brief The In-Line functions of the <i>cuChemInterface.hpp</i> file.
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

inline OpenHurricane::cuChem::cuChemInterface::
    speciesTableInterface::speciesTableInterface()
    : WPtr_(nullptr), TCommonPtr_(nullptr), TLowPtr_(nullptr),
      THighPtr_(nullptr), highCpCoeffsPtr_(nullptr), lowCpCoeffsPtr_(nullptr) {}

inline OpenHurricane::cuChem::cuChemInterface::
    speciesTableInterface::~speciesTableInterface() noexcept {
    clear();
}

inline void OpenHurricane::cuChem::cuChemInterface::
    speciesTableInterface::transfer(speciesTableInterface &st) {
    clear();

    WPtr_ = st.WPtr_;
    st.WPtr_ = nullptr;

    TCommonPtr_ = st.TCommonPtr_;
    st.TCommonPtr_ = nullptr;

    TLowPtr_ = st.TLowPtr_;
    st.TLowPtr_ = nullptr;

    THighPtr_ = st.THighPtr_;
    st.THighPtr_ = nullptr;

    highCpCoeffsPtr_ = st.highCpCoeffsPtr_;
    st.highCpCoeffsPtr_ = nullptr;

    lowCpCoeffsPtr_ = st.lowCpCoeffsPtr_;
    st.lowCpCoeffsPtr_ = nullptr;
}

inline void OpenHurricane::cuChem::cuChemInterface::
    speciesTableInterface::clear() noexcept {
    HurDeleteDynArray(WPtr_);
    HurDeleteDynArray(TCommonPtr_);
    HurDeleteDynArray(TLowPtr_);
    HurDeleteDynArray(THighPtr_);
    HurDeleteDynArray(highCpCoeffsPtr_);
    HurDeleteDynArray(lowCpCoeffsPtr_);
}

inline void OpenHurricane::cuChem::cuChemInterface::
    speciesTableInterface::alloc(const cu_ushort nsp) {
    clear();
    WPtr_ = new cu_real[nsp]();
    TCommonPtr_ = new cu_real[nsp]();
    TLowPtr_ = new cu_real[nsp]();
    THighPtr_ = new cu_real[nsp]();
    highCpCoeffsPtr_ = new cu_real[7 * nsp]();
    lowCpCoeffsPtr_ = new cu_real[7 * nsp]();
}

inline OpenHurricane::cuChem::cuChemInterface::
    reactionCoeffsInterface::reactionCoeffsInterface()
    : maxSpcSizeInRec_(0), sizePtr_(nullptr), reactionCoeffsIndexPtr_(nullptr),
      stoichCoeffPtr_(nullptr), orderPtr_(nullptr) {}

inline OpenHurricane::cuChem::cuChemInterface::
    reactionCoeffsInterface::~reactionCoeffsInterface() noexcept {
    clear();
}

inline void OpenHurricane::cuChem::cuChemInterface::
    reactionCoeffsInterface::transfer(reactionCoeffsInterface &rc) {
    clear();

    maxSpcSizeInRec_ = rc.maxSpcSizeInRec_;
    rc.maxSpcSizeInRec_ = 0;

    sizePtr_ = rc.sizePtr_;
    rc.sizePtr_ = nullptr;

    reactionCoeffsIndexPtr_ = rc.reactionCoeffsIndexPtr_;
    rc.reactionCoeffsIndexPtr_ = nullptr;

    stoichCoeffPtr_ = rc.stoichCoeffPtr_;
    rc.stoichCoeffPtr_ = nullptr;

    orderPtr_ = rc.orderPtr_;
    rc.orderPtr_ = nullptr;
}

inline void OpenHurricane::cuChem::cuChemInterface::
    reactionCoeffsInterface::clear() noexcept {
    maxSpcSizeInRec_ = 0;
    HurDeleteDynArray(sizePtr_);
    HurDeleteDynArray(reactionCoeffsIndexPtr_);
    HurDeleteDynArray(stoichCoeffPtr_);
    HurDeleteDynArray(orderPtr_);
}

inline void OpenHurricane::cuChem::cuChemInterface::
    reactionCoeffsInterface::alloc(const cu_ushort nrc,
                                   const cu_ushort maxSpcSizeInReac) {
    clear();

    maxSpcSizeInRec_ = maxSpcSizeInReac;

    sizePtr_ = new cu_ushort[nrc]();

    reactionCoeffsIndexPtr_ = new cu_ushort[nrc * maxSpcSizeInReac]();

    stoichCoeffPtr_ = new cu_real[nrc * maxSpcSizeInReac]();
    orderPtr_ = new cu_real[nrc * maxSpcSizeInReac]();
}

inline OpenHurricane::cuChem::cuChemInterface::
    reactionRateCoeffsInterface::reactionRateCoeffsInterface()
    : aPtr_(nullptr), bPtr_(nullptr), TaPtr_(nullptr) {}

inline OpenHurricane::cuChem::cuChemInterface::
    reactionRateCoeffsInterface::~reactionRateCoeffsInterface() noexcept {
    clear();
}

inline void OpenHurricane::cuChem::cuChemInterface::
    reactionRateCoeffsInterface::transfer(reactionRateCoeffsInterface &rc) {
    clear();

    aPtr_ = rc.aPtr_;
    rc.aPtr_ = nullptr;

    bPtr_ = rc.bPtr_;
    rc.bPtr_ = nullptr;

    TaPtr_ = rc.TaPtr_;
    rc.TaPtr_ = nullptr;
}

inline void OpenHurricane::cuChem::cuChemInterface::
    reactionRateCoeffsInterface::clear() noexcept {
    HurDeleteDynArray(aPtr_);
    HurDeleteDynArray(bPtr_);
    HurDeleteDynArray(TaPtr_);
}

inline void OpenHurricane::cuChem::cuChemInterface::
    reactionRateCoeffsInterface::alloc(const cu_integer nrc) {
    clear();
    aPtr_ = new cu_real[nrc]();
    bPtr_ = new cu_real[nrc]();
    TaPtr_ = new cu_real[nrc]();
}

inline OpenHurricane::cuChem::cuChemInterface::
    thirdBodyCoeffInterface::thirdBodyCoeffInterface()
    : nrcThird_(0), thidrBodyIndexPtr_(nullptr), coefThirdPtr_(nullptr) {}

inline OpenHurricane::cuChem::cuChemInterface::
    thirdBodyCoeffInterface::~thirdBodyCoeffInterface() noexcept {
    clear();
}

inline void OpenHurricane::cuChem::cuChemInterface::
    thirdBodyCoeffInterface::transfer(thirdBodyCoeffInterface &rc) {
    clear();

    nrcThird_ = rc.nrcThird_;
    rc.nrcThird_ = 0;

    thidrBodyIndexPtr_ = rc.thidrBodyIndexPtr_;
    rc.thidrBodyIndexPtr_ = nullptr;

    coefThirdPtr_ = rc.coefThirdPtr_;
    rc.coefThirdPtr_ = nullptr;
}

inline void OpenHurricane::cuChem::cuChemInterface::
    thirdBodyCoeffInterface::clear() noexcept {
    nrcThird_ = 0;
    HurDeleteDynArray(thidrBodyIndexPtr_);
    HurDeleteDynArray(coefThirdPtr_);
}

inline void OpenHurricane::cuChem::cuChemInterface::
    thirdBodyCoeffInterface::alloc(const cu_ushort nsp, const cu_ushort nrc,
                                   const cu_ushort nrcThird) {
    clear();
    nrcThird_ = nrcThird;

    thidrBodyIndexPtr_ = new cu_short[nrc]();
    coefThirdPtr_ = new cu_real[nsp * nrcThird]();

    for (cu_integer i = 0; i < nrc; ++i) {
        thidrBodyIndexPtr_[i] = cu_short(-1);
    }
}

inline OpenHurricane::cuChem::cuChemInterface::
    pressureDepCoeffInterface::pressureDepCoeffInterface()
    : nrcPD_(0), indexPtr_(nullptr), aPtr_(nullptr), bPtr_(nullptr),
      TaPtr_(nullptr), fallOffTypePtr_(nullptr), fallOffCoeffPtr_(nullptr) {}

inline OpenHurricane::cuChem::cuChemInterface::
    pressureDepCoeffInterface::~pressureDepCoeffInterface() noexcept {
    clear();
}

inline void OpenHurricane::cuChem::cuChemInterface::
    pressureDepCoeffInterface::transfer(pressureDepCoeffInterface &rc) {
    clear();
    nrcPD_ = rc.nrcPD_;
    rc.nrcPD_ = 0;

    indexPtr_ = rc.indexPtr_;
    rc.indexPtr_ = nullptr;

    aPtr_ = rc.aPtr_;
    rc.aPtr_ = nullptr;

    bPtr_ = rc.bPtr_;
    rc.bPtr_ = nullptr;

    TaPtr_ = rc.TaPtr_;
    rc.TaPtr_ = nullptr;

    fallOffTypePtr_ = rc.fallOffTypePtr_;
    rc.fallOffTypePtr_ = nullptr;

    fallOffCoeffPtr_ = rc.fallOffCoeffPtr_;
    rc.fallOffCoeffPtr_ = nullptr;
}

inline void OpenHurricane::cuChem::cuChemInterface::
    pressureDepCoeffInterface::clear() noexcept {
    nrcPD_ = 0;
    HurDeleteDynArray(indexPtr_);
    HurDeleteDynArray(aPtr_);
    HurDeleteDynArray(bPtr_);
    HurDeleteDynArray(TaPtr_);
    HurDeleteDynArray(fallOffTypePtr_);
    HurDeleteDynArray(fallOffCoeffPtr_);
}

inline void OpenHurricane::cuChem::cuChemInterface::
    pressureDepCoeffInterface::alloc(const cu_ushort nrc,
                                     const cu_ushort nrcPD) {
    clear();

    nrcPD_ = nrcPD;
    indexPtr_ = new cu_short[nrc]();
    aPtr_ = new cu_real[nrcPD]();
    bPtr_ = new cu_real[nrcPD]();
    TaPtr_ = new cu_real[nrcPD]();
    fallOffTypePtr_ = new cu_ushort[nrcPD]();
    fallOffCoeffPtr_ = new cu_real[nrcPD * 5]();

    for (cu_ushort i = 0; i < nrc; ++i) {
        indexPtr_[i] = cu_short(-1);
    }
}

inline OpenHurricane::cuChem::cuChemInterface::
    cuChemInterface()
    : nsp_(0), nrc_(0), sptInt_(), reacInt_(), revReacInt_(), thirdBInt_(),
      pressDepInt_(), sumStoiIntPtr_(nullptr), reactionTypeIntPtr_(nullptr),
      thirdBodyTypeIntPtr_(nullptr) {}

inline OpenHurricane::cuChem::cuChemInterface::
    ~cuChemInterface() noexcept {
    clear();
}

inline void OpenHurricane::cuChem::cuChemInterface::transfer(
    cuChemInterface &rc) {
    clear();

    nsp_ = rc.nsp_;
    rc.nsp_ = 0;

    nrc_ = rc.nrc_;
    rc.nrc_ = 0;

    sptInt_.transfer(rc.sptInt_);
    reacInt_.transfer(rc.reacInt_);
    revReacInt_.transfer(rc.revReacInt_);
    thirdBInt_.transfer(rc.thirdBInt_);
    pressDepInt_.transfer(rc.pressDepInt_);
    kfInt_.transfer(rc.kfInt_);
    rfInt_.transfer(rc.rfInt_);

    sumStoiIntPtr_ = rc.sumStoiIntPtr_;
    rc.sumStoiIntPtr_ = nullptr;

    reactionTypeIntPtr_ = rc.reactionTypeIntPtr_;
    rc.reactionTypeIntPtr_ = nullptr;

    thirdBodyTypeIntPtr_ = rc.thirdBodyTypeIntPtr_;
    rc.thirdBodyTypeIntPtr_ = nullptr;
}

inline void
OpenHurricane::cuChem::cuChemInterface::clear() noexcept {
    sptInt_.clear();
    reacInt_.clear();
    revReacInt_.clear();
    thirdBInt_.clear();
    pressDepInt_.clear();

    kfInt_.clear();
    rfInt_.clear();

    HurDeleteDynArray(sumStoiIntPtr_);
    HurDeleteDynArray(reactionTypeIntPtr_);
    HurDeleteDynArray(thirdBodyTypeIntPtr_);
}

inline void OpenHurricane::cuChem::cuChemInterface::allocOnlyThis(
    const cu_ushort nsp, const cu_ushort nrc) {
    HurDeleteDynArray(sumStoiIntPtr_);
    HurDeleteDynArray(reactionTypeIntPtr_);
    HurDeleteDynArray(thirdBodyTypeIntPtr_);

    nsp_ = nsp;
    nrc_ = nrc;

    sumStoiIntPtr_ = new cu_real[nrc]();
    reactionTypeIntPtr_ = new cu_ushort[nrc]();
    thirdBodyTypeIntPtr_ = new cu_ushort[nrc]();
}

#endif // CUDA_PARALLEL