/*!
 * \file transportList.cpp
 * \brief Main subroutines for transport properties table.
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
#include "transportList.hpp"

hur_nodiscard const OpenHurricane::realSquareMatrix &OpenHurricane::transportList::phiPre() const {
    if (phiPrePtr_ == nullptr) {
        phiPrePtr_ = new realSquareMatrix(species().size());
        const real c8 = 1.0 / sqrt(real(8.0));
        for (integer i = 0; i < species().size(); ++i) {
            for (integer j = 0; j < species().size(); ++j) {
                (*phiPrePtr_)(i, j) = c8 / sqrt((real(1.0) + species_[i].W() / species_[j].W()));
            }
        }
    }

    return *phiPrePtr_;
}

hur_nodiscard const OpenHurricane::realSquareMatrix &OpenHurricane::transportList::WiWj() const {
    if (WiWjPtr_ == nullptr) {
        WiWjPtr_ = new realSquareMatrix(species().size());
        for (integer i = 0; i < species().size(); ++i) {
            for (integer j = 0; j < species().size(); ++j) {
                (*WiWjPtr_)(i, j) = sqrt(species_[i].W() / species_[j].W());
            }
        }
    }

    return *WiWjPtr_;
}

OpenHurricane::transportList::transportList()
    : species_(speciesList::nullObject()), tranList_(), txi_(), tmuf_(), tkppaf_(), tphi_(),
      phiPrePtr_(nullptr), WiWjPtr_(nullptr) {}

OpenHurricane::transportList::transportList(const speciesList &species)
    : species_(species), tranList_(), txi_(), tmuf_(), tkppaf_(), tphi_(), phiPrePtr_(nullptr),
      WiWjPtr_(nullptr) {}

OpenHurricane::transportList::transportList(const speciesList &species, const controller &cont)
    : species_(species), tranList_(species.size()), txi_(), tmuf_(), tkppaf_(), tphi_(),
      phiPrePtr_(nullptr), WiWjPtr_(nullptr) {
    for (integer i = 0; i < species.size(); ++i) {
        if (cont.found(transport::className_)) {
            if (cont.subController(transport::className_).found(species_[i].name())) {
                tranList_.set(
                    i,
                    transport::creator(
                        species, i,
                        cont.subController(transport::className_).subController(species[i].name()))
                        .release());
            } else {
                tranList_.set(
                    i, transport::creator(species, i, cont.subController(transport::className_))
                           .release());
            }
        } else {
            tranList_.set(
                i,
                transport::creator(species, i, cont.subController(species_[i].name())).release());
        }
    }
}

OpenHurricane::transportList::transportList(const transportList &tt)
    : species_(tt.species_), tranList_(tt.tranList_), txi_(), tmuf_(), tkppaf_(), tphi_(),
      phiPrePtr_(nullptr), WiWjPtr_(nullptr) {}

OpenHurricane::transportList::transportList(transportList &&tt) noexcept
    : species_(tt.species_), tranList_(std::move(tt.tranList_)), txi_(), tmuf_(), tkppaf_(),
      tphi_(), phiPrePtr_(nullptr), WiWjPtr_(nullptr) {}

OpenHurricane::transportList::~transportList() noexcept {
    HurDelete(phiPrePtr_);
    HurDelete(WiWjPtr_);
}

hur_nodiscard OpenHurricane::real OpenHurricane::transportList::mu(const real p, const real T,
                                                            const integer i) const {
    return tranList_[i].mu(p, T);
}

hur_nodiscard OpenHurricane::real OpenHurricane::transportList::mu(const real p, const real T,
                                                            const realArray &yi) const {
    const integer nsp = species_.size();
    if (nsp == 1) {
        return tranList_[0].mu(p, T);
    }
    const auto xi = species_.Yi2Xi(yi);
    auto tmpMu = muf(p, T);
    // Wilke expressions

    realArray phi(nsp, Zero);
    const auto &phiP = phiPre();
    const auto &wiwj = WiWj();
    for (integer i = 0; i < nsp; ++i) {
        for (integer j = 0; j < nsp; ++j) {
            //real temp = (1.0 + sqrt(tmpMu[i] / tmpMu[j] * sqrt(species_[j].W() / species_[i].W())));
            real temp = (1.0 + sqrt(tmpMu[i] / tmpMu[j] * wiwj(j, i)));
            //phi[i] += xi[j] * temp * temp / sqrt(8.0 * (1.0 + species_[i].W() / species_[j].W()));
            phi[i] += phiP(i, j) * xi[j] * sqr(temp);
        }
    }

    real mum = Zero;
    for (integer i = 0; i < nsp; ++i) {
        mum += tmpMu[i] * xi[i] / phi[i];
    }

    return mum;
}

hur_nodiscard OpenHurricane::real OpenHurricane::transportList::mu(const real p, const real T,
                                                            const PtrList<cellRealArray> &yi,
                                                            const integer cellI) const {
    const integer nsp = species_.size();
    if (nsp == 1) {
        return tranList_[0].mu(p, T);
    }
    const auto xi = species_.Yi2Xi(yi, cellI);
    auto tmpMu = muf(p, T);
    // Wilke expressions

    realArray phi(nsp, Zero);
    const auto &phiP = phiPre();
    const auto &wiwj = WiWj();
    for (integer i = 0; i < nsp; ++i) {
        for (integer j = 0; j < nsp; ++j) {
            //real temp = (1.0 + sqrt(tmpMu[i] / tmpMu[j] * sqrt(species_[j].W() / species_[i].W())));
            real temp = (1.0 + sqrt(tmpMu[i] / tmpMu[j] * wiwj(j, i)));
            //phi[i] += xi[j] * temp * temp / sqrt(8.0 * (1.0 + species_[i].W() / species_[j].W()));
            phi[i] += phiP(i, j) * xi[j] * sqr(temp);
        }
    }

    real mum = Zero;
    for (integer i = 0; i < nsp; ++i) {
        mum += tmpMu[i] * xi[i] / phi[i];
    }

    return mum;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::transportList::muf(const real p,
                                                                  const real T) const {
    realArray tmpMu(species_.size());
    for (integer i = 0; i < species_.size(); ++i) {
        tmpMu[i] = mu(p, T, i);
    }
    return tmpMu;
}

hur_nodiscard OpenHurricane::real OpenHurricane::transportList::nu(const integer i, const real p,
                                                            const real T, const real rho) const {
    return tranList_[i].mu(p, T) / rho;
}

hur_nodiscard OpenHurricane::real OpenHurricane::transportList::kappa(const real p, const real T,
                                                               const real cpi,
                                                               const integer i) const {
    return tranList_[i].kappa(p, T, cpi);
}

hur_nodiscard OpenHurricane::real OpenHurricane::transportList::kappa(const real p, const real T,
                                                               const real mui, const real cpi,
                                                               const integer i) const {
    return tranList_[i].kappa(p, T, mui, cpi);
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::transportList::kappaf(const real p, const real T,
                                                                     const realArray &mui,
                                                                     const realArray &cpi) const {
    realArray tmpKappa(species_.size());
    for (integer i = 0; i < species_.size(); ++i) {
        tmpKappa[i] = kappa(p, T, mui[i], cpi[i], i);
    }
    return tmpKappa;
}

hur_nodiscard OpenHurricane::real OpenHurricane::transportList::kappa(const real p, const real T,
                                                               const realArray &mui,
                                                               const realArray &cpi,
                                                               const realArray &yi) const {
    const integer nsp = species_.size();
    auto xi = species_.Yi2Xi(yi);
    auto tmpKappa = kappaf(p, T, mui, cpi);
    if (nsp == 1) {
        return tranList_[0].kappa(p, T, mui[0], cpi[0]);
    }
    // Wassilewa expressions

    realArray phi(nsp, Zero);

    const auto &phiP = phiPre();
    const auto &wiwj = WiWj();
    for (integer i = 0; i < nsp; ++i) {
        for (integer j = 0; j < nsp; ++j) {
            //real temp = (1.0 + sqrt(mui[i] / mui[j] * sqrt(species_[j].W() / species_[i].W())));
            real temp = (1.0 + sqrt(mui[i] / mui[j] * wiwj(j, i)));
            phi[i] += phiP(i, j) * xi[j] * sqr(temp);
        }
    }

    real kappam = Zero;
    for (integer i = 0; i < nsp; ++i) {
        kappam += tmpKappa[i] * xi[i] / phi[i];
    }

    return kappam;
}

void OpenHurricane::transportList::muKappa(const real p, const real T, const realArray &cpi,
                                        const realArray &yi, real &mum, real &kappam) const {
    const integer nsp = species_.size();
    if (nsp == 1) {
        mum = mu(p, T, 0);
        kappam = tranList_[0].kappa(p, T, mum, cpi[0]);
        return;
    }
    if (txi_.size() == 0) {
        txi_.resize(nsp, Zero);
    }
    auto &xi = txi_;
    species_.Yi2Xi(yi, xi);
    if (tmuf_.size() == 0) {
        tmuf_.resize(nsp, Zero);
    }
    auto &mui = tmuf_;
    muf(p, T, mui);
    if (tkppaf_.size() == 0) {
        tkppaf_.resize(nsp, Zero);
    }
    auto &kappai = tkppaf_;
    kappaf(p, T, mui, cpi, kappai);
    if (tphi_.size() == 0) {
        tphi_.resize(nsp, Zero);
    }
    auto &phi = tphi_;
    phi = Zero;
    const auto &phiP = phiPre();
    const auto &wiwj = WiWj();
    for (integer i = 0; i < nsp; ++i) {
        for (integer j = 0; j < nsp; ++j) {
            //real temp = (1.0 + sqrt(mui[i] / mui[j] * sqrt(species_[j].W() / species_[i].W())));
            real temp = (1.0 + sqrt(mui[i] / mui[j] * wiwj(j, i)));
            //phi[i] += xi[j] * temp * temp / sqrt(8.0 * (1.0 + species_[i].W() / species_[j].W()));
            phi[i] += phiP(i, j) * xi[j] * sqr(temp);
        }
    }

    mum = Zero;
    kappam = Zero;
    for (integer i = 0; i < nsp; ++i) {
        mum += mui[i] * xi[i] / phi[i];
        //kappam += kappai[i] * xi[i] / (1.0650 * phi[i]);
        kappam += kappai[i] * xi[i] / phi[i];
    }
}

void OpenHurricane::transportList::muKappa(const real p, const real T, const realArray &cpi,
                                        const PtrList<cellRealArray> &yi, const integer cellI,
                                        real &mum, real &kappam) const {
    const integer nsp = species_.size();
    if (nsp == 1) {
        mum = mu(p, T, 0);
        kappam = kappa(p, T, mum, cpi[0], 0);
        return;
    }
    if (txi_.size() == 0) {
        txi_.resize(nsp, Zero);
    }
    auto &xi = txi_;
    species_.Yi2Xi(yi, cellI, xi);
    if (tmuf_.size() == 0) {
        tmuf_.resize(nsp, Zero);
    }
    auto &mui = tmuf_;
    muf(p, T, mui);

    if (tkppaf_.size() == 0) {
        tkppaf_.resize(nsp, Zero);
    }
    auto &kappai = tkppaf_;
    kappaf(p, T, mui, cpi, kappai);
    if (tphi_.size() == 0) {
        tphi_.resize(nsp, Zero);
    }
    auto &phi = tphi_;
    phi = Zero;
    const auto &phiP = phiPre();
    const auto &wiwj = WiWj();
    for (integer i = 0; i < nsp; ++i) {
        for (integer j = 0; j < nsp; ++j) {
            //real temp = (1.0 + sqrt(mui[i] / mui[j] * sqrt(species_[j].W() / species_[i].W())));
            real temp = (1.0 + sqrt(mui[i] / mui[j] * wiwj(j, i)));
            //phi[i] += xi[j] * temp * temp / sqrt(8.0 * (1.0 + species_[i].W() / species_[j].W()));
            phi[i] += phiP(i, j) * xi[j] * sqr(temp);
        }
    }

    mum = Zero;
    kappam = Zero;
    for (integer i = 0; i < nsp; ++i) {
        mum += mui[i] * xi[i] / phi[i];
        //kappam += kappai[i] * xi[i] / (1.0650 * phi[i]);
        kappam += kappai[i] * xi[i] / phi[i];
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::transportList::DiffMix(const real p, const real T,
                                                                 const PtrList<cellRealArray> &yi,
                                                                 const PtrList<cellRealArray> &Dim,
                                                                 const integer cellI) const {
    const integer nsp = species_.size();

    if (nsp == 1) {
        return Dim[0][cellI];
    }
    if (txi_.size() == 0) {
        txi_.resize(nsp, Zero);
    }
    auto &xi = txi_;
    species_.Yi2Xi(yi, cellI, xi);
    if (tmuf_.size() == 0) {
        tmuf_.resize(nsp, Zero);
    }
    auto &mui = tmuf_;
    muf(p, T, mui);

    if (tphi_.size() == 0) {
        tphi_.resize(nsp, Zero);
    }
    auto &phi = tphi_;
    phi = Zero;
    const real c8 = 1.0 / sqrt(real(8.0));
    for (integer i = 0; i < nsp; ++i) {
        for (integer j = 0; j < nsp; ++j) {
            real temp = (1.0 + sqrt(mui[i] / mui[j] * sqrt(species_[j].W() / species_[i].W())));
            phi[i] += c8 * xi[j] * temp * temp / sqrt((1.0 + species_[i].W() / species_[j].W()));
        }
    }

    real Dmm = Zero;
    for (integer i = 0; i < nsp; ++i) {
        Dmm += Dim[i][cellI] * xi[i] / phi[i];
    }
    return Dmm;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::transportList::Df(const real p, const real T,
                                                                 const realArray &yi) const {
    realArray tmpD(species_.size());
    for (integer i = 0; i < species_.size(); ++i) {
        tmpD[i] = D(p, T, yi, i);
    }
    return tmpD;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::transportList::Df(const real p, const real T,
                                                                 const PtrList<cellRealArray> &yi,
                                                                 const integer cellI) const {
    realArray tmpD(species_.size());
    for (integer i = 0; i < species_.size(); ++i) {
        tmpD[i] = D(p, T, yi, cellI, i);
    }
    return tmpD;
}

OpenHurricane::realArray OpenHurricane::transportList::Df(const real p, const real T,
                                                   const PtrList<cellRealArray> &yi,
                                                   const integer cellI, realArray &Dij) const {
    realArray tmpD(species_.size());
    for (integer i = 0; i < species_.size(); ++i) {
        tmpD[i] = D(p, T, yi, cellI, i, Dij);
    }
    return tmpD;
}

void OpenHurricane::transportList::Df(realArray &Dff, const real p, const real T,
                                   const PtrList<cellRealArray> &yi, const integer cellI) const {
    if (txi_.size() == 0) {
        const integer nsp = species_.size();
        txi_.resize(nsp, Zero);
    }
    auto &xi = txi_;
    species_.Yi2Xi(yi, cellI, xi);
    const real Tv1p5 = T * sqrt(T);
    for (integer i = 0; i < species_.size(); ++i) {
        Dff[i] = D_Xi(p, T, Tv1p5, xi, i);
    }
}

void OpenHurricane::transportList::Df(realArray &Dff, const real p, const real T,
                                   const PtrList<cellRealArray> &yi, const integer cellI,
                                   realArray &Dij) const {
    auto xi = species_.Yi2Xi(yi, cellI);
    const real Tv1p5 = T * sqrt(T);
    for (integer i = 0; i < species_.size(); ++i) {
        Dff[i] = D_Xi(p, T, Tv1p5, xi, i, Dij);
    }
}
