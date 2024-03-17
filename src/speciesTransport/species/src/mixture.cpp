/*!
 * \file mixture.cpp
 * \brief Main subroutines for mixture.
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
#include "mixture.hpp"
#include "faceInterpolation.hpp"
#include "perfectGas.hpp"
#include "thermoList.hpp"

#include "chemkinFileRead.hpp"

OpenHurricane::basicMixture::basicMixture(const runtimeMesh &mesh, const integer speciesSize)
    : isSingular_(false), mesh_(mesh), species_(speciesSize), Yi_(speciesSize), Dim_(speciesSize),
      hi_(speciesSize) {}

OpenHurricane::basicMixture::basicMixture(const runtimeMesh &mesh, const speciesList &species)
    : isSingular_(false), mesh_(mesh), species_(species), Yi_(species.size()), Dim_(species.size()),
      hi_(species.size()) {
    if (species.size() == 1 || species.size() == 0) {
        isSingular_ = true;
    } else {
        isSingular_ = false;
    }
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::basicMixture::Yic(const integer cellId) const {
    realArray Yii(species_.size(), Zero);
    if (species_.size() == 1) {
        Yii = 1.0;
        return Yii;
    }
    for (integer i = 0; i < species_.size(); ++i) {
        Yii[i] = Yi_[i][cellId];
    }

    return Yii;
}

hur_nodiscard OpenHurricane::realArrayArray
OpenHurricane::basicMixture::Yif(const integer faceZoneId) const {
    const auto &fZL = mesh_.faceZones();
    const auto &fs = mesh_.faces();
    const auto &fw = mesh_.faceWgt();
    realArrayArray tmpYi(fZL[faceZoneId].size());
    if (isSingular_) {
        for (integer fi = 0; fi < fZL[faceZoneId].size(); ++fi) {
            tmpYi[fi].resize(1, real(1.0));
        }
    } else {
        for (integer fi = 0; fi < fZL[faceZoneId].size(); ++fi) {
            const integer index = fi + fZL[faceZoneId].firstIndex();
            const integer cl = fs[index].leftCell();
            const integer cr = fs[index].rightCell();
            tmpYi[fi].resize(species_.size());
            const real wl = fw[index];
            const real wr = 1.0 - wl;
            real sumYi = 0.0;
            for (integer isp = 0; isp < species_.size(); ++isp) {
                tmpYi[fi][isp] = wl * Yi_[isp][cl] + wr * Yi_[isp][cr];
                sumYi += tmpYi[fi][isp];
            }
            if (sumYi == real(0)) {
                LFatal("The summation of all species mass fraction is zero in face zone: %s",
                       fZL[faceZoneId].name().c_str());
            }
            tmpYi[fi] /= sumYi;
        }
    }
    return tmpYi;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::basicMixture::Wi() const {
    const integer nsp = species_.size();
    realArray tWi(nsp, Zero);
    for (integer i = 0; i < nsp; ++i) {
        tWi[i] = species_[i].W();
    }
    return tWi;
}

void OpenHurricane::basicMixture::setYi(const runtimeMesh &mesh) {
    // No need to set mass fraction: yi for singular species
    if (isSingular_) {
        return;
    }
    if (Yi_.size() != species_.size()) {
        Yi_.resize(species_.size());
    }
    for (integer i = 0; i < species_.size(); ++i) {
        Yi_.set(i, new cellRealArray(object(species_[i].name(), mesh, object::WRITE_RELAY_OUTPUT,
                                            object::PRIMITIVE),
                                     mesh));
    }

    if (Dim_.size() != species_.size()) {
        Dim_.resize(species_.size());
    }

    for (integer i = 0; i < species_.size(); ++i) {
        Dim_.set(i, new cellRealArray(object(species_[i].name() + "_Dim", mesh, object::NOT_WRITE),
                                      mesh));
    }

    if (hi_.size() != species_.size()) {
        hi_.resize(species_.size());
    }

    for (integer i = 0; i < species_.size(); ++i) {
        hi_.set(i, new cellRealArray(object(species_[i].name() + "_hi", mesh, object::NOT_WRITE),
                                     mesh));
    }
}

void OpenHurricane::basicMixture::setYi(const runtimeMesh &mesh, const realArray &Yi0) {
    // No need to set mass fraction: yi for singular species
    if (isSingular_) {
        return;
    }

    for (integer i = 0; i < species_.size(); ++i) {
        Yi_.set(i, new cellRealArray(object(species_[i].name(), mesh, object::WRITE_RELAY_OUTPUT,
                                            object::PRIMITIVE),
                                     mesh, Yi0[i]));
    }

    if (hi_.size() != species_.size()) {
        hi_.resize(species_.size());
    }

    for (integer i = 0; i < species_.size(); ++i) {
        hi_.set(i, new cellRealArray(object(species_[i].name() + "_hi", mesh, object::NOT_WRITE),
                                     mesh));
    }
}

void OpenHurricane::basicMixture::lastSpeAndNormalized() {
    for (integer cellI = 0; cellI < mesh_.nCells(); ++cellI) {
        real yisum = Zero;
        for (integer i = 0; i < species_.size() - 1; i++) {
            Yi_[i][cellI] = min(real(1.0), max(real(0.0), Yi_[i][cellI]));
            yisum += Yi_[i][cellI];
        }
        Yi_[species_.size() - 1][cellI] = max(real(0.0), real(1.0) - yisum);
        yisum += Yi_[species_.size() - 1][cellI];
        for (integer i = 0; i < species_.size(); i++) {
            Yi_[i][cellI] /= yisum;
        }
    }
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::basicMixture::GFO(const stringList &fuelName, const stringList &oxygenName) const {
    realArray gfo(mesh_.nTotalCells(), Zero);
    if (isSingular_) {
        return gfo;
    }

    integerList fuelId;
    for (integer i = 0; i < fuelName.size(); ++i) {
        integer id;
        if (!species_.contains(fuelName[i], id)) {
            LFatal("Cannot found fuel: %s in species list", fuelName[i].c_str());
        }
        fuelId.append(id);
    }

    if (fuelId.size() == 0) {
        LFatal("Fuel is not given for computing Takeno flame index");
    }
    integerList oxyId;
    for (integer i = 0; i < oxygenName.size(); ++i) {
        integer id;
        if (!species_.contains(oxygenName[i], id)) {
            LFatal("Cannot found oxygen: %s in species list", oxygenName[i].c_str());
        }
        oxyId.append(id);
    }
    if (oxyId.size() == 0) {
        LFatal("Oxygen is not given for computing Takeno flame index");
    }
    if (fuelId.size() == 1 && oxyId.size() == 1) {
        const auto fid = fuelId[0];
        const auto oid = oxyId[0];

        const auto &fyi = Yi(fid);
        const auto &oyi = Yi(oid);

        const auto &fg = fyi.grad();
        const auto &og = oyi.grad();

        for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
            gfo[celli] = (fg[celli] * og[celli]) /
                         max((fg[celli].magnitude() * og[celli].magnitude()), veryTiny);
        }
    } else if (fuelId.size() == 1) {
        const auto fid = fuelId[0];
        const auto &fyi = Yi(fid);
        const auto &fg = fyi.grad();

        vectorArray og(mesh_.nCells(), Zero);
        for (integer oi = 0; oi < oxyId.size(); ++oi) {
            const auto oid = oxyId[oi];
            const auto &oyi = Yi(oid);
            const auto &ogi = oyi.grad();

            for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
                og[celli] += ogi[celli];
            }
        }
        for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
            gfo[celli] = (fg[celli] * og[celli]) /
                         max((fg[celli].magnitude() * og[celli].magnitude()), veryTiny);
        }
    } else if (oxyId.size() == 1) {
        const auto oid = oxyId[0];
        const auto &oyi = Yi(oid);
        const auto &og = oyi.grad();

        vectorArray fg(mesh_.nCells(), Zero);
        for (integer fi = 0; fi < fuelId.size(); ++fi) {
            const auto fid = fuelId[fi];
            const auto &fyi = Yi(fid);
            const auto &fgi = fyi.grad();

            for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
                fg[celli] += fgi[celli];
            }
        }
        for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
            gfo[celli] = (fg[celli] * og[celli]) /
                         max((fg[celli].magnitude() * og[celli].magnitude()), veryTiny);
        }
    } else {
        vectorArray fg(mesh_.nCells(), Zero);
        for (integer fi = 0; fi < fuelId.size(); ++fi) {
            const auto fid = fuelId[fi];
            const auto &fyi = Yi(fid);
            const auto &fgi = fyi.grad();

            for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
                fg[celli] += fgi[celli];
            }
        }
        vectorArray og(mesh_.nCells(), Zero);
        for (integer oi = 0; oi < oxyId.size(); ++oi) {
            const auto oid = oxyId[oi];
            const auto &oyi = Yi(oid);
            const auto &ogi = oyi.grad();

            for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
                og[celli] += ogi[celli];
            }
        }
        for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
            gfo[celli] = (fg[celli] * og[celli]) /
                         max((fg[celli].magnitude() * og[celli].magnitude()), veryTiny);
        }
    }

    return gfo;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::basicMixture::nGFO(const stringList &fuelName, const stringList &oxygenName) const {
    realArray gfo(mesh_.nTotalCells(), Zero);
    if (isSingular_) {
        return gfo;
    }

    integerList fuelId;
    for (integer i = 0; i < fuelName.size(); ++i) {
        integer id;
        if (!species_.contains(fuelName[i], id)) {
            LFatal("Cannot found fuel: %s in species list", fuelName[i].c_str());
        }
        fuelId.append(id);
    }

    if (fuelId.size() == 0) {
        LFatal("Fuel is not given for computing Takeno flame index");
    }
    integerList oxyId;
    for (integer i = 0; i < oxygenName.size(); ++i) {
        integer id;
        if (!species_.contains(oxygenName[i], id)) {
            LFatal("Cannot found oxygen: %s in species list", oxygenName[i].c_str());
        }
        oxyId.append(id);
    }
    if (oxyId.size() == 0) {
        LFatal("Oxygen is not given for computing Takeno flame index");
    }
    if (fuelId.size() == 1 && oxyId.size() == 1) {
        const auto fid = fuelId[0];
        const auto oid = oxyId[0];

        const auto &fyi = Yi(fid);
        const auto &oyi = Yi(oid);

        const auto &fg = fyi.grad();
        const auto &og = oyi.grad();

        for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
            gfo[celli] = (fg[celli] * og[celli]) / max(mag(fg[celli] * og[celli]), veryTiny);
        }
    } else if (fuelId.size() == 1) {
        const auto fid = fuelId[0];
        const auto &fyi = Yi(fid);
        const auto &fg = fyi.grad();

        vectorArray og(mesh_.nCells(), Zero);
        for (integer oi = 0; oi < oxyId.size(); ++oi) {
            const auto oid = oxyId[oi];
            const auto &oyi = Yi(oid);
            const auto &ogi = oyi.grad();

            for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
                og[celli] += ogi[celli];
            }
        }
        for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
            gfo[celli] = (fg[celli] * og[celli]) / max(mag(fg[celli] * og[celli]), veryTiny);
        }
    } else if (oxyId.size() == 1) {
        const auto oid = oxyId[0];
        const auto &oyi = Yi(oid);
        const auto &og = oyi.grad();

        vectorArray fg(mesh_.nCells(), Zero);
        for (integer fi = 0; fi < fuelId.size(); ++fi) {
            const auto fid = fuelId[fi];
            const auto &fyi = Yi(fid);
            const auto &fgi = fyi.grad();

            for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
                fg[celli] += fgi[celli];
            }
        }
        for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
            gfo[celli] = (fg[celli] * og[celli]) / max(mag(fg[celli] * og[celli]), veryTiny);
        }
    } else {
        vectorArray fg(mesh_.nCells(), Zero);
        for (integer fi = 0; fi < fuelId.size(); ++fi) {
            const auto fid = fuelId[fi];
            const auto &fyi = Yi(fid);
            const auto &fgi = fyi.grad();

            for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
                fg[celli] += fgi[celli];
            }
        }
        vectorArray og(mesh_.nCells(), Zero);
        for (integer oi = 0; oi < oxyId.size(); ++oi) {
            const auto oid = oxyId[oi];
            const auto &oyi = Yi(oid);
            const auto &ogi = oyi.grad();

            for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
                og[celli] += ogi[celli];
            }
        }
        for (integer celli = 0; celli < mesh_.nCells(); ++celli) {
            gfo[celli] = (fg[celli] * og[celli]) / max(mag(fg[celli] * og[celli]), veryTiny);
        }
    }

    return gfo;
}

#ifdef CUDA_PARALLEL
#include "JANAF.hpp"
#include "kineticTheory.hpp"

void OpenHurricane::mixture::cuAllocReaction() {
    const auto nsp = species_.size();
    const auto nrc = reactionPtr_->size();
    const auto &react = *reactionPtr_;

    cuChemInterfacePtr_->allocOnlyThis(nsp, nrc);
    cuChemInterfacePtr_->thirdBInt_.alloc(nsp, nrc, nrc);
    cuChemInterfacePtr_->pressDepInt_.alloc(nrc, nrc);
    integer maxSize = 0;
    integer maxRSize = 0;
    for (integer ir = 0; ir < nrc; ++ir) {
        maxSize = max(react[ir].forwardCoeffs().size(), maxSize);
        maxRSize = max(react[ir].backwardCoeffs().size(), maxRSize);
    }

    cuChemInterfacePtr_->reacInt_.alloc(nrc, maxSize);
    cuChemInterfacePtr_->revReacInt_.alloc(nrc, maxRSize);
    cuChemInterfacePtr_->kfInt_.alloc(nrc);
    cuChemInterfacePtr_->rfInt_.alloc(nrc);
}

void OpenHurricane::mixture::cuParsingReaction() {
    cuAllocReaction();
    const auto &react = *reactionPtr_;
    const auto nrc = react.size();

    for (integer i = 0; i < nrc; ++i) {
        react[i].cuSetReactionType(i, nrc, *cuChemInterfacePtr_);
    }
    cuCorrectReaction();
}

void OpenHurricane::mixture::cuCorrectReaction() {
    const auto &react = *reactionPtr_;
    auto &cuChemInt = *cuChemInterfacePtr_;
    for (integer ir = 0; ir < react.size(); ++ir) {
        for (integer i = 0; i < react[ir].forwardCoeffs().size(); ++i) {
            cuChemInt.sumStoiIntPtr_[ir] -= react[ir].forwardCoeffs()[i].stoichCoeff_;
        }
        for (integer i = 0; i < react[ir].backwardCoeffs().size(); ++i) {
            cuChemInt.sumStoiIntPtr_[ir] += react[ir].backwardCoeffs()[i].stoichCoeff_;
        }
    }

    integer countThird = 0;
    for (integer ir = 0; ir < react.size(); ++ir) {
        if (cuChemInt.thirdBInt_.thidrBodyIndexPtr_[ir] != -1) {
            countThird++;
        }
    }
    if (countThird == 0) {
        cuChemInt.thirdBInt_.clear();
        cuChemInt.thirdBInt_.thidrBodyIndexPtr_ = new cu_short[react.size()]();
        for (integer ir = 0; ir < react.size(); ++ir) {
            cuChemInt.thirdBInt_.thidrBodyIndexPtr_[ir] = cu_short(-1);
        }
    } else {
        cuChem::cuChemInterface::thirdBodyCoeffInterface thi;
        thi.alloc(species_.size(), react.size(), countThird);

        countThird = 0;
        for (integer ir = 0; ir < react.size(); ++ir) {
            if (cuChemInt.thirdBInt_.thidrBodyIndexPtr_[ir] != -1) {
                const auto id = cuChemInt.thirdBInt_.thidrBodyIndexPtr_[ir];
                thi.thidrBodyIndexPtr_[ir] = countThird;
                for (integer is = 0; is < species_.size(); ++is) {
                    thi.coefThirdPtr_[is * thi.nrcThird_ + countThird] =
                        cuChemInt.thirdBInt_.coefThirdPtr_[is * react.size() + id];
                }
                countThird++;
            }
        }
        cuChemInt.thirdBInt_.transfer(thi);
    }
    integer countPD = 0;
    for (integer ir = 0; ir < react.size(); ++ir) {
        if (cuChemInt.pressDepInt_.indexPtr_[ir] != -1) {
            countPD++;
        }
    }
    if (countPD == 0) {
        cuChemInt.pressDepInt_.clear();
        cuChemInt.pressDepInt_.indexPtr_ = new cu_short[react.size()]();
        for (integer ir = 0; ir < react.size(); ++ir) {
            cuChemInt.pressDepInt_.indexPtr_[ir] = cu_short(-1);
        }
    } else {
        cuChem::cuChemInterface::pressureDepCoeffInterface pdd;
        pdd.alloc(react.size(), countPD);
        countPD = 0;
        for (integer ir = 0; ir < react.size(); ++ir) {
            if (cuChemInt.pressDepInt_.indexPtr_[ir] != -1) {
                const cu_integer id = cuChemInt.pressDepInt_.indexPtr_[ir];
                pdd.indexPtr_[ir] = countPD;
                pdd.aPtr_[countPD] = cuChemInt.pressDepInt_.aPtr_[id];
                pdd.bPtr_[countPD] = cuChemInt.pressDepInt_.bPtr_[id];
                pdd.TaPtr_[countPD] = cuChemInt.pressDepInt_.TaPtr_[id];
                pdd.fallOffTypePtr_[countPD] = cuChemInt.pressDepInt_.fallOffTypePtr_[id];
                for (integer is = 0; is < 5; ++is) {
                    pdd.fallOffCoeffPtr_[is * pdd.nrcPD_ + countPD] =
                        cuChemInt.pressDepInt_.fallOffCoeffPtr_[is * react.size() + id];
                }
                countPD++;
            }
        }
        cuChemInt.pressDepInt_.transfer(pdd);
    }
}

void OpenHurricane::mixture::cuAddThermoToCUDA() {
    const integer nsp = species_.size();
    auto &cuChemInt = *cuChemInterfacePtr_;
    cuChemInt.sptInt_.alloc(nsp);
    const auto &tht = *thermoPtr_;
    auto &sptInt = cuChemInt.sptInt_;
    for (integer i = 0; i < tht.thTable().size(); ++i) {
        sptInt.WPtr_[i] = species_[i].W();

        const auto &ji = dynamic_cast<const JANAF &>(tht[i]);
        sptInt.TCommonPtr_[i] = ji.Tcommon();
        sptInt.THighPtr_[i] = ji.Thigh();
        sptInt.TLowPtr_[i] = ji.Tlow();

        for (integer j = 0; j < 7; ++j) {
            sptInt.lowCpCoeffsPtr_[j * nsp + i] = ji.lowCpCoeffs()[j] / species_[i].Ri();
            sptInt.highCpCoeffsPtr_[j * nsp + i] = ji.highCpCoeffs()[j] / species_[i].Ri();
        }
    }
}

void OpenHurricane::mixture::makeTransportTableCUDA() const {
    if (!isSingular()) {
        if (!transportCUDAPtr_) {
            /*CUDAReactions::speciesTable* speciesCUDAPtr = new CUDAReactions::speciesTable(species().size(), reactionCUDA().sptInt_);*/
            cu_real *ekb = new cu_real[species().size()]();
            cu_real *sigma = new cu_real[species().size()]();

            for (integer i = 0; i < species().size(); ++i) {
                const auto &kty =
                    static_cast<const kineticTheory &>((*transportPtr_).tranTable()[i]);
                ekb[i] = kty.ekb();
                sigma[i] = kty.sigma();
            }
            transportCUDAPtr_.reset(
                new CUDATransport::transportListCUDA(species().size(), ekb, sigma, speciesCUDA()));

            /*speciesCUDAPtr->destroy();
            HurDelete(speciesCUDAPtr);*/
        }
    } else {
        LFatal("Attempt to make transport table in platform for singular species");
    }
}

void OpenHurricane::mixture::makeTransportTableCUDA(const cudaStreams &streams) const {
    if (!isSingular()) {
        if (!transportCUDAPtr_) {
            /*CUDAReactions::speciesTable* speciesCUDAPtr = new CUDAReactions::speciesTable(species().size(), reactionCUDA().sptInt_);*/
            cu_real *ekb = new cu_real[species().size()]();
            cu_real *sigma = new cu_real[species().size()]();

            for (integer i = 0; i < species().size(); ++i) {
                const auto &kty =
                    static_cast<const kineticTheory &>((*transportPtr_).tranTable()[i]);
                ekb[i] = kty.ekb();
                sigma[i] = kty.sigma();
            }
            transportCUDAPtr_.reset(new CUDATransport::transportListCUDA(
                species().size(), ekb, sigma, speciesCUDA(streams), streams));

            /*speciesCUDAPtr->destroy();
            HurDelete(speciesCUDAPtr);*/
        }
    } else {
        LFatal("Attempt to make transport table in platform for singular species");
    }
}

void OpenHurricane::mixture::makeSpeciesCUDA() const {
    if (!isSingular()) {
        if (!speciesCUDAPtr_) {
            speciesCUDAPtr_.reset(
                new cuChem::speciesTable(species().size(), reactionCUDA().sptInt_));
        }
    } else {
        LFatal("Attempt to make species list in platform for singular species");
    }
}

void OpenHurricane::mixture::makeSpeciesCUDA(const cudaStreams &streams) const {
    if (!isSingular()) {
        if (!speciesCUDAPtr_) {
            speciesCUDAPtr_.reset(
                new cuChem::speciesTable(species().size(), reactionCUDA().sptInt_, streams));
        }
    } else {
        LFatal("Attempt to make species list in platform for singular species");
    }
}

#endif // CUDA_PARALLEL

OpenHurricane::mixture::mixture(const runtimeMesh &mesh, const integer speciesSize,
                                const bool noReaction)
    : basicMixture(mesh, speciesSize), thermoPtr_(nullptr), transportPtr_(nullptr),
      reactionPtr_(nullptr),
#ifdef CUDA_PARALLEL
      cuChemInterfacePtr_(nullptr), transportCUDAPtr_(nullptr), speciesCUDAPtr_(nullptr),
#endif // CUDA_PARALLEL
      noReaction_(noReaction), mixtureTabulation_(mixtureTabulations::none) {
    thermoPtr_.reset(new thermoList(species_, mesh));
    transportPtr_.reset(new transportList(species_));
    if (!noReaction_) {
        reactionPtr_.reset(new reactionList(species_, *thermoPtr_));
#ifdef CUDA_PARALLEL
        cuChemInterfacePtr_.reset(new cuChem::cuChemInterface());
#endif // CUDA_PARALLEL
    }
}

OpenHurricane::mixture::mixture(const runtimeMesh &mesh, const speciesList &species,
                                const bool noReaction)
    : basicMixture(mesh, species), thermoPtr_(nullptr), transportPtr_(nullptr),
      reactionPtr_(nullptr),
#ifdef CUDA_PARALLEL
      cuChemInterfacePtr_(nullptr), transportCUDAPtr_(nullptr), speciesCUDAPtr_(nullptr),
#endif // CUDA_PARALLEL
      noReaction_(noReaction), mixtureTabulation_(mixtureTabulations::none) {
    controller cont("mix");
    controller eosCont("eos", cont);
    controller therCont("thermo", cont);
    controller tranCont("tran", cont);

    eosCont.add(equationOfState::className_, string("perfectGas"));
    therCont.add(thermo::className_, string("constCp"));
    tranCont.add(transport::className_, string("sutherlandTwo"));
    cont.add(equationOfState::className_, eosCont);
    cont.add(thermo::className_, therCont);
    cont.add(transport::className_, tranCont);
    thermoPtr_.reset(new thermoList(species, cont, mesh));
    transportPtr_.reset(new transportList(species, cont));
}

OpenHurricane::mixture::mixture(const runtimeMesh &mesh, const speciesList &species,
                                const controller &cont, const bool inviscous)
    : basicMixture(mesh, species), thermoPtr_(nullptr), transportPtr_(nullptr),
      reactionPtr_(nullptr),
#ifdef CUDA_PARALLEL
      cuChemInterfacePtr_(nullptr), transportCUDAPtr_(nullptr), speciesCUDAPtr_(nullptr),
#endif // CUDA_PARALLEL
      noReaction_(false), mixtureTabulation_(mixtureTabulations::none) {
    if (cont.found("chemical")) {
        string chemicalType = cont.findWord("chemical");
        stringToLowerCase(chemicalType);
        if (chemicalType == "mixing") {
            noReaction_ = true;
            thermoPtr_.reset(new thermoList(species_, mesh));
            if (!inviscous) {
                transportPtr_.reset(new transportList(species_));
            }
        } else if (chemicalType == "reaction") {
            noReaction_ = false;
            thermoPtr_.reset(new thermoList(species_, mesh));
            if (!inviscous) {
                transportPtr_.reset(new transportList(species_));
            }
        } else if (chemicalType == "singleSp") {
            species_.resize(1);
            species_[0].changeW(cont.findOrDefault<real>("MW", 28.97));
            isSingular_ = true;
            noReaction_ = true;
            thermoPtr_.reset(new thermoList(species_, mesh));
            thermoPtr_->thTable().append(
                thermo::creator(cont.subController("thermo"), thermoPtr_->eos(), 0));
            if (!inviscous) {
                transportPtr_.reset(new transportList(species_));
                transportPtr_->tranTable().append(
                    transport::creator(species_, 0, cont.subController("transport")));
            }
            if (cont.found(equationOfState::className_)) {
                if (cont.isController(equationOfState::className_)) {
                    thermoPtr_->setEos(cont.subController(equationOfState::className_));
                } else {
                    thermoPtr_->setEos(cont);
                }
            } else {
                controller eqoCont("eos");
                eqoCont.add(equationOfState::className_, string(perfectGas::className_));
                thermoPtr_->setEos(eqoCont);
            }
            return;
        } else {
            errorAbortStr(("Unknown choice for \"chemical\": " + chemicalType));
        }

    } else {
        species_.resize(1);
        species_[0].changeW(cont.findOrDefault<real>("MW", 28.97));
        noReaction_ = true;
        isSingular_ = true;
        thermoPtr_.reset(new thermoList(species_, mesh));
        thermoPtr_->thTable().append(thermo::creator(cont, thermoPtr_->eos(), 0));
        if (!inviscous) {
            transportPtr_.reset(new transportList(species_));
            transportPtr_->tranTable().append(
                transport::creator(species_, 0, cont.subController("transport")));
        }
        if (cont.found(equationOfState::className_)) {
            if (cont.isController(equationOfState::className_)) {
                thermoPtr_->setEos(cont.subController(equationOfState::className_));
            } else {
                thermoPtr_->setEos(cont);
            }
        } else {
            controller eqoCont("eos");
            eqoCont.add(equationOfState::className_, string(perfectGas::className_));
            thermoPtr_->setEos(eqoCont);
        }
        return;
    }
    if (!cont.found("mixtureFiles")) {
        errorAbortStr(("Cannot fing mixtureFiles in " + cont.name()));
    }
    const auto &mixCont = cont.subController("mixtureFiles");
    std::string chemStr;
    if (mixCont.found("chemFile", true)) {
        fileName chemFile = mixCont.findParameter("chemFile", true);
        if (!chemFile.isAbsolute()) {
            fileName inputPath = mixCont.findParameter("inputPath", true);
            chemFile = inputPath / chemFile;
        }
        chemStr = readFileToString(chemFile);
    }

    std::string thermoStr;

    if (mixCont.found("thermoFile")) {
        fileName thermoFile = mixCont.findParameter("thermoFile");
        if (!thermoFile.isAbsolute()) {
            fileName inputPath = mixCont.findParameter("inputPath", true);
            thermoFile = inputPath / thermoFile;
        }
        thermoStr = readFileToString(thermoFile);
    }

    std::string tansportStr;
    if (mixCont.found("tansportFile")) {
        fileName tansportFile = mixCont.findParameter("tansportFile");
        if (!tansportFile.isAbsolute()) {
            fileName inputPath = mixCont.findParameter("inputPath", true);
            tansportFile = inputPath / tansportFile;
        }
        tansportStr = readFileToString(tansportFile);
    }

    if (cont.found("mechanism", true)) {
        std::string mechStr;
        mechStr = cont.findParameter("mechanism", true);
        if (tansportStr.empty()) {
            auto pos = mechStr.find("tran");
            if (pos == std::string::npos) {
                pos = mechStr.find("TRAN");
            }
            auto endPos = mechStr.find("END", pos);
            if (endPos == std::string::npos) {
                endPos = mechStr.find("end");
            }
            if (pos == std::string::npos || endPos == std::string::npos) {
                LFatal("Thermo section missing in chemkin file.");
            }
            tansportStr = mechStr.substr(pos, endPos - pos + 3);
            mechStr.erase(pos, endPos - pos + 3);
        }
        if (chemStr.empty()) {
            if (thermoStr.empty()) {
                auto pos = mechStr.find("THER");
                if (pos == std::string::npos) {
                    pos = mechStr.find("ther");
                }
                auto endPos = mechStr.find("END", pos);
                if (endPos == std::string::npos) {
                    endPos = mechStr.find("end");
                }
                if (pos == std::string::npos || endPos == std::string::npos) {
                    LFatal("Thermo section missing in chemkin file.");
                }
                thermoStr = mechStr.substr(pos, endPos - pos + 3);
                mechStr.erase(pos, endPos - pos + 3);
            }
            chemStr = std::move(mechStr);
        }

    } else {
        std::string mechStr;
        mechStr = chemStr + "\n";
        mechStr += thermoStr + "\ntransport\n";
        mechStr += tansportStr + "\nend\n";
        const_cast<controller &>(cont).add("mechanism", mechStr);
        const_cast<controller &>(mixCont).remove("chemFile");
        const_cast<controller &>(mixCont).remove("thermoFile");
        const_cast<controller &>(mixCont).remove("tansportFile");
    }

    if (thermoStr.empty()) {
        auto pos = chemStr.find("THER");
        if (pos == std::string::npos) {
            pos = chemStr.find("ther");
        }
        auto endPos = chemStr.find("END", pos);
        if (endPos == std::string::npos) {
            endPos = chemStr.find("end");
        }
        if (pos == std::string::npos || endPos == std::string::npos) {
            LFatal("Thermo section missing in chemkin file.");
        }
        thermoStr = chemStr.substr(pos, endPos - pos + 3);
        chemStr.erase(pos, endPos - pos + 3);
    }

    if (!noReaction_) {
        if (thermoPtr_) {
            reactionPtr_.reset(new reactionList(species_, *thermoPtr_));
        }
    }
#ifdef CUDA_PARALLEL
    if (HurGPU::useGPU()) {
        cuChemInterfacePtr_.reset(new cuChem::cuChemInterface());
    }
#endif // CUDA_PARALLEL

    chemkinFileRead myParsings(cont, *reactionPtr_, species_, *thermoPtr_, *transportPtr_, chemStr,
                               thermoStr, tansportStr, noReaction_, inviscous);
    myParsings.parsing();
#ifdef CUDA_PARALLEL
    if (HurGPU::useGPU()) {
        cuAddThermoToCUDA();
        if (!noReaction_) {
            cuParsingReaction();
        }
    }
#endif // CUDA_PARALLEL

    this->setYi(mesh);

    if (cont.found(equationOfState::className_)) {
        if (cont.isController(equationOfState::className_)) {
            thermoPtr_->setEos(cont.subController(equationOfState::className_));
        } else {
            thermoPtr_->setEos(cont);
        }
    } else {
        controller eqoCont("eos");
        eqoCont.add(equationOfState::className_, string(perfectGas::className_));
        thermoPtr_->setEos(eqoCont);
    }

    if (!inviscous) {
        /*if (cont.subController("transport").found("thermalDiffusion"))
        {
            transportPtr_->setStdPtr(speciesThermoDiffusivity::creator(cont.subController("transport").subController("thermalDiffusion"), species_));
        }*/
    }
}

OpenHurricane::mixture::mixture(const runtimeMesh &mesh, const controller &cont,
                                const bool inviscous)
    : basicMixture(mesh, 0), thermoPtr_(nullptr), transportPtr_(nullptr), reactionPtr_(nullptr),
#ifdef CUDA_PARALLEL
      cuChemInterfacePtr_(nullptr), transportCUDAPtr_(nullptr), speciesCUDAPtr_(nullptr),
#endif // CUDA_PARALLEL
      noReaction_(false), mixtureTabulation_(mixtureTabulations::none) {
    if (cont.found("chemical")) {
        string chemicalType = cont.findWord("chemical");
        stringToLowerCase(chemicalType);
        if (chemicalType == "mixing") {
            noReaction_ = true;
            thermoPtr_.reset(new thermoList(species_, mesh));
            if (!inviscous) {
                transportPtr_.reset(new transportList(species_));
            }
        } else if (chemicalType == "reaction") {
            noReaction_ = false;
            thermoPtr_.reset(new thermoList(species_, mesh));
            if (!inviscous) {
                transportPtr_.reset(new transportList(species_));
            }
        } else if (chemicalType == "singlesp") {
            species_.resize(1);
            isSingular_ = true;
            species_[0].changeW(cont.findOrDefault<real>("MW", 28.97));
            noReaction_ = true;
            thermoPtr_.reset(new thermoList(species_, mesh));
            thermoPtr_->thTable().append(
                thermo::creator(cont.subController("thermo"), thermoPtr_->eos(), 0));
            if (!inviscous) {
                transportPtr_.reset(new transportList(species_));
                transportPtr_->tranTable().append(
                    transport::creator(species_, 0, cont.subController("transport")));
            }
            if (cont.found(equationOfState::className_)) {
                if (cont.isController(equationOfState::className_)) {
                    thermoPtr_->setEos(cont.subController(equationOfState::className_));
                } else {
                    thermoPtr_->setEos(cont);
                }
            } else {
                controller eqoCont("eos");
                eqoCont.add(equationOfState::className_, string(perfectGas::className_));
                thermoPtr_->setEos(eqoCont);
            }
            return;
        } else {
            errorAbortStr(("Unknown choice for \"chemical\": " + chemicalType));
        }
        if (cont.found(equationOfState::className_)) {
            if (cont.isController(equationOfState::className_)) {
                thermoPtr_->setEos(cont.subController(equationOfState::className_));
            } else {
                thermoPtr_->setEos(cont);
            }
        } else {
            controller eqoCont("eos");
            eqoCont.add(equationOfState::className_, string(perfectGas::className_));
            thermoPtr_->setEos(eqoCont);
        }
    } else {
        species_.resize(1);
        species_[0].changeW(cont.findOrDefault<real>("MW", 28.97));
        noReaction_ = true;
        isSingular_ = true;
        thermoPtr_.reset(new thermoList(species_, mesh));
        thermoPtr_->thTable().append(
            thermo::creator(cont.subController("thermo"), thermoPtr_->eos(), 0));
        if (!inviscous) {
            transportPtr_.reset(new transportList(species_));
            transportPtr_->tranTable().append(
                transport::creator(species_, 0, cont.subController("transport")));
        }
        if (cont.found(equationOfState::className_)) {
            if (cont.isController(equationOfState::className_)) {
                thermoPtr_->setEos(cont.subController(equationOfState::className_));
            } else {
                thermoPtr_->setEos(cont);
            }
        } else {
            controller eqoCont("eos");
            eqoCont.add(equationOfState::className_, string(perfectGas::className_));
            thermoPtr_->setEos(eqoCont);
        }
        return;
    }
    if (!cont.found("mixtureFiles")) {
        errorAbortStr(("Cannot fing mixtureFiles in " + cont.name()));
    }
    const auto &mixCont = cont.subController("mixtureFiles");
    std::string chemStr;
    if (mixCont.found("chemFile")) {
        fileName chemFile = mixCont.findWord("chemFile");
        if (!chemFile.isAbsolute()) {
            fileName inputPath = mixCont.findParameter("inputPath", true);
            chemFile = inputPath / chemFile;
        }
        chemStr = readFileToString(chemFile);
    }

    std::string thermoStr;
    if (mixCont.found("thermoFile")) {
        fileName thermoFile = mixCont.findWord("thermoFile");
        if (!thermoFile.isAbsolute()) {
            fileName inputPath = mixCont.findParameter("inputPath", true);
            thermoFile = inputPath / thermoFile;
        }
        thermoStr = readFileToString(thermoFile);
    }

    std::string tansportStr;
    if (mixCont.found("tansportFile")) {
        fileName tansportFile = mixCont.findWord("tansportFile");
        if (!tansportFile.isAbsolute()) {
            fileName inputPath = mixCont.findParameter("inputPath", true);
            tansportFile = inputPath / tansportFile;
        }
        tansportStr = readFileToString(tansportFile);
    }

    if (cont.found("mechanism", true)) {
        std::string mechStr;
        mechStr = cont.findParameter("mechanism", true);
        if (tansportStr.empty()) {
            auto pos = mechStr.find("tran");
            if (pos == std::string::npos) {
                pos = mechStr.find("TRAN");
            }
            auto endPos = mechStr.find("END", pos);
            if (endPos == std::string::npos) {
                endPos = mechStr.find("end");
            }
            if (pos == std::string::npos || endPos == std::string::npos) {
                LFatal("Thermo section missing in chemkin file.");
            }
            tansportStr = mechStr.substr(pos, endPos - pos + 3);
            mechStr.erase(pos, endPos - pos + 3);
        }
        if (chemStr.empty()) {
            if (thermoStr.empty()) {
                auto pos = mechStr.find("THER");
                if (pos == std::string::npos) {
                    pos = mechStr.find("ther");
                }
                auto endPos = mechStr.find("END", pos);
                if (endPos == std::string::npos) {
                    endPos = mechStr.find("end");
                }
                if (pos == std::string::npos || endPos == std::string::npos) {
                    LFatal("Thermo section missing in chemkin file.");
                }
                thermoStr = mechStr.substr(pos, endPos - pos + 3);
                mechStr.erase(pos, endPos - pos + 3);
            }
            chemStr = std::move(mechStr);
        }

    } else {
        std::string mechStr;
        mechStr = chemStr + "\n";
        mechStr += thermoStr + "\ntransport\n";
        mechStr += tansportStr + "\nend\n";
        const_cast<controller &>(cont).add("mechanism", mechStr);
        const_cast<controller &>(mixCont).remove("chemFile");
        const_cast<controller &>(mixCont).remove("thermoFile");
        const_cast<controller &>(mixCont).remove("tansportFile");
    }

    if (thermoStr.empty()) {
        auto pos = chemStr.find("THER");
        if (pos == std::string::npos) {
            pos = chemStr.find("ther");
        }
        auto endPos = chemStr.find("END", pos);
        if (endPos == std::string::npos) {
            endPos = chemStr.find("end");
        }
        if (pos == std::string::npos || endPos == std::string::npos) {
            LFatal("Thermo section missing in chemkin file.");
        }
        thermoStr = chemStr.substr(pos, endPos - pos + 3);
        chemStr.erase(pos, endPos - pos + 3);
    }

    if (!noReaction_) {
        if (thermoPtr_) {
            reactionPtr_.reset(new reactionList(species_, *thermoPtr_));
        }
    }
#ifdef CUDA_PARALLEL
    if (HurGPU::useGPU()) {
        cuChemInterfacePtr_.reset(new cuChem::cuChemInterface());
    }
#endif // CUDA_PARALLEL

    chemkinFileRead myParsings(cont, *reactionPtr_, species_, *thermoPtr_, *transportPtr_, chemStr,
                               thermoStr, tansportStr, noReaction_, inviscous);

    myParsings.parsing();
#ifdef CUDA_PARALLEL
    if (HurGPU::useGPU()) {
        cuAddThermoToCUDA();
        if (!noReaction_) {
            cuParsingReaction();
        }
    }
#endif // CUDA_PARALLEL

    this->setYi(mesh);

    if (!inviscous) {
        /*if (cont.subController("transport").found("thermalDiffusion"))
        {
            transportPtr_->setStdPtr(speciesThermoDiffusivity::creator(cont.subController("transport").subController("thermalDiffusion"), species_));
        }*/
    }
}

OpenHurricane::mixture::~mixture() noexcept {
    thermoPtr_.clear();
    transportPtr_.clear();
    reactionPtr_.clear();

#ifdef CUDA_PARALLEL
    if (transportCUDAPtr_) {
        transportCUDAPtr_->destroy();
    }
    transportCUDAPtr_.clear();

    if (speciesCUDAPtr_) {
        speciesCUDAPtr_->destroy();
    }
    speciesCUDAPtr_.clear();

#endif // CUDA_PARALLEL
}

hur_nodiscard OpenHurricane::cellRealArray OpenHurricane::mixture::W() const {
    const PtrList<cellRealArray> &Y = basicMixture::Yi();

    cellRealArray tempW(object("Wm", Y[0].mesh(), object::NOT_WRITE), Y[0].mesh(), real(0.0));

    for (integer i = 0; i < Y.size(); ++i) {
        tempW += Y[i] / Wi(i);
    }
    tempW = 1.0 / tempW;
    return tempW;
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::W(const realArray &_yi) const {
    real tempW = 0.0;
    for (integer i = 0; i < _yi.size(); ++i) {
        tempW += _yi[i] / Wi(i);
    }
    tempW = 1.0 / tempW;
    return tempW;
}

hur_nodiscard OpenHurricane::cellRealArray
OpenHurricane::mixture::rhoi(const cellRealArray &p, const cellRealArray &T, const integer i,
                             const bool onlyInternal) const {
    std::string str;
    str = "rho_";
    str += toString(i);
    cellRealArray rhom(object(str.c_str(), p.mesh(), object::NOT_WRITE), p.mesh());
    integer nSize;
    if (onlyInternal) {
        nSize = p.mesh().nCells();
    } else {
        nSize = p.mesh().nTotalCells();
    }
    for (integer n = 0; n < nSize; ++n) {
        rhom[n] = thermoPtr_->eos().rhoi(p[n], T[n], i);
    }

    return rhom;
}

hur_nodiscard OpenHurricane::cellRealArray
OpenHurricane::mixture::rho(const cellRealArray &p, const cellRealArray &T,
                            const bool onlyInternal) const {
    cellRealArray rhom(object("rhom", p.mesh(), object::NOT_WRITE), p.mesh());
    integer nSize;
    if (onlyInternal) {
        nSize = p.mesh().nCells();
    } else {
        nSize = p.mesh().nTotalCells();
    }
    for (integer n = 0; n < nSize; ++n) {
        rhom[n] = thermoPtr_->eos().rhom(p[n], T[n], Yic(n));
    }
    return rhom;
}

hur_nodiscard OpenHurricane::cellRealArray
OpenHurricane::mixture::rho(const cellRealArray &p, const cellRealArray &T,
                            const PtrList<cellRealArray> &_Yi, const bool onlyInternal) const {
    cellRealArray rhom(object("rhom", p.mesh(), object::NOT_WRITE), p.mesh());
    integer nSize;
    if (onlyInternal) {
        nSize = p.mesh().nCells();
    } else {
        nSize = p.mesh().nTotalCells();
    }
    for (integer n = 0; n < nSize; ++n) {
        realArray Yii(species_.size(), Zero);

        for (integer i = 0; i < species_.size(); ++i) {
            Yii[i] = _Yi[i][n];
        }
        rhom[n] = thermoPtr_->eos().rhom(p[n], T[n], Yii);
    }
    return rhom;
}

void OpenHurricane::mixture::rho(cellRealArray &rhom, const cellRealArray &p,
                                 const cellRealArray &T, const PtrList<cellRealArray> &_Yi,
                                 const bool onlyInternal) const {
    integer nSize;
    if (onlyInternal) {
        nSize = p.mesh().nCells();
    } else {
        nSize = p.mesh().nTotalCells();
    }
    for (integer n = 0; n < nSize; ++n) {
        realArray Yii(species_.size(), Zero);

        for (integer i = 0; i < species_.size(); ++i) {
            Yii[i] = _Yi[i][n];
        }
        rhom[n] = thermoPtr_->eos().rhom(p[n], T[n], Yii);
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::rhoi(const real pi, const real Ti,
                                                               const integer i) const {
    return thermoPtr_->eos().rhoi(pi, Ti, i);
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::rho(const real p, const real T,
                                                              const realArray &yi) const {
    return thermoPtr_->eos().rhom(p, T, yi);
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::rho(const real p, const real T,
                                                              const integer cellI) const {
    return rho(p, T, Yic(cellI));
}

hur_nodiscard OpenHurricane::cellRealArray
OpenHurricane::mixture::mui(const cellRealArray &p, const cellRealArray &T, const integer i,
                            const bool onlyInternal) const {

    std::string str;
    str = "mu_";
    str += toString(i);
    cellRealArray mum(object(str.c_str(), p.mesh(), object::NOT_WRITE), p.mesh());
    integer nSize;
    if (onlyInternal) {
        nSize = p.mesh().nCells();
    } else {
        nSize = p.mesh().nTotalCells();
    }
    for (integer n = 0; n < nSize; ++n) {
        mum[n] = transportPtr_->mu(p[n], T[n], i);
    }

    return mum;
}

OpenHurricane::cellRealArray OpenHurricane::mixture::mu(const cellRealArray &p,
                                                        const cellRealArray &T,
                                                        const bool onlyInternal) const {
    cellRealArray mum(object("mum", p.mesh(), object::NOT_WRITE), p.mesh());
    integer nSize;
    if (onlyInternal) {
        nSize = p.mesh().nCells();
    } else {
        nSize = p.mesh().nTotalCells();
    }
    for (integer n = 0; n < nSize; ++n) {
        mum[n] = transportPtr_->mu(p[n], T[n], Yic(n));
    }

    return mum;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::mixture::mu(const cellRealArray &p,
                                                                  const cellRealArray &T,
                                                                  const integer faceZoneId) const {
    const auto &fZL = mesh_.faceZones();
    const integer si = fZL[faceZoneId].size();
    realArray tmpG(si, Zero);
    if (isSingular()) {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        for (integer fi = 0; fi < si; ++fi) {
            tmpG[fi] = transportPtr_->mu(pf[fi], tf[fi], 0);
        }
    } else {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        const auto yif = Yif(faceZoneId);
        for (integer fi = 0; fi < si; ++fi) {
            tmpG[fi] = transportPtr_->mu(pf[fi], tf[fi], yif[fi]);
        }
    }
    return tmpG;
}

void OpenHurricane::mixture::mu(const cellRealArray &p, const cellRealArray &T, cellRealArray &mum,
                                const bool onlyInternal, const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer fieldSize = p.size();
    if (onlyInternal || ((!onlyInternal) && isExtrapolate)) {
        fieldSize = mesh().nCells();
    }

    for (integer cellI = 0; cellI < fieldSize; ++cellI) {
        mum[cellI] = transportPtr_->mu(p[cellI], T[cellI], Yi(), cellI);
        if (std::isnan(mum[cellI])) {
            if (cellI >= nCells) {
                flag = false;
                countFlag = 1;
                break;
            }
        }
    }
    if ((!onlyInternal) && !isExtrapolate) {
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }
    if (((!onlyInternal) && isExtrapolate) || !flag) {
        realTransfer myTransfer(mesh_, mum, false, true);
        myTransfer.transferInit();
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                const auto mumf = mu(p, T, fzi);
                fv::extrapolate(mumf, mum, fzi);
            }
        }
        myTransfer.transferring();
        //fv::transfer(mum);
    }
}

hur_nodiscard OpenHurricane::cellRealArray
OpenHurricane::mixture::cpi(const cellRealArray &rhoi, const cellRealArray &Ti, const integer i,
                            const bool onlyInternal) const {
    std::string str;
    str = "cp_";
    str += toString(i);
    cellRealArray cpi(object(str.c_str(), rhoi.mesh(), object::NOT_WRITE), rhoi.mesh());
    integer nSize;
    if (onlyInternal) {
        nSize = rhoi.mesh().nCells();
    } else {
        nSize = rhoi.mesh().nTotalCells();
    }
    for (integer n = 0; n < nSize; ++n) {
        cpi[n] = (*thermoPtr_)[i].cp_rho(rhoi[n], Ti[n]);
    }

    return cpi;
}

hur_nodiscard OpenHurricane::cellRealArray
OpenHurricane::mixture::cp(const cellRealArray &rho, const cellRealArray &T,
                           const bool onlyInternal) const {
    if (species_.size() == 1) {
        return cpi(rho, T, 0);
    }
    cellRealArray cpm(object("cpm", rho.mesh(), object::NOT_WRITE), rho.mesh());
    integer nSize;
    if (onlyInternal) {
        nSize = rho.mesh().nCells();
    } else {
        nSize = rho.mesh().nTotalCells();
    }
    for (integer n = 0; n < nSize; ++n) {
        //cpm[n] = thermoPtr_->cp_rho(rho[n] * rho.refValue(), T[n] * T.refValue(), Yic(n));
        cpm[n] = thermoPtr_->cp_rho(rho[n], T[n], Yi(), n);
    }

    return cpm;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::mixture::cp(const cellRealArray &p,
                                                                  const cellRealArray &T,
                                                                  const integer faceZoneId) const {
    const auto &fZL = mesh_.faceZones();
    const integer si = fZL[faceZoneId].size();
    realArray tmpG(si, Zero);
    if (isSingular()) {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        for (integer fi = 0; fi < si; ++fi) {
            tmpG[fi] = thermoPtr_->cpi_p(pf[fi], tf[fi], 0);
        }
    } else {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        const auto yif = Yif(faceZoneId);
        for (integer fi = 0; fi < si; ++fi) {
            tmpG[fi] = thermoPtr_->cp_p(pf[fi], tf[fi], yif[fi]);
        }
    }
    return tmpG;
}

void OpenHurricane::mixture::cp(const cellRealArray &p, const cellRealArray &T, cellRealArray &cpm,
                                const bool onlyInternal, const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer fieldSize = p.size();
    if (onlyInternal || ((!onlyInternal) && isExtrapolate)) {
        fieldSize = mesh().nCells();
    }
    for (integer n = 0; n < fieldSize; ++n) {
        //cpm[n] = thermoPtr_->cp_rho(rho[n] * rho.refValue(), T[n] * T.refValue(), Yic(n));
        cpm[n] = thermoPtr_->cp_p(p[n], T[n], Yi(), n);
        if (std::isnan(cpm[n]) && n >= nCells) {
            flag = false;
            countFlag = 1;
            break;
        }
    }
    if ((!onlyInternal) && !isExtrapolate) {
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }
    if (((!onlyInternal) && isExtrapolate) || !flag) {
        realTransfer myTransfer(mesh_, cpm, false, true);
        myTransfer.transferInit();
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                const auto cpf = cp(p, T, fzi);
                fv::extrapolate(cpf, cpm, fzi);
            }
        }
        myTransfer.transferring();
        //fv::transfer(cpm);
    }
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::mixture::hai(const cellRealArray &p,
                                                                   const cellRealArray &T,
                                                                   const integer faceZoneId,
                                                                   const integer isp) const {
    const auto &fZL = mesh_.faceZones();
    const integer si = fZL[faceZoneId].size();
    realArray tmpG(si, Zero);

    const auto pf = fv::interpolate(p, faceZoneId);
    const auto tf = fv::interpolate(T, faceZoneId);
    for (integer fi = 0; fi < si; ++fi) {
        tmpG[fi] = thermoPtr_->ha_p(pf[fi], tf[fi], isp);
    }

    return tmpG;
}

hur_nodiscard OpenHurricane::cellRealArray
OpenHurricane::mixture::hai(const cellRealArray &p, const cellRealArray &T, const integer i,
                            const bool onlyInternal, const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    cellRealArray haiT(
        object(species()[i].name() + "ha", p.mesh(), object::NOT_WRITE, object::TEMPORARY),
        p.mesh());
    bool flag = true;
    integer countFlag = 0;
    integer nSize;
    nSize = p.size();

    if (onlyInternal || ((!onlyInternal) && isExtrapolate)) {
        nSize = mesh().nCells();
    }

    for (integer n = 0; n < nSize; ++n) {
        haiT[n] = thermoPtr_->ha_p(p[n], T[n], i);
        if (std::isnan(haiT[n]) && n >= nCells) {
            flag = false;
            countFlag = 1;
            break;
        }
    }
    if ((!onlyInternal) && !isExtrapolate) {
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }
    if (((!onlyInternal) && isExtrapolate) || !flag) {
        realTransfer myTransfer(mesh_, haiT, false, true);
        myTransfer.transferInit();
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                const auto haif = hai(p, T, fzi, i);
                fv::extrapolate(haif, haiT, fzi);
            }
        }
        myTransfer.transferring();
        //fv::transfer(haiT);
    }
    return haiT;
}

void OpenHurricane::mixture::hai(const cellRealArray &p, const cellRealArray &T,
                                 cellRealArray &haiT, const integer i, const bool onlyInternal,
                                 const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer nSize;
    nSize = p.size();

    if (onlyInternal || ((!onlyInternal) && isExtrapolate)) {
        nSize = mesh().nCells();
    }

    for (integer n = 0; n < nSize; ++n) {
        haiT[n] = thermoPtr_->ha_p(p[n], T[n], i);
        if (std::isnan(haiT[n]) && n >= nCells) {
            flag = false;
            countFlag = 1;
            break;
        }
    }
    if ((!onlyInternal) && !isExtrapolate) {
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }
    if (((!onlyInternal) && isExtrapolate) || !flag) {
        realTransfer myTransfer(mesh_, haiT, false, true);
        myTransfer.transferInit();
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                const auto haif = hai(p, T, fzi, i);
                fv::extrapolate(haif, haiT, fzi);
            }
        }
        myTransfer.transferring();
        //fv::transfer(haiT);
    }
}

void OpenHurricane::mixture::updateHai(const cellRealArray &p, const cellRealArray &T,
                                       const bool onlyInternal, const bool isExtrapolate) {
    for (integer i = 0; i < species().size(); ++i) {
        hai(p, T, hi(i), i, onlyInternal, isExtrapolate);
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::gamma(const real p, const real T,
                                                                const realArray &yi) const {
    return thermoPtr_->gamma(p, T, yi);
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::gamma(const real p, const real T,
                                                                const integer cellI) const {
    if (isSingular()) {
        return gammai(p, T, 0);
    }
    return thermoPtr_->gamma(p, T, Yi(), cellI);
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::mixture::gamma(const cellRealArray &p, const cellRealArray &T,
                              const integer faceZoneId) const {
    const auto &fZL = mesh_.faceZones();
    const integer si = fZL[faceZoneId].size();
    realArray tmpG(si, Zero);
    if (isSingular()) {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        for (integer fi = 0; fi < si; ++fi) {
            tmpG[fi] = gammai(pf[fi], tf[fi], 0);
        }
    } else {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        const auto yif = Yif(faceZoneId);
        for (integer fi = 0; fi < si; ++fi) {
            tmpG[fi] = gamma(pf[fi], tf[fi], yif[fi]);
        }
    }
    return tmpG;
}

void OpenHurricane::mixture::gamma(const cellRealArray &p, const cellRealArray &T,
                                   cellRealArray &gm, const bool onlyInternal,
                                   const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer fieldSize = p.size();
    if (onlyInternal || ((!onlyInternal) && isExtrapolate)) {
        fieldSize = mesh().nCells();
    }

    for (integer cellI = 0; cellI < fieldSize; ++cellI) {
        //gm[cellI] = gamma(p[cellI] * p.refValue(), T[cellI] * T.refValue(), Yic(cellI));
        gm[cellI] = thermoPtr_->gamma(p[cellI], T[cellI], Yi(), cellI);
        if (std::isnan(gm[cellI]) && cellI >= nCells) {
            flag = false;
            countFlag = 1;
            break;
        }
    }
    if ((!onlyInternal) && !isExtrapolate) {
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }

    if (((!onlyInternal) && isExtrapolate) || !flag) {
        realTransfer myTransfer(mesh_, gm, false, false);
        myTransfer.transferInit(0);
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                const auto gamaf = gamma(p, T, fzi);
                fv::extrapolate(gamaf, gm, fzi);
            }
        }
        myTransfer.transferring(0);
        for (integer layerI = 1; layerI < mesh_.ghostCellsType(); layerI++) {
            myTransfer.transferInit(layerI);
            myTransfer.transferring(layerI);
        }
        //fv::transfer(gm);
    }
}

void OpenHurricane::mixture::gammaBoundary(const cellRealArray &p, const cellRealArray &T,
                                           cellRealArray &gm, const bool isExtrapolate) const {
    bool flag = true;
    integer countFlag = 0;
    const integer nCells = mesh().nCells();

    const integer fieldSize = mesh().nTotalCells();

    if (!isExtrapolate) {
        for (integer cellI = nCells; cellI < fieldSize; ++cellI) {
            //gm[cellI] = gamma(p[cellI] * p.refValue(), T[cellI] * T.refValue(), Yic(cellI));
            gm[cellI] = thermoPtr_->gamma(p[cellI], T[cellI], Yi(), cellI);
            if (std::isnan(gm[cellI])) {
                flag = false;
                countFlag = 1;
                break;
            }
        }
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }

    if (isExtrapolate || !flag) {
        realTransfer myTransfer(mesh_, gm, false, false);
        myTransfer.transferInit(0);
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                const auto gamaf = gamma(p, T, fzi);
                fv::extrapolate(gamaf, gm, fzi);
            }
        }
        myTransfer.transferring(0);
        for (integer layerI = 1; layerI < mesh_.ghostCellsType(); layerI++) {
            myTransfer.transferInit(layerI);
            myTransfer.transferring(layerI);
        }
        //fv::transfer(gm);
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::gammai(const real p, const real T,
                                                                 const integer i) const {
    return thermoPtr_->gammai(p, T, i);
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::gammaTh(const real p, const real rho,
                                                                  const real T,
                                                                  const realArray &yi) const {
    return thermoPtr_->gammaTh(p, rho, T, yi);
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::gammaTh(const real p, const real rho,
                                                                  const real T,
                                                                  const integer cellI) const {
    if (isSingular()) {
        return gammaThi(p, rho, T, 0);
    }
    return gammaTh(p, rho, T, Yic(cellI));
}

hur_nodiscard OpenHurricane::real OpenHurricane::mixture::gammaThi(const real p, const real rho,
                                                                   const real T,
                                                                   const integer i) const {
    return thermoPtr_->gammaThi(p, rho, T, i);
}

hur_nodiscard OpenHurricane::PtrList<OpenHurricane::cellRealArray>
OpenHurricane::mixture::Diff(const cellRealArray &p, const cellRealArray &T,
                             const bool onlyInternal) const {
    PtrList<cellRealArray> Dm(species().size());
    for (integer i = 0; i < species().size(); ++i) {
        Dm.set(i, new cellRealArray(object(species_[i].name() + "Dm-temp", mesh(),
                                           object::NOT_WRITE, object::NOT_PRIMITIVE),
                                    mesh()));
    }

    integer fieldSize = p.size();
    if (onlyInternal) {
        fieldSize = mesh().nCells();
    }
    for (integer cellI = 0; cellI < fieldSize; ++cellI) {
        realArray Df = Diffi(p, T, Yi(), cellI);
        for (integer i = 0; i < species().size(); ++i) {
            Dm[i][cellI] = Df[i];
            Dm[i][cellI] = max(real(0.0), Dm[i][cellI]);
        }
    }

    return Dm;
}

hur_nodiscard OpenHurricane::realArrayArray
OpenHurricane::mixture::Diff(const cellRealArray &p, const cellRealArray &T,
                             const integer faceZoneId) const {
    const auto &fZL = mesh_.faceZones();
    const integer si = fZL[faceZoneId].size();
    realArrayArray tmpD(species_.size());
    realArray tmpDf(species_.size());
    const auto pf = fv::interpolate(p, faceZoneId);
    const auto tf = fv::interpolate(T, faceZoneId);
    const auto yif = Yif(faceZoneId);
    for (integer i = 0; i < species_.size(); ++i) {
        tmpD[i].resize(si);
    }
    for (integer fi = 0; fi < si; ++fi) {
        Diffi(tmpDf, pf[fi], tf[fi], yif[fi]);
        for (integer isp = 0; isp < species_.size(); ++isp) {
            tmpD[isp][fi] = tmpDf[isp];
        }
        /*auto xi = species_.Yi2Xi(yif[fi]);
        const real Tv1p5 = tf[fi] * sqrt(tf[fi]);
        for (integer isp = 0; isp < species_.size(); ++isp)
        {
            tmpD[isp][fi] = transportPtr_->D_Xi(pf[fi], tf[fi], Tv1p5, xi, isp, tmpDij);
        }*/
    }
    return tmpD;
}

void OpenHurricane::mixture::Diff(const cellRealArray &p, const cellRealArray &T,
                                  PtrList<cellRealArray> &Di, const bool onlyInternal,
                                  const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer fieldSize = p.size();
    if (onlyInternal || ((!onlyInternal) && isExtrapolate)) {
        fieldSize = mesh().nCells();
    }
    realArray Df(species_.size());
    //realArray Dij(species_.size());
    for (integer cellI = 0; cellI < fieldSize; ++cellI) {
        Diffi(Df, p, T, Yi(), cellI);
        //transportPtr_->Df
        //(
        //    Df,
        //    p[cellI],
        //    T[cellI],
        //    Yi(),
        //    cellI//, Dij
        //);
        for (integer i = 0; i < species().size(); ++i) {
            Di[i][cellI] = Df[i];
            if (std::isnan(Df[i]) && cellI >= nCells) {
                flag = false;
                countFlag = 1;
                break;
            }
        }
    }
    if ((!onlyInternal) && !isExtrapolate) {
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }
    if (((!onlyInternal) && isExtrapolate) || !flag) {
        PtrList<realTransfer> myTransfer(species_.size());
        for (integer isp = 0; isp < species_.size(); ++isp) {
            myTransfer.set(isp, new realTransfer(mesh_, Di[isp], false, true));
            myTransfer[isp].transferInit();
        }
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                const auto dif = Diff(p, T, fzi);
                for (integer isp = 0; isp < species_.size(); ++isp) {
                    fv::extrapolate(dif[isp], Di[isp], fzi);
                }
            }
        }
        for (integer isp = 0; isp < species_.size(); ++isp) {
            //fv::transfer(Di[isp]);
            myTransfer[isp].transferring();
        }
    }
}

void OpenHurricane::mixture::DiffBoundary(const cellRealArray &p, const cellRealArray &T,
                                          PtrList<cellRealArray> &Di,
                                          const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer fieldSize = p.size();

    if (!isExtrapolate) {
        realArray Df(species_.size());
        realArray Dij(species_.size());
        for (integer cellI = nCells; cellI < fieldSize; ++cellI) {
            Diffi(Df, p, T, Yi(), cellI);

            /*transportPtr_->Df
            (
                Df,
                p[cellI],
                T[cellI],
                Yi(),
                cellI,
                Dij
            );*/
            for (integer i = 0; i < species().size(); ++i) {
                Di[i][cellI] = Df[i];
                if (std::isnan(Df[i]) && cellI >= nCells) {
                    flag = false;
                    countFlag = 1;
                    break;
                }
            }
        }

        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }

    if (isExtrapolate || !flag) {
        PtrList<realTransfer> myTransfer(species_.size());
        for (integer isp = 0; isp < species_.size(); ++isp) {
            myTransfer.set(isp, new realTransfer(mesh_, Di[isp], false, true));
            myTransfer[isp].transferInit();
        }
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                const auto dif = Diff(p, T, fzi);
                for (integer isp = 0; isp < species_.size(); ++isp) {
                    fv::extrapolate(dif[isp], Di[isp], fzi);
                }
            }
        }
        for (integer isp = 0; isp < species_.size(); ++isp) {
            //fv::transfer(Di[isp]);
            myTransfer[isp].transferring();
        }
    }
}

void OpenHurricane::mixture::Diffi(realArray &Df, const cellRealArray &p, const cellRealArray &T,
                                   const PtrList<cellRealArray> &yi, const integer cellI) const {
    transportPtr_->Df(Df, p[cellI], T[cellI], yi, cellI);
}

void OpenHurricane::mixture::Diffi(realArray &Df, const real p, const real T,
                                   const realArray &yif) const {
    auto xi = species_.Yi2Xi(yif);
    const real Tv1p5 = T * sqrt(T);
    for (integer isp = 0; isp < species_.size(); ++isp) {
        Df[isp] = transportPtr_->D_Xi(p, T, Tv1p5, xi, isp);
    }
}

void OpenHurricane::mixture::muKappa(const cellRealArray &rho, const cellRealArray &p,
                                     const cellRealArray &T, cellRealArray &mum,
                                     cellRealArray &kappam, const bool onlyInternal,
                                     const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer fieldSize = p.size();
    if (onlyInternal || ((!onlyInternal) && isExtrapolate)) {
        fieldSize = mesh().nCells();
    }
    for (integer cellI = 0; cellI < fieldSize; ++cellI) {
        realArray cpi(species().size());
        for (integer i = 0; i < species().size(); ++i) {
            cpi[i] = (*thermoPtr_)[i].cp_rho(rho[cellI], T[cellI]);
        }

        transportPtr_->muKappa(p[cellI], T[cellI], cpi, Yi(), cellI, mum[cellI], kappam[cellI]);

        if ((std::isnan(mum[cellI]) || std::isnan(kappam[cellI])) && cellI >= nCells) {
            flag = false;
            countFlag = 1;
            break;
        }
    }
    if ((!onlyInternal) && !isExtrapolate) {
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }

    if (((!onlyInternal) && isExtrapolate) || !flag) {
        realTransfer myTransfer1(mesh_, mum, false, true);
        realTransfer myTransfer2(mesh_, kappam, false, true);
        myTransfer1.transferInit();
        myTransfer2.transferInit();
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                realArray mumf(fZL[fzi].size(), Zero);
                realArray kappamf(fZL[fzi].size(), Zero);
                muKappa(rho, p, T, mumf, kappamf, fzi);
                fv::extrapolate(mumf, mum, fzi);
                fv::extrapolate(kappamf, kappam, fzi);
            }
        }
        //fv::transfer(mum);
        //fv::transfer(kappam);
        myTransfer1.transferring();
        myTransfer2.transferring();
    }
}

void OpenHurricane::mixture::muKappaBoundary(const cellRealArray &rho, const cellRealArray &p,
                                             const cellRealArray &T, cellRealArray &mum,
                                             cellRealArray &kappam,
                                             const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer fieldSize = p.size();
    if (!isExtrapolate) {
        for (integer cellI = nCells; cellI < fieldSize; ++cellI) {
            realArray cpi(species().size());
            for (integer i = 0; i < species().size(); ++i) {
                cpi[i] = (*thermoPtr_)[i].cp_rho(rho[cellI], T[cellI]);
            }

            transportPtr_->muKappa(p[cellI], T[cellI], cpi, Yi(), cellI, mum[cellI], kappam[cellI]);

            if ((std::isnan(mum[cellI]) || std::isnan(kappam[cellI])) && cellI >= nCells) {
                flag = false;
                countFlag = 1;
                break;
            }
        }
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }

    if (isExtrapolate || !flag) {
        realTransfer myTransfer1(mesh_, mum, false, true);
        realTransfer myTransfer2(mesh_, kappam, false, true);
        myTransfer1.transferInit();
        myTransfer2.transferInit();
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                realArray mumf(fZL[fzi].size(), Zero);
                realArray kappamf(fZL[fzi].size(), Zero);
                muKappa(rho, p, T, mumf, kappamf, fzi);
                fv::extrapolate(mumf, mum, fzi);
                fv::extrapolate(kappamf, kappam, fzi);
            }
        }
        //fv::transfer(mum);
        //fv::transfer(kappam);
        myTransfer1.transferring();
        myTransfer2.transferring();
    }
}

void OpenHurricane::mixture::muKappa(const cellRealArray &rho, const cellRealArray &p,
                                     const cellRealArray &T, realArray &mum, realArray &kappam,
                                     const integer faceZoneId) const {
    const auto &fZL = mesh_.faceZones();
    const integer si = fZL[faceZoneId].size();
    if (isSingular()) {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        realArray yif(1, real(1.0));
        realArray cpi(1);
        for (integer fi = 0; fi < si; ++fi) {
            cpi[0] = (*thermoPtr_)[0].cp_p(pf[fi], tf[fi]);
            transportPtr_->muKappa(pf[fi], tf[fi], cpi, yif, mum[fi], kappam[fi]);
        }
    } else {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        const auto yif = Yif(faceZoneId);
        realArray cpi(species().size());
        for (integer fi = 0; fi < si; ++fi) {
            for (integer i = 0; i < species().size(); ++i) {
                cpi[i] = (*thermoPtr_)[i].cp_p(pf[fi], tf[fi]);
            }
            transportPtr_->muKappa(pf[fi], tf[fi], cpi, yif[fi], mum[fi], kappam[fi]);
        }
    }
}

void OpenHurricane::mixture::muKappaCp(const cellRealArray &rho, const cellRealArray &p,
                                       const cellRealArray &T, cellRealArray &mum,
                                       cellRealArray &kappam, cellRealArray &cpm,
                                       const bool onlyInternal, const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer fieldSize = p.size();
    if (onlyInternal || ((!onlyInternal) && isExtrapolate)) {
        fieldSize = mesh().nCells();
    }
    realArray cpi(species().size());
    for (integer cellI = 0; cellI < fieldSize; ++cellI) {
        cpm[cellI] = 0.0;
        if (isSingular_) {
            cpi[0] = (*thermoPtr_)[0].cp_rho(rho[cellI], T[cellI]);
            cpm[cellI] = cpi[0];
        } else {
            for (integer i = 0; i < species().size(); ++i) {
                cpi[i] = (*thermoPtr_)[i].cp_rho(rho[cellI], T[cellI]);
                cpm[cellI] += cpi[i] * Yi()[i][cellI];
            }
        }
        transportPtr_->muKappa(p[cellI], T[cellI], cpi, Yi(), cellI, mum[cellI], kappam[cellI]);

        if ((std::isnan(mum[cellI]) || std::isnan(kappam[cellI])) && cellI >= nCells) {
            flag = false;
            countFlag = 1;
            break;
        }
    }
    if ((!onlyInternal) && !isExtrapolate) {
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }

    if (((!onlyInternal) && isExtrapolate) || !flag) {
        realTransfer myTransfer1(mesh_, mum, false, true);
        realTransfer myTransfer2(mesh_, kappam, false, true);
        realTransfer myTransfer3(mesh_, cpm, false, true);
        myTransfer1.transferInit();
        myTransfer2.transferInit();
        myTransfer3.transferInit();
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                realArray mumf(fZL[fzi].size(), Zero);
                realArray kappamf(fZL[fzi].size(), Zero);
                realArray cpmf(fZL[fzi].size(), Zero);
                muKappaCp(rho, p, T, mumf, kappamf, cpmf, fzi);
                fv::extrapolate(mumf, mum, fzi);
                fv::extrapolate(kappamf, kappam, fzi);
                fv::extrapolate(cpmf, cpm, fzi);
            }
        }

        //fv::transfer(mum);
        //fv::transfer(kappam);
        //fv::transfer(cpm);
        myTransfer1.transferring();
        myTransfer2.transferring();
        myTransfer3.transferring();
    }
}

void OpenHurricane::mixture::muKappaCp(const cellRealArray &rho, const cellRealArray &p,
                                       const cellRealArray &T, realArray &mum, realArray &kappam,
                                       realArray &cpm, const integer faceZoneId) const {
    const auto &fZL = mesh_.faceZones();
    const integer si = fZL[faceZoneId].size();
    if (isSingular()) {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        realArray yif(1, real(1.0));
        realArray cpi(1);
        for (integer fi = 0; fi < si; ++fi) {
            cpi[0] = (*thermoPtr_)[0].cp_p(pf[fi], tf[fi]);
            transportPtr_->muKappa(pf[fi], tf[fi], cpi, yif, mum[fi], kappam[fi]);
            cpm[fi] = cpi[0];
        }
    } else {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        const auto yif = Yif(faceZoneId);
        realArray cpi(species().size());
        for (integer fi = 0; fi < si; ++fi) {
            cpm[fi] = 0.0;
            for (integer i = 0; i < species().size(); ++i) {
                cpi[i] = (*thermoPtr_)[i].cp_p(pf[fi], tf[fi]);
                cpm[fi] += cpi[i] * yif[fi][i];
            }
            transportPtr_->muKappa(pf[fi], tf[fi], cpi, yif[fi], mum[fi], kappam[fi]);
        }
    }
}

void OpenHurricane::mixture::gasProperties(const cellRealArray &rho, const cellRealArray &p,
                                           const cellRealArray &T, cellRealArray &mum,
                                           cellRealArray &kappam, cellRealArray &cpm,
                                           PtrList<cellRealArray> &Di, PtrList<cellRealArray> &hhai,
                                           const bool onlyInternal,
                                           const bool isExtrapolate) const {
    muKappaCp(rho, p, T, mum, kappam, cpm, onlyInternal, isExtrapolate);
    Diff(p, T, Di, onlyInternal, isExtrapolate);

    for (integer i = 0; i < species().size(); ++i) {
        hai(p, T, hhai[i], i, onlyInternal, isExtrapolate);
    }
}

void OpenHurricane::mixture::gasProperties(const cellRealArray &rho, const cellRealArray &p,
                                           const cellRealArray &T, realArray &mumf,
                                           realArray &kappamf, realArray &cpmf, realArrayArray &Dif,
                                           realArrayArray &haif, const integer faceZoneId) const {
    const auto &fZL = mesh_.faceZones();
    mumf.resize(fZL[faceZoneId].size(), Zero);
    kappamf.resize(fZL[faceZoneId].size(), Zero);
    cpmf.resize(fZL[faceZoneId].size(), Zero);
    const auto cStep = mesh().Iteration().cStep();

    muKappaCp(rho, p, T, mumf, kappamf, cpmf, faceZoneId);

    Dif = Diff(p, T, faceZoneId);
    haif.resize(species().size());
    for (integer i = 0; i < species().size(); ++i) {
        haif[i] = hai(p, T, faceZoneId, i);
    }
}

void OpenHurricane::mixture::E(const cellRealArray &p, const cellRealArray &T,
                               const cellVectorArray &v, cellRealArray &e,
                               const bool internalOrBoundary, const bool isExtrapolate) const {
    const integer nCells = mesh().nCells();
    bool flag = true;
    integer countFlag = 0;
    integer fieldSize = p.size();
    if (internalOrBoundary) {
        for (integer cellI = 0; cellI < mesh().nCells(); cellI++) {
            e[cellI] =
                thermoPtr_->ea_p(p[cellI], T[cellI], Yi(), cellI) + real(0.5) * v[cellI] * v[cellI];
        }
        return;
    }
    if (isExtrapolate) {
        fieldSize = mesh().nCells();
    }

    for (integer cellI = mesh().nCells(); cellI < fieldSize; cellI++) {
        e[cellI] =
            thermoPtr_->ea_p(p[cellI], T[cellI], Yi(), cellI) + real(0.5) * v[cellI] * v[cellI];
        if (std::isnan(e[cellI]) && cellI >= nCells) {
            flag = false;
            countFlag = 1;
            break;
        }
    }
    if ((!internalOrBoundary) && !isExtrapolate) {
        HurMPI::reduce(countFlag, MPI_SUM);
        HurMPI::bcast(&countFlag, 1, feature<integer>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        if (countFlag != 0) {
            flag = false;
        }
    }

    if (((!internalOrBoundary) && isExtrapolate) || !flag) {
        realTransfer myTransfer(mesh_, e, false, false);
        myTransfer.transferInit(0);
        const auto &fZL = mesh_.faceZones();
        for (integer fzi = 0; fzi < fZL.size(); ++fzi) {
            if (!fZL[fzi].isInterior() && !fZL[fzi].isCutFace() && !fZL[fzi].isPeriodic() &&
                !fZL[fzi].isPeriodicShadow()) {
                realArray Ef = E(p, T, v, fzi);
                fv::extrapolate(Ef, e, fzi);
            }
        }
        //fv::transfer(e);
        myTransfer.transferring(0);

        for (integer layerI = 1; layerI < mesh_.ghostCellsType(); layerI++) {
            myTransfer.transferInit(layerI);
            myTransfer.transferring(layerI);
        }
    }
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::mixture::E(const cellRealArray &p,
                                                                 const cellRealArray &T,
                                                                 const cellVectorArray &v,
                                                                 const integer faceZoneId) const {
    const auto &fZL = mesh_.faceZones();
    const integer si = fZL[faceZoneId].size();
    realArray tmpE(si, Zero);
    if (isSingular()) {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        const auto vf = fv::interpolate(v, faceZoneId);

        realArray yi(1, Zero);
        for (integer fi = 0; fi < si; ++fi) {
            tmpE[fi] = thermoPtr_->ea_p(pf[fi], tf[fi], yi) + real(0.5) * vf[fi] * vf[fi];
        }
    } else {
        const auto pf = fv::interpolate(p, faceZoneId);
        const auto tf = fv::interpolate(T, faceZoneId);
        const auto vf = fv::interpolate(v, faceZoneId);
        const auto yif = Yif(faceZoneId);
        for (integer fi = 0; fi < si; ++fi) {
            tmpE[fi] = thermoPtr_->ea_p(pf[fi], tf[fi], yif[fi]) + real(0.5) * vf[fi] * vf[fi];
        }
    }
    return tmpE;
}