/*!
 * \file spongeRegion.cpp
 * \brief Main subroutines for spongeRegion.
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

#include "spongeRegion.hpp"
#include "searchProcedures.hpp"

 namespace OpenHurricane{
	createClassNameStr(spongeRegion,"spongeRegion");
}

namespace OpenHurricane {
     registerObjFty(sourceTerm,spongeRegion,controller);
}

void OpenHurricane::spongeRegion::makeDist() const {
    if (dPtr_ == nullptr) {
        HurMPI::barrier();
        dPtr_ = new cellRealArray(object("dist_" + faceZoneName_, mesh(), object::NOT_WRITE),
                                  mesh(), Zero);
        integerList fzl(1);
        fzl[0] = fzid_;
        distanceMethod *distMPtr = new searchProcedures(mesh(), fzl);
        Pout << std::endl
             << "    Info: getting minium distance from face zone: " + faceZoneName_ << std::endl;
        distMPtr->getDistance(*dPtr_);
    }
}

hur_nodiscard const OpenHurricane::realArray &OpenHurricane::spongeRegion::sigma() const {
    if (mesh().moving()) {
        HurDelete(dPtr_);
        HurDelete(sigmaPtr_);
    }
    makeDist();

    if (sigmaPtr_ == nullptr) {
        sigmaPtr_ = new realArray(mesh().nCells(), Zero);
        const auto &cV = mesh().cellVol();

        auto &sg = *sigmaPtr_;

        const auto &sgl = spongeLayer();
        for (integer i = 0; i < sgl.size(); ++i) {
            const auto id = sgl[i];
            sg[id] = sigma0_ * pow((Delta_ - d()[id]) / max(Delta_, tiny), m_) * cV[id];
        }
    }

    return *sigmaPtr_;
}

OpenHurricane::spongeRegion::spongeRegion(const flowModel &flows, const iteration &iter,
                                      const controller &cont)
    : regionSourceTerms(flows, iter, cont), faceZoneName_(), fzid_(1), Delta_(0),
      sigma0_(cont.findOrDefault<real>("sigma0", 20)), m_(cont.findOrDefault<real>("m", 2)),
      dPtr_(nullptr), sigmaPtr_(nullptr) {
    if (cont.found("faceZoneName")) {
        faceZoneName_ = cont.findWord("faceZoneName");
        trim(faceZoneName_);
        const auto &fz = mesh().faceZones();
        bool found = false;
        for (integer fzi = 0; fzi < fz.size(); ++fzi) {
            if (fz[fzi].name() == faceZoneName_) {
                fzid_ = fzi;
                found = true;
                break;
            }
        }

        if (!found) {
            errorAbortStr(("Unknown face zone name: " + faceZoneName_ + " in: " + cont.name()));
        }
    } else {
        errorAbortStr(("The face zone name \"faceZoneName\" must be specified in " + cont.name()));
    }

    if (cont.found("Delta")) {
        Delta_ = cont.findType<real>("Delta", Delta_);
        if (Delta_ < tiny) {
            errorAbortStr(("The value of \"Delta\" is too small: " + toString(Delta_) + " in " +
                           cont.name()));
        }
    } else {
        errorAbortStr(("The layer width \"Delta\" must be specified in " + cont.name()));
    }

    if (regionList().size() != 1) {
        errorAbortStr(("The number of region must be 1 for each sponge region, "
                       "while it is " +
                       toString(regionList().size()) + " in " + faceZoneName_));
    }

    Pout << "        The sponge region: " << cont.name().c_str()
         << " is connected to face zone: " << faceZoneName_.c_str() << std::endl;
}

void OpenHurricane::spongeRegion::addSourceTerms(cellRealArray &rho) const {
    const auto &sg = sigma();
    const auto &sgl = spongeLayer();
    const auto rhoRef = getRefState<real>(rho.name());
    for (integer i = 0; i < sgl.size(); ++i) {
        const auto id = sgl[i];
        rho.rhs()[id] -= (rho[id] - rhoRef) * sg[id];
    }
}

void OpenHurricane::spongeRegion::addSourceTerms(const cellRealArray &rho, cellVectorArray &U) const {
    const auto &sg = sigma();
    const auto &sgl = spongeLayer();
    const auto rhoRef = getRefState<real>(rho.name());
    const auto URef = getRefState<vector>(U.name());
    for (integer i = 0; i < sgl.size(); ++i) {
        const auto id = sgl[i];
        U.rhs()[id] -= sg[id] * (rho[id] * U[id] - rhoRef * URef);
    }
}

void OpenHurricane::spongeRegion::addSourceTerms(const cellRealArray &rho, cellRealArray &phi) const {
    const auto &sg = sigma();
    const auto &sgl = spongeLayer();
    const auto rhoRef = getRefState<real>(rho.name());
    const auto phiRef = getRefState<real>(phi.name());
    for (integer i = 0; i < sgl.size(); ++i) {
        const auto id = sgl[i];
        phi.rhs()[id] -= sg[id] * (rho[id] * phi[id] - rhoRef * phiRef);
    }
}
