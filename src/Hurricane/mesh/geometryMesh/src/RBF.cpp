/*!
 * \file RBF.cpp
 * \brief The subroutines and functions of RBF.
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

#include "RBF.hpp"
#include "Eigen/Eigen"
#include "LUDecompose.hpp"
#include "geometryMesh.hpp"

 namespace OpenHurricane{
	createClassNameStr(RBF,"RBF");
}

OpenHurricane::RBF::RBF(const controller &cont, const geometryMesh &mesh, const faceZone &movingZone,
                    const integer fzid)
    : mesh_(mesh), movingZone_(movingZone), basePoints_(), oldBasePoint_(), CSFPtr_(nullptr),
      M_(3, Zero), alphax_(3, Zero), alphay_(3, Zero), alphaz_(3, Zero),
      selectedPointList_(3, Zero), firstGuess_(true),
      epsilon_(cont.findOrDefault<real>("tol", 1e-10)) {
    integer sss = 0;
    if (HurMPI::master()) {
        basePoints_ = mesh_.globalFaceZoneInfo(fzid).facePoints();
        sss = basePoints_.size();
    }
    HurMPI::bcast(&sss, 1, feature<integer>::MPIType, HurMPI::masterNo());
    if (!HurMPI::master()) {
        basePoints_.resize(sss);
    }
    HurMPI::bcastVectorList(basePoints_, HurMPI::masterNo(), HurMPI::getComm());

    if (cont.found("compactSupportFunction")) {
        CSFPtr_ = compactSupportFunction::creator(cont);
    } else {
        LFatal("Please specify the compact support function");
    }
}

hur_nodiscard OpenHurricane::vector OpenHurricane::RBF::g(const vector &pi) const {
    vector gg(Zero);
    for (integer i = 0; i < selectedPointList_.size(); ++i) {
        const integer j = selectedPointList_[i];
        const auto phij = CSFPtr_->phi(dist(pi, basePoints_[j]));
        gg.x() += alphax_[j] * phij;
        gg.y() += alphay_[j] * phij;
        gg.z() += alphaz_[j] * phij;
    }
    return gg;
}

void OpenHurricane::RBF::g(const vector &pi, vector &gg) const {
    gg = Zero;
    for (integer i = 0; i < selectedPointList_.size(); ++i) {
        const integer j = selectedPointList_[i];
        const auto phij = CSFPtr_->phi(dist(pi, basePoints_[j]));
        gg.x() += alphax_[j] * phij;
        gg.y() += alphay_[j] * phij;
        gg.z() += alphaz_[j] * phij;
    }
}

hur_nodiscard OpenHurricane::vector OpenHurricane::RBF::newPoint(const vector &oldPi) const {
    return oldPi + g(oldPi);
}

void OpenHurricane::RBF::newPoint(const vector &oldPi, vector &np) const {
    g(oldPi, np);
    np += oldPi;
}

void OpenHurricane::RBF::getAlpha() {
    integer istart = 0;
    if (firstGuess_) {
        firstGuess_ = false;
    } else {
        istart = selectedPointList_.size() - 1;
    }
    for (integer i = istart; i < selectedPointList_.size(); ++i) {
        const integer pi = selectedPointList_[i];
        alphax_[i] = oldBasePoint_[pi].x() - basePoints_[pi].x();
        alphay_[i] = oldBasePoint_[pi].y() - basePoints_[pi].y();
        alphaz_[i] = oldBasePoint_[pi].z() - basePoints_[pi].z();
        for (integer j = i; j < selectedPointList_.size(); ++j) {
            const integer pj = selectedPointList_[j];

            // It is also phi_ji
            const auto phi_ij = CSFPtr_->phi(dist(basePoints_[pi], basePoints_[pj]));
            M_[i][j] = phi_ij;
        }
    }

    integerList pivotIndices_(selectedPointList_.size());
    LUDecomposeDPivoting(M_, pivotIndices_);
    LUBacksubstituteDPivoting(M_, alphax_, pivotIndices_);
    LUBacksubstituteDPivoting(M_, alphay_, pivotIndices_);
    LUBacksubstituteDPivoting(M_, alphaz_, pivotIndices_);
}

void OpenHurricane::RBF::initialSelected() {
    firstGuess_ = true;
    alphax_.resize(3, Zero);
    alphay_.resize(3, Zero);
    alphaz_.resize(3, Zero);
    M_.resize(3);

    selectedPointList_.resize(3);
    selectedPointList_[0] = 0;
    selectedPointList_[1] = integer(basePoints_.size() / 2);
    selectedPointList_[2] = basePoints_.size() - 1;
}

void OpenHurricane::RBF::solving() {
    initialSelected();

    while (true) {
        getAlpha();
        real maxError = Zero;
        integer flag = -1;
        for (integer pi = 0; pi < basePoints_.size(); ++pi) {
            const auto np = newPoint(oldBasePoint_[pi]);
            const auto err = Derror(np, basePoints_[pi]);
            if (err > maxError) {
                maxError = err;
                flag = pi;
            }
        }

        if (isSatisfied(maxError)) {
            break;
        } else {
            const integer oldSize = selectedPointList_.size();
            alphax_.resize(oldSize + 1);
            alphay_.resize(oldSize + 1);
            alphaz_.resize(oldSize + 1);
            selectedPointList_.push_back(flag);
            M_.resize(oldSize + 1);
        }
    }
}

hur_nodiscard bool OpenHurricane::RBF::isSatisfied(const real err) const {
    return err < epsilon_;
}