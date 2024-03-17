/*!
 * \file kNearestNeighbour.cpp
 * \brief The subroutines and functions of CFD time advance iteration
 * \author Yang Hongzhen
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

#include "kNNGrid.hpp"
#include "KDTree.hpp"

OpenHurricane::kNearestNeighbour::kNearestNeighbour(const controller &cont, const vectorArray &sor,
                                                const vectorArray &tar)
    : sor_(sor), tar_(tar), nbrPtr_(nullptr) {
    if (cont.found("nearNodes")) {
        k_ = cont.findType<integer>("nearNodes", k_);
    } else {
        k_ = 1;
    }

    if (k_ >= sor_.size()) {
        std::string errMsg =
            "Value of k_ is larger than the number of source nodes: " + toString(sor_.size());
        errMsg += "\n  Current k_ = ";
        errMsg += toString(k_);
        errMsg += "\n  Please check if the new domain is larger than the "
                  "source domain?";
        errorAbortStr(errMsg);
    }
}

OpenHurricane::kNearestNeighbour::kNearestNeighbour(const integer k, const vectorArray &sor,
                                                const vectorArray &tar)
    : k_(k), sor_(sor), tar_(tar), nbrPtr_(nullptr) {
    if (k_ >= sor_.size()) {
        std::string errMsg =
            "Value of k_ is larger than the number of source nodes: " + toString(sor_.size());
        errMsg += "\nCurrent k_ = ";
        errMsg += toString(k_);
        errMsg += "\n  Please check if the new domain is larger than the "
                  "source domain?";
        errorAbortStr(errMsg);
    }
}

OpenHurricane::kNNGrid::kNNGrid(const controller &cont, const vectorArray &sor, const vectorArray &tar)
    : kNearestNeighbour(cont, sor, tar), l_(), ref_(), contain_() {}

OpenHurricane::kNNGrid::kNNGrid(const integer k, const vectorArray &sor, const vectorArray &tar)
    : kNearestNeighbour(k, sor, tar), l_(), ref_(), contain_() {}

void OpenHurricane::kNNGrid::getNearestNeigbhour() {
    creatBoundBox();
    vector bb = bb_.span();

    integer nx = (integer)std::ceil(bb.x() / max(veryTiny, l_.x())) + 1;
    integer ny = (integer)std::ceil(bb.y() / max(veryTiny, l_.y())) + 1;
    integer nz = (integer)std::ceil(bb.z() / max(veryTiny, l_.z())) + 1;

    integerListListList nt(nx);
    contain_.resize(nx);

    for (integer i = 0; i < nx; i++) {
        nt[i].resize(ny);
        contain_[i].resize(ny);
        for (integer j = 0; j < ny; j++) {
            nt[i][j].resize(nz, Zero);
            contain_[i][j].resize(nz);
        }
    }

    //pre-treatment of searching

    integer q = 0;
    integer w = 0;

    for (integer i = 0; i < sor_.size(); i++) {
        integerVector ijk = getIndex(sor_[i]);

        nt[ijk.x()][ijk.y()][ijk.z()]++;
    }

    for (integer i = 0; i < nx; i++) {
        for (integer j = 0; j < ny; j++) {
            for (integer k = 0; k < nz; k++) {
                if (nt[i][j][k] == 0) {
                    w++;
                } else if (nt[i][j][k] >= 20) {
                    q++;
                }
            }
        }
    }

    for (integer i = 0; i < nx; i++) {
        for (integer j = 0; j < ny; j++) {
            for (integer k = 0; k < nz; k++) {
                contain_[i][j][k].resize(nt[i][j][k]);
                nt[i][j][k] = 0;
            }
        }
    }

    for (integer i = 0; i < sor_.size(); i++) {
        integerVector ijk = getIndex(sor_[i]);
        contain_[ijk.x()][ijk.y()][ijk.z()][nt[ijk.x()][ijk.y()][ijk.z()]++] = i;
    }

    //searching
    nbrPtr_.reset(new integerListList(tar_.size()));
    for (integer i = 0; i < tar_.size(); i++) {
        nbrPtr_->operator[](i).resize(k_);
    }

    for (integer n = 0; n < tar_.size(); n++) {
        integerVector ijk1 = getIndex(tar_[n]);
        integerVector ijk2 = ijk1;

        realArray d(k_, 1e50);
        bool flag = true;

        while (flag) {
            for (integer i = ijk1.x(); i <= ijk2.x(); i++) {
                for (integer j = ijk1.y(); j <= ijk2.y(); j++) {
                    for (integer k = ijk1.z(); k <= ijk2.z(); k++) {
                        for (integer m = 0; m < nt[i][j][k]; m++) {
                            integer nn = contain_[i][j][k][m];
                            real dist = (sor_[nn] - tar_[n]).magnitude() / ref_;

                            for (integer o = 0; o < k_; o++) {
                                if (dist < d[o]) {
                                    for (integer o1 = k_ - 1; o1 > o; o1--) {
                                        d[o1] = d[o1 - 1];
                                        (*nbrPtr_)[n][o1] = (*nbrPtr_)[n][o1 - 1];
                                    }
                                    d[o] = dist;
                                    (*nbrPtr_)[n][o] = nn;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            realArray uu(2 * feature<vector>::nElements_);
            for (integer dim = 0; dim < feature<vector>::nElements_; dim++) {
                uu[2 * dim] = tar_[n][dim] / ref_ - bb_.min()[dim] - l_[dim] * ijk1[dim];
                uu[2 * dim + 1] = bb_.min()[dim] + l_[dim] * (ijk2[dim] + 1) - tar_[n][dim] / ref_;
            }

            integerList ib(2 * feature<vector>::nElements_, Zero);
            expandIJK(ijk1, ijk2, nx, ny, nz, uu, d[k_ - 1], flag, ib);
        }
    }

#ifdef HUR_DEBUG

    std::cout << HurMPI::getProcRank() << " done " << std::endl;
#endif // HUR_DEBUG

    HurMPI::barrier();
}

hur_nodiscard OpenHurricane::integerVector OpenHurricane::kNNGrid::getIndex(const vector point) const {
    integerVector index;
    for (integer i = 0; i < feature<integerVector>::nElements_; i++) {
        if (l_[i] <= veryTiny) {
            index[i] = 0;
        } else {
            index[i] = integer((point[i] / ref_ - bb_.min()[i]) / l_[i]);
        }
    }
    return index;
}

void OpenHurricane::kNNGrid::creatBoundBox() {
    boundBox sorBB(sor_);
    boundBox tarBB(tar_);

    vector min = sorBB.min();
    vector max = sorBB.max();

    if (tarBB.min().x() < min.x()) {
        min.x() = tarBB.min().x();
    }
    if (tarBB.min().y() < min.y()) {
        min.y() = tarBB.min().y();
    }
    if (tarBB.min().z() < min.z()) {
        min.z() = tarBB.min().z();
    }

    if (tarBB.max().x() > max.x()) {
        max.x() = tarBB.max().x();
    }
    if (tarBB.max().y() > max.y()) {
        max.y() = tarBB.max().y();
    }
    if (tarBB.max().z() > max.z()) {
        max.z() = tarBB.max().z();
    }

    real *length = new real[feature<vector>::nElements_];
    for (integer i = 0; i < feature<vector>::nElements_; i++) {
        length[i] = max[i] - min[i];
    }
    ref_ = *std::max_element(length, length + feature<vector>::nElements_);

    max /= ref_;
    min /= ref_;

    max += veryTiny;
    min -= veryTiny;

    bb_.min() = min;
    bb_.max() = max;

    real num = real(feature<vector>::nElements_);
    real volume = bb_.volume();

    real beta = std::pow(5.0 * real(k_), 1.0 / num);
    l_ = beta * std::pow(volume / real(sor_.size()), 1.0 / num);

    vector span = bb_.span();

    vector qwe = l_ * ref_;

    for (integer i = 0; i < feature<vector>::nElements_; i++) {
        if (l_[i] < 0.01 * span[i]) {
            l_[i] = 0.01 * span[i];
        }
    }
}

void OpenHurricane::kNNGrid::expandIJK(integer &i1, integer &i2, integer &j1, integer &j2, integer &k1,
                                   integer &k2, const integer nx, const integer ny,
                                   const integer nz, const realArray &uu, const real dm, bool &flag,
                                   integerList &ib) {
    //integerList ib(2 * feature<vector>::nElements_, Zero);
    integer nn = -1;
    real dmax = 1e50;
    for (integer i = 0; i < 2 * feature<vector>::nElements_; i++) {
        if (ib[i] == 1) {
            continue;
        }
        if (uu[i] < dmax) {
            dmax = uu[i];
            nn = i;
        }
    }
    if (dmax >= dm) {
        flag = false;
        return;
    }
    switch (nn) {
    case 0:
        if (i1 == 0) {
            ib[0] = 1;
            expandIJK(i1, i2, j1, j2, k1, k2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        i1 -= 1;
        break;
    case 1:
        if (i2 == nx - 1) {
            ib[1] = 1;
            expandIJK(i1, i2, j1, j2, k1, k2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        i2 += 1;
        break;
    case 2:
        if (j1 == 0) {
            ib[2] = 1;
            expandIJK(i1, i2, j1, j2, k1, k2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        j1 -= 1;
        break;
    case 3:
        if (j2 == ny - 1) {
            ib[3] = 1;
            expandIJK(i1, i2, j1, j2, k1, k2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        j2 += 1;
        break;
    case 4:
        if (k1 == 0) {
            ib[4] = 1;
            expandIJK(i1, i2, j1, j2, k1, k2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        k1 -= 1;
        break;
    case 5:
        if (k2 == nz - 1) {
            ib[5] = 1;
            expandIJK(i1, i2, j1, j2, k1, k2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        k2 += 1;
        break;
    }
}

void OpenHurricane::kNNGrid::expandIJK(integerVector &ijk1, integerVector &ijk2, const integer nx,
                                   const integer ny, const integer nz, const realArray &uu,
                                   const real dm, bool &flag, integerList &ib) {
    //integerList ib(2 * feature<vector>::nElements_, Zero);
    integer nn = -1;
    real dmax = 1e50;
    for (integer i = 0; i < 2 * feature<vector>::nElements_; i++) {
        if (ib[i] == 1) {
            continue;
        }
        if (uu[i] < dmax) {
            dmax = uu[i];
            nn = i;
        }
    }
    if (dmax >= dm) {
        flag = false;
        return;
    }
    switch (nn) {
    case 0:
        if (ijk1.x() == 0) {
            ib[0] = 1;
            expandIJK(ijk1, ijk2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        ijk1.x() -= 1;
        break;
    case 1:
        if (ijk2.x() == nx - 1) {
            ib[1] = 1;
            expandIJK(ijk1, ijk2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        ijk2.x() += 1;
        break;
    case 2:
        if (ijk1.y() == 0) {
            ib[2] = 1;
            expandIJK(ijk1, ijk2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        ijk1.y() -= 1;
        break;
    case 3:
        if (ijk2.y() == ny - 1) {
            ib[3] = 1;
            expandIJK(ijk1, ijk2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        ijk2.y() += 1;
        break;
    case 4:
        if (ijk1.z() == 0) {
            ib[4] = 1;
            expandIJK(ijk1, ijk2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        ijk1.z() -= 1;
        break;
    case 5:
        if (ijk2.z() == nz - 1) {
            ib[5] = 1;
            expandIJK(ijk1, ijk2, nx, ny, nz, uu, dm, flag, ib);
            return;
        }
        ijk2.z() += 1;
        break;
    }
}