/*!
 * \file CRSMatrixAddressing.cpp
 * \brief Main subroutines of CRSMatrixAddressing.
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
#include "CRSMatrixAddressing.hpp"

void OpenHurricane::CRSMatrixAddressing::allocate(const size_type nCells) {
    if (leftCell_.size() != rightCell_.size()) {
        LFatal("The interior face size of left cell and right cell are not equal!");
    }
    nRows_ = nCells;
    nColumns_ = nCells;
    rowPtr_.resize(nRows_ + 1);
    rowPtr_[0] = 0;

    List<size_type> rowElem(nCells, 1);

    // Diagonal part
    NNZ_ = nCells;

    for (integer fi = 0; fi < leftCell_.size(); ++fi) {
        const auto cl = leftCell_[fi];
        const auto cr = rightCell_[fi];

        rowElem[cl]++;
        NNZ_++;
        rowElem[cr]++;
        NNZ_++;
    }

    col_.resize(NNZ_);

    List<List<size_type>> cellNeber(nCells);
    diagIndex_.resize(nCells);
    cellNeberFace_.resize(nCells);

    for (integer i = 0; i < nCells; ++i) {
        rowPtr_[i + 1] = rowElem[i] + rowPtr_[i];
        cellNeber[i].resize(rowElem[i]);
        cellNeberFace_[i].resize(rowElem[i] - 1);
        cellNeber[i][0] = i;
    }
    rowElem = Zero;

    for (integer fi = 0; fi < leftCell_.size(); ++fi) {
        const auto cl = leftCell_[fi];
        const auto cr = rightCell_[fi];

        cellNeberFace_[cl][rowElem[cl]] = fi;
        cellNeberFace_[cr][rowElem[cr]] = fi;
        rowElem[cl]++;
        rowElem[cr]++;
        cellNeber[cl][rowElem[cl]] = cr;
        cellNeber[cr][rowElem[cr]] = cl;
    }
    for (integer i = 0; i < nCells; ++i) {
        std::sort(cellNeber[i].begin(), cellNeber[i].end());
        for (integer j = 0; j < cellNeber[i].size(); ++j) {
            col_[rowPtr_[i] + j] = cellNeber[i][j];
        }
    }

    List<std::map<size_type, size_type>> cellCellIdMap(nCells);
    for (integer i = 0; i < nCells; ++i) {
        for (integer j = 0; j < rowPtr_[i + 1] - rowPtr_[i]; ++j) {
            cellCellIdMap[i].emplace(col_[rowPtr_[i] + j], rowPtr_[i] + j);
        }
        diagIndex_[i] = cellCellIdMap[i].at(i);
    }
    leftCellIndex_.resize(leftCell_.size());
    rightCellIndex_.resize(leftCell_.size());

    for (integer fi = 0; fi < leftCell_.size(); ++fi) {
        const auto cl = leftCell_[fi];
        const auto cr = rightCell_[fi];
        leftCellIndex_[fi] = cellCellIdMap[cl].at(cr);
        rightCellIndex_[fi] = cellCellIdMap[cr].at(cl);
    }
}

OpenHurricane::CRSMatrixAddressing::CRSMatrixAddressing(const size_type nCells,
                                                    const List<size_type> &leftCell,
                                                    const List<size_type> &rightCell,
                                                    const CRSMatrixAddrCutInterface &proceIntf)
    : leftCell_(leftCell), rightCell_(rightCell), proceIntf_(proceIntf) {
    allocate(nCells);
}

OpenHurricane::CRSMatrixAddressing::CRSMatrixAddressing(const size_type nCells,
                                                    List<size_type> &leftCell,
                                                    List<size_type> &rightCell,
                                                    CRSMatrixAddrCutInterface &proceIntf,
                                                    const bool istransfer) {
    if (istransfer) {
        leftCell_.transfer(leftCell);
        rightCell_.transfer(rightCell);
        proceIntf_.transfer(proceIntf);
    } else {
        leftCell_ = leftCell;
        rightCell_ = rightCell;
        proceIntf_ = proceIntf;
    }
    allocate(nCells);
}

void OpenHurricane::CRSMatrixAddressing::clear() noexcept {
    nRows_ = 0;
    nColumns_ = 0;
    NNZ_ = 0;
    rowPtr_.clear();
    col_.clear();
    diagIndex_.clear();
    cellNeberFace_.clear();
    leftCell_.clear();
    rightCell_.clear();
    leftCellIndex_.clear();
    rightCellIndex_.clear();
    proceIntf_.clear();
}

void OpenHurricane::CRSMatrixAddressing::transfer(CRSMatrixAddressing &other)noexcept {
    nRows_ = other.nRows_;
    nColumns_ = other.nColumns_;
    NNZ_ = other.NNZ_;
    rowPtr_.transfer(other.rowPtr_);
    col_.transfer(other.col_);
    diagIndex_.transfer(other.diagIndex_);
    cellNeberFace_.transfer(other.cellNeberFace_);
    leftCell_.transfer(other.leftCell_);
    rightCell_.transfer(other.rightCell_);
    leftCellIndex_.transfer(other.leftCellIndex_);
    rightCellIndex_.transfer(other.rightCellIndex_);

    proceIntf_.transfer(other.proceIntf_);
}

OpenHurricane::CRSMatrixAddressing &
OpenHurricane::CRSMatrixAddressing::operator=(const CRSMatrixAddressing &other) {
    if (this == std::addressof(other)) {
        return *this;
    }
    nRows_ = other.nRows_;
    nColumns_ = other.nColumns_;
    NNZ_ = other.NNZ_;
    rowPtr_ = other.rowPtr_;
    col_ = other.col_;
    diagIndex_ = other.diagIndex_;
    cellNeberFace_ = other.cellNeberFace_;
    leftCell_ = other.leftCell_;
    rightCell_ = other.rightCell_;
    leftCellIndex_ = other.leftCellIndex_;
    rightCellIndex_ = other.rightCellIndex_;
    proceIntf_ = other.proceIntf_;
    return *this;
}
