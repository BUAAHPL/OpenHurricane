/*!
 * \file CRSMatrixAddressing.hpp
 * \brief Header of sparse square matrices CRS addressing
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

#include "CRSMatrixProcessTransfer.hpp"
#include "List.hpp"
#include "dataStructure.hpp"

namespace OpenHurricane {

    /**
     * \brief sparse square matrices CRS addressing .
     */
    class CRSMatrixAddressing {
    public:
        using size_type = integer;

    private:
        /** \brief Number of rows and columns*/
        size_type nRows_, nColumns_;

        /** \brief The number of non-zero elements. */
        size_type NNZ_;

        /**
         * \brief The index offset of non-zero elements in every row.
         * The size must be nRows_+1.
         */
        List<size_type> rowPtr_;

        /**
         * \brief The column index of non-zero elements.
         * The size must be NNZ_.
         */
        List<size_type> col_;

        List<size_type> diagIndex_;

        List<List<size_type>> cellNeberFace_;

        /**
         * \brief Left cell origin index.
         */
        List<size_type> leftCell_;

        /**
         * \brief Right cell origin index.
         */
        List<size_type> rightCell_;

        /**
         * \brief Left cell index in the CRS data.
         */
        List<size_type> leftCellIndex_;

        /**
         * \brief Right cell index in the CRS data.
         */
        List<size_type> rightCellIndex_;

        void allocate(const size_type nCells);

        /** \brief The interface for data transfer. */
        CRSMatrixAddrCutInterface proceIntf_;

    public:
        // Constructors

        /** \brierf Construct null .*/
        inline CRSMatrixAddressing()
            : nRows_(), nColumns_(), NNZ_(), rowPtr_(), col_(), diagIndex_(), cellNeberFace_(),
              leftCell_(), rightCell_(), leftCellIndex_(), rightCellIndex_(), proceIntf_() {}

        CRSMatrixAddressing(const size_type nCells, const List<size_type> &leftCell,
                            const List<size_type> &rightCell,
                            const CRSMatrixAddrCutInterface &proceIntf);

        CRSMatrixAddressing(const size_type nCells, List<size_type> &leftCell,
                            List<size_type> &rightCell, CRSMatrixAddrCutInterface &proceIntf,
                            const bool istransfer);

        /** \brierf Copy constructor.*/
        inline CRSMatrixAddressing(const CRSMatrixAddressing &other)
            : nRows_(other.nRows_), nColumns_(other.nColumns_), NNZ_(other.NNZ_),
              rowPtr_(other.rowPtr_), col_(other.col_), diagIndex_(other.diagIndex_),
              cellNeberFace_(other.cellNeberFace_), leftCell_(other.leftCell_),
              rightCell_(other.rightCell_), leftCellIndex_(other.leftCellIndex_),
              rightCellIndex_(other.rightCellIndex_), proceIntf_(other.proceIntf_) {}

        /** \brierf Copy constructor.*/
        inline CRSMatrixAddressing(CRSMatrixAddressing &&other) noexcept
            : nRows_(std::move(other.nRows_)), nColumns_(std::move(other.nColumns_)),
              NNZ_(std::move(other.NNZ_)), rowPtr_(std::move(other.rowPtr_)),
              col_(std::move(other.col_)), diagIndex_(std::move(other.diagIndex_)),
              cellNeberFace_(std::move(other.cellNeberFace_)),
              leftCell_(std::move(other.leftCell_)), rightCell_(std::move(other.rightCell_)),
              leftCellIndex_(std::move(other.leftCellIndex_)),
              rightCellIndex_(std::move(other.rightCellIndex_)),
              proceIntf_(std::move(other.proceIntf_)) {}

        /** \brierf Destructor.*/
        inline ~CRSMatrixAddressing() noexcept { clear(); }

        /**
         * \brief Clear.
         */
        void clear() noexcept;

        /** \brief Number of rows and columns*/
        hur_nodiscard inline size_type nRows() const noexcept { return nRows_; }

        hur_nodiscard inline size_type nColumns() const noexcept { return nColumns_; }

        /** \brief The number of non-zero elements. */
        hur_nodiscard inline size_type NNZ() const noexcept { return NNZ_; }

        /**
         * \brief The index offset of non-zero elements in every row.
         * The size must be nRows_+1.
         */
        hur_nodiscard inline const List<size_type> &rowPtr() const noexcept { return rowPtr_; }

        /**
         * \brief The column index of non-zero elements.
         * The size must be NNZ_.
         */
        hur_nodiscard inline const List<size_type> &col() const noexcept { return col_; }

        hur_nodiscard inline const List<size_type> &diagIndex() const noexcept {
            return diagIndex_;
        }

        hur_nodiscard inline const List<List<size_type>> &cellNeberFace() const noexcept {
            return cellNeberFace_;
        }

        /**
         * \brief Left cell origin index.
         */
        hur_nodiscard inline const List<size_type> &leftCell() const noexcept { return leftCell_; }

        /**
         * \brief Right cell origin index.
         */
        hur_nodiscard inline const List<size_type> &rightCell() const noexcept {
            return rightCell_;
        }

        /**
         * \brief Left cell index in the CRS data.
         */
        hur_nodiscard inline const List<size_type> &leftCellIndex() const noexcept {
            return leftCellIndex_;
        }

        /**
         * \brief Right cell index in the CRS data.
         */
        hur_nodiscard inline const List<size_type> &rightCellIndex() const noexcept {
            return rightCellIndex_;
        }

        /** \brierf Transfer the contents of the argument Matrix into this Matrix
         *  and annul the argument Matrix.
         */
        void transfer(CRSMatrixAddressing &other) noexcept;

        /** \brierf Assignment operator.*/
        CRSMatrixAddressing &operator=(const CRSMatrixAddressing &other);

        /** \brierf Assignment operator.*/
        inline CRSMatrixAddressing &operator=(CRSMatrixAddressing &&other) noexcept {
            transfer(other);
            return *this;
        }

        /** \brief The number of interface elements. */
        hur_nodiscard inline size_type NIntf() const noexcept {
            return proceIntf_.intfCellId().size();
        }

        /** \brief The interface for data transfer. */
        hur_nodiscard inline const CRSMatrixAddrCutInterface &proceIntf() const noexcept {
            return proceIntf_;
        }
    };

} //  namespace OpenHurricane
