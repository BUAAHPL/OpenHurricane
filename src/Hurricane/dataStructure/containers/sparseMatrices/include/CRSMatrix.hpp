/*!
 * \file CRSMatrix.hpp
 * \brief Header of CRS (Compressed Sparse Row) sparse CRSMatrix
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
#include "Array.hpp"
#include "CRSMatrixAddressing.hpp"
#include "HurMPIBase.hpp"
#include "List.hpp"
#include "dataStructure.hpp"
#include "matrices.hpp"

namespace OpenHurricane {
    template <class Type> class CRSMatrixSolType {
    public:
        using solution_type = Type;
    };

    template <> class CRSMatrixSolType<real> {
    public:
        using solution_type = real;
    };

    template <> class CRSMatrixSolType<realSquareMatrix> {
    public:
        using solution_type = realArray;
    };

    /**
     * \brief Compressed Row Storage (CRS) for n * n sparse matrix for CFD.
     *      |  10   1   0   0   -2  0  |
     *      |  2    9   8   0   6   1  |
     * A =  |  0    7   8   8   0   0  |
     *      |  0    0   6   4   7   0  |
     *      |  1    7   0   9   1   1  |
     *      |  0    5   0   0   2   3  |
     * then, n = 6, NNZ (Number of non-zero elements) = 21
     * v      = (10,1,-2,2,9,8,6,1,7,8,8,6,4,7,1,7,9,1,1,5,2,3) of size NNZ
     * col    = ( 0,1, 4,0,1,2,4,5,1,2,3,2,3,4,0,1,3,4,5,1,4,5) of size NNZ
     * rowPtr = (0,3,8,11,14,19,22) of size n+1: if v[k] = A(i,j), then j = col[k] and rowPtr[i] <= k < rowPtr[i+1]
     * diagIndex = (0,4,9,12,17,21): index in v and col for diagonal elements, i.e, A(i,i) = v[k], k = diagIndex[i]
     *
     * The CFD sparse matrix is usually non-symmetric.
     */
    template <class Type> class CRSMatrix {
    public:
        using value_type = Type;
        using reference = Type &;
        using const_reference = const Type &;
        using difference_type = integer;
        using size_type = integer;

        using mType = CRSMatrix<Type>;

        using elementType = Type;

        using solution_type = typename CRSMatrixSolType<Type>::solution_type;

    protected:
        CRSMatrixAddressing addressing_;

        /**
         * \brief Value array pointer.
         * The size must be NNZ_;
         */
        Type *hur_restrict v_;

        bool ownData_;

        List<value_type> invAii_;

        List<value_type> interfaceA_;

    public:
        /** \brierf Return a null matrix. */
        inline static const CRSMatrix &nullObject() {
            return NullRefObj::nullRef<CRSMatrix<Type>>();
        }

        CRSMatrix() = delete;

        inline CRSMatrix(const CRSMatrixAddressing &addressing)
            : addressing_(addressing), v_(nullptr), ownData_(true), invAii_(),
              interfaceA_(addressing.NIntf()) {
            if (addressing_.NNZ() != 0) {
                v_ = new elementType[NNZ()];
            }
        }

        /** \brierf Copy constructor. */
        inline CRSMatrix(const CRSMatrix &other)
            : addressing_(other.addressing_), v_(nullptr), ownData_(true), invAii_(other.invAii_),
              interfaceA_(other.interfaceA_) {
            if (rowPtr().size() != 0 && col().size() != 0 && other.v_ != nullptr) {
                v_ = new elementType[NNZ()];

                for (integer i = 0; i < nRows(); ++i) {
                    for (integer j = rowPtr()[i]; j < rowPtr()[i + 1]; ++j) {
                        v_[j] = other.v_[j];
                    }
                }
            }
        }
        inline CRSMatrix &operator=(const CRSMatrix &other) {
            if (this == std::addressof(other)) {
                return *this;
            }
            clear();
            addressing_ = other.addressing_;
            invAii_ = other.invAii_;
            ownData_ = true;
            interfaceA_ = other.interfaceA_;
            if (rowPtr().size() != 0 && col().size() != 0 && other.v_ != nullptr) {
                v_ = new elementType[NNZ()];

                for (integer i = 0; i < nRows(); ++i) {
                    for (integer j = rowPtr()[i]; j < rowPtr()[i + 1]; ++j) {
                        v_[j] = other.v_[j];
                    }
                }
            }
            return *this;
        }

        /** \brierf Copy constructor. */
        inline CRSMatrix(CRSMatrix &&other) noexcept
            : addressing_(std::move(other.addressing_)), v_(std::move(other.v_)), ownData_(true),
              invAii_(std::move(other.invAii_)), interfaceA_(std::move(other.interfaceA_)) {
            other.v_ = nullptr;
            other.ownData_ = false;
        }
        inline CRSMatrix &operator=(CRSMatrix &&other) noexcept {
            transfer(other);
            return *this;
        }

        /** \brierf Clone. */
        hur_nodiscard inline uniquePtr<CRSMatrix> clone() const {
            return uniquePtr<CRSMatrix>(new CRSMatrix(*this));
        }

        /** \brierf Destructor. */
        inline ~CRSMatrix() noexcept { clear(); }

        hur_nodiscard inline const CRSMatrixAddressing &addressing() const noexcept {
            return addressing_;
        }

        /** \brierf Return the number of rows. */
        hur_nodiscard inline size_type nRows() const noexcept { return addressing_.nRows(); }

        /** \brierf Return the number of columns. */
        hur_nodiscard inline size_type nColumns() const noexcept { return addressing_.nColumns(); }

        /** \brierf Return the number of elements in matrix. */
        hur_nodiscard inline size_type size() const noexcept { return nRows() * nColumns(); }

        /** \brierf Return element vector of the constant Matrix. */
        hur_nodiscard inline const Type *v() const noexcept { return v_; }

        /** \brierf Return element vector of the Matrix. */
        hur_nodiscard inline Type *v() noexcept { return v_; }

        /** \brierf Return element vector of the constant Matrix. */
        hur_nodiscard inline const Type *data() const noexcept { return v_; }

        /** \brierf Return element vector of the Matrix. */
        hur_nodiscard inline Type *data() noexcept { return v_; }

        /** \brief The number of non-zero elements. */
        hur_nodiscard inline size_type NNZ() const noexcept { return addressing_.NNZ(); }

        /**
         * \brief The index offset of non-zero elements in every row.
         * The size must be nRows_+1.
         */
        hur_nodiscard inline const List<size_type> &rowPtr() const noexcept {
            return addressing_.rowPtr();
        }

        /**
         * \brief The column index of non-zero elements.
         * The size must be NNZ_.
         */
        hur_nodiscard inline const List<size_type> &col() const noexcept {
            return addressing_.col();
        }

        hur_nodiscard inline bool ownData() const noexcept { return ownData_; }

        hur_nodiscard inline const List<size_type> &diagIndex() const noexcept {
            return addressing_.diagIndex();
        }

        hur_nodiscard inline size_type diagIndex(const size_type i) const noexcept {
            return addressing_.diagIndex()[i];
        }

        hur_nodiscard inline const List<size_type> &leftCellIndex() const noexcept {
            return addressing_.leftCellIndex();
        }

        hur_nodiscard inline const List<size_type> &rightCellIndex() const noexcept {
            return addressing_.rightCellIndex();
        }

        hur_nodiscard inline value_type &Aii(const size_type i) noexcept {
            return this->v_[diagIndex(i)];
        }

        hur_nodiscard inline const value_type &Aii(const size_type i) const noexcept {
            return this->v_[diagIndex(i)];
        }

        hur_nodiscard inline value_type &invAii(const size_type i) {
            if (!hasInversedDiag()) {
                LFatal("Attempt to access null inversed diagonal parts");
            }
            return invAii_[i];
        }

        hur_nodiscard inline const value_type &invAii(const size_type i) const {
            if (!hasInversedDiag()) {
                LFatal("Attempt to access null inversed diagonal parts");
            }
            return invAii_[i];
        }
        void setInvAii() { invAii_.resize(nCells()); }
        hur_nodiscard inline bool hasInversedDiag() const noexcept { return !invAii_.empty(); }

        hur_nodiscard inline List<value_type> &interfaceA() noexcept { return interfaceA_; }
        hur_nodiscard inline const List<value_type> &interfaceA() const noexcept {
            return interfaceA_;
        }

        hur_nodiscard inline size_type nCells() const noexcept { return diagIndex().size(); }

        /** \brierf Clear the Matrix. */
        inline void clear() noexcept {
            addressing_.clear();
            if (ownData_) {
                HurDeleteDynArray(v_);
            } else {
                v_ = nullptr;
            }
            invAii_.clear();
            interfaceA_.clear();
        }

        /** \brierf Transfer the contents of the argument Matrix into this Matrix
          and annul the argument Matrix. */
        void transfer(CRSMatrix &other) noexcept {
            if (this == std::addressof(other)) {
                return;
            }
            clear();
            addressing_.transfer(other.addressing_);
            v_ = other.v_;
            other.v_ = nullptr;
            other.ownData_ = false;
            ownData_ = true;
            invAii_.transfer(other.invAii_);
            interfaceA_.transfer(other.interfaceA_);
        }

        void multiple(const Array<solution_type> &x, Array<solution_type> &Ax) const {
            const size_type nCells = nRows();
            Ax = Zero;
            size_type uEnd = rowPtr()[0];
            size_type uStart;
            for (size_type i = 0; i < nCells; i++) {
                uStart = uEnd;
                uEnd = rowPtr()[i + 1];
                for (integer j0 = uStart; j0 < uEnd; ++j0) {
                    size_type j = col()[j0];
                    Ax[i] += v_[j0] * x[j];
                }
            }
        }

        void interfaceTransfer(const Array<solution_type> &x, Array<solution_type> &b) const {
            if (interfaceA_.size() == 0) {
                return;
            }

            CRSMatrixProcessTransfer<solution_type> procT(addressing_.proceIntf());
            Array<solution_type> recvX(interfaceA_.size());
            procT.transferData(x, recvX);

            const auto &intfCI = addressing_.proceIntf().intfCellId();
            for (integer i = 0; i < intfCI.size(); ++i) {
                b[intfCI[i]] -= interfaceA_[i] * recvX[i];
            }
        }
    };

    template <>
    void CRSMatrix<realSquareMatrix>::interfaceTransfer(const Array<solution_type> &x,
                                                        Array<solution_type> &b) const;

} // namespace OpenHurricane
