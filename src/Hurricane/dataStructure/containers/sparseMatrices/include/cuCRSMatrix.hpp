/*!
 * \file cuCRSMatrix.hpp
 * \brief Header of sparse square matrices CRS matrix for CUDA
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

#include "cuArrayBase.hpp"

#ifdef CUDA_PARALLEL

namespace OpenHurricane {

    /**
     * \brief sparse square matrices CRS matrix addressing for CUDA managed.
     */
    class cuMCRSMatrixAddressing {
    public:
        using size_type = cu_integer;

    private:
        /** \brief Number of rows and columns*/
        size_type nRows_, nColumns_;

        /** \brief The number of non-zero elements. */
        size_type NNZ_;

        /**
         * \brief The index offset of non-zero elements in every row.
         * The size must be nRows_+1.
         */
        cuMArray<size_type, size_type> rowPtr_;

        /**
         * \brief The column index of non-zero elements.
         * The size must be NNZ_.
         */
        cuMArray<size_type, size_type> col_;

        cuMArray<size_type, size_type> diagIndex_;

    public:
        /** \brierf Construct null .*/
        cu_dual inline cuMCRSMatrixAddressing()
            : nRows_(0), nColumns_(0), NNZ_(0), rowPtr_(), col_(), diagIndex_() {}

        cu_host cuMCRSMatrixAddressing(const size_type nRows, const size_type nColumns,
                                       const size_type NNZ, const size_type *__restrict__ rowPtr,
                                       const size_type *__restrict__ col);

        /** \brierf Copy constructor.*/
        cu_dual inline cuMCRSMatrixAddressing(const cuMCRSMatrixAddressing &other)
            : nRows_(other.nRows_), nColumns_(other.nColumns_), NNZ_(other.NNZ_),
              rowPtr_(other.rowPtr_), col_(other.col_), diagIndex_(other.diagIndex_) {}

        /** \brierf Assignment operator.*/
        cuMCRSMatrixAddressing &operator=(const cuMCRSMatrixAddressing &other) = delete;

        /** \brierf Destructor.*/
        cu_dual inline ~cuMCRSMatrixAddressing() noexcept {}

        /**
         * \brief Clear.
         */
        cu_host void clear() noexcept;

        /** \brief Number of rows and columns*/
        hur_nodiscard cu_dual inline size_type nRows() const noexcept { return nRows_; }

        hur_nodiscard cu_dual inline size_type nColumns() const noexcept { return nColumns_; }

        /** \brief The number of non-zero elements. */
        hur_nodiscard cu_dual inline size_type NNZ() const noexcept { return NNZ_; }

        /**
         * \brief The index offset of non-zero elements in every row.
         * The size must be nRows_+1.
         */
        hur_nodiscard cu_dual inline const auto &rowPtr() const noexcept { return rowPtr_; }

        /**
         * \brief The column index of non-zero elements.
         * The size must be NNZ_.
         */
        hur_nodiscard cu_dual inline const auto &col() const noexcept { return col_; }

        hur_nodiscard cu_dual inline const auto &diagIndex() const noexcept { return diagIndex_; }
        cu_host inline void memPrefetchAsync(const int devId, cudaStream_t s = 0) const {
            rowPtr_.memPrefetchAsync(devId, s);
            col_.memPrefetchAsync(devId, s);
            diagIndex_.memPrefetchAsync(devId, s);
        }
    };

    /**
     * \brief sparse square matrices CRS matrix for CUDA managed.
     */
    template <class Type> class cuMCRSMatrix : public cuMCRSMatrixAddressing {
    public:
        using value_type = Type;
        using size_type = cuMCRSMatrixAddressing::size_type;

        using mType = cuMCRSMatrix<Type>;

        using solution_type = cuMArray<Type, size_type>;

    private:
        /**
         * \brief Value array pointer.
         * The size must be NNZ_;
         */
        cuMArray<Type, size_type> value_;

    public:
        cu_dual inline cuMCRSMatrix() : cuMCRSMatrixAddressing(), value_() {}

        cu_host inline cuMCRSMatrix(const size_type nRows, const size_type nColumns,
                                    const size_type NNZ, const size_type *__restrict__ rowPtr,
                                    const size_type *__restrict__ col)
            : cuMCRSMatrixAddressing(nRows, nColumns, NNZ, rowPtr, col), value_(NNZ) {}

        cu_dual inline cuMCRSMatrix(const cuMCRSMatrix &other)
            : cuMCRSMatrixAddressing(other), value_(other.value_) {}

        cuMCRSMatrix &operator=(const cuMCRSMatrix &other) = delete;

        cu_dual inline ~cuMCRSMatrix() noexcept {}

        cu_host inline void clear() noexcept {
            cuMCRSMatrixAddressing::clear();
            value_.clear();
        }

        /** \brierf Return element vector of the constant Matrix. */
        hur_nodiscard cu_dual inline const Type *data() const noexcept { return value_.data(); }

        /** \brierf Return element vector of the Matrix. */
        hur_nodiscard cu_dual inline Type *data() noexcept { return value_.data(); }

        cu_host inline void memPrefetchAsync(const int devId, cudaStream_t s = 0) const {
            cuMCRSMatrixAddressing::memPrefetchAsync(devId, s);
            value_.memPrefetchAsync(devId, s);
        }
    };

    /**
     * \brief sparse square matrices CRS matrix for CUDA managed.
     */
    template <class Type> class cuMCRSBlockMatrix : public cuMCRSMatrixAddressing {
    public:
        using value_type = Type;
        using size_type = cuMCRSMatrixAddressing::size_type;

        using mType = cuMCRSBlockMatrix<Type>;
        using solution_type = cuMBlockArray<Type, size_type>;

    private:
        /**
         * \brief Value array pointer.
         * The size must be NNZ_;
         */
        cuMBlockArray<Type, size_type> value_;

    public:
        cu_dual inline cuMCRSBlockMatrix() : cuMCRSMatrixAddressing(), value_() {}

        cu_host inline cuMCRSBlockMatrix(const size_type nRows, const size_type nColumns,
                                         const size_type NNZ, const size_type *__restrict__ rowPtr,
                                         const size_type *__restrict__ col,
                                         const cu_ushort blockDim)
            : cuMCRSMatrixAddressing(nRows, nColumns, NNZ, rowPtr, col), value_(NNZ, blockDim) {}

        cu_host inline cuMCRSBlockMatrix(const size_type nRows, const size_type nColumns,
                                         const size_type NNZ, const size_type *__restrict__ rowPtr,
                                         const size_type *__restrict__ col,
                                         const cu_ushort blockDimx, const cu_ushort blockDimy)
            : cuMCRSMatrixAddressing(nRows, nColumns, NNZ, rowPtr, col),
              value_(NNZ, blockDimx, blockDimy) {}

        cu_dual inline cuMCRSBlockMatrix(const cuMCRSBlockMatrix &other)
            : cuMCRSMatrixAddressing(other), value_(other.value_) {}

        cuMCRSBlockMatrix &operator=(const cuMCRSBlockMatrix &other) = delete;

        cu_dual inline ~cuMCRSBlockMatrix() noexcept {}

        cu_host inline void clear() noexcept {
            cuMCRSMatrixAddressing::clear();
            value_.clear();
        }

        /** \brierf Return element vector of the constant Matrix. */
        hur_nodiscard cu_dual inline const Type *data() const noexcept { return value_.data(); }

        /** \brierf Return element vector of the Matrix. */
        hur_nodiscard cu_dual inline Type *data() noexcept { return value_.data(); }

        /**
         * \brief Value array pointer.
         * The size must be NNZ_;
         */
        hur_nodiscard cu_dual inline cuMBlockArray<Type, size_type> &value() noexcept {
            return value_;
        }

        /**
         * \brief Value array pointer.
         * The size must be NNZ_;
         */
        hur_nodiscard cu_dual inline const cuMBlockArray<Type, size_type> &value() const noexcept {
            return value_;
        }

        cu_host inline void memPrefetchAsync(const int devId, cudaStream_t s = 0) const {
            cuMCRSMatrixAddressing::memPrefetchAsync(devId, s);
            value_.memPrefetchAsync(devId, s);
        }
    };

    /**
     * \brief sparse square matrices CRS matrix addressing for CUDA.
     */
    class cuCRSMatrixAddressing {
    public:
        using size_type = cu_integer;

    private:
        /** \brief Number of rows and columns*/
        size_type nRows_, nColumns_;

        /** \brief The number of non-zero elements. */
        size_type NNZ_;

        /**
         * \brief The index offset of non-zero elements in every row.
         * The size must be nRows_+1.
         */
        cuArray<size_type, size_type> rowPtr_;

        /**
         * \brief The column index of non-zero elements.
         * The size must be NNZ_.
         */
        cuArray<size_type, size_type> col_;

        cuArray<size_type, size_type> diagIndex_;

    public:
        /** \brierf Construct null .*/
        cu_dual inline cuCRSMatrixAddressing()
            : nRows_(0), nColumns_(0), NNZ_(0), rowPtr_(), col_(), diagIndex_() {}

        cu_host cuCRSMatrixAddressing(const size_type nRows, const size_type nColumns,
                                      const size_type NNZ, const size_type *__restrict__ rowPtr,
                                      const size_type *__restrict__ col);

        /** \brierf Copy constructor.*/
        cu_dual inline cuCRSMatrixAddressing(const cuCRSMatrixAddressing &other)
            : nRows_(other.nRows_), nColumns_(other.nColumns_), NNZ_(other.NNZ_),
              rowPtr_(other.rowPtr_), col_(other.col_), diagIndex_(other.diagIndex_) {}

        /** \brierf Assignment operator.*/
        cuCRSMatrixAddressing &operator=(const cuCRSMatrixAddressing &other) = delete;

        /** \brierf Destructor.*/
        cu_dual inline ~cuCRSMatrixAddressing() noexcept {}

        /**
         * \brief Clear.
         */
        cu_host void clear() noexcept;

        /** \brief Number of rows and columns*/
        hur_nodiscard cu_dual inline size_type nRows() const noexcept { return nRows_; }

        hur_nodiscard cu_dual inline size_type nColumns() const noexcept { return nColumns_; }

        /** \brief The number of non-zero elements. */
        hur_nodiscard cu_dual inline size_type NNZ() const noexcept { return NNZ_; }

        /**
         * \brief The index offset of non-zero elements in every row.
         * The size must be nRows_+1.
         */
        hur_nodiscard cu_dual inline const auto &rowPtr() const noexcept { return rowPtr_; }

        /**
         * \brief The column index of non-zero elements.
         * The size must be NNZ_.
         */
        hur_nodiscard cu_dual inline const auto &col() const noexcept { return col_; }

        hur_nodiscard cu_dual inline const auto &diagIndex() const noexcept { return diagIndex_; }

        cu_host inline void buildAddressingInDevice() {
            rowPtr_.copyFromHost();
            col_.copyFromHost();
            diagIndex_.copyFromHost();
        }
        cu_host inline void buildAddressingInDeviceAsync(cudaStream_t s) {
            rowPtr_.copyFromHostAsync(s);
            col_.copyFromHostAsync(s);
            diagIndex_.copyFromHostAsync(s);
        }
    };

    /**
     * \brief sparse square matrices CRS matrix for CUDA.
     */
    template <class Type> class cuCRSBlockMatrix : public cuCRSMatrixAddressing {
    public:
        using value_type = Type;
        using size_type = cuCRSMatrixAddressing::size_type;

        using mType = cuCRSBlockMatrix<Type>;
        
        using solution_type = cuBlockArray1<Type, size_type>;
    private:
        /**
         * \brief Value array pointer.
         * The size must be NNZ_;
         */
        cuBlockArray1<Type, size_type> value_;

    public:
        cu_dual inline cuCRSBlockMatrix() : cuCRSMatrixAddressing(), value_() {}

        cu_host inline cuCRSBlockMatrix(const size_type nRows, const size_type nColumns,
                                        const size_type NNZ, const size_type *__restrict__ rowPtr,
                                        const size_type *__restrict__ col, const cu_ushort blockDim)
            : cuCRSMatrixAddressing(nRows, nColumns, NNZ, rowPtr, col), value_(NNZ, blockDim) {}

        cu_host inline cuCRSBlockMatrix(const size_type nRows, const size_type nColumns,
                                        const size_type NNZ, const size_type *__restrict__ rowPtr,
                                        const size_type *__restrict__ col,
                                        const cu_ushort blockDimx, const cu_ushort blockDimy)
            : cuCRSMatrixAddressing(nRows, nColumns, NNZ, rowPtr, col),
              value_(NNZ, blockDimx, blockDimy) {}

        cu_dual inline cuCRSBlockMatrix(const cuCRSBlockMatrix &other)
            : cuCRSMatrixAddressing(other), value_(other.value_) {}

        cuCRSBlockMatrix &operator=(const cuCRSBlockMatrix &other) = delete;

        cu_dual inline ~cuCRSBlockMatrix() noexcept {}

        cu_host inline void clear() noexcept {
            cuCRSMatrixAddressing::clear();
            value_.clear();
        }

        /** \brierf Return element vector of the constant Matrix. */
        hur_nodiscard cu_device inline const Type *ddata() const noexcept { return value_.ddata(); }

        /** \brierf Return element vector of the Matrix. */
        hur_nodiscard cu_device inline Type *ddata() noexcept { return value_.ddata(); }

        /** \brierf Return element vector of the constant Matrix. */
        hur_nodiscard cu_host inline const Type *hdata() const noexcept { return value_.hdata(); }

        /** \brierf Return element vector of the Matrix. */
        hur_nodiscard cu_host inline Type *hdata() noexcept { return value_.hdata(); }

        /**
         * \brief Value array pointer.
         * The size must be NNZ_;
         */
        hur_nodiscard cu_dual inline cuBlockArray1<Type, size_type> &value() noexcept {
            return value_;
        }

        /**
         * \brief Value array pointer.
         * The size must be NNZ_;
         */
        hur_nodiscard cu_dual inline const cuBlockArray1<Type, size_type> &value() const noexcept {
            return value_;
        }

        cu_host inline void copyFromHost() { value_.copyFromHost(); }
        cu_host inline void copyFromHostAsync(cudaStream_t stream) {
            value_.copyFromHostAsync(stream);
        }

        cu_host inline void copyToHost() const { value_.copyToHost(); }
        cu_host inline void copyToHostAsync(cudaStream_t stream) const {
            value_.copyToHostAsync(stream);
        }
    };
} //  namespace OpenHurricane

#endif // CUDA_PARALLEL