/*!
 * \file Matrix.hpp
 * \brief Header of Matrix
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
#include "HurMPIBase.hpp"
#include "List.hpp"
#include "dataStructure.hpp"
#include <Eigen/Dense>

namespace OpenHurricane {

    /**
     * \brief The base class of matrix.
     * \tparam Matrices - The specific matrix type.
     * \tparam Type - The type of elements.
     */
    template <class Matrices, class Type>
    class Matrix : public Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> {
    public:
        using Base = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using mType = Matrix<Matrices, Type>;
        using value_type = Type;
        using size_type = integer;
        using elementType = Type;

    protected:
        // Private data

        /** \brief Number of rows and columns. */
        size_type nRows_, nColumns_;

        inline void allocate() {
            if (nRows_ > 0 && nColumns_ > 0) {
                Base::resize(nRows_, nColumns_);
            }
        }

    public:
#ifdef MPI_PARALLEL
        /** \brief MPI datatype. */
        static constexpr MPI_Datatype MPIType = feature<elementType>::MPIType;
#endif // MPI_PARALLEL

        /** \brief Return a null matrix. */
        hur_nodiscard inline static const mType &nullObject() {
            return NullRefObj::nullRef<Matrix<Matrices, Type>>();
        }

    public:
        /** \brief Construct null . */
        inline Matrix() : Base(), nRows_(0), nColumns_(0) {}

        /** \brief Construct with m rows and n columns. */
        inline Matrix(const integer m, const integer n) : Base(), nRows_(m), nColumns_(n) {
#ifdef HUR_DEBUG
            if (m < 0 || n < 0) {
                LFatal("invalid size");
            }
#endif // HUR_DEBUG
            allocate();
        }

        /** \brief Construct with m rows and n columns
         *  initializing all elements to zero.
         */
        inline Matrix(const integer m, const integer n, const zero)
            : Base(), nRows_(m), nColumns_(n) {
#ifdef HUR_DEBUG
            if (m < 0 || n < 0) {
                LFatal("invalid size");
            }
#endif // HUR_DEBUG
            allocate();
            Base::setZero();
        }

        /** \brief Construct with m rows and n columns
         * initializing all elements to the given value.
         */
        inline Matrix(const integer m, const integer n, const Type &v)
            : Base(), nRows_(m), nColumns_(n) {
#ifdef HUR_DEBUG
            if (m < 0 || n < 0) {
                LFatal("invalid size");
            }
#endif // HUR_DEBUG
            allocate();
            Base::setConstant(v);
        }

        /**
         * \brief Copy constructor..
         */
        inline Matrix(const mType &Mat) : Base(Mat), nRows_(Mat.nRows_), nColumns_(Mat.nColumns_) {}

        /** \brief Copy constructor. */
        inline Matrix(mType &&Mat) noexcept
            : Base(std::move(Mat)), nRows_(std::move(Mat.nRows_)),
              nColumns_(std::move(Mat.nColumns_)) {}

        /** \brief Copy constructor from matrix of different matrices. */
        template <class Matrices2>
        inline explicit Matrix(const Matrix<Matrices2, Type> &Mat)
            : Base(), nRows_(Mat.nRows_), nColumns_(Mat.nColumns_) {
            allocate();
            Base::operator=(Mat);
        }

        /** \brief Clone. */
        hur_nodiscard inline uniquePtr<mType> clone() const {
            return uniquePtr<Matrix<Matrices, Type>>(new Matrix<Matrices, Type>(*this));
        }

        /** \brief Destructor. */
        inline ~Matrix() noexcept { clear(); }

        // Member Functions

        // Access

        /** \brief Return the number of rows. */
        hur_nodiscard inline integer m() const noexcept { return nRows_; }

        /** \brief Return the number of rows. */
        hur_nodiscard inline integer nRows() const noexcept { return nRows_; }

        /** \brief Return the number of rows. */
        hur_nodiscard inline integer rows() const noexcept { return nRows_; }

        /** \brief Return the number of columns. */
        hur_nodiscard inline integer n() const noexcept { return nColumns_; }

        /** \brief Return the number of columns. */
        hur_nodiscard inline integer nColumns() const noexcept { return nColumns_; }

        /** \brief Return the number of columns. */
        hur_nodiscard inline integer cols() const noexcept { return nColumns_; }

        /** \brief Return the number of columns. */
        hur_nodiscard inline integer nCols() const noexcept { return nColumns_; }

        /** \brief Return the number of elements in matrix. */
        hur_nodiscard inline integer size() const noexcept { return nRows_ * nColumns_; }

        /** \brief Return element vector of the constant Matrix. */
        hur_nodiscard inline const Type *v() const { return Base::data(); }

        /** \brief Return element vector of the Matrix. */
        hur_nodiscard inline Type *v() { return Base::data(); }

        /** \brief Return element vector of the constant Matrix. */
        hur_nodiscard inline const Type *data() const { return Base::data(); }

        /** \brief Return element vector of the Matrix. */
        hur_nodiscard inline Type *data() { return Base::data(); }

        /** \brief Return i_th element of the matrix element vector. */
        hur_nodiscard inline Type &element(const integer i) {
#ifdef HUR_DEBUG
            if (i < 0 || i > size()) {
                LFatal("Index %d out of range 0,...%d", i, size() - 1);
            }
#endif // HUR_DEBUG
            return v()[i];
        }

        /** \brief Return i_th element of the matrix element vector. */
        hur_nodiscard inline const Type &element(const integer i) const {
#ifdef HUR_DEBUG
            if (i < 0 || i > size()) {
                LFatal("Index %d out of range 0,...%d", i, size() - 1);
            }
#endif // HUR_DEBUG
            return v()[i];
        }

        /** \brief Check index i is within valid range (0,..., m-1). */
        inline void checki(const integer i) const {
#ifdef HUR_DEBUG
            if (nRows_ <= 0 || nColumns_ <= 0) {
                LFatal("Attempt to access element from empty matrix");
            } else if (i < 0 || i >= nRows_) {
                LFatal("Index: %d out of range 0 ... %d", i, nRows_ - 1);
            }
#endif // HUR_DEBUG
        }

        /** \brief Check index j is within valid range (0,..., n-1). */
        inline void checkj(const integer j) const {
#ifdef HUR_DEBUG
            if (nRows_ <= 0 || nColumns_ <= 0) {
                LFatal("Attempt to access element from empty matrix");
            } else if (j < 0 || j >= nColumns_) {
                LFatal("Index: %d out of range 0 ... %d", 1, nColumns_ - 1);
            }
#endif // HUR_DEBUG
        }

        /** \brief Clear the Matrix. */
        inline void clear() noexcept {
            Base::resize(0, 0);
            nRows_ = 0;
            nColumns_ = 0;
        }

        /** \brief Transfer the contents of the argument Matrix into this Matrix
         *  and annul the argument Matrix.
         */
        inline void transfer(mType &mt) {
            if (this == std::addressof(mt)) {
                return;
            }
            clear();
            nRows_ = mt.nRows_;
            mt.nRows_ = 0;
            nColumns_ = mt.nColumns_;
            mt.nColumns_ = 0;

            mt.swap(*this);
            mt.clear();
        }

        /**
         * \brief Set the size of the matrix.
         * \param[in] m - The number of rows
         * \param[in] n - The number of columns
         */
        inline void resize(const integer m, const integer n) {
            mType newM(m, n, Zero);
            integer minM = min(m, nRows_);
            integer minN = min(n, nColumns_);

            for (integer i = 0; i < minM; ++i) {
                for (integer j = 0; j < minN; ++j) {
                    newM(i, j) = (*this)(i, j);
                }
            }

            transfer(newM);
        }

        /**
         * \brief Reset the size of the matrix.
         * \param[in] m - The number of rows and columns (Square matrix)
         */
        inline void resize(const integer m) { resize(m, m); }

        /**
         * \brief Shallow reset the size of the matrix.
         * \param[in] m - The number of rows
         * \param[in] n - The number of columns
         */
        inline void shallowResize(const integer m, const integer n) {
            nRows_ = m;
            nColumns_ = n;
        }

        /** \brief Return the transpose of the matrix. */
        hur_nodiscard inline Matrices transpose() const {
            const Matrix<Matrices, Type> &A = *this;
            Matrices At(n(), m());
            for (integer i = 0; i < m(); ++i) {
                for (integer j = 0; j < n(); ++j) {
                    At(j, i) = A(i, j);
                }
            }
            return At;
        }

        // Member Operators

        /** \brief Return subscript-checked row of Matrix. */
        hur_nodiscard inline Type *operator[](const integer i) {
            checki(i);
            return v() + i * nColumns_;
        }

        /** \brief Return subscript-checked row of Matrix. */
        hur_nodiscard inline const Type *operator[](const integer i) const {
            checki(i);
            return v() + i * nColumns_;
        }

        /** \brief (i, j) const element access operator. */
        hur_nodiscard inline const Type &operator()(const integer i, const integer j) const {
            checki(i);
            checkj(j);

            return v()[i * nColumns_ + j];
        }

        /** \brief (i, j) element access operator. */
        hur_nodiscard inline Type &operator()(const integer i, const integer j) {
            checki(i);
            checkj(j);

            return v()[i * nColumns_ + j];
        }

        /** \brief Assignment operator. */
        inline Matrix &operator=(const mType &M) {
            if (this != std::addressof(M)) {
                if (nRows_ != M.nRows_ || nColumns_ != M.nColumns_) {
                    clear();
                    nRows_ = M.nRows_;
                    nColumns_ = M.nColumns_;
                    allocate();
                }

                if (v() != nullptr) {
                    const integer mn = size();
                    for (integer i = 0; i < mn; ++i) {
                        v()[i] = M.v()[i];
                    }
                }
            }
            return *this;
        }

        /** \brief Assignment operator. */
        inline Matrix &operator=(mType &&M) noexcept {
            nRows_ = M.nRows_;
            nColumns_ = M.nColumns_;
            Base::operator=(std::move(static_cast<Base &&>(M)));
            M.nRows_ = 0;
            M.nColumns_ = 0;
            return *this;
        }

        /** \brief Assignment of all elements to zero. */
        inline Matrix &operator=(const zero) {
            if (size() != 0) {
                Base::setZero();
            }
            return *this;
        }

        /** \brief Assignment of all elements to the given value. */
        inline Matrix &operator=(const Type &a) {
            if (size() != 0) {
                Base::setConstant(a);
            }
            return *this;
        }
    };

    // Global Functions and operators

    /**
     * \brief Get the maximum element of the matrix M.
     */
    template <class Matrices, class Type> inline const Type &max(const Matrix<Matrices, Type> &M) {
        const integer mn = M.size();
        if (mn) {
            integer curMaxI = 0;
            const Type *Mv = M.v();
            for (integer i = 1; i < mn; ++i) {
                if (Mv[i] > Mv[curMaxI]) {
                    curMaxI = i;
                }
            }

            return Mv[curMaxI];
        } else {
            LFatal("Matrix is empty");
            return M(0, 0);
        }
    }

    /**
     * \brief Get the minimum element of the matrix M.
     */
    template <class Matrices, class Type> inline const Type &min(const Matrix<Matrices, Type> &M) {
        const integer mn = M.size();
        if (mn) {
            integer curMinI = 0;
            const Type *Mv = M.v();
            for (integer i = 1; i < mn; ++i) {
                if (Mv[i] < Mv[curMinI]) {
                    curMinI = i;
                }
            }

            return Mv[curMinI];
        } else {
            LFatal("Matrix is empty");
            return M(0, 0);
        }
    }

    /**
     * \brief Negative the matrix.
     */
    template <class Matrices, class Type>
    inline Matrices operator-(const Matrix<Matrices, Type> &M) {
        Matrices nM(M.m(), M.n());
        if (M.m() && M.n()) {
            Type *nMv = nM.v();
            const Type *Mv = M.v();

            const integer mn = M.size();
            for (integer i = 0; i < mn; ++i) {
                nMv[i] = -Mv[i];
            }
        }
        return nM;
    }

    /**
     * \brief Get (M1+M2).
     */
    template <class Matrices, class Type>
    inline Matrices operator+(const Matrix<Matrices, Type> &M1, const Matrix<Matrices, Type> &M2) {
        if (M1.m() != M2.m()) {
            std::string errMsg;
            errMsg = "Attempt to add matrices with different numbers of rows: ";
            errMsg += toString(M1.m());
            errMsg += ", ";
            errMsg += toString(M2.m());
            errorAbortStr(errMsg);
        } else if (M1.n() != M2.n()) {
            std::string errMsg;
            errMsg = "Attempt to add matrices with different numbers of columns: ";
            errMsg += toString(M1.n());
            errMsg += ", ";
            errMsg += toString(M2.n());
            errorAbortStr(errMsg);
        }

        Matrices M1M2(M1.m(), M1.n());

        Type *M1M2v = M1M2.v();
        const Type *M1v = M1.v();
        const Type *M2v = M2.v();

        const integer mn = M1.size();
        for (integer i = 0; i < mn; ++i) {
            M1M2v[i] = M1v[i] + M2v[i];
        }

        return M1M2;
    }

    /**
     * \brief Get (M1-M2).
     */
    template <class Matrices, class Type>
    inline Matrices operator-(const Matrix<Matrices, Type> &M1, const Matrix<Matrices, Type> &M2) {
        if (M1.m() != M2.m()) {
            std::string errMsg;
            errMsg = "Attempt to subtract matrices with different numbers of rows: ";
            errMsg += toString(M1.m());
            errMsg += ", ";
            errMsg += toString(M2.m());
            errorAbortStr(errMsg);
        } else if (M1.n() != M2.n()) {
            std::string errMsg;
            errMsg = "Attempt to subtract matrices with different numbers of columns: ";
            errMsg += toString(M1.n());
            errMsg += ", ";
            errMsg += toString(M2.n());
            errorAbortStr(errMsg);
        }

        Matrices M1M2(M1.m(), M1.n());

        Type *M1M2v = M1M2.v();
        const Type *M1v = M1.v();
        const Type *M2v = M2.v();

        const integer mn = M1.size();
        for (integer i = 0; i < mn; ++i) {
            M1M2v[i] = M1v[i] - M2v[i];
        }

        return M1M2;
    }

    template <class Matrices, class Type>
    inline Matrices operator*(const real s, const Matrix<Matrices, Type> &M) {
        Matrices sM(M.m(), M.n());
        if (M.m() && M.n()) {
            Type *sMv = sM.v();
            const Type *Mv = M.v();

            const integer mn = M.size();
            for (integer i = 0; i < mn; ++i) {
                sMv[i] = s * Mv[i];
            }
        }
        return sM;
    }

    template <class Matrices, class Type>
    inline Matrices operator*(const Matrix<Matrices, Type> &M, const real s) {
        Matrices sM(M.m(), M.n());
        if (M.m() && M.n()) {
            Type *sMv = sM.v();
            const Type *Mv = M.v();

            const integer mn = M.size();
            for (integer i = 0; i < mn; ++i) {
                sMv[i] = Mv[i] * s;
            }
        }
        return sM;
    }

    template <class Matrices, class Type>
    inline Matrices operator/(const Matrix<Matrices, Type> &M, const real s) {
        Matrices sM(M.m(), M.n());
        if (M.m() && M.n()) {
            Type *sMv = sM.v();
            const Type *Mv = M.v();

            const integer mn = M.size();
            for (integer i = 0; i < mn; ++i) {
                sMv[i] = Mv[i] / s;
            }
        }
        return sM;
    }

    template <class Matrices, class Type>
    inline Array<Type> operator*(const Matrix<Matrices, Type> &M, const Array<Type> &f) {
        if (M.n() != f.size()) {
            LFatal("Attempt to multiply a matrix and a field: Matrix(%d, %d), Array(%d)", M.m(),
                   M.n(), f.size());
        }

        Array<Type> MF(M.m(), Zero);
        for (integer i = 0; i < M.m(); ++i) {
            for (integer j = 0; j < M.n(); ++j) {
                MF[i] += M[i][j] * f[j];
            }
        }
        return MF;
    }

    template <class Matrices1, class Matrices2, class Type>
    hur_nodiscard inline typename innerProductType<Type, Matrices1, Matrices2>::type
    operator*(const Matrix<Matrices1, Type> &M1, const Matrix<Matrices2, Type> &M2) {
        if (M1.n() != M2.m()) {
            LFatal("Attempt to multiply two matrix: Matrix1(%d, %d), Matrix2(%d, %d)", M1.m(),
                   M1.n(), M2.m(), M2.n());
        }

        typename innerProductType<Type, Matrices1, Matrices2>::type M12(M1.m(), M2.n(), Zero);
        for (integer i = 0; i < M12.m(); ++i) {
            for (integer j = 0; j < M12.n(); ++j) {
                for (integer k = 0; k < M2.m(); ++k) {
                    M12(i, j) += M1(i, k) * M2(k, j);
                }
            }
        }
        return M12;
    }

} // namespace OpenHurricane
