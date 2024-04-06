/*!
 * \file SquareMatrix.hpp
 * \brief Header of  square matrix
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
#include "Matrix.hpp"

namespace OpenHurricane {
    template <class Type> class SquareMatrix : public Matrix<SquareMatrix<Type>, Type> {
    public:
        using Base = Matrix<SquareMatrix<Type>, Type>;
        using value_type = Type;
        using size_type = integer;
        using elementType = Type;

    protected:
        /**
         * \brief The flag to specify that the square matrix is actually a diagonal matrix.
         * It must be set explicitly.
         */
        bool onlyDiagFillNonZero_;

    public:
        hur_nodiscard inline static const SquareMatrix &nullObject() {
            return NullRefObj::nullRef<SquareMatrix<Type>>();
        }

        /** \brief Construct null. */
        inline SquareMatrix() : Base(), onlyDiagFillNonZero_(false) {}

        /** \brief Construct with n rows and n columns. */
        inline SquareMatrix(const size_type n) : Base(n, n), onlyDiagFillNonZero_(false) {}
        inline SquareMatrix(const size_type m, const size_type n, const zero z = Zero)
            : Base(m, n, Zero), onlyDiagFillNonZero_(false) {
            if (n != m) {
                LFatal("Attempt to construct a square matrix %d x %d", m, n);
            }
        }
        inline SquareMatrix(const size_type n, const zero)
            : Base(n, n, Zero), onlyDiagFillNonZero_(false) {}
        inline SquareMatrix &operator=(const zero) {
            Base::operator=(Zero);
            return *this;
        }
        inline SquareMatrix(const size_type m, const size_type n, const Type &t)
            : Base(m, n, t), onlyDiagFillNonZero_(false) {
            if (n != m) {
                LFatal("Attempt to construct a square matrix %d x %d", m, n);
            }
        }
        inline SquareMatrix(const size_type n, const value_type &t)
            : Base(n, n, t), onlyDiagFillNonZero_(false) {}
        inline SquareMatrix &operator=(const value_type &t) {
            Base::operator=(t);
            return *this;
        }

        inline SquareMatrix(const size_type n, const Identity<value_type>)
            : Base(n, n, Zero), onlyDiagFillNonZero_(false) {
            for (size_type i = 0; i < n; i++) {
                (*this)(i, i) = value_type(I);
            }
        }

        inline SquareMatrix &operator=(const Identity<Type>) {
            Base::operator=(Zero);
            for (size_type i = 0; i < this->n(); ++i) {
                this->operator()(i, i) = Type(I);
            }
            return *this;
        }

        inline SquareMatrix(const SquareMatrix &SM) : Base(SM), onlyDiagFillNonZero_(false) {}
        inline SquareMatrix &operator=(const SquareMatrix &sm) {
            if (this != std::addressof(sm)) {
                Base::operator=(sm);
            }
            return *this;
        }

        inline SquareMatrix(SquareMatrix &&SM) noexcept
            : Base(std::move(SM)), onlyDiagFillNonZero_(false) {}

        inline SquareMatrix &operator=(SquareMatrix &&sm) noexcept {
            Base::operator=(std::move(sm));
            return *this;
        }

        hur_nodiscard inline uniquePtr<SquareMatrix<Type>> clone() const {
            return uniquePtr<SquareMatrix<Type>>(new SquareMatrix<Type>(*this));
        }

        /**
         * \brief Destructor.
         */
        inline virtual ~SquareMatrix() noexcept {}

        /** \brief Resize the square matrix. */
        inline void resize(const size_type m) { Base::resize(m, m); }

        /**
         * \brief The flag to specify that the square matrix is actually a diagonal matrix.
         * It must be set explicitly.
         */
        hur_nodiscard inline bool onlyDiagFillNonZero() const noexcept {
            return onlyDiagFillNonZero_;
        }

        inline void setOnlyDiagFillNonZero() noexcept { onlyDiagFillNonZero_ = true; }
        inline void unsetOnlyDiagFillNonZero() noexcept { onlyDiagFillNonZero_ = false; }

        hur_nodiscard inline SquareMatrix inv() const {
            if (onlyDiagFillNonZero_) {
                SquareMatrix<Type> iiv(*this);
                for (size_type i = 0; i < this->n(); ++i) {
                    iiv(i, i) = inv(this->operator()(i, i));
                }
                iiv.setOnlyDiagFillNonZero();
                return iiv;
            } else {
                auto iv = Base::Base::inverse();
                size_type nC = (size_type)iv.matrix().cols();
                size_type nR = (size_type)iv.matrix().rows();
                SquareMatrix<Type> iiv(nR);
                for (size_type i = 0; i < nR; ++i) {
                    for (size_type j = 0; j < nC; ++j) {
                        iiv(i, j) = iv(i, j);
                    }
                }
                return iiv;
            }
        }

        inline void inverseAndStore() {
            if (onlyDiagFillNonZero_) {
                for (size_type i = 0; i < this->n(); ++i) {
                    this->operator()(i, i) = inv(this->operator()(i, i));
                }
            } else {
                auto iv = Base::Base::colPivHouseholderQr().inverse();

                for (size_type i = 0; i < this->nRows(); ++i) {
                    for (size_type j = 0; j < this->nCols(); ++j) {
                        (*this)(i, j) = iv(i, j);
                    }
                }
            }
        }
    };

    /**
     * \brief Doolittle LU decompose with pivoting.
     *  \f[A = \left[ {\begin{array}{*{20}{c}}
        {{u_{11}}}&{{u_{12}}}& \cdots &{{u_{1n}}} \\
        {{l_{21}}}&{{u_{22}}}& \cdots &{{u_{2n}}} \\
         \vdots & \ddots & \ddots & \vdots  \\
        {{l_{n1}}}& \cdots &{{l_{n,n - 1}}}&{{u_{nn}}}
        \end{array}} \right]\f]
     */
    template <class Type>
    void LUDecomposeDPivoting(SquareMatrix<Type> &A, integerList &M, integer &sign) {
        integer n = A.n();
        sign = 1;
        if (n != M.size()) {
            M.resize(n);
        }

        List<Type> s(n);
        for (integer k = 0; k < n; ++k) {
            Type maxs = A(k, k);
            integer ik = k;
            for (integer i = k; i < n; ++i) {
                Type sumlu = Zero;
                for (integer t = 0; t < k; ++t) {
                    sumlu += A(i, t) * A(t, k);
                }
                s[i] = A(i, k) - sumlu;
                if (i == k) {
                    maxs = s[i];
                }
                if (mag(s[i]) > mag(maxs)) {
                    ik = i;
                    maxs = s[i];
                }
            }
            M[k] = ik;
            if (ik != k) {
                for (integer t = 0; t <= k - 1; ++t) {
                    Swap(A(k, t), A(ik, t));
                }
                for (integer t = k; t < n; ++t) {
                    Swap(A(k, t), A(ik, t));
                }
                Swap(s[k], s[ik]);
                sign *= -1;
            }

            A(k, k) = s[k];
            for (integer j = k + 1; j < n; ++j) {
                if (k < n) {
                    Type sumlu = Zero;
                    for (integer t = 0; t < k; ++t) {
                        sumlu += A(k, t) * A(t, j);
                    }
                    A(k, j) -= sumlu;
                }
            }
            for (integer i = k + 1; i < n; ++i) {
                if (k < n - 1) {
                    A(i, k) = s[i] / A(k, k);
                }
            }
        }
    }

    /** \brief Return the LU decomposed SquareMatrix det. */
    template <class Type>
    hur_nodiscard real detDecomposed(const SquareMatrix<Type> &M, const integer sign) {
        real d = 1.0;
        for (integer i = 0; i < M.m(); ++i) {
            d *= M(i, i);
        }
        d *= sign;
        return d;
    }

    /** \brief Returnthe det for SquareMatrix. */
    template <class Type> hur_nodiscard inline real det(const SquareMatrix<Type> &M) {
        return M.determinant();
    }

    template <class Type> class innerProductType<Type, SquareMatrix<Type>, SquareMatrix<Type>> {
    public:
        using type = SquareMatrix<Type>;
    };

    template <class Type>
    inline Array<Type> operator*(const SquareMatrix<Type> &M, const Array<Type> &f) {
        if (M.n() != f.size()) {
            LFatal("Attempt to multiply a matrix and a field: Matrix(%d, %d), Array(%d)", M.m(),
                   M.n(), f.size());
        }

        Array<Type> MF(M.m(), Zero);
        if (M.onlyDiagFillNonZero()) {
            for (integer i = 0; i < M.m(); ++i) {
                MF[i] += M(i, i) * f[i];
            }
        } else {
            for (integer i = 0; i < M.m(); ++i) {
                for (integer j = 0; j < M.n(); ++j) {
                    MF[i] += M[i][j] * f[j];
                }
            }
        }
        return MF;
    }

    template <class Type>
    hur_nodiscard inline SquareMatrix<Type> operator*(const SquareMatrix<Type> &M1,
                                                      const SquareMatrix<Type> &M2) {
        if (M1.n() != M2.m()) {
            LFatal("Attempt to multiply two matrix: Matrix1(%d, %d), Matrix2(%d, %d)", M1.m(),
                   M1.n(), M2.m(), M2.n());
        }
        SquareMatrix<Type> M12(M1.m(), Zero);
        if (M1.onlyDiagFillNonZero() && M2.onlyDiagFillNonZero()) {
            for (integer i = 0; i < M12.m(); ++i) {
                M12(i, i) = M1(i, i) * M2(i, i);
            }
            M12.setOnlyDiagFillNonZero();
        } else if (M1.onlyDiagFillNonZero()) {
            for (integer i = 0; i < M12.m(); ++i) {
                for (integer j = 0; j < M12.n(); ++j) {
                    M12(i, j) += M1(i, i) * M2(i, j);
                }
            }
        } else if (M2.onlyDiagFillNonZero()) {
            for (integer i = 0; i < M12.m(); ++i) {
                for (integer j = 0; j < M12.n(); ++j) {
                    M12(i, j) += M1(i, j) * M2(j, j);
                }
            }
        } else {
            for (integer i = 0; i < M12.m(); ++i) {
                for (integer j = 0; j < M12.n(); ++j) {
                    for (integer k = 0; k < M2.m(); ++k) {
                        M12(i, j) += M1(i, k) * M2(k, j);
                    }
                }
            }
        }
        return M12;
    }
} //  namespace OpenHurricane
