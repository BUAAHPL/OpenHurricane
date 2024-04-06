/*!
 * \file SymmMatrix.hpp
 * \brief Header of symmetric matrix
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
#include "SquareMatrix.hpp"

namespace OpenHurricane {
    template <class Type> class SymmMatrix : public SquareMatrix<Type> {
    public:
        using Base = SquareMatrix<Type>;
        using value_type = typename Base::value_type;
        using size_type = typename Base::size_type;
        using elementType = typename Base::elementType;

    public:
        /*!\brief Construct null.*/
        inline SymmMatrix() : Base() {}

        /*!\brief Construct with m rows and m columns.*/
        inline SymmMatrix(const size_type m) : Base(m) {}

        /*!\brief Construct with m rows and m columns.*/
        inline SymmMatrix(const size_type m, const zero) : Base(m, Zero) {}
        inline SymmMatrix &operator=(const zero) {
            Base::operator=(Zero);
            return *this;
        }

        /*!\brief Construct with m rows and m columns.*/
        inline SymmMatrix(const size_type m, const value_type &t) : Base(m, t) {}
        inline SymmMatrix &operator=(const Type &a) {
            Base::operator=(a);
            return *this;
        }

        /*!\brief Construct as copy.*/
        inline SymmMatrix(const SymmMatrix &other) : Base(other) {}

        inline SymmMatrix &operator=(const SymmMatrix &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        inline SymmMatrix(SymmMatrix &&other) noexcept : Base(std::move(other)) {}
        inline SymmMatrix &operator=(SymmMatrix &&other) noexcept {
            Base::operator=(std::move(other));
            return *this;
        }
        /*!\brief Construct with n rows and n columns.*/
        inline SymmMatrix(const size_type n, const Identity<value_type>) : Base(n, Zero) {
            for (size_type i = 0; i < n; i++) {
                (*this)(i, i) = value_type(I);
            }
        }
        inline SymmMatrix &operator=(const Identity<Type>) {
            Base::operator=(Zero);
            for (size_type i = 0; i < this->n(); ++i) {
                this->operator()(i, i) = Type(I);
            }
            return *this;
        }

        hur_nodiscard inline uniquePtr<SymmMatrix<Type>> clone() const {
            return uniquePtr<SymmMatrix<Type>>(new SymmMatrix<Type>(*this));
        }
    };

    template <class Type> class innerProductType<Type, SymmMatrix<Type>, SymmMatrix<Type>> {
    public:
        using type = SquareMatrix<Type>;
    };

    template <class Type> class innerProductType<Type, SymmMatrix<Type>, SquareMatrix<Type>> {
    public:
        using type = SquareMatrix<Type>;
    };

    template <class Type> class innerProductType<Type, SquareMatrix<Type>, SymmMatrix<Type>> {
    public:
        using type = SquareMatrix<Type>;
    };

    /*!\brief Return the det for SymmMatrix.*/
    template <class Type>hur_nodiscard inline real det(const SymmMatrix<Type> &M) {
        return M.determinant();
    }

} // namespace OpenHurricane

