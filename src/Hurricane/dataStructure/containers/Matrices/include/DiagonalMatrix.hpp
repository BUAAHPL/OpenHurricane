/*!
 * \file DiagonalMatrix.hpp
 * \brief Header of  diagonal matrix.
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
#include "Matrix.hpp"
#include "RectangularMatrix.hpp"

namespace OpenHurricane {
    template <class Type> class DiagonalMatrix : public List<Type> {
    public:
        using Base = List<Type>;
        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;

    public:

        /*!\brief Construct null.*/
        inline DiagonalMatrix() : Base() {}

        /*!\brief Construct from Matrix.*/
        template <class Matrices>
        inline DiagonalMatrix(const Matrix<Matrices, Type> &m) : Base(min(m.m(), m.n())) {
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) = m(i, i);
            }
        }

        /*!\brief Construct with size.*/
        inline DiagonalMatrix(const size_type size) : Base(size) {}

        /*!\brief Construct with size and value.*/
        inline DiagonalMatrix(const size_type size, const value_type &t) : Base(size, t) {}
        inline DiagonalMatrix(const size_type size, const zero) : Base(size, Zero) {}

        inline DiagonalMatrix(const DiagonalMatrix &other) : Base(other) {}

        inline DiagonalMatrix(DiagonalMatrix &&other) noexcept : Base(std::move(other)) {}

        /*!\brief Destructor.*/
        inline ~DiagonalMatrix() noexcept {}

        using Base::operator=;

        inline DiagonalMatrix &operator=(const zero) {
            Base::operator=(Zero);
            return *this;
        }

        inline DiagonalMatrix &operator=(DiagonalMatrix &&other) noexcept {
            Base::operator=(std::move(other));
            return *this;
        }

        /*!\brief Invert and return itself.*/
        DiagonalMatrix &invert() {
            for (size_type i = 0; i < this->size(); ++i) {
                value_type x = this->operator[](i);
                if (mag(x) < veryTiny) {
                    this->operator[](i) = value_type(0);
                } else {
                    this->operator[](i) = value_type(1) / x;
                }
            }
            return *this;
        }
    };

    // Global Functions

    template <class Type> hur_nodiscard DiagonalMatrix<Type> inv(const DiagonalMatrix<Type> &A) {
        DiagonalMatrix<Type> invA(A.size());
        for (typename DiagonalMatrix<Type>::size_type i = 0; i < A.size(); ++i) {
            Type x = A[i];
            if (mag(x) < veryTiny) {
                invA[i] = Type(0);
            } else {
                invA[i] = Type(1) / x;
            }
        }
        return invA;
    }

} // namespace OpenHurricane

