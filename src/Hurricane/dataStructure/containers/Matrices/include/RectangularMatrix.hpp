/*!
 * \file RectangularMatrix.hpp
 * \brief Header of rectangular matrix
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
#include "SquareMatrix.hpp"

namespace OpenHurricane {
    template <class Type> class RectangularMatrix : public Matrix<RectangularMatrix<Type>, Type> {
    public:
        using Base = Matrix<RectangularMatrix<Type>, Type>;
        using value_type = Type;
        using size_type = integer;
        using elementType = Type;

    public:

        /*!\brief Construct null.*/
        inline RectangularMatrix() : Base() {}

        /*!\brief Construct with m rows and n columns.*/
        inline RectangularMatrix(const integer m, const integer n) : Base(m, n) {}
        inline RectangularMatrix(const integer m, const integer n, const zero) : Base(m, n, Zero) {}
        inline RectangularMatrix &operator=(const zero) {
            Base::operator=(Zero);
            return *this;
        }
        inline RectangularMatrix(const integer m, const integer n, const Type &a) : Base(m, n, a) {}
        inline RectangularMatrix &operator=(const value_type &a) {
            Base::operator=(a);
            return *this;
        }

        /*!\brief Construct as copy.*/
        inline RectangularMatrix(const RectangularMatrix &other) : Base(other) {}
        inline RectangularMatrix &operator=(const RectangularMatrix &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        inline RectangularMatrix(const RectangularMatrix &&other) noexcept
            : Base(std::move(other)) {}
        inline RectangularMatrix &operator=(RectangularMatrix &&other) noexcept {
            Base::operator=(std::move(other));
            return *this;
        }

        hur_nodiscard inline uniquePtr<RectangularMatrix<Type>> clone() const {
            return uniquePtr<RectangularMatrix<Type>>(new RectangularMatrix<Type>(*this));
        }

        inline virtual ~RectangularMatrix() noexcept {}
    };

    template <class Type>
    class innerProductType<Type, RectangularMatrix<Type>, RectangularMatrix<Type>> {
    public:
        using type = RectangularMatrix<Type>;
    };

    template <class Type>
    class innerProductType<Type, RectangularMatrix<Type>, SquareMatrix<Type>> {
    public:
        using type = RectangularMatrix<Type>;
    };

    template <class Type>
    class innerProductType<Type, SquareMatrix<Type>, RectangularMatrix<Type>> {
    public:
        using type = RectangularMatrix<Type>;
    };

} // namespace OpenHurricane