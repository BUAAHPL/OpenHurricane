/*!
 * \file tensorArray.hpp
 * \brief Headers of the tensor Array.
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
#include "vectorArray.hpp"

namespace OpenHurricane {
    using tensorArray = Array<tensor>;
    using symmTensorArray = Array<symmTensor>;
    using diagTensorArray = Array<diagTensor>;
    using sphericalTensorArray = Array<sphericalTensor>;

    template <> fileOsstream &Array<tensor>::writeToStream(fileOsstream &fos) const;
    template <> fileOsstream &Array<symmTensor>::writeToStream(fileOsstream &fos) const;
    template <> fileOsstream &Array<diagTensor>::writeToStream(fileOsstream &fos) const;
    template <> fileOsstream &Array<sphericalTensor>::writeToStream(fileOsstream &fos) const;

    template <>
    fileOsstream &Array<tensor>::writeToStream(fileOsstream &fos, const size_type outSize) const;
    template <>
    fileOsstream &Array<symmTensor>::writeToStream(fileOsstream &fos,
                                                   const size_type outSize) const;
    template <>
    fileOsstream &Array<diagTensor>::writeToStream(fileOsstream &fos,
                                                   const size_type outSize) const;
    template <>
    fileOsstream &Array<sphericalTensor>::writeToStream(fileOsstream &fos,
                                                        const size_type outSize) const;

    template <> fileOsstream &Array<tensor>::writeToStreamWithFactor(fileOsstream &fos) const;
    template <> fileOsstream &Array<symmTensor>::writeToStreamWithFactor(fileOsstream &fos) const;
    template <> fileOsstream &Array<diagTensor>::writeToStreamWithFactor(fileOsstream &fos) const;
    template <>
    fileOsstream &Array<sphericalTensor>::writeToStreamWithFactor(fileOsstream &fos) const;

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &Array<tensor>::writeToStreamWithFactor(fileOsstream &fos,
                                                         const size_type outSize) const;
    template <>
    fileOsstream &Array<symmTensor>::writeToStreamWithFactor(fileOsstream &fos,
                                                             const size_type outSize) const;
    template <>
    fileOsstream &Array<diagTensor>::writeToStreamWithFactor(fileOsstream &fos,
                                                             const size_type outSize) const;
    template <>
    fileOsstream &Array<sphericalTensor>::writeToStreamWithFactor(fileOsstream &fos,
                                                                  const size_type outSize) const;

    template <>
    void Array<tensor>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                      const size_type outSize) const;
    template <>
    void Array<symmTensor>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                          const size_type outSize) const;
    template <>
    void Array<diagTensor>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                          const size_type outSize) const;
    template <>
    void Array<sphericalTensor>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                               const size_type outSize) const;

    template <>
    void Array<tensor>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                              const size_type outSize) const;
    template <>
    void Array<symmTensor>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                                  const size_type outSize) const;
    template <>
    void Array<diagTensor>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                                  const size_type outSize) const;
    template <>
    void
    Array<sphericalTensor>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                                  const size_type outSize) const;

    template <>
    void Array<tensor>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const;
    template <>
    void Array<symmTensor>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const;
    template <>
    void Array<diagTensor>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const;
    template <>
    void Array<sphericalTensor>::writeMinMaxToStream(fileOsstream &fos,
                                                     const size_type outSize) const;

    /*!\brief Hodge Dual operator (tensor -> vector).*/
    hur_nodiscard vectorArray operator*(const tensorArray &tf);
    hur_nodiscard tensorArray operator*(const vectorArray &vf);

    /*!\brief Outer product*/
    hur_nodiscard tensorArray operator&(const vectorArray &f1, const vectorArray &f2);

    /*!\brief Inner product*/
    hur_nodiscard tensorArray operator*(const tensorArray &f1, const tensorArray &f2);
    hur_nodiscard tensorArray operator*(const tensorArray &tf, const Identity<real> &i);
    hur_nodiscard tensorArray operator*(tensorArray &&f1, const tensorArray &f2);
    hur_nodiscard tensorArray operator*(const tensorArray &f1, tensorArray &&f2);
    hur_nodiscard tensorArray operator*(tensorArray &&f1, tensorArray &&f2);
    hur_nodiscard vectorArray operator*(const tensorArray &f1, const vectorArray &f2);
    hur_nodiscard vectorArray operator*(const tensorArray &f1, Array<vector> &&f2);
    hur_nodiscard vectorArray operator*(const vectorArray &f1, const tensorArray &f2);
    hur_nodiscard vectorArray operator*(Array<vector> &&f1, const tensorArray &f2);
    hur_nodiscard vectorArray operator/(const vectorArray &vf, const tensorArray &tf);
    hur_nodiscard vectorArray operator/(vectorArray &&f1, const tensorArray &f2);

    hur_nodiscard vectorArray operator*(const symmTensorArray &stf);
    hur_nodiscard tensorArray operator*(const symmTensorArray &st1f, const symmTensorArray &st2f);
    hur_nodiscard realArray operator&&(const symmTensorArray &st1f, const symmTensorArray &st2f);
    hur_nodiscard vectorArray operator*(const symmTensorArray &stf, const vectorArray &vf);
    hur_nodiscard vectorArray operator*(const symmTensorArray &stf, vectorArray &&vf);
    hur_nodiscard vectorArray operator*(const vectorArray &vf, const symmTensorArray &stf);
    hur_nodiscard vectorArray operator*(vectorArray &&vf, const symmTensorArray &stf);

    hur_nodiscard diagTensorArray operator*(const Array<diagTensor> &st1f,
                                            const Array<diagTensor> &st2f);
    hur_nodiscard vectorArray operator*(const Array<diagTensor> &stf, const Array<vector> &vf);
    hur_nodiscard vectorArray operator*(const Array<diagTensor> &stf, vectorArray &&vf);
    hur_nodiscard vectorArray operator*(const Array<vector> &vf, const Array<diagTensor> &stf);
    hur_nodiscard vectorArray operator*(vectorArray &&vf, const Array<diagTensor> &stf);

    hur_nodiscard sphericalTensorArray operator*(const sphericalTensorArray &st1f,
                                                 const sphericalTensorArray &st2f);
    hur_nodiscard vectorArray operator*(const sphericalTensorArray &stf, const vectorArray &vf);
    hur_nodiscard vectorArray operator*(const sphericalTensorArray &stf, vectorArray &&vf);
    hur_nodiscard vectorArray operator*(const vectorArray &vf, const sphericalTensorArray &stf);
    hur_nodiscard vectorArray operator*(vectorArray &&vf, const sphericalTensorArray &stf);

    /*!\brief Transpose for tensor.*/
    inline void transpose(tensorArray &tran, const tensorArray &lf) {
        checkArraysSize(tran, lf, "Transpose");
        for (typename tensorArray::size_type i = 0; i < tran.size(); ++i) {
            tran[i] = lf[i].transpose();
        }
    }

    /*!\brief Transpose for tensor.*/
    hur_nodiscard inline tensorArray transpose(const tensorArray &lf) {
        tensorArray tran(lf.size());
        transpose(tran, lf);
        return tran;
    }

    /*!\brief Transpose for tensor.*/
    hur_nodiscard inline tensorArray transpose(tensorArray &&lf) {
        tensorArray tf(std::move(lf));
        transpose(tf, tf);
        return tf;
    }

    /*!\brief Return the trace of a tensor.*/
    hur_nodiscard realArray tr(const tensorArray &tf);
    hur_nodiscard realArray tr(const symmTensorArray &stf);
    hur_nodiscard realArray tr(const Array<diagTensor> &stf);
    hur_nodiscard realArray tr(const sphericalTensorArray &stf);

    /*!\brief Return the skew-symmetric part of a tensor.*/
    hur_nodiscard tensorArray skew(const tensorArray &tf);

    /*!\brief Return the skew-symmetric part of a tensor.*/
    hur_nodiscard tensorArray skew(tensorArray &&tf);

    /*!\brief Return the determinant of a tensor.*/
    hur_nodiscard realArray det(const tensorArray &tf);
    hur_nodiscard realArray det(const symmTensorArray &stf);
    hur_nodiscard realArray det(const Array<diagTensor> &stf);
    hur_nodiscard realArray det(const sphericalTensorArray &stf);

    /*!\brief Return the cofactor tensor of a tensor.*/
    hur_nodiscard tensorArray cof(const tensorArray &tf);
    hur_nodiscard symmTensorArray cof(const symmTensorArray &stf);

    /*!\brief Return the cofactor tensor of a tensor.*/
    hur_nodiscard tensorArray cof(tensorArray &&tf);
    hur_nodiscard symmTensorArray cof(symmTensorArray &&stf);

    /*!\brief Return the inverse of a tensor.*/
    hur_nodiscard tensorArray inv(const tensorArray &tf);
    hur_nodiscard symmTensorArray inv(const symmTensorArray &stf);
    hur_nodiscard diagTensorArray inv(const Array<diagTensor> &stf);

    /*!\brief Return the inverse of a tensor.*/
    hur_nodiscard tensorArray inv(tensorArray &&tf);
    hur_nodiscard symmTensorArray inv(symmTensorArray &&stf);
    hur_nodiscard diagTensorArray inv(diagTensorArray &&stf);

    /*!\brief Return the 1st invariant of a tensor.*/
    hur_nodiscard realArray invariantI(const tensorArray &tf);
    hur_nodiscard realArray invariantI(const symmTensorArray &stf);
    hur_nodiscard realArray invariantI(const Array<diagTensor> &stf);
    hur_nodiscard realArray invariantI(const sphericalTensorArray &stf);

    /*!\brief Return the 2nd invariant of a tensor.*/
    hur_nodiscard realArray invariantII(const tensorArray &tf);
    hur_nodiscard realArray invariantII(const symmTensorArray &stf);

    /*!\brief Return the 3rd invariant of a tensor.*/
    hur_nodiscard realArray invariantIII(const tensorArray &tf);
    hur_nodiscard realArray invariantIII(const symmTensorArray &stf);
    hur_nodiscard realArray invariantIII(const Array<diagTensor> &stf);
    hur_nodiscard realArray invariantIII(const sphericalTensorArray &stf);

    /*!\brief Return the cofactor symmetric tensor of a symmetric tensor.*/
    hur_nodiscard diagTensorArray diag(const Array<tensor> &stf);

    /*!\brief Return the symmetric part of a tensor array.*/
    hur_nodiscard Array<symmTensor> symm(const tensorArray &tb);

    /*!\brief Return twice the symmetric part of a tensor array.*/
    hur_nodiscard Array<symmTensor> twoSymm(const tensorArray &tb);

    hur_nodiscard vectorArray diagToVector(const tensorArray &tb);

    /*!\brief Return magnitude of a skew-symmetric tensor array*/
    hur_nodiscard realArray skewMagnitude(const tensorArray &tb);

    hur_nodiscard tensorArray operator+(const Array<tensor> &tf, const Array<diagTensor> &df);
    hur_nodiscard tensorArray operator+(tensorArray &&tf, const Array<diagTensor> &df);
    hur_nodiscard tensorArray operator-(const Array<tensor> &tf, const Array<diagTensor> &df);
    hur_nodiscard tensorArray operator-(tensorArray &&tf, const Array<diagTensor> &df);
    hur_nodiscard Array<symmTensor> operator+(const Array<symmTensor> &tf,
                                              const Array<diagTensor> &df);
    hur_nodiscard Array<symmTensor> operator+(Array<symmTensor> &&tf, const Array<diagTensor> &df);
    hur_nodiscard Array<symmTensor> operator-(const Array<symmTensor> &tf,
                                              const Array<diagTensor> &df);
    hur_nodiscard Array<symmTensor> operator-(Array<symmTensor> &&tf, const Array<diagTensor> &df);

    hur_nodiscard tensorArray operator+(const tensorArray &tf, const sphericalTensorArray &df);
    hur_nodiscard tensorArray operator+(tensorArray &&tf, const sphericalTensorArray &df);
    hur_nodiscard tensorArray operator-(const tensorArray &tf, const sphericalTensorArray &df);
    hur_nodiscard tensorArray operator-(tensorArray &&tf, const sphericalTensorArray &df);
    hur_nodiscard Array<symmTensor> operator+(const Array<symmTensor> &tf,
                                              const sphericalTensorArray &df);
    hur_nodiscard Array<symmTensor> operator+(Array<symmTensor> &&tf,
                                              const sphericalTensorArray &df);
    hur_nodiscard Array<symmTensor> operator-(const Array<symmTensor> &tf,
                                              const sphericalTensorArray &df);
    hur_nodiscard Array<symmTensor> operator-(Array<symmTensor> &&tf,
                                              const sphericalTensorArray &df);
    hur_nodiscard sphericalTensorArray operator*(const realArray &sf,
                                                 const sphericalTensorArray &df);
    hur_nodiscard sphericalTensorArray operator*(const realArray &sf, sphericalTensorArray &&df);
    hur_nodiscard sphericalTensorArray operator*(const realArray &sf, const Identity<real> &ii);

} // namespace OpenHurricane