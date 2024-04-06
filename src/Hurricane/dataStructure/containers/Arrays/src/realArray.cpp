/*!
 * \file realArray.cpp
 * \brief Main subroutines for real array.
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

#include "realArray.hpp"
#include "HurMPI.hpp"
#include "Pout.hpp"
#include <cmath>
#include <iomanip>

namespace OpenHurricane {
    template <> void Array<real>::replace(const int i, const realArray &l) {
        this->operator=(l);
    }

    template <> void realArrayArray::replace(const int i, const realArray &l) {

        for (size_type j = 0; j < this->size(); ++j) {
            this->operator[](j)[i] = l[j];
        }
    }

    template <> void Array<real>::replace(const int i, const real &c) {

        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c;
        }
    }

    template <> void realArrayArray::replace(const int i, const real &c) {

        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j)[i] = c;
        }
    }

    template <> void Array<real>::setComponent(const int d, const real &c) {

        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c;
        }
    }

    template <> void Array<real>::setComponent(const int d, const Array<real> &c) {

        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c[j];
        }
    }

    template <> void realArrayArray::setComponent(const int d, const real &c) {

        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j)[d] = c;
        }
    }

    template <> void Array<real>::setComponent(const real &c) {

        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c;
        }
    }

    template <> void realArrayArray::setComponent(const real &c) {

        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c;
        }
    }

    template <>
    void component(Array<typename Array<real>::elementType> &comp, const realArray &lf,
                   const int d) {
        comp = lf;
    }

    template <>
    void component(Array<typename realArrayArray::elementType> &comp, const Array<realArray> &lf,
                   const int d) {

        for (Array<realArray>::size_type i = 0; i < lf.size(); i++) {
            comp[i] = lf[i][d];
        }
    }

    template <> Array<typename Array<real>::elementType> Array<real>::component(const int d) const {
        return *this;
    }

    template <>
    Array<typename realArrayArray::elementType> realArrayArray::component(const int d) const {
        Array<typename realArrayArray::elementType> comp(this->size());

        for (size_type i = 0; i < this->size(); i++) {
            comp[i] = (*this)[i][d];
        }
        return comp;
    }

    template <> fileOsstream &Array<real>::writeToStream(fileOsstream &fos) const {

        return writeToStreamWithFactor(fos);
    }

    template <>
    fileOsstream &Array<real>::writeToStream(fileOsstream &fos, const size_type outSize) const {
        return writeToStreamWithFactor(fos, outSize);
    }

    template <> fileOsstream &Array<real>::writeToStreamWithFactor(fileOsstream &fos) const {
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            size_type i = 0;
            while (i < this->size()) {
                sstr << std::setprecision(feature<real>::precision) << this->operator[](i++) << " ";
                if (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << this->operator[](i++)
                         << " ";
                }
                if (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << this->operator[](i++)
                         << " ";
                }
                sstr << "\n";
            }
            fos.os() << sstr.str().c_str();
        } else {
            if (this->size()) {
                realArray newF;
                newF.transfer(*(this->clone()));
                fos.write(reinterpret_cast<const char *>(&newF[0]), this->byteSize());
            }
        }

        fos.unsetRealPrecision();
        return fos;
    }

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &Array<real>::writeToStreamWithFactor(fileOsstream &fos,
                                                       const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            size_type i = 0;
            while (i < outSize) {
                sstr << std::setprecision(feature<real>::precision) << this->operator[](i++) << " ";
                if (i < outSize) {
                    sstr << std::setprecision(feature<real>::precision) << this->operator[](i++)
                         << " ";
                }
                if (i < outSize) {
                    sstr << std::setprecision(feature<real>::precision) << this->operator[](i++)
                         << " ";
                }
                sstr << "\n";
            }
            fos.os() << sstr.str().c_str();
        } else {
            if (outSize) {
                realArray newF;
                newF.transfer(*(this->clone()));
                fos.write(reinterpret_cast<const char *>(&newF[0]), outSize * sizeof(real));
            }
        }

        fos.unsetRealPrecision();
        return fos;
    }

    template <>
    void Array<real>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        real minV = veryLarge;
        real maxV = -veryLarge;

        for (size_type i = 0; i < minSize; ++i) {
            minV = min(this->operator[](i), minV);
            maxV = max(this->operator[](i), maxV);
        }

        fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
    }
    template <>
    void Array<real>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                    const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        real minV = veryLarge;
        real maxV = -veryLarge;

        for (size_type i = 0; i < minSize; ++i) {
            minV = min(this->operator[](i), minV);
            maxV = max(this->operator[](i), maxV);
        }

        fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
    }

    template <>
    void Array<real>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                            const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        real minV = veryLarge;
        real maxV = -veryLarge;

        for (size_type i = 0; i < minSize; ++i) {
            minV = min(this->operator[](i), minV);
            maxV = max(this->operator[](i), maxV);
        }
        HurMPI::reduce(minV, MPI_MIN);
        HurMPI::reduce(maxV, MPI_MAX);
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
        }
    }
    template <> fileOsstream &realArrayArray::writeToStream(fileOsstream &fos) const {
        if (this->size()) {
            for (size_type d = 0; d < this->operator[](0).size(); ++d) {
                this->component(d).writeToStream(fos);
            }
        }
        return fos;
    }

    template <>
    fileOsstream &realArrayArray::writeToStream(fileOsstream &fos, const size_type outSize) const {
        if (this->size()) {
            for (size_type d = 0; d < this->operator[](0).size(); ++d) {
                this->component(d).writeToStream(fos, outSize);
            }
        }
        return fos;
    }

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &realArrayArray::writeToStreamWithFactor(fileOsstream &fos,
                                                          const size_type outSize) const {
        if (this->size()) {
            for (size_type d = 0; d < this->operator[](0).size(); ++d) {
                this->component(d).writeToStreamWithFactor(fos, outSize);
            }
        }
        return fos;
    }

    template <> fileOsstream &realArrayArray::writeToStreamWithFactor(fileOsstream &fos) const {
        if (this->size()) {
            for (size_type d = 0; d < this->operator[](0).size(); ++d) {
                this->component(d).writeToStreamWithFactor(fos);
            }
        }
        return fos;
    }

    hur_nodiscard realArray operator*(const realArray &f1, const realArray &f2) {
        checkArraysSize(f1, f2, "*");
        realArray f(f1.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] * f2[i];
        }

        return f;
    }

    hur_nodiscard realArray operator*(Array<real> &&f1, Array<real> &&f2) noexcept {
        checkArraysSize(f1, f2, "*");
        Array<real> tf(std::move(f1));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] *= f2[i];
        }
        return tf;
    }

    hur_nodiscard realArray operator*(Array<real> &&tft, const real t) noexcept {
        Array<real> tf(std::move(tft));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] *= t;
        }
        return tf;
    }

    hur_nodiscard realArray operator*(const real t, const realArray &f) {
        realArray tf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            tf[i] = t * f[i];
        }
        return tf;
    }

    hur_nodiscard realArray operator*(const realArray &f, const real &t) {
        realArray tf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            tf[i] = f[i] * t;
        }
        return tf;
    }

    hur_nodiscard realArray operator*(const real t, Array<real> &&tft) noexcept {
        Array<real> tf(std::move(tft));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] *= t;
        }
        return tf;
    }

    hur_nodiscard realArray operator/(const realArray &f1, realArray &&f2) noexcept {
        checkArraysSize(f1, f2, "/");
        Array<real> tf(std::move(f2));

        for (realArray::size_type i = 0; i < f1.size(); ++i) {
            tf[i] = f1[i] / tf[i];
        }
        return tf;
    }

    hur_nodiscard realArray operator/(realArray &&f1, realArray &&f2) noexcept {
        checkArraysSize(f1, f2, "/");
        Array<real> tf(std::move(f1));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] /= f2[i];
        }
        return tf;
    }

    hur_nodiscard realArray operator/(const real t, const realArray &f) {
        realArray tf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            tf[i] = t / f[i];
        }
        return tf;
    }

    hur_nodiscard realArray operator/(const real t, realArray &&f) noexcept {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = t / tf[i];
        }
        return tf;
    }

    hur_nodiscard realArray exp(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::exp(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray exp(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::sin(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray sin(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::sin(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray sin(realArray &&f) {
        realArray tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::sin(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray sinh(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::sinh(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray sinh(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::sinh(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray cos(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::cos(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray cos(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::cos(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray cosh(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::cosh(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray cosh(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::cosh(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray tan(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::tan(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray tan(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::tan(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray tanh(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::tanh(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray tanh(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::tanh(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray asin(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::asin(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray asin(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::asin(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray asinh(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::asinh(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray asinh(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::asinh(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray acos(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::acos(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray acos(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::acos(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray acosh(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::acosh(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray acosh(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::acosh(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray atan(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::atan(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray atan(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::atan(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray atan2(const realArray &y, const realArray &x) {
        realArray sf(x.size());

        for (realArray::size_type i = 0; i < x.size(); ++i) {
            sf[i] = std::atan2(y[i], x[i]);
        }
        return sf;
    }

    hur_nodiscard realArray atan2(realArray &&y, const realArray &x) {
        Array<real> ty(std::move(y));

        for (realArray::size_type i = 0; i < x.size(); ++i) {
            ty[i] = std::atan2(ty[i], x[i]);
        }
        return ty;
    }

    hur_nodiscard realArray atan2(const realArray &y, realArray &&x) {
        Array<real> tx(std::move(x));

        for (realArray::size_type i = 0; i < tx.size(); ++i) {
            tx[i] = std::atan2(y[i], tx[i]);
        }
        return tx;
    }

    hur_nodiscard realArray atan2(realArray &&y, realArray &&x) {
        Array<real> ty(std::move(y));

        for (realArray::size_type i = 0; i < x.size(); ++i) {
            ty[i] = std::atan2(ty[i], x[i]);
        }
        return ty;
    }

    hur_nodiscard realArray atanh(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::atanh(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray atanh(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::atanh(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray log(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::log(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray log(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::log(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray log10(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::log10(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray log10(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::log10(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray pow(const realArray &f, const real p) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::pow(f[i], p);
        }
        return sf;
    }

    hur_nodiscard realArray pow(realArray &&f, const real p) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::pow(tf[i], p);
        }
        return tf;
    }

    hur_nodiscard realArray pow(const realArray &f, const realArray &p) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::pow(f[i], p[i]);
        }
        return sf;
    }

    hur_nodiscard realArray pow(realArray &&f, const realArray &p) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::pow(tf[i], p[i]);
        }
        return tf;
    }

    hur_nodiscard realArray pow(const realArray &f, realArray &&p) {
        Array<real> tp(std::move(p));

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            tp[i] = std::pow(f[i], tp[i]);
        }
        return tp;
    }

    hur_nodiscard realArray pow(realArray &&f, realArray &&p) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::pow(tf[i], p[i]);
        }
        return tf;
    }

    hur_nodiscard realArray sqr(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = sqr(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray sqr(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = sqr(tf[i]);
        }
        return tf;
    }

    hur_nodiscard realArray sqrt(const realArray &f) {
        realArray sf(f.size());

        for (realArray::size_type i = 0; i < f.size(); ++i) {
            sf[i] = std::sqrt(f[i]);
        }
        return sf;
    }

    hur_nodiscard realArray sqrt(realArray &&f) {
        Array<real> tf(std::move(f));

        for (realArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = std::sqrt(tf[i]);
        }
        return tf;
    }

    template <>
    void Array<real>::writeAveToPout(fileOsstream &fos, const Array<real> &rhs,
                                     const Array<real> &cV, const size_type n, const size_type allN,
                                     const real &rhs0, const bool calRhs0,
                                     const bool modifyRhs0) const {
        real rhsAve = 0.0;
        real zRhs = 0.0;

        for (size_type i = 0; i < n; ++i) {
            zRhs = fabs(rhs[i]) / cV[i];
            rhsAve += sqr(zRhs);
        }
        realArray rhsAveL(HurMPI::getProcSize(), Zero);
        HurMPI::gather(&rhsAve, 1, feature<real>::MPIType, rhsAveL.data(), 1,
                       feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

        if (HurMPI::master()) {
            rhsAve = 0.0;
            for (size_type i = 0; i < rhsAveL.size(); i++) {
                rhsAve += rhsAveL[i];
            }
            rhsAve = sqrt(rhsAve / real(allN));

            if (calRhs0 || rhs0 == 0) {
                if (rhsAve < veryTiny) {
                    const_cast<real &>(rhs0) = 0.0;
                } else {
                    const_cast<real &>(rhs0) = rhsAve;
                    rhsAve = 1.0;
                }
            } else {
                if (modifyRhs0) {
                    real temRhs = rhsAve / rhs0;
                    if (rhs0 < veryTiny && rhsAve < veryTiny) {
                        temRhs = 0.0;
                    } else if (temRhs > 1e3) {
                        const_cast<real &>(rhs0) = rhsAve;
                        rhsAve = 1.0;
                    } else {
                        rhsAve = temRhs;
                    }
                } else {
                    if (rhs0 < veryTiny && rhsAve < veryTiny) {
                        // doing nothing
                    } else if (rhs0 < veryTiny && rhsAve > veryTiny) {
                        const_cast<real &>(rhs0) = rhsAve;
                        rhsAve = 1.0;
                    } else {
                        rhsAve /= rhs0;
                    }
                }
            }
        }
        if (HurMPI::master()) {
            std::cout.setf(std::ios::showpoint);
            std::cout << std::left << std::setfill(' ') << std::setw(12) << std::setprecision(5)
                      << rhsAve;
            std::cout.unsetf(std::ios::showpoint);

            if (fos.opened()) {
                fos.os() << '\t' << std::setprecision(5) << fabs(rhsAve);
            }
        }
        HurMPI::bcast(&rhsAve, 1, feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        setCurRhs(rhsAve);
    }
} // namespace OpenHurricane