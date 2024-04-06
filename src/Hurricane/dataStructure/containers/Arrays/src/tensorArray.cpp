/*!
 * \file tensorArray.cpp
 * \brief Main subroutines for tensor array.
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

#include "tensorArray.hpp"
#include "HurMPI.hpp"

namespace OpenHurricane {

    template <>
    fileOsstream &Array<tensor>::writeToStream(fileOsstream &fos) const {
        return writeToStreamWithFactor(fos);
    }

    template <>
    fileOsstream &Array<symmTensor>::writeToStream(fileOsstream &fos) const {

        return writeToStreamWithFactor(fos);
    }

    template <>
    fileOsstream &Array<diagTensor>::writeToStream(fileOsstream &fos) const {
        return writeToStreamWithFactor(fos);
    }

    template <>
    fileOsstream &
    sphericalTensorArray::writeToStream(fileOsstream &fos) const {

        return writeToStreamWithFactor(fos);
    }

    template <>
    fileOsstream &
    Array<tensor>::writeToStream(fileOsstream &fos, const size_type outSize) const {
        return writeToStreamWithFactor(fos, outSize);
    }

    template <>
    fileOsstream &
    Array<symmTensor>::writeToStream(fileOsstream &fos, const size_type outSize) const {
        return writeToStreamWithFactor(fos, outSize);
    }

    template <>
    fileOsstream &
    Array<diagTensor>::writeToStream(fileOsstream &fos, const size_type outSize) const {
        return writeToStreamWithFactor(fos, outSize);
    }
    template <>
    fileOsstream &
    sphericalTensorArray::writeToStream(fileOsstream &fos,
                                                   const size_type outSize) const {
        return writeToStreamWithFactor(fos, outSize);
    }

    template <>
    fileOsstream &
    Array<tensor>::writeToStreamWithFactor(fileOsstream &fos) const {
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            for (int j = 0; j < tensor::nElements_; j++) {
                sstr.clear();
                sstr.str("");
                size_type i = 0;
                while (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << this->operator[](i++)[j]
                         << " ";
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision)
                             << this->operator[](i++)[j] << " ";
                    }
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision)
                             << this->operator[](i++)[j] << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            }
        } else {
            if (this->size()) {
                for (int j = 0; j < tensor::nElements_; j++) {
                    Array<real> components = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&components[0]),
                              components.byteSize());
                }
            }
        }
        fos.unsetRealPrecision();
        return fos;
    }

    template <>
    fileOsstream &
    Array<symmTensor>::writeToStreamWithFactor(fileOsstream &fos) const {
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            for (int j = 0; j < symmTensor::nElements_; j++) {
                size_type i = 0;
                while (i < this->size()) {
                    fos.os() << this->operator[](i++)[j] << " ";
                    if (i < this->size()) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    if (i < this->size()) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    fos.os() << std::endl;
                }
            }
        } else {
            if (this->size()) {
                for (int j = 0; j < symmTensor::nElements_; j++) {
                    Array<real> components = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&components[0]),
                              components.byteSize());
                }
            }
        }
        fos.unsetRealPrecision();
        return fos;
    }

    template <>
    fileOsstream &
    Array<diagTensor>::writeToStreamWithFactor(fileOsstream &fos) const {
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            for (int j = 0; j < diagTensor::nElements_; j++) {
                size_type i = 0;
                while (i < this->size()) {
                    fos.os() << this->operator[](i++)[j] << " ";
                    if (i < this->size()) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    if (i < this->size()) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    fos.os() << std::endl;
                }
            }
        } else {
            if (this->size()) {
                for (int j = 0; j < diagTensor::nElements_; j++) {
                    Array<real> components = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&components[0]),
                              components.byteSize());
                }
            }
        }
        fos.unsetRealPrecision();
        return fos;
    }

    template <>
    fileOsstream &
    sphericalTensorArray::writeToStreamWithFactor(fileOsstream &fos) const {
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            for (int j = 0; j < sphericalTensor::nElements_; j++) {
                size_type i = 0;
                while (i < this->size()) {
                    fos.os() << this->operator[](i++)[j] << " ";
                    if (i < this->size()) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    if (i < this->size()) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    fos.os() << std::endl;
                }
            }
        } else {
            if (this->size()) {
                for (int j = 0; j < sphericalTensor::nElements_; j++) {
                    Array<real> components = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&components[0]),
                              components.byteSize());
                }
            }
        }
        fos.unsetRealPrecision();
        return fos;
    }

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &
    Array<tensor>::writeToStreamWithFactor(fileOsstream &fos,
                                                      const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            for (int j = 0; j < tensor::nElements_; j++) {
                sstr.clear();
                sstr.str("");
                size_type i = 0;
                while (i < minSize) {
                    sstr << std::setprecision(feature<real>::precision) << this->operator[](i++)[j]
                         << " ";
                    if (i < minSize) {
                        sstr << std::setprecision(feature<real>::precision)
                             << this->operator[](i++)[j] << " ";
                    }
                    if (i < minSize) {
                        sstr << std::setprecision(feature<real>::precision)
                             << this->operator[](i++)[j] << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            }
        } else {
            if (minSize) {
                for (int j = 0; j < tensor::nElements_; j++) {
                    Array<real> components = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&components[0]),
                              minSize * sizeof(real));
                }
            }
        }
        fos.unsetRealPrecision();
        return fos;
    }

    template <>
    fileOsstream &
    Array<symmTensor>::writeToStreamWithFactor(fileOsstream &fos,
                                                          const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            for (int j = 0; j < symmTensor::nElements_; j++) {
                size_type i = 0;
                while (i < minSize) {
                    fos.os() << this->operator[](i++)[j] << " ";
                    if (i < minSize) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    if (i < minSize) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    fos.os() << std::endl;
                }
            }
        } else {
            if (minSize) {
                for (int j = 0; j < symmTensor::nElements_; j++) {
                    Array<real> components = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&components[0]),
                              minSize * sizeof(real));
                }
            }
        }
        fos.unsetRealPrecision();
        return fos;
    }

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &
    Array<diagTensor>::writeToStreamWithFactor(fileOsstream &fos,
                                                          const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            for (int j = 0; j < diagTensor::nElements_; j++) {
                size_type i = 0;
                while (i < minSize) {
                    fos.os() << this->operator[](i++)[j] << " ";
                    if (i < minSize) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    if (i < minSize) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    fos.os() << std::endl;
                }
            }
        } else {
            if (minSize) {
                for (int j = 0; j < diagTensor::nElements_; j++) {
                    Array<real> components = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&components[0]),
                              minSize * sizeof(real));
                }
            }
        }
        fos.unsetRealPrecision();
        return fos;
    }

    template <>
    fileOsstream &
    sphericalTensorArray::writeToStreamWithFactor(fileOsstream &fos,
                                                             const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            for (int j = 0; j < sphericalTensor::nElements_; j++) {
                size_type i = 0;
                while (i < minSize) {
                    fos.os() << this->operator[](i++)[j] << " ";
                    if (i < minSize) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    if (i < minSize) {
                        fos.os() << this->operator[](i++)[j] << " ";
                    }
                    fos.os() << std::endl;
                }
            }
        } else {
            if (minSize) {
                for (int j = 0; j < sphericalTensor::nElements_; j++) {
                    Array<real> components = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&components[0]),
                              minSize * sizeof(real));
                }
            }
        }
        fos.unsetRealPrecision();
        return fos;
    }

    template <>
    void Array<tensor>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                                 const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < tensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
        }
    }

    template <>
    void
    Array<symmTensor>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                                const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < symmTensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
        }
    }

    template <>
    void
    Array<diagTensor>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                                const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < diagTensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
        }
    }

    template <>
    void
    sphericalTensorArray::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                                   const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < sphericalTensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
        }
    }

    template <>
    void
    Array<tensor>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                                    const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < tensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            HurMPI::reduce(minV, MPI_MIN);
            HurMPI::reduce(maxV, MPI_MAX);
            if (HurMPI::master()) {
                fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
                fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
            }
        }
    }

    template <>
    void Array<symmTensor>::writeMinMaxToStreamWithFactorByMaster(
        fileOsstream &fos, const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < symmTensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            HurMPI::reduce(minV, MPI_MIN);
            HurMPI::reduce(maxV, MPI_MAX);
            if (HurMPI::master()) {
                fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
                fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
            }
        }
    }

    template <>
    void Array<diagTensor>::writeMinMaxToStreamWithFactorByMaster(
        fileOsstream &fos, const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < diagTensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            HurMPI::reduce(minV, MPI_MIN);
            HurMPI::reduce(maxV, MPI_MAX);
            if (HurMPI::master()) {
                fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
                fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
            }
        }
    }

    template <>
    void sphericalTensorArray::writeMinMaxToStreamWithFactorByMaster(
        fileOsstream &fos, const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < sphericalTensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            HurMPI::reduce(minV, MPI_MIN);
            HurMPI::reduce(maxV, MPI_MAX);
            if (HurMPI::master()) {
                fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
                fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
            }
        }
    }

    template <>
    void Array<tensor>::writeMinMaxToStream(fileOsstream &fos,
                                                       const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < tensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
        }
    }

    template <>
    void Array<symmTensor>::writeMinMaxToStream(fileOsstream &fos,
                                                           const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < symmTensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
        }
    }

    template <>
    void Array<diagTensor>::writeMinMaxToStream(fileOsstream &fos,
                                                           const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < diagTensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
        }
    }

    template <>
    void sphericalTensorArray::writeMinMaxToStream(fileOsstream &fos,
                                                              const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < sphericalTensor::nElements_; j++) {
            real minV = veryLarge;
            real maxV = -veryLarge;
            for (size_type i = 0; i < minSize; ++i) {
                minV = min(this->operator[](i)[j], minV);
                maxV = max(this->operator[](i)[j], maxV);
            }
            fos.write(reinterpret_cast<const char *>(&minV), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxV), sizeof(double));
        }
    }

    hur_nodiscard vectorArray operator*(const tensorArray &tf) {
        vectorArray f(tf.size());
        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = *tf[i];
        }
        return f;
    }

    hur_nodiscard tensorArray operator*(const vectorArray &vf) {
        tensorArray f(vf.size());
        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = *vf[i];
        }
        return f;
    }

    hur_nodiscard tensorArray operator&(const vectorArray &f1, const vectorArray &f2) {
        checkArraysSize(f1, f2, "&");
        tensorArray f(f1.size());
        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] & f2[i];
        }

        return f;
    }

    hur_nodiscard tensorArray operator*(const tensorArray &f1, const tensorArray &f2) {
        checkArraysSize(f1, f2, "*");
        tensorArray f(f1.size());

        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] * f2[i];
        }

        return f;
    }

    hur_nodiscard tensorArray operator*(const tensorArray &tf, const Identity<real> &i) {
        tensorArray f(tf.size());

        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = tf[i];
        }

        return f;
    }

    hur_nodiscard tensorArray operator*(tensorArray &&f1, const tensorArray &f2) {
        checkArraysSize(f1, f2, "*");
        tensorArray tff(std::move(f1));

        for (tensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i] *= f2[i];
        }

        return tff;
    }

    hur_nodiscard tensorArray operator*(const tensorArray &f1, tensorArray &&f2) {
        checkArraysSize(f1, f2, "*");
        tensorArray tff(std::move(f2));

        for (tensorArray::size_type i = 0; i < f1.size(); ++i) {
            tff[i] = f1[i] * tff[i];
        }

        return tff;
    }

    hur_nodiscard tensorArray operator*(tensorArray &&f1, tensorArray &&f2) {
        checkArraysSize(f1, f2, "*");
        tensorArray tff(std::move(f1));

        for (tensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i] *= f2[i];
        }

        return tff;
    }

    hur_nodiscard vectorArray operator*(const tensorArray &f1, const vectorArray &f2) {
        checkArraysSize(f1, f2, "*");
        vectorArray f(f1.size());

        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] * f2[i];
        }

        return f;
    }

    hur_nodiscard vectorArray operator*(const tensorArray &f1, Array<vector> &&f2) {
        checkArraysSize(f1, f2, "*");
        Array<vector> tf(std::move(f2));

        for (tensorArray::size_type i = 0; i < f1.size(); ++i) {
            tf[i] = f1[i] * tf[i];
        }
        return tf;
    }

    hur_nodiscard vectorArray operator*(const vectorArray &f1, const tensorArray &f2) {
        checkArraysSize(f1, f2, "*");
        vectorArray f(f1.size());

        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] * f2[i];
        }

        return f;
    }

    hur_nodiscard vectorArray operator*(Array<vector> &&f1, const tensorArray &f2) {
        checkArraysSize(f1, f2, "*");
        vectorArray tf(std::move(f1));

        for (tensorArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = tf[i] * f2[i];
        }
        return tf;
    }

    hur_nodiscard vectorArray operator/(const vectorArray &f1, const tensorArray &f2) {
        checkArraysSize(f1, f2, "*");
        vectorArray f(f1.size());

        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] / f2[i];
        }

        return f;
    }

    hur_nodiscard vectorArray operator/(vectorArray &&f1, const tensorArray &f2) {
        checkArraysSize(f1, f2, "*");
        vectorArray tf(std::move(f1));

        for (tensorArray::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = tf[i] / f2[i];
        }
        return tf;
    }

    hur_nodiscard vectorArray operator*(const symmTensorArray &stf) {
        vectorArray f(stf.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = *stf[i];
        }
        return f;
    }

    hur_nodiscard tensorArray operator*(const symmTensorArray &st1f, const symmTensorArray &st2f) {
        checkArraysSize(st1f, st2f, "*");
        tensorArray f(st1f.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = st1f[i] * st2f[i];
        }
        return f;
    }

    hur_nodiscard realArray operator&&(const symmTensorArray &st1f, const symmTensorArray &st2f) {
        checkArraysSize(st1f, st2f, "*");
        realArray f(st1f.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = st1f[i] && st2f[i];
        }
        return f;
    }

    hur_nodiscard vectorArray operator*(const symmTensorArray &stf, const vectorArray &vf) {
        checkArraysSize(stf, vf, "*");
        vectorArray f(stf.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = stf[i] * vf[i];
        }
        return f;
    }

    hur_nodiscard vectorArray operator*(const symmTensorArray &stf, vectorArray &&vf) {
        checkArraysSize(stf, vf, "*");
        vectorArray tvf(std::move(vf));

        for (symmTensorArray::size_type i = 0; i < tvf.size(); ++i) {
            tvf[i] = stf[i] * tvf[i];
        }
        return tvf;
    }

    hur_nodiscard vectorArray operator*(const vectorArray &vf, const symmTensorArray &stf) {
        checkArraysSize(vf, stf, "*");
        vectorArray f(stf.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = vf[i] * stf[i];
        }
        return f;
    }

    hur_nodiscard vectorArray operator*(vectorArray &&vf, const symmTensorArray &stf) {
        checkArraysSize(vf, stf, "*");
        vectorArray tvf(std::move(vf));

        for (symmTensorArray::size_type i = 0; i < tvf.size(); ++i) {
            tvf[i] = tvf[i] * stf[i];
        }
        return tvf;
    }

    hur_nodiscard diagTensorArray operator*(const Array<diagTensor> &st1f,
                                            const Array<diagTensor> &st2f) {
        checkArraysSize(st1f, st2f, "*");
        diagTensorArray f(st1f.size());

        for (Array<diagTensor>::size_type i = 0; i < f.size(); ++i) {
            f[i] = st1f[i] * st2f[i];
        }
        return f;
    }

    hur_nodiscard vectorArray operator*(const Array<diagTensor> &stf, const Array<vector> &vf) {
        checkArraysSize(stf, vf, "*");
        vectorArray f(stf.size());

        for (Array<diagTensor>::size_type i = 0; i < f.size(); ++i) {
            f[i] = stf[i] * vf[i];
        }
        return f;
    }

    hur_nodiscard vectorArray operator*(const Array<diagTensor> &stf, vectorArray &&vf) {
        checkArraysSize(stf, vf, "*");
        vectorArray tvf(std::move(vf));

        for (Array<diagTensor>::size_type i = 0; i < tvf.size(); ++i) {
            tvf[i] = stf[i] * tvf[i];
        }
        return tvf;
    }

    hur_nodiscard vectorArray operator*(const Array<vector> &vf, const Array<diagTensor> &stf) {
        checkArraysSize(vf, stf, "*");
        vectorArray f(stf.size());
        for (Array<diagTensor>::size_type i = 0; i < f.size(); ++i) {
            f[i] = vf[i] * stf[i];
        }
        return f;
    }

    hur_nodiscard vectorArray operator*(vectorArray &&vf, const Array<diagTensor> &stf) {
        checkArraysSize(vf, stf, "*");
        vectorArray tvf(std::move(vf));
        for (Array<diagTensor>::size_type i = 0; i < tvf.size(); ++i) {
            tvf[i] = tvf[i] * stf[i];
        }
        return tvf;
    }

    hur_nodiscard sphericalTensorArray operator*(const sphericalTensorArray &st1f,
                                                 const sphericalTensorArray &st2f) {
        checkArraysSize(st1f, st2f, "*");
        sphericalTensorArray f(st1f.size());

        for (sphericalTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = st1f[i] * st2f[i];
        }
        return f;
    }

    hur_nodiscard vectorArray operator*(const sphericalTensorArray &stf, const vectorArray &vf) {
        checkArraysSize(stf, vf, "*");
        vectorArray f(stf.size());

        for (sphericalTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = stf[i] * vf[i];
        }
        return f;
    }

    hur_nodiscard vectorArray operator*(const sphericalTensorArray &stf, vectorArray &&vf) {
        checkArraysSize(stf, vf, "*");
        vectorArray tvf(std::move(vf));

        for (sphericalTensorArray::size_type i = 0; i < tvf.size(); ++i) {
            tvf[i] = stf[i] * tvf[i];
        }
        return tvf;
    }

    hur_nodiscard vectorArray operator*(const vectorArray &vf, const sphericalTensorArray &stf) {
        checkArraysSize(vf, stf, "*");
        vectorArray f(stf.size());

        for (sphericalTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = vf[i] * stf[i];
        }
        return f;
    }

    hur_nodiscard vectorArray operator*(vectorArray &&vf, const sphericalTensorArray &stf) {
        checkArraysSize(vf, stf, "*");
        vectorArray tvf(std::move(vf));

        for (sphericalTensorArray::size_type i = 0; i < tvf.size(); ++i) {
            tvf[i] = tvf[i] * stf[i];
        }
        return tvf;
    }

    hur_nodiscard realArray tr(const tensorArray &tf) {
        realArray f(tf.size());
        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = tr(tf[i]);
        }
        return f;
    }

    hur_nodiscard realArray tr(const symmTensorArray &stf) {
        realArray f(stf.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = tr(stf[i]);
        }
        return f;
    }

    hur_nodiscard realArray tr(const Array<diagTensor> &stf) {
        realArray f(stf.size());
        for (Array<diagTensor>::size_type i = 0; i < f.size(); ++i) {
            f[i] = tr(stf[i]);
        }
        return f;
    }

    hur_nodiscard realArray tr(const sphericalTensorArray &stf) {
        realArray f(stf.size());
        for (sphericalTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = tr(stf[i]);
        }
        return f;
    }

    hur_nodiscard tensorArray skew(const tensorArray &tf) {
        tensorArray f(tf.size());

        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = skew(tf[i]);
        }
        return f;
    }

    hur_nodiscard tensorArray skew(tensorArray &&tf) {
        tensorArray tff(std::move(tf));

        for (tensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = skew(tff[i]);
        }
        return tff;
    }

    hur_nodiscard realArray det(const tensorArray &tf) {
        realArray f(tf.size());
        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = det(tf[i]);
        }
        return f;
    }

    hur_nodiscard realArray det(const symmTensorArray &stf) {
        realArray f(stf.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = det(stf[i]);
        }
        return f;
    }

    hur_nodiscard realArray det(const Array<diagTensor> &stf) {
        realArray f(stf.size());
        for (Array<diagTensor>::size_type i = 0; i < f.size(); ++i) {
            f[i] = det(stf[i]);
        }
        return f;
    }

    hur_nodiscard realArray det(const sphericalTensorArray &stf) {
        realArray f(stf.size());
        for (sphericalTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = det(stf[i]);
        }
        return f;
    }

    hur_nodiscard tensorArray cof(const tensorArray &tf) {
        tensorArray f(tf.size());

        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = cof(tf[i]);
        }
        return f;
    }
    hur_nodiscard symmTensorArray cof(const symmTensorArray &stf) {
        symmTensorArray f(stf.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = cof(stf[i]);
        }
        return f;
    }

    hur_nodiscard tensorArray cof(tensorArray &&tf) {
        tensorArray tff(std::move(tf));

        for (tensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = cof(tff[i]);
        }
        return tff;
    }

    hur_nodiscard symmTensorArray cof(symmTensorArray &&stf) {
        symmTensorArray sttf(std::move(stf));

        for (symmTensorArray::size_type i = 0; i < sttf.size(); ++i) {
            sttf[i] = cof(sttf[i]);
        }
        return sttf;
    }

    hur_nodiscard tensorArray inv(const tensorArray &tf) {
        tensorArray f(tf.size());

        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = inv(tf[i]);
        }
        return f;
    }
    hur_nodiscard symmTensorArray inv(const symmTensorArray &stf) {
        symmTensorArray f(stf.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = inv(stf[i]);
        }
        return f;
    }

    hur_nodiscard diagTensorArray inv(const Array<diagTensor> &stf) {
        diagTensorArray f(stf.size());
        for (Array<diagTensor>::size_type i = 0; i < f.size(); ++i) {
            f[i] = inv(stf[i]);
        }
        return f;
    }

    hur_nodiscard tensorArray inv(tensorArray &&tf) {
        tensorArray tff(std::move(tf));
        for (tensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = inv(tff[i]);
        }
        return tff;
    }

    hur_nodiscard symmTensorArray inv(symmTensorArray &&stf) {
        symmTensorArray sttf(std::move(stf));

        for (symmTensorArray::size_type i = 0; i < sttf.size(); ++i) {
            sttf[i] = inv(sttf[i]);
        }
        return sttf;
    }

    hur_nodiscard diagTensorArray inv(diagTensorArray &&stf) {
        diagTensorArray sttf(std::move(stf));
        for (Array<diagTensor>::size_type i = 0; i < sttf.size(); ++i) {
            sttf[i] = inv(sttf[i]);
        }
        return sttf;
    }

    hur_nodiscard realArray invariantI(const tensorArray &tf) {
        return tr(tf);
    }

    hur_nodiscard realArray invariantI(const symmTensorArray &stf) {
        return tr(stf);
    }

    hur_nodiscard realArray invariantI(const Array<diagTensor> &stf) {
        return tr(stf);
    }

    hur_nodiscard realArray invariantI(const sphericalTensorArray &stf) {
        return tr(stf);
    }

    hur_nodiscard realArray invariantII(const tensorArray &tf) {
        realArray f(tf.size());
        for (tensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = invariantII(tf[i]);
        }
        return f;
    }

    hur_nodiscard realArray invariantII(const symmTensorArray &stf) {
        realArray f(stf.size());

        for (symmTensorArray::size_type i = 0; i < f.size(); ++i) {
            f[i] = invariantII(stf[i]);
        }
        return f;
    }

    hur_nodiscard realArray invariantIII(const tensorArray &tf) {
        return det(tf);
    }

    hur_nodiscard realArray invariantIII(const symmTensorArray &stf) {
        return det(stf);
    }

    hur_nodiscard realArray invariantIII(const Array<diagTensor> &stf) {
        return det(stf);
    }

    hur_nodiscard realArray invariantIII(const sphericalTensorArray &stf) {
        return det(stf);
    }

    hur_nodiscard diagTensorArray diag(const Array<tensor> &stf) {
        diagTensorArray df(stf.size());

        for (Array<diagTensor>::size_type j = 0; j < df.size(); j++) {
            df[j].xx() = stf[j].xx();
            df[j].yy() = stf[j].yy();
            df[j].zz() = stf[j].zz();
        }
        return df;
    }

    hur_nodiscard Array<symmTensor> symm(const tensorArray &tb) {
        Array<symmTensor> sy(tb.size());

        for (tensorArray::size_type i = 0; i < sy.size(); ++i) {
            sy[i] = symm(tb[i]);
        }
        return sy;
    }

    hur_nodiscard Array<symmTensor> twoSymm(const tensorArray &tb) {
        Array<symmTensor> sy(tb.size());

        for (tensorArray::size_type i = 0; i < sy.size(); ++i) {
            sy[i] = twoSymm(tb[i]);
        }
        return sy;
    }

    hur_nodiscard vectorArray diagToVector(const tensorArray &tb) {
        vectorArray sy(tb.size());

        for (tensorArray::size_type i = 0; i < sy.size(); ++i) {
            sy[i] = diagToVector(tb[i]);
        }
        return sy;
    }

    hur_nodiscard realArray skewMagnitude(const tensorArray &tb) {
        realArray sy(tb.size());
        for (tensorArray::size_type i = 0; i < sy.size(); ++i) {
            sy[i] = skewMagnitude(tb[i]);
        }
        return sy;
    }

    hur_nodiscard tensorArray operator+(const Array<tensor> &tf, const Array<diagTensor> &df) {
        checkArraysSize(tf, df, "+");
        tensorArray rf(tf.size());

        for (Array<diagTensor>::size_type i = 0; i < rf.size(); ++i) {
            rf[i] = tf[i] + df[i];
        }
        return rf;
    }

    hur_nodiscard tensorArray operator+(tensorArray &&tf, const Array<diagTensor> &df) {
        checkArraysSize(tf, df, "+");
        tensorArray tff(std::move(tf));

        for (Array<diagTensor>::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = tff[i] + df[i];
        }
        return tff;
    }

    hur_nodiscard tensorArray operator-(const Array<tensor> &tf, const Array<diagTensor> &df) {
        checkArraysSize(tf, df, "-");
        tensorArray rf(tf.size());

        for (Array<diagTensor>::size_type i = 0; i < rf.size(); ++i) {
            rf[i] = tf[i] - df[i];
        }
        return rf;
    }
    hur_nodiscard tensorArray operator-(tensorArray &&tf, const Array<diagTensor> &df) {
        checkArraysSize(tf, df, "-");
        tensorArray tff(std::move(tf));

        for (Array<diagTensor>::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = tff[i] - df[i];
        }
        return tff;
    }

    hur_nodiscard Array<symmTensor> operator+(const Array<symmTensor> &tf,
                                              const Array<diagTensor> &df) {
        checkArraysSize(tf, df, "+");
        Array<symmTensor> rf(tf.size());

        for (Array<diagTensor>::size_type i = 0; i < rf.size(); ++i) {
            rf[i] = tf[i] + df[i];
        }
        return rf;
    }
    hur_nodiscard Array<symmTensor> operator+(Array<symmTensor> &&tf, const Array<diagTensor> &df) {
        checkArraysSize(tf, df, "+");
        Array<symmTensor> tff(std::move(tf));

        for (Array<diagTensor>::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = tff[i] + df[i];
        }
        return tff;
    }

    hur_nodiscard Array<symmTensor> operator-(const Array<symmTensor> &tf,
                                              const Array<diagTensor> &df) {
        checkArraysSize(tf, df, "-");
        Array<symmTensor> rf(tf.size());

        for (Array<diagTensor>::size_type i = 0; i < rf.size(); ++i) {
            rf[i] = tf[i] - df[i];
        }
        return rf;
    }
    hur_nodiscard Array<symmTensor> operator-(Array<symmTensor> &&tf, const Array<diagTensor> &df) {
        checkArraysSize(tf, df, "-");
        Array<symmTensor> tff(std::move(tf));

        for (Array<diagTensor>::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = tff[i] - df[i];
        }
        return tff;
    }

    hur_nodiscard tensorArray operator+(const tensorArray &tf, const sphericalTensorArray &df) {
        checkArraysSize(tf, df, "+");
        tensorArray rf(tf.size());

        for (sphericalTensorArray::size_type i = 0; i < rf.size(); ++i) {
            rf[i] = tf[i] + df[i];
        }
        return rf;
    }

    hur_nodiscard tensorArray operator+(tensorArray &&tf, const sphericalTensorArray &df) {
        checkArraysSize(tf, df, "+");
        tensorArray tff(std::move(tf));

        for (sphericalTensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = tff[i] + df[i];
        }
        return tff;
    }

    hur_nodiscard tensorArray operator-(const tensorArray &tf, const sphericalTensorArray &df) {
        checkArraysSize(tf, df, "-");
        tensorArray rf(tf.size());

        for (sphericalTensorArray::size_type i = 0; i < rf.size(); ++i) {
            rf[i] = tf[i] - df[i];
        }
        return rf;
    }

    hur_nodiscard tensorArray operator-(tensorArray &&tf, const sphericalTensorArray &df) {
        checkArraysSize(tf, df, "-");
        tensorArray tff(std::move(tf));

        for (sphericalTensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = tff[i] - df[i];
        }
        return tff;
    }

    hur_nodiscard Array<symmTensor> operator+(const Array<symmTensor> &tf,
                                              const sphericalTensorArray &df) {
        checkArraysSize(tf, df, "+");
        Array<symmTensor> rf(tf.size());

        for (sphericalTensorArray::size_type i = 0; i < rf.size(); ++i) {
            rf[i] = tf[i] + df[i];
        }
        return rf;
    }

    hur_nodiscard Array<symmTensor> operator+(Array<symmTensor> &&tf,
                                              const sphericalTensorArray &df) {
        checkArraysSize(tf, df, "+");
        Array<symmTensor> tff(std::move(tf));

        for (sphericalTensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = tff[i] + df[i];
        }
        return tff;
    }

    hur_nodiscard Array<symmTensor> operator-(const Array<symmTensor> &tf,
                                              const sphericalTensorArray &df) {
        checkArraysSize(tf, df, "-");
        Array<symmTensor> rf(tf.size());

        for (sphericalTensorArray::size_type i = 0; i < rf.size(); ++i) {
            rf[i] = tf[i] - df[i];
        }
        return rf;
    }

    hur_nodiscard Array<symmTensor> operator-(Array<symmTensor> &&tf,
                                              const sphericalTensorArray &df) {
        checkArraysSize(tf, df, "-");
        Array<symmTensor> tff(std::move(tf));

        for (sphericalTensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i] = tff[i] - df[i];
        }
        return tff;
    }

    hur_nodiscard sphericalTensorArray operator*(const realArray &sf,
                                                 const sphericalTensorArray &df) {
        checkArraysSize(sf, df, "*");
        sphericalTensorArray tff(sf.size());

        for (sphericalTensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i].ii() = sf[i] * df[i].ii();
        }
        return tff;
    }

    hur_nodiscard sphericalTensorArray operator*(const realArray &sf, sphericalTensorArray &&df) {
        checkArraysSize(sf, df, "*");
        sphericalTensorArray tff(std::move(df));

        for (sphericalTensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i].ii() *= sf[i];
        }
        return tff;
    }

    hur_nodiscard sphericalTensorArray operator*(const realArray &sf, const Identity<real> &ii) {
        sphericalTensorArray tff(sf.size());

        for (sphericalTensorArray::size_type i = 0; i < tff.size(); ++i) {
            tff[i].ii() = sf[i];
        }
        return tff;
    }
} // namespace OpenHurricane