/*!
 * \file vectorArray.cpp
 * \brief Main subroutines for vector array.
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

#include "vectorArray.hpp"
#include "HurMPI.hpp"

namespace OpenHurricane {
    template <> fileOsstream &Array<vector>::writeToStream(fileOsstream &fos) const {
        return writeToStreamWithFactor(fos);
    }

    template <>
    fileOsstream &Array<vector>::writeToStream(fileOsstream &fos, const size_type outSize) const {
        return writeToStreamWithFactor(fos, outSize);
    }

    template <> fileOsstream &Array<vector>::writeToStreamWithFactor(fileOsstream &fos) const {
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            for (int j = 0; j < vector::nElements_; j++) {
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
                for (int j = 0; j < vector::nElements_; j++) {
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
    fileOsstream &Array<vector>::writeToStreamWithFactor(fileOsstream &fos,
                                                         const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            for (int j = 0; j < vector::nElements_; j++) {
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
                for (int j = 0; j < vector::nElements_; j++) {
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
    void Array<vector>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                      const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < vector::nElements_; j++) {
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
    void Array<vector>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                              const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < vector::nElements_; j++) {
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
    void Array<vector>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < vector::nElements_; j++) {
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
    void Array<vector>::writeAveToPout(fileOsstream &fos, const Array<vector> &rhs,
                                       const Array<real> &cV, const size_type n,
                                       const size_type allN, const vector &rhs0, const bool calRhs0,
                                       const bool modifyRhs0) const {
        vector rhsAve = 0.0;
        vector zRhs = 0.0;

        for (size_type i = 0; i < n; ++i) {
            zRhs[0] = fabs(rhs[i][0]) / cV[i];
            rhsAve[0] += sqr(zRhs[0]);
            zRhs[1] = fabs(rhs[i][1]) / cV[i];
            rhsAve[1] += sqr(zRhs[1]);
            zRhs[2] = fabs(rhs[i][2]) / cV[i];
            rhsAve[2] += sqr(zRhs[2]);
        }

        realArray rhsAveL0(HurMPI::getProcSize(), Zero);
        HurMPI::gather(&rhsAve[0], 1, feature<real>::MPIType, rhsAveL0.data(), 1,
                       feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        realArray rhsAveL1(HurMPI::getProcSize(), Zero);
        HurMPI::gather(&rhsAve[1], 1, feature<real>::MPIType, rhsAveL1.data(), 1,
                       feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

        realArray rhsAveL2(HurMPI::getProcSize(), Zero);
        HurMPI::gather(&rhsAve[2], 1, feature<real>::MPIType, rhsAveL2.data(), 1,
                       feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

        if (HurMPI::master()) {
            //rhsAve[0] = 0.0;
            rhsAve = 0.0;

            for (size_type i = 0; i < rhsAveL0.size(); i++) {
                rhsAve[0] += rhsAveL0[i];
                rhsAve[1] += rhsAveL1[i];
                rhsAve[2] += rhsAveL2[i];
            }
            rhsAve /= real(allN);
            rhsAve.x() = sqrt(rhsAve.x());
            rhsAve.y() = sqrt(rhsAve.y());
            rhsAve.z() = sqrt(rhsAve.z());
            if (calRhs0 || (rhs0.x() == 0 && rhs0.y() == 0 && rhs0.z() == 0)) {
                if (rhsAve[0] < veryTiny) {
                    const_cast<vector &>(rhs0)[0] = 1.0;
                } else {
                    const_cast<vector &>(rhs0)[0] = rhsAve[0];
                }
                if (rhsAve[1] < veryTiny) {
                    const_cast<vector &>(rhs0)[1] = 1.0;
                } else {
                    const_cast<vector &>(rhs0)[1] = rhsAve[1];
                }
                if (rhsAve[2] < veryTiny) {
                    const_cast<vector &>(rhs0)[2] = 1.0;
                } else {
                    const_cast<vector &>(rhs0)[2] = rhsAve[2];
                }
                rhsAve = 1.0;
            } else {
                if (modifyRhs0) {
                    real tem1 = rhsAve[0] / rhs0[0];
                    if (tem1 > 1e3) {
                        const_cast<vector &>(rhs0)[0] = rhsAve[0];
                        rhsAve[0] = 1.0;
                    } else {
                        rhsAve[0] = tem1;
                    }

                    real tem2 = rhsAve[1] / rhs0[1];
                    if (tem2 > 1e3) {
                        const_cast<vector &>(rhs0)[1] = rhsAve[1];
                        rhsAve[1] = 1.0;
                    } else {
                        rhsAve[1] = tem2;
                    }

                    real tem3 = rhsAve[2] / rhs0[2];
                    if (tem3 > 1e3) {
                        const_cast<vector &>(rhs0)[2] = rhsAve[2];
                        rhsAve[2] = 1.0;
                    } else {
                        rhsAve[2] = tem3;
                    }
                } else {
                    rhsAve[0] /= rhs0[0];
                    rhsAve[1] /= rhs0[1];
                    rhsAve[2] /= rhs0[2];
                }
            }
        }
        if (HurMPI::master()) {
            std::cout.setf(std::ios::showpoint);
            std::cout << std::left << std::setfill(' ') << std::setw(12) << std::setprecision(5)
                      << rhsAve[0];
            std::cout << std::left << std::setfill(' ') << std::setw(12) << std::setprecision(5)
                      << rhsAve[1];
            std::cout << std::left << std::setfill(' ') << std::setw(12) << std::setprecision(5)
                      << rhsAve[2];
            std::cout.unsetf(std::ios::showpoint);

            if (fos.opened()) {
                fos.os() << '\t' << std::setprecision(5) << fabs(rhsAve[0]);
                fos.os() << '\t' << std::setprecision(5) << fabs(rhsAve[1]);
                fos.os() << '\t' << std::setprecision(5) << fabs(rhsAve[2]);
            }
        }
        HurMPI::bcast(&rhsAve[0], 3, feature<vector>::MPIType, HurMPI::masterNo(),
                      HurMPI::getComm());
        setCurRhs(rhsAve);
    }

    hur_nodiscard realArray operator*(const Array<vector> &f1, const Array<vector> &f2) {
        checkArraysSize(f1, f2, "*");
        realArray f(f1.size());

        for (Array<vector>::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] * f2[i];
        }

        return f;
    }
    hur_nodiscard vectorArray operator^(const Array<vector> &f1, const Array<vector> &f2) {
        checkArraysSize(f1, f2, "^");
        vectorArray f(f1.size());

        for (Array<vector>::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] ^ f2[i];
        }

        return f;
    }

    hur_nodiscard realArray cos(const Array<vector> &vf1, const Array<vector> &vf2) {
        checkArraysSize(vf1, vf2, "cos");
        realArray f(vf1.size());

        for (Array<vector>::size_type i = 0; i < f.size(); ++i) {
            f[i] = cos(vf1[i], vf2[i]);
        }
        return f;
    }

    hur_nodiscard realArray cos(const Array<vector> &vf1, const vector &v2) {
        realArray f(vf1.size());

        for (Array<vector>::size_type i = 0; i < f.size(); ++i) {
            f[i] = cos(vf1[i], v2);
        }
        return f;
    }

    hur_nodiscard realArray cos(const vector &v1, const Array<vector> &vf2) {
        realArray f(vf2.size());

        for (Array<vector>::size_type i = 0; i < f.size(); ++i) {
            f[i] = cos(v1, vf2[i]);
        }
        return f;
    }

    hur_nodiscard realArray dist(const Array<vector> &vf1, const Array<vector> &vf2) {
        checkArraysSize(vf1, vf2, "dist");
        realArray f(vf1.size());

        for (Array<vector>::size_type i = 0; i < f.size(); ++i) {
            f[i] = dist(vf1[i], vf2[i]);
        }
        return f;
    }

    hur_nodiscard realArray div(const Array<vector> &vf1) {
        realArray f(vf1.size());

        for (Array<vector>::size_type i = 0; i < f.size(); ++i) {
            f[i] = vf1[i].x() + vf1[i].y() + vf1[i].z();
        }
        return f;
    }

    template <> fileOsstream &Array<vector2D>::writeToStream(fileOsstream &fos) const {
        return writeToStreamWithFactor(fos);
    }

    template <>
    fileOsstream &Array<vector2D>::writeToStream(fileOsstream &fos, const size_type outSize) const {
        return writeToStreamWithFactor(fos, outSize);
    }

    template <> fileOsstream &Array<vector2D>::writeToStreamWithFactor(fileOsstream &fos) const {
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            for (int j = 0; j < vector2D::nElements_; j++) {
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
                for (int j = 0; j < vector2D::nElements_; j++) {
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
    fileOsstream &Array<vector2D>::writeToStreamWithFactor(fileOsstream &fos,
                                                           const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        fos.setRealPrecision();
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            for (int j = 0; j < vector2D::nElements_; j++) {
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
                for (int j = 0; j < vector2D::nElements_; j++) {
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
    void Array<vector2D>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                        const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < vector2D::nElements_; j++) {
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
    void Array<vector2D>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                                const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < vector2D::nElements_; j++) {
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
    void Array<vector2D>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        for (int j = 0; j < vector2D::nElements_; j++) {
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
} // namespace OpenHurricane