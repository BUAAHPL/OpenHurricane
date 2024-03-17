/*!
 * \file integerArray.cpp
 * \brief Main subroutines for integer Array.
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

#include "integerArray.hpp"
#include "HurMPI.hpp"

namespace OpenHurricane {
    template <> void Array<integer>::replace(const int i, const Array<integer> &l) {
        this->operator=(l);
    }

    template <> void Array<integer>::replace(const int i, const integer &c) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c;
        }
    }

    template <> void Array<integer>::setComponent(const int d, const integer &c) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c;
        }
    }

    template <> void Array<integer>::setComponent(const integer &c) {
        for (size_type j = 0; j < this->size(); j++) {
            this->operator[](j) = c;
        }
    }

    template <>
    void component(Array<typename Array<integer>::value_type> &comp, const Array<integer> &lf,
                   const int d) {
        comp = lf;
    }

    template <>
    Array<typename Array<integer>::value_type> Array<integer>::component(const int d) const {
        return *this;
    }

    template <> fileOsstream &Array<integer>::writeToStream(fileOsstream &fos) const {
        return writeToStreamWithFactor(fos);
    }

    template <>
    fileOsstream &Array<integer>::writeToStream(fileOsstream &fos, const integer outSize) const {
        return writeToStreamWithFactor(fos, outSize);
    }

    template <> fileOsstream &Array<integer>::writeToStreamWithFactor(fileOsstream &fos) const {
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            size_type i = 0;
            while (i < this->size()) {
                fos.os() << this->operator[](i++) << " ";
                if (i < this->size()) {
                    fos.os() << this->operator[](i++) << " ";
                }
                if (i < this->size()) {
                    fos.os() << this->operator[](i++) << " ";
                }
                if (i < this->size()) {
                    fos.os() << this->operator[](i++) << " ";
                }
                if (i < this->size()) {
                    fos.os() << this->operator[](i++) << " ";
                }
                fos.os() << std::endl;
            }
        } else {
            if (this->size()) {
                Array<integer> newF;
                newF.transfer(*(this->clone()));

                fos.write(reinterpret_cast<const char *>(&newF[0]), this->byteSize());
            }
        }

        return fos;
    }

    /*!\brief Start from v_[0].*/
    template <>
    fileOsstream &Array<integer>::writeToStreamWithFactor(fileOsstream &fos,
                                                          const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            size_type i = 0;
            while (i < minSize) {
                fos.os() << this->operator[](i++) << " ";
                if (i < minSize) {
                    fos.os() << this->operator[](i++) << " ";
                }
                if (i < minSize) {
                    fos.os() << this->operator[](i++) << " ";
                }
                if (i < minSize) {
                    fos.os() << this->operator[](i++) << " ";
                }
                if (i < minSize) {
                    fos.os() << this->operator[](i++) << " ";
                }
                fos.os() << std::endl;
            }
        } else {
            if (minSize) {
                Array<integer> newF;
                newF.transfer(*(this->clone()));

                fos.write(reinterpret_cast<const char *>(&newF[0]), minSize * sizeof(integer));
            }
        }
        return fos;
    }

    template <>
    void Array<integer>::writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                       const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        integer minV = feature<integer>::max;
        integer maxV = feature<integer>::min;
        for (size_type i = 0; i < minSize; ++i) {
            minV = min(this->operator[](i), minV);
            maxV = max(this->operator[](i), maxV);
        }
        double minx = (double)minV;
        double maxx = (double)maxV;
        fos.write(reinterpret_cast<const char *>(&minx), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&maxx), sizeof(double));
    }

    template <>
    void Array<integer>::writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                               const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        integer minV = feature<integer>::max;
        integer maxV = feature<integer>::min;
        for (size_type i = 0; i < minSize; ++i) {
            minV = min(this->operator[](i), minV);
            maxV = max(this->operator[](i), maxV);
        }
        HurMPI::reduce(minV, MPI_MIN);
        HurMPI::reduce(maxV, MPI_MAX);
        double minx = (double)minV;
        double maxx = (double)maxV;
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minx), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxx), sizeof(double));
        }
    }

    template <>
    void Array<integer>::writeMinMaxToStream(fileOsstream &fos, const size_type outSize) const {
        size_type minSize = min(outSize, this->size());
        integer minV = feature<integer>::max;
        integer maxV = feature<integer>::min;
        for (size_type i = 0; i < minSize; ++i) {
            minV = min(this->operator[](i), minV);
            maxV = max(this->operator[](i), maxV);
        }
        double minx = (double)minV;
        fos.write(reinterpret_cast<const char *>(&minx), sizeof(double));
        double maxx = (double)maxV;
        fos.write(reinterpret_cast<const char *>(&maxx), sizeof(double));
    }

} // namespace OpenHurricane
