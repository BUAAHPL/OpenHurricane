/*!
 * \file faceZoneArray.cpp
 * \brief The subroutines and functions of face zone arrays
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

#include "faceZoneArray.hpp"
#include "HurMPI.hpp"
#include "Pout.hpp"
#include <cmath>
#include <iomanip>

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::integer>::writeOutput(fileOsstream &fos) const {
    if (fos.format() == IOsstream::ASCII_FORMAT) {
        if (!HurMPI::parRun()) {
            std::stringstream sstr;
            integer i = 0;
            while (i < this->size()) {
                sstr << std::setprecision(feature<real>::precision) << (*this)[i++] << " ";
                if (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << (*this)[i++] << " ";
                }
                if (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << (*this)[i++] << " ";
                }
                sstr << "\n";
            }
            fos.os() << sstr.str().c_str();
        } else {
            if (this->size()) {
                fos.write(reinterpret_cast<const char *>(&(*this)[0]), this->size() * sizeof(real));
            }
        }
    } else {
        integerList recvcnt(HurMPI::getProcSize(), Zero);
        recvcnt[HurMPI::getProcRank()] = this->size();
        HurMPI::gatherList(recvcnt, HurMPI::masterNo(), HurMPI::getComm());
        integerList displs;
        integer nTotal = 0;
        realArray writeArray;
        if (HurMPI::master()) {
            for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                nTotal += recvcnt[ip];
            }
            displs.resize(HurMPI::getProcSize(), Zero);
            for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                displs[ip] = displs[ip - 1] + recvcnt[ip - 1];
            }
            writeArray.resize(nTotal, Zero);
        }

        HurMPI::Request request;
        HurMPI::igatherv(&const_cast<faceZoneArray<integer> &>(*this)[0], this->size(),
                         feature<integer>::MPIType, writeArray.data(), recvcnt.data(),
                         displs.data(), feature<integer>::MPIType, HurMPI::masterNo(),
                         HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::master()) {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                std::stringstream sstr;
                integer i = 0;
                while (i < writeArray.size()) {
                    sstr << std::setprecision(feature<real>::precision) << writeArray[i++] << " ";
                    if (i < writeArray.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                             << " ";
                    }
                    if (i < writeArray.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                             << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            } else {
                if (writeArray.size()) {
                    fos.write(reinterpret_cast<const char *>(&writeArray[0]),
                              writeArray.size() * sizeof(real));
                }
            }
        }
    }
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::integer>::writeOutput(fileOsstream &fos,
                                                               const integer fzid) const {
    writeOutput(fos);
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::integer>::writeMinMaxOutput(fileOsstream &fos) const {
    integer minX = feature<integer>::max;
    integer maxX = feature<integer>::min;
    for (integer i = 0; i < this->size(); ++i) {
        minX = min(minX, this->operator[](i));
        maxX = max(maxX, this->operator[](i));
    }

    if (HurMPI::parRun()) {
        HurMPI::reduce(minX, MPI_MIN);
        HurMPI::reduce(maxX, MPI_MAX);
    }
    if (HurMPI::master()) {
        double minx = (double)minX;
        fos.write(reinterpret_cast<const char *>(&minx), sizeof(double));
        double maxx = (double)maxX;
        fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
    }
}
template <>
void OpenHurricane::faceZoneArray<OpenHurricane::integer>::writeMinMaxOutput(fileOsstream &fos,
                                                                     const integer fzid) const {
    writeMinMaxOutput(fos);
}


template <> void OpenHurricane::faceZoneArray<OpenHurricane::real>::writeOutput(fileOsstream &fos) const {
    if (!HurMPI::parRun()) {
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            integer i = 0;
            while (i < this->size()) {
                sstr << std::setprecision(feature<real>::precision) << (*this)[i++] << " ";
                if (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << (*this)[i++] << " ";
                }
                if (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << (*this)[i++] << " ";
                }
                sstr << "\n";
            }
            fos.os() << sstr.str().c_str();
        } else {
            if (this->size()) {
                fos.write(reinterpret_cast<const char *>(&(*this)[0]), this->size() * sizeof(real));
            }
        }
    } else {
        integerList recvcnt(HurMPI::getProcSize(), Zero);
        recvcnt[HurMPI::getProcRank()] = this->size();
        HurMPI::gatherList(recvcnt, HurMPI::masterNo(), HurMPI::getComm());
        integerList displs;
        integer nTotal = 0;
        realArray writeArray;
        if (HurMPI::master()) {
            for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                nTotal += recvcnt[ip];
            }
            displs.resize(HurMPI::getProcSize(), Zero);
            for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                displs[ip] = displs[ip - 1] + recvcnt[ip - 1];
            }
            writeArray.resize(nTotal, Zero);
        }

        HurMPI::Request request;
        HurMPI::igatherv(&const_cast<faceZoneArray<real> &>(*this)[0], this->size(),
                         feature<real>::MPIType, writeArray.data(), recvcnt.data(), displs.data(),
                         feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::master()) {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                std::stringstream sstr;
                integer i = 0;
                while (i < writeArray.size()) {
                    sstr << std::setprecision(feature<real>::precision) << writeArray[i++] << " ";
                    if (i < writeArray.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                             << " ";
                    }
                    if (i < writeArray.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                             << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            } else {
                if (writeArray.size()) {
                    fos.write(reinterpret_cast<const char *>(&writeArray[0]),
                              writeArray.size() * sizeof(real));
                }
            }
        }
    }
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::real>::writeOutput(fileOsstream &fos,
                                                            const integer fzid) const {
    writeOutput(fos);
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::real>::writeMinMaxOutput(fileOsstream &fos) const {
    real minX = veryLarge;
    real maxX = -veryLarge;
    for (integer i = 0; i < this->size(); ++i) {
        minX = min(minX, this->operator[](i));
        maxX = max(maxX, this->operator[](i));
    }

    if (HurMPI::parRun()) {
        HurMPI::reduce(minX, MPI_MIN);
        HurMPI::reduce(maxX, MPI_MAX);
    }
    if (HurMPI::master()) {
        fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
    }
}
template <>
void OpenHurricane::faceZoneArray<OpenHurricane::real>::writeMinMaxOutput(fileOsstream &fos,
                                                                  const integer fzid) const {
    writeMinMaxOutput(fos);
}


template <>
void OpenHurricane::faceZoneArray<OpenHurricane::symmTensor>::writeOutput(fileOsstream &fos) const {
    for (int j = 0; j < symmTensor::nElements_; ++j) {
        if (!HurMPI::parRun()) {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                std::stringstream sstr;
                integer i = 0;
                while (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j] << " ";
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j]
                             << " ";
                    }
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j]
                             << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            } else {
                if (this->size()) {
                    realArray cnf = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&cnf[0]), this->size() * sizeof(real));
                }
            }
        } else {
            realArray cnf = this->component(j);
            integerList recvcnt(HurMPI::getProcSize(), Zero);
            recvcnt[HurMPI::getProcRank()] = this->size();
            HurMPI::gatherList(recvcnt, HurMPI::masterNo(), HurMPI::getComm());
            integerList displs;
            integer nTotal = 0;
            realArray writeArray;
            if (HurMPI::master()) {
                for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                    nTotal += recvcnt[ip];
                }
                displs.resize(HurMPI::getProcSize(), Zero);
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + recvcnt[ip - 1];
                }
                writeArray.resize(nTotal, Zero);
            }

            HurMPI::Request request;
            HurMPI::igatherv(&cnf[0], this->size(), feature<real>::MPIType, writeArray.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeArray.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                             << " ";
                        if (i < writeArray.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                                 << " ";
                        }
                        if (i < writeArray.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeArray.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeArray[0]),
                                  writeArray.size() * sizeof(real));
                    }
                }
            }
        }
    }
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::symmTensor>::writeOutput(fileOsstream &fos,
                                                                  const integer fzid) const {
    writeOutput(fos);
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::symmTensor>::writeMinMaxOutput(fileOsstream &fos) const {
    for (integer j = 0; j < symmTensor::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer i = 0; i < this->size(); ++i) {
            minX = min(minX, this->operator[](i)[j]);
            maxX = max(maxX, this->operator[](i)[j]);
        }

        if (HurMPI::parRun()) {
            HurMPI::reduce(minX, MPI_MIN);
            HurMPI::reduce(maxX, MPI_MAX);
        }
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
        }
    }
}
template <>
void OpenHurricane::faceZoneArray<OpenHurricane::symmTensor>::writeMinMaxOutput(fileOsstream &fos,
                                                                        const integer fzid) const {
    writeMinMaxOutput(fos);
}


template <> void OpenHurricane::faceZoneArray<OpenHurricane::tensor>::writeOutput(fileOsstream &fos) const {
    for (int j = 0; j < tensor::nElements_; ++j) {
        if (!HurMPI::parRun()) {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                std::stringstream sstr;
                integer i = 0;
                while (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j] << " ";
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j]
                             << " ";
                    }
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j]
                             << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            } else {
                if (this->size()) {
                    realArray cnf = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&cnf[0]), this->size() * sizeof(real));
                }
            }
        } else {
            realArray cnf = this->component(j);
            integerList recvcnt(HurMPI::getProcSize(), Zero);
            recvcnt[HurMPI::getProcRank()] = this->size();
            HurMPI::gatherList(recvcnt, HurMPI::masterNo(), HurMPI::getComm());
            integerList displs;
            integer nTotal = 0;
            realArray writeArray;
            if (HurMPI::master()) {
                for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                    nTotal += recvcnt[ip];
                }
                displs.resize(HurMPI::getProcSize(), Zero);
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + recvcnt[ip - 1];
                }
                writeArray.resize(nTotal, Zero);
            }

            HurMPI::Request request;
            HurMPI::igatherv(&cnf[0], this->size(), feature<real>::MPIType, writeArray.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeArray.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                             << " ";
                        if (i < writeArray.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                                 << " ";
                        }
                        if (i < writeArray.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeArray.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeArray[0]),
                                  writeArray.size() * sizeof(real));
                    }
                }
            }
        }
    }
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::tensor>::writeOutput(fileOsstream &fos,
                                                              const integer fzid) const {
    writeOutput(fos);
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::tensor>::writeMinMaxOutput(fileOsstream &fos) const {
    for (integer j = 0; j < tensor::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer i = 0; i < this->size(); ++i) {
            minX = min(minX, this->operator[](i)[j]);
            maxX = max(maxX, this->operator[](i)[j]);
        }

        if (HurMPI::parRun()) {
            HurMPI::reduce(minX, MPI_MIN);
            HurMPI::reduce(maxX, MPI_MAX);
        }
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
        }
    }
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::tensor>::writeMinMaxOutput(fileOsstream &fos,
                                                                    const integer fzid) const {
    writeMinMaxOutput(fos);
}

template <> void OpenHurricane::faceZoneArray<OpenHurricane::vector>::writeOutput(fileOsstream &fos) const {
    for (int j = 0; j < vector::nElements_; ++j) {
        if (!HurMPI::parRun()) {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                std::stringstream sstr;
                integer i = 0;
                while (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j] << " ";
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j]
                             << " ";
                    }
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j]
                             << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            } else {
                if (this->size()) {
                    realArray cnf = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&cnf[0]), this->size() * sizeof(real));
                }
            }
        } else {
            realArray cnf = this->component(j);
            integerList recvcnt(HurMPI::getProcSize(), Zero);
            recvcnt[HurMPI::getProcRank()] = this->size();
            HurMPI::gatherList(recvcnt, HurMPI::masterNo(), HurMPI::getComm());
            integerList displs;
            integer nTotal = 0;
            realArray writeArray;
            if (HurMPI::master()) {
                for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                    nTotal += recvcnt[ip];
                }
                displs.resize(HurMPI::getProcSize(), Zero);
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + recvcnt[ip - 1];
                }
                writeArray.resize(nTotal, Zero);
            }

            HurMPI::Request request;
            HurMPI::igatherv(&cnf[0], this->size(), feature<real>::MPIType, writeArray.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeArray.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                             << " ";
                        if (i < writeArray.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                                 << " ";
                        }
                        if (i < writeArray.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeArray.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeArray[0]),
                                  writeArray.size() * sizeof(real));
                    }
                }
            }
        }
    }
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::vector>::writeOutput(fileOsstream &fos,
                                                              const integer fzid) const {
    writeOutput(fos);
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::vector>::writeMinMaxOutput(fileOsstream &fos) const {
    for (integer j = 0; j < vector::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer i = 0; i < this->size(); ++i) {
            minX = min(minX, this->operator[](i)[j]);
            maxX = max(maxX, this->operator[](i)[j]);
        }

        if (HurMPI::parRun()) {
            HurMPI::reduce(minX, MPI_MIN);
            HurMPI::reduce(maxX, MPI_MAX);
        }
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
        }
    }
}
template <>
void OpenHurricane::faceZoneArray<OpenHurricane::vector>::writeMinMaxOutput(fileOsstream &fos,
                                                                    const integer fzid) const {
    writeMinMaxOutput(fos);
}
template <>
void OpenHurricane::faceZoneArray<OpenHurricane::vector2D>::writeOutput(fileOsstream &fos) const {
    for (int j = 0; j < vector2D::nElements_; ++j) {
        if (!HurMPI::parRun()) {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                std::stringstream sstr;
                integer i = 0;
                while (i < this->size()) {
                    sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j] << " ";
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j]
                             << " ";
                    }
                    if (i < this->size()) {
                        sstr << std::setprecision(feature<real>::precision) << (*this)[i++][j]
                             << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            } else {
                if (this->size()) {
                    realArray cnf = this->component(j);
                    fos.write(reinterpret_cast<const char *>(&cnf[0]), this->size() * sizeof(real));
                }
            }
        } else {
            realArray cnf = this->component(j);
            integerList recvcnt(HurMPI::getProcSize(), Zero);
            recvcnt[HurMPI::getProcRank()] = this->size();
            HurMPI::gatherList(recvcnt, HurMPI::masterNo(), HurMPI::getComm());
            integerList displs;
            integer nTotal = 0;
            realArray writeArray;
            if (HurMPI::master()) {
                for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                    nTotal += recvcnt[ip];
                }
                displs.resize(HurMPI::getProcSize(), Zero);
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + recvcnt[ip - 1];
                }
                writeArray.resize(nTotal, Zero);
            }

            HurMPI::Request request;
            HurMPI::igatherv(&cnf[0], this->size(), feature<real>::MPIType, writeArray.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeArray.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                             << " ";
                        if (i < writeArray.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                                 << " ";
                        }
                        if (i < writeArray.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeArray[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeArray.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeArray[0]),
                                  writeArray.size() * sizeof(real));
                    }
                }
            }
        }
    }
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::vector2D>::writeOutput(fileOsstream &fos,
                                                                const integer fzid) const {
    writeOutput(fos);
}

template <>
void OpenHurricane::faceZoneArray<OpenHurricane::vector2D>::writeMinMaxOutput(fileOsstream &fos) const {
    for (integer j = 0; j < vector2D::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer i = 0; i < this->size(); ++i) {
            minX = min(minX, this->operator[](i)[j]);
            maxX = max(maxX, this->operator[](i)[j]);
        }

        if (HurMPI::parRun()) {
            HurMPI::reduce(minX, MPI_MIN);
            HurMPI::reduce(maxX, MPI_MAX);
        }
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
        }
    }
}
template <>
void OpenHurricane::faceZoneArray<OpenHurricane::vector2D>::writeMinMaxOutput(fileOsstream &fos,
                                                                      const integer fzid) const {
    writeMinMaxOutput(fos);
}