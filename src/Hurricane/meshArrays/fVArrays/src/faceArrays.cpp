/*!
 * \file faceArrays.cpp
 * \brief Main subroutines for faceArrays.
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
#include "faceArrays.hpp"
#include "cellMesh.hpp"
#include "faceMesh.hpp"

namespace OpenHurricane {
    template <>
    void OpenHurricane::geometryArray<real, faceMesh>::writeOutput(fileOsstream &fos,
                                                               const integer fzid) const {
        const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
        const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
        const auto &faces = mesh().faces();
        //const auto& fw = mesh().fW();
        const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
        const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
        const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

        realArray faceValue(fsize, Zero);

        integer count = 0;
        for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
            faceValue[count] = this->operator[](fi);
        }
        if (!HurMPI::parRun()) {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                std::stringstream sstr;
                integer i = 0;
                while (i < fsize) {
                    sstr << std::setprecision(feature<real>::precision) << faceValue[i++] << " ";
                    if (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                    }
                    if (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            } else {
                if (fsize) {
                    fos.write(reinterpret_cast<const char *>(&faceValue[0]), fsize * sizeof(real));
                }
            }
            return;
        }
        realArray writeField;
        if (HurMPI::master()) {
            writeField.resize(nTotalFaces, Zero);
        }

        HurMPI::Request request;
        HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                         recvcnt.data(), displs.data(), feature<real>::MPIType, HurMPI::masterNo(),
                         HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::master()) {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                std::stringstream sstr;
                integer i = 0;
                while (i < writeField.size()) {
                    sstr << std::setprecision(feature<real>::precision) << writeField[i++] << " ";
                    if (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                    }
                    if (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                    }
                    sstr << "\n";
                }
                fos.os() << sstr.str().c_str();
            } else {
                if (writeField.size()) {
                    fos.write(reinterpret_cast<const char *>(&writeField[0]),
                              writeField.size() * sizeof(real));
                }
            }
        }
    }

    template <>
    void OpenHurricane::geometryArray<vector, faceMesh>::writeOutput(fileOsstream &fos,
                                                                 const integer fzid) const {
        const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
        const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
        const auto &faces = mesh().faces();
        //const auto& fw = mesh().fW();
        const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
        const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
        const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

        realArray faceValue(fsize, Zero);
        realArray writeField;
        if (!HurMPI::parRun()) {
            if (HurMPI::master()) {
                writeField.resize(nTotalFaces, Zero);
            }
        }
        for (int j = 0; j < vector::nElements_; ++j) {
            integer count = 0;
            for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
                faceValue[count] = this->operator[](fi)[j];
            }
            if (!HurMPI::parRun()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (fsize) {
                        fos.write(reinterpret_cast<const char *>(&faceValue[0]),
                                  fsize * sizeof(real));
                    }
                }
                continue;
            }

            HurMPI::Request request;
            HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeField.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeField[0]),
                                  writeField.size() * sizeof(real));
                    }
                }
            }
        }
    }

    template <>
    void OpenHurricane::geometryArray<vector2D, faceMesh>::writeOutput(fileOsstream &fos,
                                                                   const integer fzid) const {
        const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
        const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
        const auto &faces = mesh().faces();
        //const auto& fw = mesh().fW();
        const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
        const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
        const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

        realArray faceValue(fsize, Zero);
        realArray writeField;
        if (!HurMPI::parRun()) {
            if (HurMPI::master()) {
                writeField.resize(nTotalFaces, Zero);
            }
        }
        for (int j = 0; j < vector2D::nElements_; ++j) {
            integer count = 0;
            for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
                faceValue[count] = this->operator[](fi)[j];
            }
            if (!HurMPI::parRun()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (fsize) {
                        fos.write(reinterpret_cast<const char *>(&faceValue[0]),
                                  fsize * sizeof(real));
                    }
                }
                continue;
            }

            HurMPI::Request request;
            HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeField.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeField[0]),
                                  writeField.size() * sizeof(real));
                    }
                }
            }
        }
    }

    template <>
    void OpenHurricane::geometryArray<tensor, faceMesh>::writeOutput(fileOsstream &fos,
                                                                 const integer fzid) const {
        const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
        const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
        const auto &faces = mesh().faces();
        //const auto& fw = mesh().fW();
        const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
        const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
        const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

        realArray faceValue(fsize, Zero);
        realArray writeField;
        if (!HurMPI::parRun()) {
            if (HurMPI::master()) {
                writeField.resize(nTotalFaces, Zero);
            }
        }
        for (int j = 0; j < tensor::nElements_; ++j) {
            integer count = 0;
            for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
                faceValue[count] = this->operator[](fi)[j];
            }
            if (!HurMPI::parRun()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (fsize) {
                        fos.write(reinterpret_cast<const char *>(&faceValue[0]),
                                  fsize * sizeof(real));
                    }
                }
                continue;
            }

            HurMPI::Request request;
            HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeField.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeField[0]),
                                  writeField.size() * sizeof(real));
                    }
                }
            }
        }
    }

    template <>
    void OpenHurricane::geometryArray<symmTensor, faceMesh>::writeOutput(fileOsstream &fos,
                                                                     const integer fzid) const {
        const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
        const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
        const auto &faces = mesh().faces();
        //const auto& fw = mesh().fW();
        const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
        const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
        const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

        realArray faceValue(fsize, Zero);
        realArray writeField;
        if (!HurMPI::parRun()) {
            if (HurMPI::master()) {
                writeField.resize(nTotalFaces, Zero);
            }
        }
        for (int j = 0; j < symmTensor::nElements_; ++j) {
            integer count = 0;
            for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
                faceValue[count] = this->operator[](fi)[j];
            }
            if (!HurMPI::parRun()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (fsize) {
                        fos.write(reinterpret_cast<const char *>(&faceValue[0]),
                                  fsize * sizeof(real));
                    }
                }
                continue;
            }

            HurMPI::Request request;
            HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeField.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeField[0]),
                                  writeField.size() * sizeof(real));
                    }
                }
            }
        }
    }

} // namespace OpenHurricane

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::faceMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    real minX = veryLarge;
    real maxX = -veryLarge;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        real vv = this->operator[](fi);
        minX = min(minX, vv);
        maxX = max(maxX, vv);
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
void OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::faceMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    for (integer j = 0; j < vector::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
            real vv = this->operator[](fi)[j];
            minX = min(minX, vv);
            maxX = max(maxX, vv);
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
void OpenHurricane::geometryArray<OpenHurricane::vector2D, OpenHurricane::faceMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    for (integer j = 0; j < vector2D::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
            real vv = this->operator[](fi)[j];
            minX = min(minX, vv);
            maxX = max(maxX, vv);
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
void OpenHurricane::geometryArray<OpenHurricane::tensor, OpenHurricane::faceMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    for (integer j = 0; j < tensor::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
            real vv = this->operator[](fi)[j];
            minX = min(minX, vv);
            maxX = max(maxX, vv);
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
void OpenHurricane::geometryArray<OpenHurricane::symmTensor, OpenHurricane::faceMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    for (integer j = 0; j < symmTensor::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
            real vv = this->operator[](fi)[j];
            minX = min(minX, vv);
            maxX = max(maxX, vv);
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