#include "hdf5IO.hpp"
/*!
 * \file hdf5IO.inl
 * \brief In-Line subroutines of the <i>hdf5IO.hpp</i> file.
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

#define readDataSet(Type1, Type2, dType, size, dataset, list, memspace, dataspace) \
    Type1 *data = new Type1[size];                                                 \
    dataset.read(data, dType, memspace, dataspace);                                \
    for (integer i = 0; i < size; ++i) {                                           \
        list[i] = (typename Type2)data[i];                                         \
    }                                                                              \
    delete[] data;

#define readDataSetMultiCmpt(Type1, Type2, dType, size, nCmpts, dataset, list, memspace, \
                             dataspace)                                                  \
    Type1 *data = new Type1[size * nCmpts];                                              \
    dataset.read(data, dType, memspace, dataspace);                                      \
    for (integer i = 0; i < size; ++i) {                                                 \
        for (integer j = 0; j < nCmpts; ++j) {                                           \
            list[i][j] = (typename Type2)data[i * nCmpts + j];                           \
        }                                                                                \
    }                                                                                    \
    delete[] data;

template <class Type>
inline void OpenHurricane::hdf5IO::write(const Type *data, const integer size, const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ != 1) {
        LFatal("Cannot write multiple-components data");
    }

    const int Rank = 1;
    hsize_t dimsf[1];
    dimsf[0] = size;

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);
    H5::DataSet dataSet = filePtr_->createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    dataSet.write(data, dType);
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    if (dataTypePtr != nullptr) {
        delete dataTypePtr;
        dataTypePtr = nullptr;
    }
}

template <class Type>
inline void OpenHurricane::hdf5IO::write(const Type *data, const integer size, const string &groupName,
                                     const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ != 1) {
        LFatal("Cannot write multiple-components data");
    }

    const int Rank = 1;
    hsize_t dimsf[1];
    dimsf[0] = size;

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);
    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataSet = groupSet.createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    dataSet.write(data, dType);
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    if (dataTypePtr != nullptr) {
        delete dataTypePtr;
        dataTypePtr = nullptr;
    }
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeSingleCmpt(const List<Type> &l, const Type factor,
                                               const integer size, const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ != 1) {
        checkWarning("Only write single component data list");
        return;
    }

    const int Rank = 1;
    hsize_t dimsf[1];
    const integer minSize = min(l.size(), size);
    dimsf[0] = minSize;

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);
    H5::DataSet dataSet = filePtr_->createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    Type *data = new Type[minSize];
    for (integer i = 0; i < minSize; ++i) {
        data[i] = l[i] * factor;
    }
    dataSet.write(data, dType);
    delete[] data;

    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    if (dataTypePtr != nullptr) {
        delete dataTypePtr;
        dataTypePtr = nullptr;
    }
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeSingleCmpt(const List<Type> &l, const Type factor,
                                               const integer size, const string &groupName,
                                               const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ != 1) {
        checkWarning("Only write single component data list");
        return;
    }

    const int Rank = 1;
    hsize_t dimsf[1];
    const integer minSize = min(l.size(), size);
    dimsf[0] = minSize;

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);
    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataSet = groupSet.createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    Type *data = new Type[minSize];
    for (integer i = 0; i < minSize; ++i) {
        data[i] = l[i] * factor;
    }
    dataSet.write(data, dType);
    delete[] data;

    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    if (dataTypePtr != nullptr) {
        delete dataTypePtr;
        dataTypePtr = nullptr;
    }
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeSingleCmpt(const List<Type> &l, const integer size,
                                               const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ != 1) {
        checkWarning("Only write single component data list");
        return;
    }

    const int Rank = 1;
    hsize_t dimsf[1];
    const integer minSize = min(l.size(), size);
    dimsf[0] = minSize;

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);
    H5::DataSet dataSet = filePtr_->createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);

    dataSet.write(l.data(), dType);

    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    if (dataTypePtr != nullptr) {
        delete dataTypePtr;
        dataTypePtr = nullptr;
    }
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeSingleCmpt(const List<Type> &l, const integer size,
                                               const string &groupName, const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ != 1) {
        checkWarning("Only write single component data list");
        return;
    }

    const int Rank = 1;
    hsize_t dimsf[1];
    const integer minSize = min(l.size(), size);
    dimsf[0] = minSize;

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);
    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataSet = groupSet.createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);

    dataSet.write(l.data(), dType);

    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    if (dataTypePtr != nullptr) {
        delete dataTypePtr;
        dataTypePtr = nullptr;
    }
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeSingleCmpt(const List<Type> &l, const string &dataName) {
    writeSingleCmpt(l, l.size(), dataName);
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeSingleCmpt(const List<Type> &l, const string &groupName,
                                               const string &dataName) {
    writeSingleCmpt(l, l.size(), groupName, dataName);
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeMultipleCmpt(const List<Type> &l,
                                                 const typename feature<Type>::elementType factor,
                                                 const integer size, const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ == 1) {
        errorAbortStr(("Attempt to write 1D data: " + dataName));
    }
    if (l.size() == 0) {
        return;
    }
    const int Rank = 2;
    hsize_t dimsf[2];
    const integer minSize = min(l.size(), size);
    dimsf[0] = (hsize_t)minSize;
    dimsf[1] = (hsize_t)feature<Type>::nElements_;
    typename feature<Type>::elementType *dataPtr =
        new typename feature<Type>::elementType[minSize * feature<Type>::nElements_];

    for (integer i = 0; i < minSize; ++i) {
        for (integer j = 0; j < feature<Type>::nElements_; ++j) {
            dataPtr[i * feature<Type>::nElements_ + j] = l[i][j] * factor;
        }
    }

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);

    H5::DataSet dataSet = filePtr_->createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    dataSet.write(dataPtr, dType);
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    delete dataTypePtr;
    delete[] dataPtr;
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeMultipleCmpt(const List<Type> &l,
                                                 const typename feature<Type>::elementType factor,
                                                 const integer size, const string &groupName,
                                                 const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ == 1) {
        errorAbortStr(("Attempt to write 1D data: " + dataName));
    }
    if (l.size() == 0) {
        return;
    }
    const int Rank = 2;
    hsize_t dimsf[2];
    const integer minSize = min(l.size(), size);
    dimsf[0] = (hsize_t)minSize;
    dimsf[1] = (hsize_t)feature<Type>::nElements_;
    typename feature<Type>::elementType *dataPtr =
        new typename feature<Type>::elementType[minSize * feature<Type>::nElements_];

    for (integer i = 0; i < minSize; ++i) {
        for (integer j = 0; j < feature<Type>::nElements_; ++j) {
            dataPtr[i * feature<Type>::nElements_ + j] = l[i][j] * factor;
        }
    }

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);

    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataSet = groupSet.createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    dataSet.write(dataPtr, dType);
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    delete dataTypePtr;
    delete[] dataPtr;
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeMultipleCmpt(const List<Type> &l, const integer size,
                                                 const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ == 1) {
        errorAbortStr(("Attempt to write 1D data: " + dataName));
    }
    if (l.size() == 0) {
        return;
    }
    const int Rank = 2;
    hsize_t dimsf[2];
    const integer minSize = min(l.size(), size);
    dimsf[0] = (hsize_t)minSize;
    dimsf[1] = (hsize_t)feature<Type>::nElements_;
    typename feature<Type>::elementType *dataPtr =
        new typename feature<Type>::elementType[minSize * feature<Type>::nElements_];

    for (integer i = 0; i < minSize; ++i) {
        for (integer j = 0; j < feature<Type>::nElements_; ++j) {
            dataPtr[i * feature<Type>::nElements_ + j] = l[i][j];
        }
    }

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);

    H5::DataSet dataSet = filePtr_->createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    dataSet.write(dataPtr, dType);
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    delete dataTypePtr;
    delete[] dataPtr;
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeMultipleCmpt(const List<Type> &l, const integer size,
                                                 const string &groupName, const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    if (feature<Type>::nElements_ == 1) {
        errorAbortStr(("Attempt to write 1D data: " + dataName));
    }
    if (l.size() == 0) {
        return;
    }
    const int Rank = 2;
    hsize_t dimsf[2];
    const integer minSize = min(l.size(), size);
    dimsf[0] = (hsize_t)minSize;
    dimsf[1] = (hsize_t)feature<Type>::nElements_;
    typename feature<Type>::elementType *dataPtr =
        new typename feature<Type>::elementType[minSize * feature<Type>::nElements_];

    for (integer i = 0; i < minSize; ++i) {
        for (integer j = 0; j < feature<Type>::nElements_; ++j) {
            dataPtr[i * feature<Type>::nElements_ + j] = l[i][j];
        }
    }

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);

    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataSet = groupSet.createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    dataSet.write(dataPtr, dType);
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    delete dataTypePtr;
    delete[] dataPtr;
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeMultipleCmpt(const List<Type> &l, const string &dataName) {
    writeMultipleCmpt(l, l.size(), dataName);
}

template <class Type>
inline void OpenHurricane::hdf5IO::writeMultipleCmpt(const List<Type> &l, const string &groupName,
                                                 const string &dataName) {
    writeMultipleCmpt(l, l.size(), groupName, dataName);
}

template <class Type>
inline void OpenHurricane::hdf5IO::write(const List<Type> &l,
                                     const typename feature<Type>::elementType factor,
                                     const integer size, const string &dataName) {
    if (feature<Type>::nElements_ == 1) {
        LFatal("Cannot write single-components data");
    }
    if (l.size() == 0) {
        return;
    }
    writeMultipleCmpt(l, factor, size, dataName);
}

template <class Type>
inline void
OpenHurricane::hdf5IO::write(const List<Type> &l, const typename feature<Type>::elementType factor,
                         const integer size, const string &groupName, const string &dataName) {
    if (feature<Type>::nElements_ == 1) {
        LFatal("Cannot write single-components data");
    }
    if (l.size() == 0) {
        return;
    }
    writeMultipleCmpt(l, factor, size, groupName, dataName);
}

template <class Type>
inline void OpenHurricane::hdf5IO::write(const List<Type> &l, const integer size,
                                     const string &dataName) {
    if (feature<Type>::nElements_ == 1) {
        LFatal("Cannot write single-components data");
    }
    if (l.size() == 0) {
        return;
    }
    writeMultipleCmpt(l, size, dataName);
}

template <class Type>
inline void OpenHurricane::hdf5IO::write(const List<Type> &l, const integer size,
                                     const string &groupName, const string &dataName) {
    if (feature<Type>::nElements_ == 1) {
        LFatal("Cannot write single-components data");
    }
    if (l.size() == 0) {
        return;
    }
    writeMultipleCmpt(l, size, groupName, dataName);
}

template <class Type>
inline void OpenHurricane::hdf5IO::write(const List<Type> &l, const string &dataName) {
    if (feature<Type>::nElements_ == 1) {
        LFatal("Cannot write single-components data");
    }
    if (l.size() == 0) {
        return;
    }
    writeMultipleCmpt(l, l.size(), dataName);
}

template <class Type>
inline void OpenHurricane::hdf5IO::write(const List<Type> &l, const string &groupName,
                                     const string &dataName) {
    if (feature<Type>::nElements_ == 1) {
        LFatal("Cannot write single-components data");
    }
    if (l.size() == 0) {
        return;
    }
    writeMultipleCmpt(l, l.size(), groupName, dataName);
}

template <template <typename T> class Form, typename Type>
inline void OpenHurricane::hdf5IO::writeArrayArray(const Form<Type> &l, const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    const int Rank = 2;
    hsize_t dimsf[2];
    integer maxSize = 0;
    for (integer i = 0; i < l.size(); ++i) {
        maxSize = max(l[i].size(), maxSize);
    }
    if (maxSize == 0) {
        return;
    }
    dimsf[0] = (hsize_t)l.size();
    dimsf[1] = (hsize_t)maxSize;
    typename feature<Type>::elementType *dataPtr =
        new typename feature<Type>::elementType[maxSize * l.size()];

    for (integer i = 0; i < l.size(); ++i) {
        for (integer j = 0; j < maxSize; ++j) {
            if (j < l[i].size()) {
                dataPtr[i * maxSize + j] = l[i][j];
            } else {
                dataPtr[i * maxSize + j] = Zero;
            }
        }
    }

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);

    H5::DataSet dataSet = filePtr_->createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    dataSet.write(dataPtr, dType);
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    delete dataTypePtr;
    delete[] dataPtr;
}

template <template <typename T> class Form, typename Type>
inline void OpenHurricane::hdf5IO::writeArrayArray(const Form<Type> &l, const string &groupName,
                                               const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    const int Rank = 2;
    hsize_t dimsf[2];
    integer maxSize = 0;
    for (integer i = 0; i < l.size(); ++i) {
        maxSize = max(l[i].size(), maxSize);
    }
    if (maxSize == 0) {
        return;
    }
    dimsf[0] = (hsize_t)l.size();
    dimsf[1] = (hsize_t)maxSize;
    typename feature<Type>::elementType *dataPtr =
        new typename feature<Type>::elementType[maxSize * l.size()];

    for (integer i = 0; i < l.size(); ++i) {
        for (integer j = 0; j < maxSize; ++j) {
            if (j < l[i].size()) {
                dataPtr[i * maxSize + j] = l[i][j];
            } else {
                dataPtr[i * maxSize + j] = Zero;
            }
        }
    }

    const int dataType = feature<typename feature<Type>::elementType>::dataFormat;
    H5::IntType *dataTypePtr = nullptr;
    auto dType = H5::PredType::NATIVE_FLOAT;
    if (dataType == 1) // float
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_FLOAT);
        dType = H5::PredType::NATIVE_FLOAT;
    } else if (dataType == 2) // double
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_DOUBLE);
        dType = H5::PredType::NATIVE_DOUBLE;
    } else if (dataType == 3) // int64
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT64);
        dType = H5::PredType::NATIVE_INT64;
    } else if (dataType == 4) // int32
    {
        dataTypePtr = new H5::IntType(H5::PredType::NATIVE_INT32);
        dType = H5::PredType::NATIVE_INT32;
    } else {
        LFatal("Unknown data type");
    }
    dataTypePtr->setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);

    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataSet = groupSet.createDataSet(dataName.c_str(), *dataTypePtr, dataSpace);
    dataSet.write(dataPtr, dType);
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::Attribute atrr =
        dataSet.createAttribute(std::string("dataType"), H5::PredType::NATIVE_INT32, atrrdataSpace);
    atrr.write(H5::PredType::NATIVE_INT32, &dataType);
    delete dataTypePtr;
    delete[] dataPtr;
}

template <class Type>
inline void OpenHurricane::hdf5IO::readSingleCmpt(List<Type> &l, const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read data from closed file: " + dataName));
    }
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());
    //H5T_class_t type_class = dataset.getTypeClass();
    auto dType = dataset.getDataType();
    H5::DataSpace dataspace = dataset.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();

    if (rank != 1) {
        LFatal("The data is not in 1D");
    }

    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    hsize_t offset[1];
    offset[0] = 0;
    dataspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    H5::DataSpace memspace(rank, dims_out);
    memspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    if (l.size() < (integer)dims_out[0]) {
        l.resize((integer)dims_out[0]);
    }
    int dataType = 0;
    readTypeFromAttribute(dataset, dataType);
    if (feature<typename feature<Type>::elementType>::dataFormat == dataType) {
        dataset.read(l.data(), dType, memspace, dataspace);
    } else if (dataType == 1) // float
    {
        readDataSet(float, feature<Type>::elementType, dType, (integer)dims_out[0], dataset, l,
                    memspace, dataspace);
    } else if (dataType == 2) // double
    {
        readDataSet(double, feature<Type>::elementType, dType, (integer)dims_out[0], dataset, l,
                    memspace, dataspace);
    } else if (dataType == 3) // int64
    {
        readDataSet(int64_t, feature<Type>::elementType, dType, (integer)dims_out[0], dataset, l,
                    memspace, dataspace);
    } else if (dataType == 4) // int32
    {
        readDataSet(int32_t, feature<Type>::elementType, dType, (integer)dims_out[0], dataset, l,
                    memspace, dataspace);
    } else {
        LFatal("Unknown data type");
    }
}

template <class Type>
inline void OpenHurricane::hdf5IO::readSingleCmpt(List<Type> &l, const string &groupName,
                                              const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read data from closed file: " + dataName));
    }
    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());
    //H5T_class_t type_class = dataset.getTypeClass();
    auto dType = dataset.getDataType();
    H5::DataSpace dataspace = dataset.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();

    if (rank != 1) {
        LFatal("The data is not in 1D");
    }

    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    hsize_t offset[1];
    offset[0] = 0;
    dataspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    H5::DataSpace memspace(rank, dims_out);
    memspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    if (l.size() < (integer)dims_out[0]) {
        l.resize((integer)dims_out[0]);
    }
    int dataType = 0;
    readTypeFromAttribute(dataset, dataType);
    if (feature<typename feature<Type>::elementType>::dataFormat == dataType) {
        dataset.read(l.data(), dType, memspace, dataspace);
    } else if (dataType == 1) // float
    {
        readDataSet(float, feature<Type>::elementType, dType, (integer)dims_out[0], dataset, l,
                    memspace, dataspace);
    } else if (dataType == 2) // double
    {
        readDataSet(double, feature<Type>::elementType, dType, (integer)dims_out[0], dataset, l,
                    memspace, dataspace);
    } else if (dataType == 3) // int64
    {
        readDataSet(int64_t, feature<Type>::elementType, dType, (integer)dims_out[0], dataset, l,
                    memspace, dataspace);
    } else if (dataType == 4) // int32
    {
        readDataSet(int32_t, feature<Type>::elementType, dType, (integer)dims_out[0], dataset, l,
                    memspace, dataspace);
    } else {
        LFatal("Unknown data type");
    }
}

template <class Type>
inline void OpenHurricane::hdf5IO::readMultipleCmpt(List<Type> &l, const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read data from closed file: " + dataName));
    }
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());
    //H5T_class_t type_class = dataset.getTypeClass();
    auto dType = dataset.getDataType();
    H5::DataSpace dataspace = dataset.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();

    if (rank != 2) {
        LFatal("The data is not in 2D");
    }

    hsize_t dims_out[2];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    hsize_t *offset = new hsize_t[dims_out[1]];
    for (hsize_t i = 0; i < dims_out[1]; ++i) {
        offset[i] = 0;
    }
    dataspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    H5::DataSpace memspace(rank, dims_out);
    memspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    if (l.size() < (integer)dims_out[0]) {
        l.resize((integer)dims_out[0]);
    }
    int dataType = 0;
    readTypeFromAttribute(dataset, dataType);
    if (feature<typename feature<Type>::elementType>::dataFormat == dataType) {
        readDataSetMultiCmpt(typename feature<Type>::elementType, feature<Type>::elementType, dType,
                             (integer)dims_out[0], (integer)dims_out[1], dataset, l, memspace,
                             dataspace);
    } else if (dataType == 1) // float
    {
        readDataSetMultiCmpt(float, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 2) // double
    {
        readDataSetMultiCmpt(double, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 3) // int64
    {
        readDataSetMultiCmpt(int64_t, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 4) // int32
    {
        readDataSetMultiCmpt(int32_t, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else {
        LFatal("Unknown data type");
    }
    delete[] offset;
}

template <class Type>
inline void OpenHurricane::hdf5IO::readMultipleCmpt(List<Type> &l, const string &groupName,
                                                const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read data from closed file: " + dataName));
    }
    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());
    //H5T_class_t type_class = dataset.getTypeClass();
    auto dType = dataset.getDataType();
    H5::DataSpace dataspace = dataset.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();

    if (rank != 2) {
        LFatal("The data is not in 2D");
    }

    hsize_t dims_out[2];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    hsize_t *offset = new hsize_t[dims_out[1]];
    for (hsize_t i = 0; i < dims_out[1]; ++i) {
        offset[i] = 0;
    }
    dataspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    H5::DataSpace memspace(rank, dims_out);
    memspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    if (l.size() < (integer)dims_out[0]) {
        l.resize((integer)dims_out[0]);
    }
    int dataType = 0;
    readTypeFromAttribute(dataset, dataType);
    if (feature<typename feature<Type>::elementType>::dataFormat == dataType) {
        readDataSetMultiCmpt(typename feature<Type>::elementType, feature<Type>::elementType, dType,
                             (integer)dims_out[0], (integer)dims_out[1], dataset, l, memspace,
                             dataspace);
    } else if (dataType == 1) // float
    {
        readDataSetMultiCmpt(float, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 2) // double
    {
        readDataSetMultiCmpt(double, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 3) // int64
    {
        readDataSetMultiCmpt(int64_t, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 4) // int32
    {
        readDataSetMultiCmpt(int32_t, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else {
        LFatal("Unknown data type");
    }
    delete[] offset;
}

template <template <typename T> class Form, typename Type>
inline void OpenHurricane::hdf5IO::readArrayArray(Form<Type> &l, const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read data from closed file: " + dataName));
    }
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());
    //H5T_class_t type_class = dataset.getTypeClass();
    auto dType = dataset.getDataType();
    H5::DataSpace dataspace = dataset.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();

    if (rank != 2) {
        LFatal("The data is not in 2D");
    }

    hsize_t dims_out[2];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    hsize_t *offset = new hsize_t[dims_out[1]];
    for (integer i = 0; i < dims_out[1]; ++i) {
        offset[i] = 0;
    }
    dataspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    H5::DataSpace memspace(rank, dims_out);
    memspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    if (l.size() < dims_out[0]) {
        l.resize((integer)dims_out[0]);
    }
    for (integer i = 0; i < l.size(); ++i) {
        l[i].resize((integer)dims_out[1]);
    }
    int dataType = 0;
    readTypeFromAttribute(dataset, dataType);
    if (feature<typename feature<Type>::elementType>::dataFormat == dataType) {
        readDataSetMultiCmpt(typename feature<Type>::elementType, feature<Type>::elementType, dType,
                             (integer)dims_out[0], (integer)dims_out[1], dataset, l, memspace,
                             dataspace);
    } else if (dataType == 1) // float
    {
        readDataSetMultiCmpt(float, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 2) // double
    {
        readDataSetMultiCmpt(double, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 3) // int64
    {
        readDataSetMultiCmpt(int64_t, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 4) // int32
    {
        readDataSetMultiCmpt(int32_t, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else {
        LFatal("Unknown data type");
    }
    delete[] offset;
}

template <template <typename T> class Form, typename Type>
inline void OpenHurricane::hdf5IO::readArrayArray(Form<Type> &l, const string &groupName,
                                              const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read data from closed file: " + dataName));
    }
    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());
    //H5T_class_t type_class = dataset.getTypeClass();
    auto dType = dataset.getDataType();
    H5::DataSpace dataspace = dataset.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();

    if (rank != 2) {
        LFatal("The data is not in 2D");
    }

    hsize_t dims_out[2];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    hsize_t *offset = new hsize_t[dims_out[1]];
    for (hsize_t i = 0; i < dims_out[1]; ++i) {
        offset[i] = 0;
    }
    dataspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    H5::DataSpace memspace(rank, dims_out);
    memspace.selectHyperslab(H5S_SELECT_SET, dims_out, offset);
    if (l.size() < (integer)dims_out[0]) {
        l.resize((integer)dims_out[0]);
    }
    for (integer i = 0; i < l.size(); ++i) {
        l[i].resize((integer)dims_out[1]);
    }
    int dataType = 0;
    readTypeFromAttribute(dataset, dataType);
    if (feature<typename feature<Type>::elementType>::dataFormat == dataType) {
        readDataSetMultiCmpt(typename feature<Type>::elementType, feature<Type>::elementType, dType,
                             (integer)dims_out[0], (integer)dims_out[1], dataset, l, memspace,
                             dataspace);
    } else if (dataType == 1) // float
    {
        readDataSetMultiCmpt(float, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 2) // double
    {
        readDataSetMultiCmpt(double, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 3) // int64
    {
        readDataSetMultiCmpt(int64_t, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else if (dataType == 4) // int32
    {
        readDataSetMultiCmpt(int32_t, feature<Type>::elementType, dType, (integer)dims_out[0],
                             (integer)dims_out[1], dataset, l, memspace, dataspace);
    } else {
        LFatal("Unknown data type");
    }
    delete[] offset;
}

template <class Type>
inline void OpenHurricane::hdf5IO::read(List<Type> &l, const string &dataName) const {
    if (feature<Type>::nElements_ == 1) {
        LFatal("Cannot read single-components data");
    }
    readMultipleCmpt(l, dataName);
}

template <class Type>
inline void OpenHurricane::hdf5IO::read(List<Type> &l, const string &groupName,
                                    const string &dataName) const {
    if (feature<Type>::nElements_ == 1) {
        LFatal("Cannot read single-components data");
    }
    readMultipleCmpt(l, groupName, dataName);
}

inline void OpenHurricane::hdf5IO::setOpened() noexcept {
    openClosed_ = OPENED;
}

inline void OpenHurricane::hdf5IO::setClosed() noexcept {
    openClosed_ = CLOSED;
}

inline OpenHurricane::hdf5IO::hdf5IO()
    : filename_(), filePtr_(nullptr), flag_(H5F_ACC_DEFAULT), openClosed_(CLOSED) {}

inline OpenHurricane::hdf5IO::hdf5IO(const fileName &fN)
    : filename_(fN), filePtr_(nullptr), flag_(H5F_ACC_DEFAULT), openClosed_(CLOSED) {}

inline OpenHurricane::hdf5IO::hdf5IO(const fileName &fN, const unsigned int flg)
    : filename_(fN), filePtr_(nullptr), flag_(flg), openClosed_(CLOSED) {}

inline OpenHurricane::hdf5IO::~hdf5IO() noexcept {
    if (filePtr_) {
        if (opened()) {
            close();
        }
    }
    filePtr_.clear();
}

inline void OpenHurricane::hdf5IO::open(const fileName &fN) {
    if (opened()) {
        LFatal("Attempt to open an opened file: %s with new file: %s", filename_.c_str(),
               fN.c_str());
    }
    filename_ = fN;
    try {
        filePtr_.reset(new H5::H5File(fN.c_str(), flag_));
    } catch (H5::Exception &e) {
        auto errMsg = e.getDetailMsg();
        LFatal("Can not open file for\"%s\"", errMsg.c_str());
    }
    setOpened();
}

inline void OpenHurricane::hdf5IO::open() {
    if (opened()) {
        return;
    }
    try {
        filePtr_.reset(new H5::H5File(filename_.c_str(), flag_));
    } catch (H5::Exception &e) {
        auto errMsg = e.getDetailMsg();
        LFatal("Can not open file for\"%s\"", errMsg.c_str());
    }
    setOpened();
}

inline void OpenHurricane::hdf5IO::open(const unsigned int flg) {
    if (opened()) {
        close();
    }
    flag_ = flg;

    try {
        filePtr_.reset(new H5::H5File(filename_.c_str(), flag_));
    } catch (H5::Exception &e) {
        auto errMsg = e.getDetailMsg();
        LFatal("Can not open file for\"%s\"", errMsg.c_str());
    }

    setOpened();
}

inline void OpenHurricane::hdf5IO::close() {
    if (opened()) {
        try {
            filePtr_->close();
        } catch (H5::Exception &e) {
            auto errMsg = e.getDetailMsg();
            LError("Can not open file for\"%s\"", errMsg.c_str());
        }
    }
    setClosed();
}

hur_nodiscard inline bool OpenHurricane::hdf5IO::opened() const noexcept {
    return openClosed_ == OPENED;
}

hur_nodiscard inline bool OpenHurricane::hdf5IO::closed() const noexcept {
    return openClosed_ == CLOSED;
}

inline void OpenHurricane::hdf5IO::createGroup(const string &groupName) {
    auto groupSet = filePtr_->createGroup(groupName);
}

#undef readDataSet
#undef readDataSetMultiCmpt