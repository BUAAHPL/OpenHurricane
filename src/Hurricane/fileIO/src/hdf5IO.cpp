/*!
 * \file hdf5IO.cpp
 * \brief Main subroutines of the <i>hdf5IO.hpp</i> file.
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
#include "hdf5IO.hpp"
template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const typename feature<real>::elementType factor,
                              const integer size, const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, factor, size, dataName);
}
template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const typename feature<real>::elementType factor,
                              const integer size, const string &groupName, const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, factor, size, groupName, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const integer size, const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, size, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const integer size, const string &groupName,
                              const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, size, groupName, dataName);
}

template <> void OpenHurricane::hdf5IO::write(const List<real> &l, const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<real> &l, const string &groupName,
                              const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, groupName, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l,
                              const typename feature<integer>::elementType factor,
                              const integer size, const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, factor, size, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l,
                              const typename feature<integer>::elementType factor,
                              const integer size, const string &groupName, const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, factor, size, groupName, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l, const integer size, const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, size, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l, const integer size, const string &groupName,
                              const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, size, groupName, dataName);
}

template <> void OpenHurricane::hdf5IO::write(const List<integer> &l, const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<integer> &l, const string &groupName,
                              const string &dataName) {
    if (l.size() == 0) {
        return;
    }
    writeSingleCmpt(l, groupName, dataName);
}

void OpenHurricane::hdf5IO::writeIntegerAttributeToFile(const integer l, const string &attrName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    auto attr =
        filePtr_->createAttribute(attrName.c_str(), H5::PredType::NATIVE_INT32, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_INT32, &l);
}

void OpenHurricane::hdf5IO::writeRealAttributeToFile(const real l, const string &attrName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);

    double d[1];
    d[0] = (double)l;

    auto attr =
        filePtr_->createAttribute(attrName.c_str(), H5::PredType::NATIVE_DOUBLE, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_DOUBLE, d);
}

void OpenHurricane::hdf5IO::writeVectorAttributeToFile(const vector &l, const string &attrName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {3};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);

    double d[3];
    d[0] = (double)l[0];
    d[1] = (double)l[1];
    d[2] = (double)l[2];

    auto attr =
        filePtr_->createAttribute(attrName.c_str(), H5::PredType::NATIVE_DOUBLE, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_DOUBLE, d);
}

void OpenHurricane::hdf5IO::writeStringAttributeToFile(const std::string &l, const string &attrName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    atrrdimsf[0] = l.size();
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);

    auto attr =
        filePtr_->createAttribute(attrName.c_str(), H5::PredType::NATIVE_CHAR, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_CHAR, l.c_str());
}

void OpenHurricane::hdf5IO::writeIntegerAttributeToDataset(const integer l, const string &attrName,
                                                       const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    auto attr =
        dataset.createAttribute(attrName.c_str(), H5::PredType::NATIVE_INT32, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_INT32, &l);
}

void OpenHurricane::hdf5IO::writeIntegerAttributeToDataset(const integer l, const string &attrName,
                                                       const string &groupName,
                                                       const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    auto groupSet = filePtr_->openGroup(groupName.c_str());
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    auto attr =
        dataset.createAttribute(attrName.c_str(), H5::PredType::NATIVE_INT32, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_INT32, &l);
}

void OpenHurricane::hdf5IO::writeRealAttributeToDataset(const real l, const string &attrName,
                                                    const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);

    double d[1];
    d[0] = (double)l;
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());
    auto attr =
        dataset.createAttribute(attrName.c_str(), H5::PredType::NATIVE_DOUBLE, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_DOUBLE, d);
}

void OpenHurricane::hdf5IO::writeRealAttributeToDataset(const real l, const string &attrName,
                                                    const string &groupName,
                                                    const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);

    double d[1];
    d[0] = (double)l;
    auto groupSet = filePtr_->openGroup(groupName.c_str());
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());
    auto attr =
        dataset.createAttribute(attrName.c_str(), H5::PredType::NATIVE_DOUBLE, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_DOUBLE, d);
}

void OpenHurricane::hdf5IO::writeVectorAttributeToDataset(const vector &l, const string &attrName,
                                                      const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {3};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);

    double d[3];
    d[0] = (double)l[0];
    d[1] = (double)l[1];
    d[2] = (double)l[2];
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());
    auto attr =
        dataset.createAttribute(attrName.c_str(), H5::PredType::NATIVE_DOUBLE, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_DOUBLE, d);
}

void OpenHurricane::hdf5IO::writeVectorAttributeToDataset(const vector &l, const string &attrName,
                                                      const string &groupName,
                                                      const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {3};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);

    double d[3];
    d[0] = (double)l[0];
    d[1] = (double)l[1];
    d[2] = (double)l[2];
    auto groupSet = filePtr_->openGroup(groupName.c_str());
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());
    auto attr =
        dataset.createAttribute(attrName.c_str(), H5::PredType::NATIVE_DOUBLE, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_DOUBLE, d);
}

void OpenHurricane::hdf5IO::writeStringAttributeToDataset(const std::string &l, const string &attrName,
                                                      const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    atrrdimsf[0] = l.size();
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());

    auto attr = dataset.createAttribute(attrName.c_str(), H5::PredType::NATIVE_CHAR, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_CHAR, l.c_str());
}

void OpenHurricane::hdf5IO::writeStringAttributeToDataset(const std::string &l, const string &attrName,
                                                      const string &groupName,
                                                      const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    atrrdimsf[0] = l.size();
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    auto groupSet = filePtr_->openGroup(groupName.c_str());
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());

    auto attr = dataset.createAttribute(attrName.c_str(), H5::PredType::NATIVE_CHAR, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_CHAR, l.c_str());
}

void OpenHurricane::hdf5IO::writeIntegerAttributeToGroup(const integer l, const string &attrName,
                                                     const string &groupName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    auto groupSet = filePtr_->openGroup(groupName.c_str());
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    auto attr =
        groupSet.createAttribute(attrName.c_str(), H5::PredType::NATIVE_INT32, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_INT32, &l);
}

void OpenHurricane::hdf5IO::writeRealAttributeToGroup(const real l, const string &attrName,
                                                  const string &groupName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);

    double d[1];
    d[0] = (double)l;
    auto groupSet = filePtr_->openGroup(groupName.c_str());
    auto attr =
        groupSet.createAttribute(attrName.c_str(), H5::PredType::NATIVE_DOUBLE, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_DOUBLE, d);
}

void OpenHurricane::hdf5IO::writeVectorAttributeToGroup(const vector &l, const string &attrName,
                                                    const string &groupName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {3};
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);

    double d[3];
    d[0] = (double)l[0];
    d[1] = (double)l[1];
    d[2] = (double)l[2];
    auto groupSet = filePtr_->openGroup(groupName.c_str());
    auto attr =
        groupSet.createAttribute(attrName.c_str(), H5::PredType::NATIVE_DOUBLE, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_DOUBLE, d);
}

void OpenHurricane::hdf5IO::writeStringAttributeToGroup(const std::string &l, const string &attrName,
                                                    const string &groupName) {
    if (closed()) {
        errorAbortStr(("Attempt to write attribute to closed file: " + attrName));
    }
    const int atrrRank = 1;
    hsize_t atrrdimsf[1] = {1};
    atrrdimsf[0] = l.size();
    H5::DataSpace atrrdataSpace(atrrRank, atrrdimsf);
    auto groupSet = filePtr_->openGroup(groupName.c_str());

    auto attr =
        groupSet.createAttribute(attrName.c_str(), H5::PredType::NATIVE_CHAR, atrrdataSpace);
    attr.write(H5::PredType::NATIVE_CHAR, l.c_str());
}

void OpenHurricane::hdf5IO::readIntegerAttributeFromFile(integer &l, const string &attrName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    H5::Attribute attr = filePtr_->openAttribute(attrName.c_str());
    H5::DataSpace dataspace = attr.getSpace();

    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 1) {
        LFatal("The size of attribute for file muts be 1");
    }

    auto dType = attr.getDataType();
    int32_t d[1];
    attr.read(dType, d);
    l = (integer)d[0];
}

void OpenHurricane::hdf5IO::readRealAttributeFromFile(real &l, const string &attrName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    H5::Attribute attr = filePtr_->openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 1) {
        LFatal("The size of attribute for file muts be 1");
    }

    auto dType = attr.getDataType();
    double d[1];
    attr.read(dType, d);
    l = (real)d[0];
}

void OpenHurricane::hdf5IO::readVectorAttributeFromFile(vector &l, const string &attrName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    H5::Attribute attr = filePtr_->openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 3) {
        LFatal("The size of attribute for file muts be 3");
    }

    auto dType = attr.getDataType();
    double d[3];
    attr.read(dType, d);
    l[0] = (real)d[0];
    l[1] = (real)d[1];
    l[2] = (real)d[2];
}

void OpenHurricane::hdf5IO::readStringAttributeFromFile(std::string &l, const string &attrName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    H5::Attribute attr = filePtr_->openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    auto dType = attr.getDataType();
    char *c = new char[dims_out[0] + 1];
    attr.read(dType, c);
    c[dims_out[0]] = '\0';
    l = c;
}

void OpenHurricane::hdf5IO::readIntegerAttributeFromDataset(integer &l, const string &attrName,
                                                        const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());

    H5::Attribute attr = dataset.openAttribute(attrName.c_str());
    H5::DataSpace dataspace = attr.getSpace();

    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 1) {
        LFatal("The size of attribute for file muts be 1");
    }

    auto dType = attr.getDataType();
    int32_t d[1];
    attr.read(dType, d);
    l = (integer)d[0];
}

void OpenHurricane::hdf5IO::readIntegerAttributeFromDataset(integer &l, const string &attrName,
                                                        const string &groupName,
                                                        const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    auto groupSet = filePtr_->openGroup(groupName.c_str());
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());

    H5::Attribute attr = dataset.openAttribute(attrName.c_str());
    H5::DataSpace dataspace = attr.getSpace();

    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 1) {
        LFatal("The size of attribute for file muts be 1");
    }

    auto dType = attr.getDataType();
    int32_t d[1];
    attr.read(dType, d);
    l = (integer)d[0];
}

void OpenHurricane::hdf5IO::readRealAttributeFromDataset(real &l, const string &attrName,
                                                     const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }

    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());

    H5::Attribute attr = dataset.openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 1) {
        LFatal("The size of attribute for file muts be 1");
    }

    auto dType = attr.getDataType();
    double d[1];
    attr.read(dType, d);
    l = (real)d[0];
}

void OpenHurricane::hdf5IO::readRealAttributeFromDataset(real &l, const string &attrName,
                                                     const string &groupName,
                                                     const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }

    auto groupSet = filePtr_->openGroup(groupName.c_str());
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());

    H5::Attribute attr = dataset.openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 1) {
        LFatal("The size of attribute for file muts be 1");
    }

    auto dType = attr.getDataType();
    double d[1];
    attr.read(dType, d);
    l = (real)d[0];
}

void OpenHurricane::hdf5IO::readVectorAttributeFromDataset(vector &l, const string &attrName,
                                                       const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }

    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());

    H5::Attribute attr = dataset.openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 3) {
        LFatal("The size of attribute for file muts be 3");
    }

    auto dType = attr.getDataType();
    double d[3];
    attr.read(dType, d);
    l[0] = (real)d[0];
    l[1] = (real)d[1];
    l[2] = (real)d[2];
}

void OpenHurricane::hdf5IO::readVectorAttributeFromDataset(vector &l, const string &attrName,
                                                       const string &groupName,
                                                       const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }

    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());

    H5::Attribute attr = dataset.openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 3) {
        LFatal("The size of attribute for file muts be 3");
    }

    auto dType = attr.getDataType();
    double d[3];
    attr.read(dType, d);
    l[0] = (real)d[0];
    l[1] = (real)d[1];
    l[2] = (real)d[2];
}

void OpenHurricane::hdf5IO::readStringAttributeFromDataset(std::string &l, const string &attrName,
                                                       const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());

    H5::Attribute attr = dataset.openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    auto dType = attr.getDataType();
    char *c = new char[dims_out[0] + 1];
    attr.read(dType, c);
    c[dims_out[0]] = '\0';
    l = c;
}

void OpenHurricane::hdf5IO::readStringAttributeFromDataset(std::string &l, const string &attrName,
                                                       const string &groupName,
                                                       const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    auto groupSet = filePtr_->openGroup(groupName.c_str());
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());

    H5::Attribute attr = dataset.openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    auto dType = attr.getDataType();
    char *c = new char[dims_out[0] + 1];
    attr.read(dType, c);
    c[dims_out[0]] = '\0';
    l = c;
}

void OpenHurricane::hdf5IO::readIntegerAttributeFromGroup(integer &l, const string &attrName,
                                                      const string &groupName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    auto groupSet = filePtr_->openGroup(groupName.c_str());

    H5::Attribute attr = groupSet.openAttribute(attrName.c_str());
    H5::DataSpace dataspace = attr.getSpace();

    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 1) {
        LFatal("The size of attribute for file muts be 1");
    }

    auto dType = attr.getDataType();
    int32_t d[1];
    attr.read(dType, d);
    l = (integer)d[0];
}

void OpenHurricane::hdf5IO::readRealAttributeFromGroup(real &l, const string &attrName,
                                                   const string &groupName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }

    auto groupSet = filePtr_->openGroup(groupName.c_str());

    H5::Attribute attr = groupSet.openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 1) {
        LFatal("The size of attribute for file muts be 1");
    }

    auto dType = attr.getDataType();
    double d[1];
    attr.read(dType, d);
    l = (real)d[0];
}

void OpenHurricane::hdf5IO::readVectorAttributeFromGroup(vector &l, const string &attrName,
                                                     const string &groupName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }

    auto groupSet = filePtr_->openGroup(groupName.c_str());

    H5::Attribute attr = groupSet.openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 3) {
        LFatal("The size of attribute for file muts be 3");
    }

    auto dType = attr.getDataType();
    double d[3];
    attr.read(dType, d);
    l[0] = (real)d[0];
    l[1] = (real)d[1];
    l[2] = (real)d[2];
}

void OpenHurricane::hdf5IO::readStringAttributeFromGroup(std::string &l, const string &attrName,
                                                     const string &groupName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read attribute from closed file: " + attrName));
    }
    auto groupSet = filePtr_->openGroup(groupName.c_str());

    H5::Attribute attr = groupSet.openAttribute(attrName.c_str());

    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for file");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    auto dType = attr.getDataType();
    char *c = new char[dims_out[0] + 1];
    attr.read(dType, c);
    c[dims_out[0]] = '\0';
    l = c;
}

void OpenHurricane::hdf5IO::readTypeFromAttribute(const H5::DataSet &dataset, int &dataTypr) const {
    H5::Attribute attr = dataset.openAttribute("dataType");
    H5::DataSpace dataspace = attr.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        LFatal("The attribute must be in 1D for data");
    }
    hsize_t dims_out[1];
    const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (dims_out[0] != 1) {
        LFatal("The size of attribute of type for data muts be 1");
    }
    attr.read(H5::PredType::NATIVE_INT32, &dataTypr);
}

void OpenHurricane::hdf5IO::writeString(const std::string &str, const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    const int Rank = 1;
    hsize_t dimsf[1];

    dimsf[0] = str.size();
    H5::IntType dataType(H5::PredType::NATIVE_CHAR);
    auto dType = H5::PredType::NATIVE_CHAR;
    dataType.setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);
    H5::DataSet dataSet = filePtr_->createDataSet(dataName.c_str(), dataType, dataSpace);
    dataSet.write(str.c_str(), dType);
}

void OpenHurricane::hdf5IO::writeString(const std::string &str, const string &groupName,
                                    const string &dataName) {
    if (closed()) {
        errorAbortStr(("Attempt to write data to closed file: " + dataName));
    }
    const int Rank = 1;
    hsize_t dimsf[1];

    dimsf[0] = str.size();
    H5::IntType dataType(H5::PredType::NATIVE_CHAR);
    auto dType = H5::PredType::NATIVE_CHAR;
    dataType.setOrder(H5T_ORDER_LE);
    H5::DataSpace dataSpace(Rank, dimsf);
    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataSet = groupSet.createDataSet(dataName.c_str(), dataType, dataSpace);
    dataSet.write(str.c_str(), dType);
}

template <> void OpenHurricane::hdf5IO::write(const List<std::string> &l, const string &dataName) {
    if (l.size() == 0) {
        return;
    }

    std::string str;
    str = l[0];
    for (integer i = 1; i < l.size(); ++i) {
        str += ",";
        str += l[i];
    }
    writeString(str, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<std::string> &l, const string &groupName,
                              const string &dataName) {
    if (l.size() == 0) {
        return;
    }

    std::string str;
    str = l[0];
    for (integer i = 1; i < l.size(); ++i) {
        str += ",";
        str += l[i];
    }
    writeString(str, groupName, dataName);
}

template <> void OpenHurricane::hdf5IO::write(const List<string> &l, const string &dataName) {
    if (l.size() == 0) {
        return;
    }

    std::string str;
    str = l[0];
    for (integer i = 1; i < l.size(); ++i) {
        str += ",";
        str += l[i];
    }
    writeString(str, dataName);
}

template <>
void OpenHurricane::hdf5IO::write(const List<string> &l, const string &groupName,
                              const string &dataName) {
    if (l.size() == 0) {
        return;
    }

    std::string str;
    str = l[0];
    for (integer i = 1; i < l.size(); ++i) {
        str += ",";
        str += l[i];
    }
    writeString(str, groupName, dataName);
}

void OpenHurricane::hdf5IO::readString(std::string &str, const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read data to closed file: " + dataName));
    }
    H5::DataSet dataset = filePtr_->openDataSet(dataName.c_str());
    H5T_class_t type_class = dataset.getTypeClass();
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

    auto dType = dataset.getDataType();
    char *c = new char[dims_out[0] + 1];
    dataset.read(c, dType, memspace, dataspace);
    c[dims_out[0]] = '\0';
    str = c;
    delete[] c;
}

void OpenHurricane::hdf5IO::readString(std::string &str, const string &groupName,
                                   const string &dataName) const {
    if (closed()) {
        errorAbortStr(("Attempt to read data to closed file: " + dataName));
    }
    auto groupSet = filePtr_->openGroup(groupName);
    H5::DataSet dataset = groupSet.openDataSet(dataName.c_str());
    H5T_class_t type_class = dataset.getTypeClass();
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

    auto dType = dataset.getDataType();
    char *c = new char[dims_out[0] + 1];
    dataset.read(c, dType, memspace, dataspace);
    c[dims_out[0]] = '\0';
    str = c;
    delete[] c;
}

bool OpenHurricane::hdf5IO::exist(const string &name) const {
    return filePtr_->nameExists(name);
}

bool OpenHurricane::hdf5IO::exist(const string &groupName, const string &name) const {
    const auto groupSte = filePtr_->openGroup(groupName);
    return groupSte.nameExists(name);
}

bool OpenHurricane::hdf5IO::isHDF5File() const {
    return filePtr_->isHdf5(filename_);
}

template <> void OpenHurricane::hdf5IO::read(List<real> &l, const string &dataName) const {
    readSingleCmpt(l, dataName);
}

template <>
void OpenHurricane::hdf5IO::read(List<real> &l, const string &groupName, const string &dataName) const {
    readSingleCmpt(l, groupName, dataName);
}

template <> void OpenHurricane::hdf5IO::read(List<integer> &l, const string &dataName) const {
    readSingleCmpt(l, dataName);
}

template <>
void OpenHurricane::hdf5IO::read(List<integer> &l, const string &groupName,
                             const string &dataName) const {
    readSingleCmpt(l, groupName, dataName);
}

template <> void OpenHurricane::hdf5IO::read(List<std::string> &l, const string &dataName) const {
    std::string str;
    readString(str, dataName);
    split(str, l, ",");
}

template <>
void OpenHurricane::hdf5IO::read(List<std::string> &l, const string &groupName,
                             const string &dataName) const {
    std::string str;
    readString(str, groupName, dataName);
    split(str, l, ",");
}

template <> void OpenHurricane::hdf5IO::read(List<string> &l, const string &dataName) const {
    std::string str;
    readString(str, dataName);
    split(str, l, ",");
}

template <>
void OpenHurricane::hdf5IO::read(List<string> &l, const string &groupName,
                             const string &dataName) const {
    std::string str;
    readString(str, groupName, dataName);
    split(str, l, ",");
}