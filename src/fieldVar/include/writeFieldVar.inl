/*!
 * \file writeFieldVar.inl
 * \brief The In-Line functions of the <i>solver.hpp</i> file.
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

inline OpenHurricane::PtrList<OpenHurricane::cellRealArray> &OpenHurricane::writeFieldVar::yiAve() {
    if (yiAvePtr_ == nullptr) {
        yiAvePtr_ = new PtrList<cellRealArray>(yi().size());
    }
    return *yiAvePtr_;
}

inline bool OpenHurricane::writeFieldVar::yiSet() const {
    return yiAvePtr_ != nullptr;
}

inline const OpenHurricane::runtimeMesh &OpenHurricane::writeFieldVar::mesh() const {
    return flows_.mesh();
}

inline const OpenHurricane::cellRealArray &OpenHurricane::writeFieldVar::rho() const {
    return flows_.rho();
}

inline const OpenHurricane::cellRealArray &OpenHurricane::writeFieldVar::p() const {
    return flows_.p();
}

inline const OpenHurricane::cellRealArray &OpenHurricane::writeFieldVar::T() const {
    return flows_.T();
}

inline const OpenHurricane::cellRealArray &OpenHurricane::writeFieldVar::E() const {
    return flows_.E();
}

inline const OpenHurricane::cellRealArray &OpenHurricane::writeFieldVar::gama() const {
    return flows_.gama();
}

inline const OpenHurricane::cellVectorArray &OpenHurricane::writeFieldVar::v() const {
    return flows_.v();
}

inline const OpenHurricane::PtrList<OpenHurricane::cellRealArray> &
OpenHurricane::writeFieldVar::yi() const {
    return flows_.mixtures().Yi();
}

inline void OpenHurricane::writeFieldVar::setOutputField(const string &varName,
                                                         const realArray &value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellRealArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(), value);
    }
}

inline void OpenHurricane::writeFieldVar::setOutputField(const string &varName,
                                                         realArray &&value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellRealArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(),
                                          std::move(value));
    }
}

inline void OpenHurricane::writeFieldVar::setOutputField(const string &varName,
                                                         const vectorArray &value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second =
            new cellVectorArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(), value);
    }
}

inline void OpenHurricane::writeFieldVar::setOutputField(const string &varName,
                                                         vectorArray &&value) const {
    auto inter = outFieldVarMap_.find(varName);
    if (inter != outFieldVarMap_.end()) {
        HurDelete(inter->second);
        inter->second = new cellVectorArray(object(varName, mesh(), object::WRITE_OUTPUT), mesh(),
                                            std::move(value));
    }
}

inline const std::map<std::string, OpenHurricane::object *> &
OpenHurricane::writeFieldVar::getoutVarMap() const {
    return outFieldVarMap_;
}

inline const OpenHurricane::string &OpenHurricane::writeFieldVar::writeId() const noexcept {
    return writeId_;
}

inline void OpenHurricane::writeFieldVar::write() const {
    if (writeNow()) {
        writeToFile();
    }
}