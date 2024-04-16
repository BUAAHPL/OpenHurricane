/*!
 * \file fileName.cpp
 * \brief The subroutines and functions of file name
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
#include "fileName.hpp"
#include "Lists.hpp"
#include <filesystem>
#include <fstream>
#include <set>

void OpenHurricane::uniqueStringList(stringList &str) {
    std::set<string> ses;
    for (stringList::iterator iter = str.begin(); iter != str.end();) {
        auto retu = ses.emplace(*iter);
        if (!retu.second) {
            iter = str.erase(iter);
            continue;
        }
        iter++;
    }
}

const OpenHurricane::fileName OpenHurricane::fileName::null;

OpenHurricane::fileName::fileName(const stringList &strl) : pathText_() {
    std::filesystem::path tmp;
    for (integer elemI = 0; elemI < strl.size(); ++elemI) {
        tmp /= (static_cast<const std::string &>(strl[elemI]));
    }
    pathText_ = tmp.string();
}

hur_nodiscard bool OpenHurricane::fileName::isAbsolute() const noexcept {
    std::filesystem::path tmp(pathText_);
    return tmp.is_absolute();
}

OpenHurricane::fileName &OpenHurricane::fileName::toAbsolute() {
    std::filesystem::path tmp(pathText_);
    if (!tmp.is_absolute()) {
        tmp = std::filesystem::absolute(tmp);
        pathText_ = tmp.string();
    }
    return *this;
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::fileName::relativePath() const {
    std::filesystem::path tmp(pathText_);
    return tmp.relative_path().string();
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::fileName::parentPath() const {
    std::filesystem::path tmp(pathText_);
    return tmp.parent_path().string();
}

hur_nodiscard std::string OpenHurricane::fileName::name(const bool noExt) const {
    std::filesystem::path tmp(pathText_);
    std::string str(tmp.filename().string());
    if (!noExt) {
        return str;
    }
    size_t dot = str.rfind('.');

    if (dot == std::string::npos) {
        return str;
    } else {
        return str.substr(0, dot);
    }
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::fileName::lessExt() const {
    std::filesystem::path tmp(pathText_);
    tmp.replace_extension();
    return tmp.string();
}

hur_nodiscard std::string OpenHurricane::fileName::extension() const {
    std::filesystem::path tmp(pathText_);
    return tmp.extension().string();
}

OpenHurricane::fileName &OpenHurricane::fileName::operator/=(const fileName &Added) {
    std::filesystem::path tmp1(pathText_);
    tmp1 /= Added.pathText_;
    pathText_ = tmp1.string();
    return *this;
}

OpenHurricane::fileName &OpenHurricane::fileName::operator/=(const std::string &Added) {
    std::filesystem::path tmp1(pathText_);
    tmp1 /= Added;
    pathText_ = tmp1.string();
    return *this;
}

OpenHurricane::fileName &OpenHurricane::fileName::operator/=(const char *c) {
    std::filesystem::path tmp1(pathText_);
    tmp1 /= c;
    pathText_ = tmp1.string();
    return *this;
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::operator/(const fileName &fn1, const fileName &fn2) {
    std::filesystem::path tmp1(fn1.pathText());
    tmp1 /= fn2.pathText();
    return tmp1.string();
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::operator/(const fileName &fn1,
                                                       const std::string &fn2) {
    std::filesystem::path tmp1(fn1.pathText());
    tmp1 /= fn2;
    return tmp1.string();
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::operator/(const std::string &fn1,
                                                       const fileName &fn2) {
    std::filesystem::path tmp1(fn1);
    tmp1 /= fn2.pathText();
    return tmp1.string();
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::operator/(const fileName &fn1, const char *fn2) {
    std::filesystem::path tmp1(fn1.pathText());
    tmp1 /= fn2;
    return tmp1.string();
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::operator/(const char *fn1, const fileName &fn2) {
    std::filesystem::path tmp1(fn1);
    tmp1 /= fn2.pathText();
    return tmp1.string();
}

hur_nodiscard std::string OpenHurricane::readFileToString(const fileName &fN) {
    std::ifstream myin(fN, std::ios::binary);
    if (!myin.good()) {
        LFatal("File: \"%s\" is not exist", fN.c_str());
    }
    char *buffer;
    std::streampos lenthc = myin.seekg(0, std::ios::end).tellg();
    buffer = new char[static_cast<std::streamsize>(lenthc) + 1];
    myin.seekg(0, std::ios::beg).read(buffer, static_cast<std::streamsize>(lenthc));
    myin.close();
    buffer[static_cast<std::streamsize>(lenthc)] = '\0';
    auto fileStr = std::string(buffer);
    delete[] buffer;
    return fileStr;
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::getCurrentPath() {
    auto path = std::filesystem::current_path();
    return path.string();
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::u8path(const std::string &Source) {
#ifdef _WIN32
    auto path = std::filesystem::u8path(Source);
    return path.string();
#else
    return Source;
#endif
}
