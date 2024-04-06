/*!
 * \file fileName.hpp
 * \brief Header of file name
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
#pragma once
#include "objectFactory.hpp"
#include "string.hpp"

namespace OpenHurricane {
    template <class T> class List;
    typedef List<string> stringList;

    class fileName {
    public:
        using value_type = char;
        using string_type = std::string;

    private:
        string_type pathText_;

    public:
        static const fileName null;

        inline fileName() : pathText_() {}

        inline fileName(const fileName &other) : pathText_(other.pathText_) {}
        inline fileName &operator=(const fileName &other) {
            if (this != std::addressof(other)) {
                pathText_ = other.pathText_;
            }
            return *this;
        }

        inline fileName(fileName &&other) noexcept : pathText_(std::move(other.pathText_)) {}
        inline fileName &operator=(fileName &&other) noexcept {
            pathText_ = std::move(other.pathText_);
            return *this;
        }

        inline fileName(OpenHurricane::string &&str) noexcept
            : pathText_(std::move(static_cast<std::string &>(str))) {}
        inline fileName &operator=(OpenHurricane::string &&str) noexcept {
            pathText_ = std::move(static_cast<std::string &&>(str));
            return *this;
        }

        inline fileName(const std::string &str) : pathText_(str) {}
        inline fileName &operator=(const std::string &str) {
            pathText_ = str;
            return *this;
        }

        inline fileName(std::string &&str) noexcept : pathText_(std::move(str)) {}
        inline fileName &operator=(std::string &&str) noexcept {
            pathText_ = std::move(str);
            return *this;
        }
        inline fileName(const char *c) : pathText_(c) {}
        inline fileName &operator=(const char *c) {
            pathText_ = c;
            return *this;
        }

        explicit fileName(const stringList &strl);

        inline ~fileName() noexcept {}

        void clear() noexcept { pathText_.clear(); }

        hur_nodiscard bool isAbsolute() const noexcept;
        hur_nodiscard inline bool isRelative() const noexcept { return !isAbsolute(); }

        fileName &toAbsolute();

        inline void swap(fileName &Rhs) noexcept { pathText_.swap(Rhs.pathText_); }

        hur_nodiscard inline bool empty() const noexcept { return pathText_.empty(); }

        hur_nodiscard fileName relativePath() const;

        hur_nodiscard fileName parentPath() const;
        hur_nodiscard inline fileName path() const { return parentPath(); }

        /**\brief Return file name, optionally without extension.*/
        hur_nodiscard std::string name(const bool noExt = false) const;

        /**\brief Return file name without extension (part before last .)*/
        hur_nodiscard fileName lessExt() const;

        /**\brief Return file name extension (part after last .)*/
        hur_nodiscard std::string extension() const;
        hur_nodiscard inline std::string ext() const { return extension(); }
        hur_nodiscard inline auto size() const noexcept { return pathText_.size(); }

        hur_nodiscard inline const value_type *c_str() const noexcept { return pathText_.data(); }
        hur_nodiscard inline const string_type &pathText() const noexcept { return pathText_; }

        hur_nodiscard inline operator const std::string &() const { return pathText_; }
        hur_nodiscard inline operator std::string &() { return pathText_; }
        hur_nodiscard inline const std::string &operator()() const { return pathText_; }
        hur_nodiscard inline std::string &operator()() { return pathText_; }

        fileName &operator/=(const fileName &Added);

        fileName &operator/=(const std::string &Added);

        fileName &operator/=(const char *c);

        inline fileName &operator+=(const fileName &Added) {
            pathText_ += Added.pathText_;
            return *this;
        }

        inline fileName &operator+=(const std::string &Added) {
            pathText_ += Added;
            return *this;
        }

        inline fileName &operator+=(const char *c) {
            pathText_ += c;
            return *this;
        }

        inline fileName &operator+=(char c) {
            pathText_ += c;
            return *this;
        }

        hur_nodiscard friend inline fileName operator+(const fileName &fn1, const fileName &fn2) {
            fileName tmp;
            tmp = fn1.pathText_ + fn2.pathText_;
            return tmp;
        }

        hur_nodiscard friend inline fileName operator+(const fileName &fn1,
                                                       const std::string &fn2) {
            fileName tmp(fn1);
            tmp += fn2;
            return tmp;
        }

        hur_nodiscard friend inline fileName operator+(const std::string &fn1,
                                                       const fileName &fn2) {
            fileName tmp(fn1);
            tmp += fn2;
            return tmp;
        }

        hur_nodiscard friend inline fileName operator+(const fileName &fn1, const char *fn2) {
            fileName tmp(fn1);
            tmp += fn2;
            return tmp;
        }

        hur_nodiscard friend inline fileName operator+(const char *fn1, const fileName &fn2) {
            fileName tmp(fn1);
            tmp += fn2;
            return tmp;
        }

        hur_nodiscard friend inline fileName operator+(const fileName &fn1, char fn2) {
            fileName tmp(fn1);
            tmp += fn2;
            return tmp;
        }

        hur_nodiscard friend inline fileName operator+(char fn1, const fileName &fn2) {
            fileName tmp(fn1);
            tmp += fn2;
            return tmp;
        }
    };

    template <> hur_nodiscard inline std::string toString(const fileName &fn) {
        return fn.pathText();
    }

    hur_nodiscard fileName operator/(const fileName &fn1, const fileName &fn2);

    hur_nodiscard fileName operator/(const fileName &fn1, const std::string &fn2);
    hur_nodiscard fileName operator/(const std::string &fn1, const fileName &fn2);

    hur_nodiscard fileName operator/(const fileName &fn1, const char *fn2);
    hur_nodiscard fileName operator/(const char *fn1, const fileName &fn2);

    hur_nodiscard std::string readFileToString(const fileName &fN);

    hur_nodiscard fileName getCurrentPath();

    hur_nodiscard fileName u8path(const std::string &Source);

} // namespace OpenHurricane
