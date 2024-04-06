/*!
 * \file controlElement.hpp
 * \brief Headers of control elements.
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
#include "stdMaps.hpp"
#include "fileName.hpp"
#include "real.hpp"
#include "smartPointer.hpp"

namespace OpenHurricane {
    class controller;

    template <class Type> class Array;
    using realArray = Array<real>;

    class controlElement {
    protected:
        std::string key_;

    public:
        inline controlElement() : key_() {}

        inline controlElement(const std::string &key) : key_(key) {}
        inline controlElement(const controlElement &other) : key_(other.key_) {}
        inline controlElement &operator=(const controlElement &other) {
            if (this != std::addressof(other)) {
                key_ = other.key_;
            }
            return *this;
        }

        inline controlElement(controlElement &&other) noexcept : key_(std::move(other.key_)) {}
        inline controlElement &operator=(controlElement &&other) noexcept {
            key_ = std::move(other.key_);
            return *this;
        }

        hur_nodiscard virtual uniquePtr<controlElement> clone(const controller &perentCont) const;

        hur_nodiscard virtual uniquePtr<controlElement> clone() const;

        inline virtual ~controlElement() noexcept {}

        hur_nodiscard inline std::string &key() noexcept { return key_; }
        hur_nodiscard inline const std::string &key() const noexcept { return key_; }

        hur_nodiscard virtual const string &wordContEle() const;
        hur_nodiscard virtual string &wordContEle();

        hur_nodiscard virtual const std::string &parameterContEle() const;
        hur_nodiscard virtual std::string &parameterContEle();

        hur_nodiscard virtual IStringStream &ISStreamContEle() const;

        hur_nodiscard virtual const std::string &textContEle() const;
        hur_nodiscard virtual std::string &textContEle();

        hur_nodiscard virtual inline bool isControllerCE() const noexcept { return false; }
        hur_nodiscard virtual inline bool isWordCE() const noexcept { return false; }
        hur_nodiscard virtual inline bool isParameterCE() const noexcept { return false; }
        hur_nodiscard virtual inline bool isTextCE() const noexcept { return false; }
        hur_nodiscard virtual inline bool isRealArrayCE() const noexcept { return false; }

        hur_nodiscard virtual const controller &contContEle() const;
        hur_nodiscard virtual controller &contContEle();

        hur_nodiscard virtual const realArray &rArrayContEle() const;
        hur_nodiscard virtual realArray &rArrayContEle();

        virtual void writeToXML(std::stringstream& sstr, const integer ilayer = 1) const {
            LFatal("Attempt to call null function");
        }
    };
} // namespace OpenHurricane