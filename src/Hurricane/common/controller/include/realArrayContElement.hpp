/*!
 * \file realArrayContElement.hpp
 * \brief Headers of realArray control elements.
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
#include "controlElement.hpp"
#include "realArray.hpp"

namespace OpenHurricane {

    class realArrayContElement : public controlElement {
    private:
        realArray v_;

    public:
        inline realArrayContElement() : controlElement(), v_() {}

        inline realArrayContElement(const std::string &key, const realArray &value)
            : controlElement(key), v_(value) {}

        inline realArrayContElement(const realArrayContElement &other)
            : controlElement(other), v_(other.v_) {}

        inline realArrayContElement &operator=(const realArrayContElement &other) {
            if (this != std::addressof(other)) {
                controlElement::operator=(other);
                v_ = other.v_;
            }
            return *this;
        }

        inline realArrayContElement(realArrayContElement &&other) noexcept
            : controlElement(std::move(other)), v_(std::move(other.v_)) {}
        inline realArrayContElement &operator=(realArrayContElement &&other) noexcept {
            controlElement::operator=(std::move(other));
            v_ = std::move(other.v_);
            return *this;
        }
        hur_nodiscard virtual uniquePtr<controlElement> clone(const controller &perentCont) const;

        inline virtual ~realArrayContElement() noexcept {}

        hur_nodiscard virtual inline const realArray &rArrayContEle() const { return v_; }
        hur_nodiscard virtual inline realArray &rArrayContEle() { return v_; }

        hur_nodiscard virtual inline bool isRealArrayCE() const noexcept { return true; }

        virtual void writeToXML(std::stringstream &sstr, const integer ilayer = 1) const;
    };

} // namespace OpenHurricane
