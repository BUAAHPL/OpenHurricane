/*!
 * \file textContElement.hpp
 * \brief Headers of text control elements.
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

namespace OpenHurricane {

    class textContElement : public controlElement {
    private:
        std::string text_;

    public:
        inline textContElement() : controlElement(), text_() {}

        inline textContElement(const std::string &key, const string &value)
            : controlElement(key), text_(value) {}

        inline textContElement(const textContElement &other)
            : controlElement(other), text_(other.text_) {}

        inline textContElement &operator=(const textContElement &other) {
            if (this != std::addressof(other)) {
                controlElement::operator=(other);
                text_ = other.text_;
            }
            return *this;
        }

        inline textContElement(textContElement &&other) noexcept
            : controlElement(std::move(other)), text_(std::move(other.text_)) {}
        inline textContElement &operator=(textContElement &&other) noexcept {
            controlElement::operator=(std::move(other));
            text_ = std::move(other.text_);
            return *this;
        }
        hur_nodiscard virtual uniquePtr<controlElement> clone(const controller &perentCont) const;

        inline virtual ~textContElement() noexcept {}

        hur_nodiscard virtual inline const std::string &textContEle() const { return text_; }
        hur_nodiscard virtual inline std::string &textContEle() { return text_; }

        hur_nodiscard virtual inline bool isTextCE() const noexcept { return true; }

        virtual void writeToXML(std::stringstream &sstr, const integer ilayer = 1) const;
    };

} // namespace OpenHurricane
