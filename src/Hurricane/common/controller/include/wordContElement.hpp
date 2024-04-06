/*!
 * \file wordContElement.hpp
 * \brief Headers of string control elements.
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

    class wordContElement : public controlElement {
    private:
        string wordValue_;

    public:
        inline wordContElement() : controlElement(), wordValue_() {}

        inline wordContElement(const std::string &key, const string &value)
            : controlElement(key), wordValue_(value) {}

        inline wordContElement(const std::string &key, string &&value)
            : controlElement(key), wordValue_(std::move(value)) {}

        inline wordContElement(const wordContElement &other)
            : controlElement(other), wordValue_(other.wordValue_) {}

        inline wordContElement &operator=(const wordContElement &other) {
            if (this != std::addressof(other)) {
                controlElement::operator=(other);
                wordValue_ = other.wordValue_;
            }
            return *this;
        }

        inline wordContElement(wordContElement &&other) noexcept
            : controlElement(std::move(other)), wordValue_(std::move(other.wordValue_)) {}
        inline wordContElement &operator=(wordContElement &&other) noexcept {
            controlElement::operator=(std::move(other));
            wordValue_ = std::move(other.wordValue_);
            return *this;
        }
        hur_nodiscard virtual uniquePtr<controlElement> clone(const controller &perentCont) const;

        inline virtual ~wordContElement() noexcept {}

        hur_nodiscard virtual inline const string &wordContEle() const { return wordValue_; }
        hur_nodiscard virtual inline string &wordContEle() { return wordValue_; }

        hur_nodiscard virtual inline bool isWordCE() const noexcept { return true; }

        virtual void writeToXML(std::stringstream &sstr, const integer ilayer = 1) const;
    };

} // namespace OpenHurricane
