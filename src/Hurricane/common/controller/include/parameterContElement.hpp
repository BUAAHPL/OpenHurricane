/*!
 * \file parameterContElement.hpp
 * \brief Headers of parameter control elements.
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
#include "controlElement.hpp"

namespace OpenHurricane {

    class parameterContElement : public controlElement {
    private:
        std::string parameterValue_;

        mutable IStringStream iss_;

        bool hasSetIS_;

    public:
        inline parameterContElement()
            : controlElement(), parameterValue_(), iss_(), hasSetIS_(false) {}

        inline parameterContElement(const std::string &key, const std::string &value)
            : controlElement(key), parameterValue_(value), iss_(), hasSetIS_(true) {
            iss_.str(parameterValue_);
        }

        inline parameterContElement(const parameterContElement &other)
            : controlElement(other), parameterValue_(other.parameterValue_), iss_(),
              hasSetIS_(true) {
            iss_.str(parameterValue_);
        }

        inline parameterContElement &operator=(const parameterContElement &other) {
            if (this != std::addressof(other)) {
                controlElement::operator=(other);
                parameterValue_ = other.parameterValue_;
                hasSetIS_ = true;
                iss_.str(parameterValue_);
            }
            return *this;
        }

        inline parameterContElement(parameterContElement &&other) noexcept
            : controlElement(std::move(other)), parameterValue_(std::move(other.parameterValue_)),
              iss_(), hasSetIS_(true) {
            iss_.str(parameterValue_);
        }
        inline parameterContElement &operator=(parameterContElement &&other) noexcept {
            controlElement::operator=(std::move(other));
            parameterValue_ = std::move(other.parameterValue_);
            hasSetIS_ = true;
            iss_.str(parameterValue_);
            return *this;
        }

        template <class Type>
        parameterContElement(const string &key, const Type &t)
            : controlElement(key), parameterValue_(), iss_(), hasSetIS_(true) {
            OStringStream oss;
            const std::streamsize defaultPrecision = oss.precision();
            oss.precision(feature<Type>::precision);
            oss << t;
            parameterValue_ = oss.str();
            iss_.str(parameterValue_);
            oss.precision(defaultPrecision);
        }

        hur_nodiscard virtual uniquePtr<controlElement> clone(const controller &perentCont) const;

        inline virtual ~parameterContElement() noexcept {}

        hur_nodiscard virtual inline const std::string &parameterContEle() const {
            return parameterValue_;
        }
        hur_nodiscard virtual inline std::string &parameterContEle() { return parameterValue_; }

        hur_nodiscard virtual IStringStream &ISStreamContEle() const;

        hur_nodiscard virtual inline bool isParameterCE() const noexcept { return true; }

        virtual void writeToXML(std::stringstream &sstr, const integer ilayer = 1) const;
    };

} // namespace OpenHurricane
