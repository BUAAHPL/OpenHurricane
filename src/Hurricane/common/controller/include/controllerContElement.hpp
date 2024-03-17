/*!
 * \file controllerContElement.hpp
 * \brief Headers of controller control elements.
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
#include "controller.hpp"

namespace OpenHurricane {

    class controllerContElement : public controlElement, public controller {
    private:
    public:
        inline controllerContElement() : controlElement(), controller() {}

        inline controllerContElement(const std::string &key, const controller &parentCont,
                                     const controller &other)
            : controlElement(key), controller(parentCont, other) {}

        inline controllerContElement(const controllerContElement &other) = delete;
        inline controllerContElement(const controller &parentCont,
                                     const controllerContElement &other)
            : controlElement(other), controller(parentCont, other) {}

        inline controllerContElement &operator=(const controllerContElement &other) {
            if (this != std::addressof(other)) {
                controlElement::operator=(other);
                controller::operator=(other);
            }
            return *this;
        }

        inline controllerContElement(controllerContElement &&other) noexcept
            : controlElement(std::move(other)), controller(std::move(other)) {}
        inline controllerContElement(const controller &parentCont,
                                     controllerContElement &&other) noexcept
            : controlElement(std::move(other)), controller(parentCont, std::move(other)) {}
        inline controllerContElement &operator=(controllerContElement &&other) noexcept {
            controlElement::operator=(std::move(other));
            controller::operator=(std::move(other));
            return *this;
        }

        hur_nodiscard virtual uniquePtr<controlElement> clone(const controller &perentCont) const;

        inline virtual ~controllerContElement() noexcept {}

        hur_nodiscard virtual const controller &contContEle() const;
        hur_nodiscard virtual controller &contContEle();

        hur_nodiscard virtual inline bool isControllerCE() const noexcept { return true; }

        virtual void writeToXML(std::stringstream &sstr, const integer ilayer = 1) const;
    };

} // namespace OpenHurricane
