/*!
 * \file thirdBodyEfficiency.hpp
 * \brief Header of third body efficiencies.
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
#include "Lists.hpp"
#include "speciesList.hpp"

namespace OpenHurricane {
    /**
     * \brief The class of third-body efficiencies of three-body reaction.
     */
    class thirdBodyEfficiency : public realArray {
    private:
        // Private data

        const speciesList &species_;

    public:
        thirdBodyEfficiency() = default;

        inline thirdBodyEfficiency(const speciesList &species, const realArray &efficiencies)
            : realArray(efficiencies), species_(species) {}

        inline thirdBodyEfficiency(const thirdBodyEfficiency &tbe)
            : realArray(tbe), species_(tbe.species_) {}

        inline thirdBodyEfficiency(const speciesList &species, const controller &cont)
            : realArray(cont.findType<realList>("thirdBody", realList())), species_(species) {}

        inline ~thirdBodyEfficiency() noexcept {}

        /**\brief Calculate and return M, the concentration of the third-body.*/
        hur_nodiscard inline real M(const realArray &ci) const {
#ifdef HUR_DEBUG

            if (ci.size() < this->size()) {
                std::string errMsg;
                errMsg = " The size of the yi: ";
                errMsg += toString(ci.size());
                errMsg += " is not euqal to the size of the alphaji: ";
                errMsg += toString(this->size());
                errorAbortStr(errMsg);
            }

#endif // HUR_DEBUG

            real m = 0.0;
            for (integer isp = 0; isp < this->size(); isp++) {
                m += this->operator[](isp) * ci[isp];
            }
            return m;
        }

        hur_nodiscard inline real dMdCi(const integer i) const {
            integer ns = species_.size() - 1;
            return this->operator[](i) - species_[i].W() / species_[ns].W() * this->operator[](ns);
        }

        hur_nodiscard inline const speciesList &species() const noexcept { return species_; }
    };

} // namespace OpenHurricane