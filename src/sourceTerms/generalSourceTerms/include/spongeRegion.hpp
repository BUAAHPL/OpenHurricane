/*!
 * \file spongeRegion.hpp
 * \brief Headers of base class of sponge region source terms.
 *        The subroutines and functions are in the <i>spongeRegion.cpp</i> file.
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

#include "regionSourceTerms.hpp"

namespace OpenHurricane {

/*!\brief The class of spongeRegion.*/
class spongeRegion : public regionSourceTerms {
public:
protected:
    /**
     * \brief The name of the connected face zone.
     */
    string faceZoneName_;

    /**
     * \brief The id of the connected face zone.
     */
    integer fzid_;

    /**
     * \brief The layer width.
     */
    real Delta_;

    real sigma0_;

    real m_;

    mutable cellRealArray *dPtr_;

    void makeDist() const;
    hur_nodiscard inline const cellRealArray &d() const;

    /**
     * \brief The sponge layer profile.
     */
    mutable realArray *sigmaPtr_;

    /**
     * \brief The sponge layer profile.
     */
    hur_nodiscard const realArray &sigma() const;

    inline const integerList &spongeLayer() const noexcept;

    template <typename Type>
    hur_nodiscard inline Type getRefState(const string &UName) const;

public:
    declareClassNames;

    // Constructors

    spongeRegion() = delete;

    spongeRegion(const flowModel &flows, const iteration &iter,
                 const controller &cont);

    spongeRegion(const spongeRegion &st) = delete;

    /*!\brief Destructor.*/
    virtual ~spongeRegion() noexcept {
        HurDelete(dPtr_);
        HurDelete(sigmaPtr_);
    }

    /**
     * \brief Add to Continuity equation.
     */
    virtual void addSourceTerms(cellRealArray &rho) const;

    /**
     * \brief Add to momentum equation.
     */
    virtual void addSourceTerms(const cellRealArray &rho,
                                cellVectorArray &U) const;

    /**
     * \brief Add to scalar equation.
     */
    virtual void addSourceTerms(const cellRealArray &rho,
                                cellRealArray &phi) const;

    void operator=(const spongeRegion &st) = delete;
};
} // namespace OpenHurricane

#include "spongeRegion.inl"