/*!
 * \file reordering.hpp
 * \brief Headers of the reordering mesh.
 *        The subroutines and functions are in the <i>reordering.cpp</i> file.
 * \author Rao Sihang
 * \version V2.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of Hurricane project (Highly Universal Rocket &
 *Ramjet sImulation Codes for ANalysis and Evaluation) \copyright Copyright (C)
 *2019-2024, Prof. Xu Xu's group at Beihang University.
 *
 * License
 *		This file is part of OpenHurricane
 *
 *		OpenHurricane is free software: you can redistribute it and/or
 *		modify it under the terms of the GNU General Public License as
 *		published by the Free Software Foundation, either version 3 of the License, or
 *		(at your option) any later version.
 *
 *		OpenHurricane is distributed in the hope that it will be useful,
 *		but WITHOUT ANY WARRANTY; without even the implied warranty of
 *		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 *		Public License for more details.
 *
 *		You should have received a copy of the GNU General Public
 *		License along with OpenHurricane.  If not, see
 *		<http://www.gnu.org/licenses/>.
 *
 *
 */
#pragma once
#include "smartPointer.hpp"
#include "controller.hpp"
#include "objectFactory.hpp"
#include "decomposingByMetis.hpp"
#include "Lists.hpp"

namespace OpenHurricane {

class reordering {
private:
    /**\brief Disallow null constructor.*/
    reordering() = delete;

    /**\brief Disallow construct as copy.*/
    reordering(const reordering &) = delete;

    /**\brief Disallow default bitwise assignment.*/
    void operator=(const reordering &) = delete;

protected:
    decomposingMeshByMetis &decomposedMesh_;

    integer nvtxs_;

    integerList xadj_;

    integerList adjncy_;

    integerList perm_;

    integerList iperm_;

public:
    inline const integerList &perm() const noexcept;
    inline integerList &perm();

    inline const integerList &iperm() const noexcept;
    inline integerList &iperm();

protected:
    /*!\brief Compute the information for reorder*/
    void computeArray();

    /*!\brief Update the cell information and order face structure*/
    void update();

    /*!\brief Update cell order for facePairMap structure.*/
    void updateFacePairMap();

    /*!\brief Update cell order for perPairMap structue.*/
    void updatePerPairMap();

    integer getBandwidth() const;

    void printBandwidthReduction(const integer bw1, const integer bw2) const;

public:
    declareClassNames;
declareObjFty(reordering,controller,
                        (decomposingMeshByMetis & dcpM), (dcpM));

    // Static data
    

    /*!\brief Constructors.*/

    /*!\brief Construct from components.*/
    inline reordering(decomposingMeshByMetis &dcpM);

    // Selectors

    /*!\brief Select null constructed.*/
    static uniquePtr<reordering> creator(decomposingMeshByMetis &dcpM,
                                       const controller &cont);

    /*!\brief Destructor.*/
    inline virtual ~reordering() noexcept {};

    /*!\brief Reorder cell by one mean of RCM and
     * multilevelNestedDissectionByMetis.*/
    virtual void reorder() = 0;
};

} // namespace OpenHurricane

#include "reordering.inl"