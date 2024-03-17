/*!
 * \file speciesElements.hpp
 * \brief All the information about the definition of the species elements.
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

#include "List.hpp"
#include "atomicWeights.hpp"
#include "real.hpp"
#include "string.hpp"

namespace OpenHurricane {
    class speciesElements;
    /**
     * @brief The list of species elements.
     */
    using speciesElementsList = List<speciesElements>;

    /**
     * \brief The class of the element of the species.
     */
    class speciesElements {
    private:
        /** \brief The name of element*/
        string name_;

        /** \brief Molecular weight of element [g/mol or kg/kmol]*/
        mutable real atomicWeight_;

        /** \brief Number of atoms of this element in the specie.*/
        integer nAtoms_;

        inline void getAtomicWeight() const {
            if (!atomicWeight::search(name_.c_str(), atomicWeight_)) {
                LFatal("Unknown element: \"%s\" please check!", name_.c_str());
            }
        }

    public:
        inline speciesElements() : name_(), atomicWeight_(-1.0), nAtoms_(0) {}

        inline speciesElements(const string &name, const integer Atoms)
            : name_(name), atomicWeight_(-1.0), nAtoms_(Atoms) {
            getAtomicWeight();
        }

        inline speciesElements(const speciesElements &ele)
            : name_(ele.name_), atomicWeight_(ele.atomicWeight_), nAtoms_(ele.nAtoms_) {
            getAtomicWeight();
        }
        inline speciesElements &operator=(const speciesElements &ele) {
            if (this != std::addressof(ele)) {
                name_ = ele.name_;
                nAtoms_ = ele.nAtoms_;
                getAtomicWeight();
            }
            return *this;
        }

        inline ~speciesElements() noexcept {}

        /** \brief Return the number of atoms of this element in the specie.*/
        hur_nodiscard inline integer nAtoms() const noexcept { return nAtoms_; }

        /** \brief Return non-const access of the number of atoms of this element in the specie.*/
        hur_nodiscard inline integer &nAtoms() noexcept { return nAtoms_; }

        /** \brief Return the non-const access to the name of the element*/
        hur_nodiscard inline const string &name() const noexcept { return name_; }

        /** \brief Return the atomic weight of the element [g/mol or kg/kmol]*/
        hur_nodiscard inline real atomicWeight() const {
            if (atomicWeight_ < 0.0) {
                getAtomicWeight();
            }
            return atomicWeight_;
        }

        /** \brief Set atomic weight of the element [g/mol or kg/kmol].*/
        inline void setAtomicWeight(const real aw) noexcept { atomicWeight_ = aw; }

        /** \brief Set atomic weight of the element [g/mol or kg/kmol].*/
        inline real setAtomicWeight() { return atomicWeight(); }
    };

    /**
     * \brief Cast element to a string.
     */
    template <> hur_nodiscard inline std::string toString(const speciesElements &e) {
        return e.name();
    }
} // namespace OpenHurricane
