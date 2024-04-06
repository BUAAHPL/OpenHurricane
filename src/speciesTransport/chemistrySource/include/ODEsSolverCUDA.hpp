/*!
 * \file ODEsSolverCUDA.hpp
 * \brief Headers of ODEs solver in CUDA.
 *        The subroutines and functions are in the <i>ODEsSolverCUDA.cpp</i> file.
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
#ifdef CUDA_PARALLEL
#include "CUDAChemistrySourceODEs.h"

namespace OpenHurricane {
/**
 * \brief The base class of ODE solver.
 */
class ODEsSolverCUDA {
protected:
    /*!\brief The maximum steps for sub-iteration.*/
    unsigned int maxSteps_;

    /*!\brief The absolute tolerance.*/
    cu_real ATOL_;

    /*!\brief The relative tolerance.*/
    cu_real RTOL_;

    /*!\brief To get the maximum error of the equations.
     * \param[in] y0 - The initial state solutions.
     * \param[in] y - The solutions.
     * \param[in] e - The errors of each equations.
     * \return Return the maximum error.
     */
    cu_device cu_real maxError(CUDAChemistrySourceODEs &odes,
                                        const cu_real *__restrict__ y0,
                                        const cu_real *__restrict__ y,
                                        cu_real *__restrict__ e,
                                        unsigned int tid,
                                        const unsigned int dimX) const;

public:
    // Constructors

    /** \brief Disallow null constructor.*/
    cu_dual ODEsSolverCUDA() = delete;

    /*!
     * \brief Construct from ODEs and by given absolute and relative torelance
     *		  as well as the maximum steps.
     */
    inline cu_host ODEsSolverCUDA(const cu_real atol,
                                         const cu_real rtol,
                                         const unsigned int maxStep = 10000);

    inline cu_dual ODEsSolverCUDA(const ODEsSolverCUDA &os);

    /*!\brief Destructor.*/
    inline cu_dual virtual ~ODEsSolverCUDA() noexcept {}

    /*!\brief Return the maximum steps for iteration.*/
    inline cu_device cu_integer maxSteps() const noexcept;

    /*!\brief Return the absolute tolerance.*/
    inline cu_device cu_real ATOL() const noexcept;

    /*!\brief Return the relative tolerance.*/
    inline cu_device cu_real RTOL() const noexcept;
};

} // namespace OpenHurricane

#include "ODEsSolverCUDA.inl"

#endif // CUDA_PARALLEL