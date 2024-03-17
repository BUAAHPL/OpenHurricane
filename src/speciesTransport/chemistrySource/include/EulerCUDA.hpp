/*!
 * \file EulerCUDA.hpp
 * \brief Headers of Euler ODEs solver in CUDA.
 *        The subroutines and functions are in the <i>EulerCUDA.cpp</i> file.
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
#ifdef CUDA_PARALLEL
#include "ODEsSolverCUDA.hpp"

namespace OpenHurricane {
/**
 * \brief The class of adaptive step Euler solver in CUDA.
 */
class EulerCUDA : public ODEsSolverCUDA {

private:
public:
    cu_dual EulerCUDA() = delete;

    inline cu_host EulerCUDA(const cu_real atol, const cu_real rtol,
                                    const unsigned int maxStep = 1000);

    inline cu_dual EulerCUDA(const EulerCUDA &es);

    inline cu_dual virtual ~EulerCUDA() noexcept {}

protected:
    /**
     *\brief Solve the ODEs.
     * \param[in] odes - ODEs
     * \param[in] t0 - The initial time value
     * \param dydt0 - The temporary array
     * \param[in] y0 - The initial solutions
     * \param[in] dt - The timestep
     * \param[out] y - The solutions
     * \param y2tmp - The temporary array
     * \param error - The temporary array for error
     * \param[in] tid - The thread index
     * \param[in] dimX - The x dimension of thread block
     * \return The maximum error
     */
    cu_device cu_real
    solve(CUDAChemistrySourceODEs &odes, const cu_real t0,
          const cu_real *__restrict__ y0, const cu_real *__restrict__ dydt0,
          const cu_real dt, cu_real *__restrict__ y,
          cu_real *__restrict__ y2tmp, cu_real *__restrict__ error,
          unsigned int tid, const unsigned int dimX);

    /**
     *\brief Solve the ODEs.
     * \param[in] odes - ODEs
     * \param[in] t0 - The initial time value
     * \param dydt0 - The temporary array
     * \param[in] y0 - The initial solutions
     * \param[in] dt - The timestep
     * \param[out] y - The solutions
     * \param y2tmp - The temporary array
     * \param error - The temporary array for error
     * \param[in] tid - The thread index
     * \param[in] dimX - The x dimension of thread block
     * \return The maximum error
     */
    cu_device cu_real
    solveFactoring(CUDAChemistrySourceODEs &odes, const cu_real t0,
                   const cu_real *__restrict__ y0,
                   const cu_real *__restrict__ dydt0, const cu_real dt,
                   const cu_real odeFactor, cu_real *__restrict__ y,
                   cu_real *__restrict__ y2tmp, cu_real *__restrict__ error,
                   unsigned int tid, const unsigned int dimX);

    cu_device void solve(CUDAChemistrySourceODEs &odes,
                                const cu_real t0,
                                const cu_real *__restrict__ y0,
                                const cu_real *__restrict__ dydt0,
                                const cu_real dt, cu_real *__restrict__ y,
                                unsigned int tid, const unsigned int dimX);

    cu_device void
    solveFactoring(CUDAChemistrySourceODEs &odes, const cu_real t0,
                   const cu_real *__restrict__ y0,
                   const cu_real *__restrict__ dydt0, const cu_real dt,
                   const cu_real odeFactor, cu_real *__restrict__ y,
                   unsigned int tid, const unsigned int dimX);

    /**
     * \brief Solve the ODEs by giving the initial timestep
     *        and try to adjust the step according to the specified tolerance.
     * \param t - The time value.
     * \param dt - The initial timestep.
     * \param y - The solutions.
     */
    cu_device void
    adaptiveSolve(CUDAChemistrySourceODEs &odes, cu_real &t, cu_real &dt,
                  const cu_real *__restrict__ y0, cu_real *__restrict__ y,
                  cu_real *__restrict__ y2tmp, cu_real *__restrict__ dydt0,
                  cu_real *__restrict__ error, unsigned int tid,
                  const unsigned int dimX);

    /**
     * \brief Solve the ODEs by giving the initial timestep
     *        and try to adjust the step according to the specified tolerance.
     * \param t - The time value.
     * \param dt - The initial timestep.
     * \param y - The solutions.
     */
    cu_device void adaptiveSolveFactoring(
        CUDAChemistrySourceODEs &odes, cu_real &t, cu_real &dt,
        const cu_real odeFactor, const cu_real *__restrict__ y0,
        cu_real *__restrict__ y, cu_real *__restrict__ y2tmp,
        cu_real *__restrict__ dydt0, cu_real *__restrict__ error,
        unsigned int tid, const unsigned int dimX);

public:
    /*!
     * \brief Solve the ODEs by giving the initial timestep in a period time
     *        and try to adjust the step according to the specified tolerance.
     * \param[in] t0 - The initial time value.
     * \param[in] tn - The end time value.
     * \param[in] dt0 - The initial timestep.
     * \param[out] y - The solutions.
     */
    cu_device void
    solve(CUDAChemistrySourceODEs &odes, const cu_real t0, cu_real &tn,
          cu_real &dt0, cu_real *__restrict__ y0, cu_real *__restrict__ y,
          cu_real *__restrict__ y2tmp, cu_real *__restrict__ dydt0,
          cu_real *__restrict__ error, unsigned int tid,
          const unsigned int dimX);

    /*!
     * \brief Solve the ODEs by giving the initial timestep in a period time
     *        and try to adjust the step according to the specified tolerance.
     * \param[in] t0 - The initial time value.
     * \param[in] tn - The end time value.
     * \param[in] dt0 - The initial timestep.
     * \param[out] y - The solutions.
     */
    cu_device void
    solveFactoring(CUDAChemistrySourceODEs &odes, const cu_real t0,
                   cu_real &tn, cu_real &dt0, const cu_real odeFactor,
                   cu_real *__restrict__ y0, cu_real *__restrict__ y,
                   cu_real *__restrict__ y2tmp, cu_real *__restrict__ dydt0,
                   cu_real *__restrict__ error, unsigned int tid,
                   const unsigned int dimX);
};

} // namespace OpenHurricane

#include "EulerCUDA.inl"

#endif // CUDA_PARALLEL