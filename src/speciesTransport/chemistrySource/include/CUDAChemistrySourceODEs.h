/*!
 * \file CUDAChemistrySourceODEs.h
 * \brief Header of only computing chemistry source ODEs in CUDA platform.
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
#include "CUDAFunctions.hpp"
#include "ODEsCUDA.hpp"
#include "cudaEvents.hpp"
#include "cudaStreams.hpp"
#include "nThreadsAndBlocks.hpp"
#include "reactionTableCUDA.hpp"

namespace OpenHurricane {
class EulerCUDA;

/**
 * \brief The class of computing chemistry source ODEs in CUDA platform.
 */
class CUDAChemistrySourceODEs : public ODEsCUDA {
private:
    cuChem::reactionTable reactions_;

    /**
     * \brief 0 - is const pressure.
     * 1 - is const volume.
     */
    short isContPreOrVol_;

    /*!\brief Molar heat capacity at constant pressure [J/(kmol K)].*/
    cu_device cu_real Cp(const cuChem::speciesTable &spt,
                                  const unsigned short nsp, const cu_real Td,
                                  const cu_real *__restrict__ ci,
                                  cu_real *__restrict Cpi, unsigned int tid,
                                  const unsigned int dimX) const;

    /*!\brief Molar heat capacity at constant volume [J/(kmol K)].*/
    cu_device cu_real Cv(const cuChem::speciesTable &spt,
                                  const unsigned short nsp, const cu_real Td,
                                  const cu_real *__restrict__ ci,
                                  cu_real *__restrict Cvi, unsigned int tid,
                                  const unsigned int dimX) const;

    cu_device cu_real DTDt(const cuChem::speciesTable &spt,
                                    const unsigned short nsp,
                                    const cu_real Td,
                                    const cu_real *__restrict__ dcdt,
                                    cu_real *__restrict Hai, unsigned int tid,
                                    const unsigned int dimX) const;

    cu_device cu_real DTDtv(const cuChem::speciesTable &spt,
                                     const unsigned short nsp,
                                     const cu_real Td,
                                     const cu_real *__restrict__ dcdt,
                                     cu_real *__restrict Eai,
                                     unsigned int tid,
                                     const unsigned int dimX) const;

public:
    /**
     * \brief Disallow null constructor.
     */
    cu_dual CUDAChemistrySourceODEs() = delete;

    cu_host
    CUDAChemistrySourceODEs(const cuChem::reactionTable &reactions,
                            const short isContPreOrVol);

    /**
     * \brief Copy constructor.
     */
    cu_dual CUDAChemistrySourceODEs(const CUDAChemistrySourceODEs &cso);

    /**
     * \brief Destructor.
     */
    inline cu_dual virtual ~CUDAChemistrySourceODEs() noexcept {}

    /**\brief Return the number of equations in the system.*/
    inline cu_dual cu_integer nEqns() const noexcept;

    inline cu_dual const cuChem::reactionTable &
    reactions() const noexcept;

    cu_device void DcDtConstPressure(
        const cuChem::reactionTable &reactions, const unsigned short nsp,
        const unsigned short nrc, const cu_real *__restrict__ ci,
        cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai,
        unsigned int tid, const unsigned int dimX);

    cu_device void DcDtConstVolume(
        const cuChem::reactionTable &reactions, const unsigned short nsp,
        const unsigned short nrc, const cu_real *__restrict__ ci,
        cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai,
        unsigned int tid, const unsigned int dimX);

    cu_device void timescale(const cu_real Td,
                                    const cu_real *__restrict__ ci,
                                    cu_real *__restrict__ GdRT,
                                    cu_real *__restrict__ tci,
                                    unsigned int tid, const unsigned int dimX);

    cu_device void MTS_whichGroup(const cu_real dtBase,
                                         const cu_real *__restrict__ tci,
                                         unsigned short *__restrict__ whichGrp,
                                         unsigned int tid,
                                         const unsigned int dimX);

    cu_device unsigned short
    MTS_GroupDT(const cu_real dtBase, cu_real *__restrict__ tci,
                unsigned short *__restrict__ whichGrp, unsigned int tid,
                const unsigned int dimX);

    /**\brief Calculate the derivatives in dydt.*/
    cu_device void DyDt(const cu_real t,
                               const cu_real *__restrict__ y,
                               cu_real *__restrict__ tmp,
                               cu_real *__restrict__ dydt, unsigned int tid,
                               const unsigned int dimX);

    /**
     *\brief Calculate the Jacobian of the system y' = f(x,y)
     *  Need by the stiff-system solvers.
     */
    cu_device void jacobian(const cu_real t,
                                   const cu_real *__restrict__ y,
                                   cu_real *__restrict__ dfdt,
                                   cu_real *__restrict__ dfdy,
                                   unsigned int tid, const unsigned int dimX);
};

void calcChemODEs(const cu_integer nsp, CUDAChemistrySourceODEs &odes,
                  EulerCUDA &solver, cu_real t0, cu_real &tn,
                  cu_real &dt0, cu_real *__restrict__ Tci);

void calcChemODEsMTS(const cu_integer nsp, CUDAChemistrySourceODEs &odes,
                     EulerCUDA &solver, cu_real t0, cu_real &tn,
                     cu_real &dt0, cu_real *__restrict__ Tci);

void calcChemODEsArray(cu_real *__restrict__ hostTYiRho,
                       const cu_integer nsp, const cu_integer nCells,
                       CUDAChemistrySourceODEs &odes, EulerCUDA &solver,
                       const nThreadsAndBlocks &nTB,
                       const cu_real *__restrict__ hostDt,
                       cu_real *__restrict__ hostSubDt,
                       const bool returnCi = false);

void calcChemODEsArrayMTS(cu_real *__restrict__ hostTYiRho,
                          const cu_integer nsp, const cu_integer nCells,
                          CUDAChemistrySourceODEs &odes, EulerCUDA &solver,
                          const nThreadsAndBlocks &nTB,
                          const cu_real *__restrict__ hostDt,
                          cu_real *__restrict__ hostSubDt,
                          const bool returnCi = false);

void calcChemODEsArray(cu_real *__restrict__ hostTYiRho,
                       const cu_integer nsp, const cu_integer nCells,
                       CUDAChemistrySourceODEs &odes, EulerCUDA &solver,
                       const cu_real *__restrict__ hostDt,
                       cu_real *__restrict__ hostSubDt,
                       const cu_real *__restrict__ hostOdeFactor);

} // namespace OpenHurricane

#include "CUDAChemistrySourceODEs.inl"

#endif // CUDA_PARALLEL