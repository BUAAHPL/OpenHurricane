#include "CUDAChemistrySourceODEs.h"
/*!
 * \file CUDAChemistrySourceODEs.cu
 * \brief Main subroutines for only computing chemistry source ODEs in CUDA platform.
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

#ifdef CUDA_PARALLEL
#include "CUDAChemistrySourceODEs.h"

#include "CUDAFunctions.hpp"
#include "EulerCUDA.hpp"
#include "cudaReduction.hpp"
#include "cudaSort.hpp"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>

namespace OpenHurricane {
cu_device cu_real CUDAChemistrySourceODEs::Cp(
    const cuChem::speciesTable &spt, const unsigned short nsp,
    const cu_real Td, const cu_real *__restrict__ ci,
    cu_real *__restrict Cpi, unsigned int tid,
    const unsigned int dimX) const {
    auto j = tid;
    while (j < nsp) {
        Cpi[j] = ci[j] * spt.Cp0(Td, j);
        /*if (cellId == 0)
        {
                printf("Cpi[%d] = %.15e ci = %.15e cCp = %.15e\n", j, spt.Cp0(Td, j), ci[j], ci[j] * spt.Cp0(Td, j));
        }*/
        j += dimX;
    }
    __syncthreads();

    /*if (cellId == 0)
    {
            j = tid;
            while (j < nsp)
            {
                    printf("Cpi[%d] = %.15e\n", j, Cpi[j]);
                    j += dimX;
            }
    }
    __syncthreads();*/

    cudaReduction::reduce(Cpi, tid, nsp, dimX);

    return Cpi[0];
}

cu_device cu_real CUDAChemistrySourceODEs::Cv(
    const cuChem::speciesTable &spt, const unsigned short nsp,
    const cu_real Td, const cu_real *__restrict__ ci,
    cu_real *__restrict Cvi, unsigned int tid,
    const unsigned int dimX) const {
    auto j = tid;
    while (j < nsp) {
        Cvi[j] = ci[j] * spt.Cv0(Td, j);
        j += blockDim.x;
    }
    __syncthreads();
    cudaReduction::reduce(Cvi, tid, nsp, dimX);
    return Cvi[0];
}

cu_device cu_real CUDAChemistrySourceODEs::DTDt(
    const cuChem::speciesTable &spt, const unsigned short nsp,
    const cu_real Td, const cu_real *__restrict__ dcdt,
    cu_real *__restrict Hai, unsigned int tid,
    const unsigned int dimX) const {
    auto j = tid;
    while (j < nsp) {
        Hai[j] = dcdt[j] * spt.Ha0(Td, j);
        j += blockDim.x;
    }
    __syncthreads();
    cudaReduction::reduce(Hai, tid, nsp, dimX);
    return Hai[0];
}

cu_device cu_real CUDAChemistrySourceODEs::DTDtv(
    const cuChem::speciesTable &spt, const unsigned short nsp,
    const cu_real Td, const cu_real *__restrict__ dcdt,
    cu_real *__restrict Eai, unsigned int tid,
    const unsigned int dimX) const {
    auto j = tid;
    while (j < nsp) {
        Eai[j] = dcdt[j] * spt.Ea0(Td, j);
        j += blockDim.x;
    }
    __syncthreads();
    cudaReduction::reduce(Eai, tid, nsp, dimX);
    return Eai[0];
}

CUDAChemistrySourceODEs::CUDAChemistrySourceODEs(
    const cuChem::reactionTable &reactions, const short isContPreOrVol)
    : ODEsCUDA(), reactions_(reactions), isContPreOrVol_(isContPreOrVol) {}

CUDAChemistrySourceODEs::CUDAChemistrySourceODEs(
    const CUDAChemistrySourceODEs &cso)
    : ODEsCUDA(cso), reactions_(cso.reactions_),
      isContPreOrVol_(cso.isContPreOrVol_) {}

cu_device void CUDAChemistrySourceODEs::DcDtConstPressure(
    const cuChem::reactionTable &reactions, const unsigned short nsp,
    const unsigned short nrc, const cu_real *__restrict__ ci,
    cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai,
    unsigned int tid, const unsigned int dimX) {
    const auto Td = ci[0];
    ci += 1;

    omegai += 1;
    auto j = tid;
    while (j < nsp) {
        omegai[j] = 0;
        GdRT[j] = reactions.species().nGdRT(j, Td);
        j += dimX;
    }
    __syncthreads();

    j = tid;
    while (j < nrc) {
        reactions.omega(j, Td, ci, GdRT, omegai);
        j += dimX;
    }
    __syncthreads();

    const auto Cpp = Cp(reactions.species(), nsp, Td, ci, GdRT, tid, dimX);
    __syncthreads();
    const auto DTDtt =
        DTDt(reactions.species(), nsp, Td, omegai, GdRT, tid, dimX);
    ci -= 1;
    omegai -= 1;
    if (tid == 0) {
        omegai[0] = -DTDtt / Cpp;
    }
}

cu_device void CUDAChemistrySourceODEs::DcDtConstVolume(
    const cuChem::reactionTable &reactions, const unsigned short nsp,
    const unsigned short nrc, const cu_real *__restrict__ ci,
    cu_real *__restrict__ GdRT, cu_real *__restrict__ omegai,
    unsigned int tid, const unsigned int dimX) {
    const auto Td = ci[0];
    omegai += 1;
    ci += 1;
    auto j = tid;
    while (j < nsp) {
        omegai[j] = 0;
        GdRT[j] = reactions.species().nGdRT(j, Td);
        j += dimX;
    }
    __syncthreads();

    j = tid;
    while (j < nrc) {
        reactions.omega(j, Td, ci, GdRT, omegai);
        j += dimX;
    }
    __syncthreads();

    const auto Cvv = Cv(reactions.species(), nsp, Td, ci, GdRT, tid, dimX);
    __syncthreads();
    const auto DTDtt =
        DTDtv(reactions.species(), nsp, Td, omegai, GdRT, tid, dimX);
    __syncthreads();
    ci -= 1;
    omegai -= 1;
    if (tid == 0) {
        omegai[0] = -DTDtt / Cvv;
    }
}

cu_device void CUDAChemistrySourceODEs::timescale(
    const cu_real Td, const cu_real *__restrict__ ci,
    cu_real *__restrict__ GdRT, cu_real *__restrict__ tci, unsigned int tid,
    const unsigned int dimX) {
    const auto nsp = this->reactions_.nsp();
    const auto nrc = this->reactions_.nrc();
    auto j = tid;
    while (j < nsp) {
        GdRT[j] = reactions_.species().nGdRT(j, Td);
        j += dimX;
    }
    __syncthreads();

    j = tid;
    cu_real kff, krr;
    while (j < nrc) {
        reactions_.dWdciAllSpecies(j, Td, ci, GdRT, kff, krr, tci);
        j += dimX;
    }
    __syncthreads();

    j = tid;
    while (j < nsp) {
        tci[j] = cu_min(200, 1.0 / (fabs(tci[j]) + 1e-20));
        j += dimX;
    }
    __syncthreads();

    //cudaSort::oddEvenSort<cu_real>(tci, nsp, tid, 0);
}

cu_device void OpenHurricane::CUDAChemistrySourceODEs::MTS_whichGroup(
    const cu_real dtBase, const cu_real *__restrict__ tci,
    unsigned short *__restrict__ whichGrp, unsigned int tid,
    const unsigned int dimX) {
    const auto nsp = this->reactions_.nsp();
    auto j = tid;
    while (j < nsp) {
        whichGrp[j] = dtBase > tci[j]
                          ? cu_max(0, short(log10(dtBase / tci[j]))) + 1
                          : 1;
        j += dimX;
    }
    __syncthreads();
}

cu_device unsigned short CUDAChemistrySourceODEs::MTS_GroupDT(
    const cu_real dtBase, cu_real *__restrict__ tci,
    unsigned short *__restrict__ whichGrp, unsigned int tid,
    const unsigned int dimX) {
    const auto nsp = this->reactions_.nsp();
    cudaSort::oddEvenSort<cu_real>(tci, nsp, tid, 1);
    unsigned short Nm =
        dtBase > tci[0]
            ? cu_max(0, short(std::floor(log10(dtBase / tci[0])))) + 1
            : 1;

    auto i = tid;

    if (Nm <= dimX) {
        cu_real dtG = -1;
        if (i < Nm) {
            bool found = false;
            for (unsigned short j = 0; j < nsp; ++j) {
                const short ig =
                    Nm - (dtBase > tci[j]
                              ? cu_max(0, short(log10(dtBase / tci[j]))) + 1
                              : 1);
                if (!found && ig == i) {
                    dtG = tci[j];
                    break;
                }
            }
        }
        __syncthreads();

        if (i < Nm) {
            tci[i] = dtG;
        }
    } else {
        unsigned short ii = (Nm + dimX - 1) / dimX;

        cu_real *dtt = new cu_real[ii];
        unsigned short count = 0;
        while (i < Nm) {
            dtt[i] = -1;
            bool found = false;
            for (unsigned short j = 0; j < nsp; ++j) {
                const short ig =
                    Nm - (dtBase > tci[j]
                              ? cu_max(0, short(log10(dtBase / tci[j]))) + 1
                              : 1);
                if (!found && ig == i) {
                    dtt[count++] = tci[j];
                    break;
                }
            }
            i += dimX;
        }
        __syncthreads();
        i = tid;
        count = 0;
        while (i < Nm) {
            tci[i] = dtt[count++];
            i += dimX;
        }
        delete[] dtt;
    }

    return Nm;
}

cu_device void CUDAChemistrySourceODEs::DyDt(
    const cu_real t, const cu_real *__restrict__ y,
    cu_real *__restrict__ tmp, cu_real *__restrict__ dydt, unsigned int tid,
    const unsigned int dimX) {
    if (isContPreOrVol_ == 0) {
        DcDtConstPressure(reactions_, reactions_.nsp(), reactions_.nrc(), y,
                          tmp, dydt, tid, dimX);
    } else {
        DcDtConstVolume(reactions_, reactions_.nsp(), reactions_.nrc(), y, tmp,
                        dydt, tid, dimX);
    }
}

cu_device void CUDAChemistrySourceODEs::jacobian(
    const cu_real t, const cu_real *__restrict__ y,
    cu_real *__restrict__ dfdt, cu_real *__restrict__ dfdy,
    unsigned int tid, const unsigned int dimX) {}

__global__ void solveChemODEsTest(CUDAChemistrySourceODEs odes,
                                  EulerCUDA solver, cu1DArray<cu_real> tt,
                                  cu1DArray<cu_real> Tci0) {
    extern __shared__ cu_real myshare[];
    cu_real *Tci = myshare;
    cu_real *tmp = (cu_real *)&Tci[odes.nEqns()];
    cu_real *dcidt = (cu_real *)&tmp[odes.nEqns()];
    cu_real *y2tmp = (cu_real *)&dcidt[odes.nEqns()];

    auto tid = threadIdx.x;
    auto dimX = blockDim.x;
    auto t0 = tt(0);
    auto tn = tt(1);
    auto dt0 = tt(2);
    solver.solve(odes, t0, tn, dt0, Tci0.dDataPtr(), Tci, y2tmp, dcidt, tmp,
                 tid, dimX);
    tid = threadIdx.x;
    if (tid == 0) {
        tt(1) = tn;
        tt(2) = dt0;
    }
    /*while (tid < odes.nEqns())
    {
            Tci0(tid) = Tci[tid];
            tid += dimX;
    }*/
}

__global__ void solveChemODEsMTSTest(CUDAChemistrySourceODEs odes,
                                     EulerCUDA solver,
                                     cu1DArray<cu_real> tt,
                                     cu1DArray<cu_real> Tci0) {
    extern __shared__ cu_real myshare[];
    cu_real *tci = myshare;
    cu_real *tmp = (cu_real *)&tci[odes.nEqns()];
    cu_real *dcidt = (cu_real *)&tmp[odes.nEqns()];
    cu_real *y2tmp = (cu_real *)&dcidt[odes.nEqns()];
    unsigned short *whichGroup = (unsigned short *)&y2tmp[odes.nEqns()];

    auto t0 = tt(0);
    auto tn = tt(1);
    auto dt0 = tt(2);
    auto tid = threadIdx.x;
    auto dimX = blockDim.x;
    auto i = threadIdx.x;
    while (i < odes.nEqns()) {
        y2tmp[i] = Tci0(i);
        i += dimX;
    }

    odes.timescale(Tci0(0), Tci0.dDataPtr() + 1, tmp, tci, tid, dimX);
    odes.MTS_whichGroup(tn, tci, whichGroup, tid, dimX);
    const auto Nm = odes.MTS_GroupDT(tn, tci, whichGroup, tid, dimX);
    __syncthreads();

    cu_real dtDid = 0;

    const cu_real beta = 0.8;

    for (unsigned short im = 0; im < Nm; ++im) {
        while (dtDid < tn) {
            cu_real dt = beta * tci[im];
            if (dtDid + dt > tn) {
                dt = tn - dtDid;
            }

            odes.DyDt(0, y2tmp, tmp, dcidt, tid, dimX);
            i = threadIdx.x;
            while (i < odes.nEqns() - 1) {
                //printf("tid = %d, im = %d, whichGroup = %d\n", i, im, whichGroup[i]);
                if (Nm - whichGroup[i] < im) {
                    dcidt[i + 1] = 0;
                }
                i += dimX;
            }
            dtDid += dt;
            auto i = threadIdx.x;
            while (i < odes.nEqns()) {
                auto err = dt * dcidt[i];
                auto ay0 = fabs(y2tmp[i]);
                y2tmp[i] += err;
                auto tol = solver.ATOL() +
                           solver.RTOL() * cu_max(fabs(y2tmp[i]), ay0);
                //printf("tid = %d, im = %d, imm = %d, dcidt = %e, dt = %e, err = %e, tol = %e\n", i, im, Nm - whichGroup[i], dcidt[i], dt, err, tol);
                tmp[i] = i == 0 ? 0 : fabs(err) / tol;

                i += dimX;
            }
            if (im != Nm - 1) {
                i = threadIdx.x;
                while (i < odes.nEqns() - 1) {
                    if (Nm - whichGroup[i] != im) {
                        tmp[i + 1] = 0;
                    }
                    //printf("tid = %d, im = %d, imm = %d, maxErr = %e\n", i, im, Nm - whichGroup[i], tmp[i + 1]);
                    i += dimX;
                }
                __syncthreads();
                cudaReduction::reduceMax(tmp, tid, odes.nEqns(), dimX);
                //printf("tid = %d, im = %d, maxErr = %e\n", threadIdx.x, im, tmp[0]);
                if (tmp[0] < 1) {
                    break;
                }
            }
        }
    }

    tid = threadIdx.x;
    if (tid == 0) {
        tt(1) = dtDid;
        tt(2) = dt0;
    }
    while (tid < odes.nEqns()) {
        Tci0(tid) = y2tmp[tid];
        tid += dimX;
    }
}

void calcChemODEs(const cu_integer nsp, CUDAChemistrySourceODEs &odes,
                  EulerCUDA &solver, cu_real t0, cu_real &tn,
                  cu_real &dt0, cu_real *__restrict__ Tci) {
    cu_real *htt = nullptr;
    checkCUDAError(cudaHostAlloc((void **)&htt, 3 * sizeof(cu_real),
                                 cudaHostAllocDefault));
    htt[0] = t0;
    htt[1] = tn;
    htt[2] = dt0;

    cu1DArray<cu_real> dtt(3, htt);

    cu1DArray<cu_real> dTci(nsp + 1, Tci);

    unsigned int dimX = 0;
    if ((nsp + 1) < 64) {
        dimX = 32;
    } else if ((nsp + 1) < 128) {
        dimX = 64;
    } else {
        dimX = 128;
    }
    dim3 blockSize(dimX, 1, 1);
    dim3 gridSize(1, 1, 1);
    const auto sharedMemSize = 4 * (nsp + 1) * sizeof(cu_real);
    solveChemODEsTest<<<gridSize, blockSize, sharedMemSize, 0>>>(odes, solver,
                                                                 dtt, dTci);

    dtt.copyToHost(htt);
    dTci.copyToHost(Tci);

    tn = htt[1];
    dt0 = htt[2];

    checkCUDAError(cudaFreeHost(htt));
}

void calcChemODEsMTS(const cu_integer nsp, CUDAChemistrySourceODEs &odes,
                     EulerCUDA &solver, cu_real t0, cu_real &tn,
                     cu_real &dt0, cu_real *__restrict__ Tci) {
    cu_real *htt = nullptr;
    checkCUDAError(cudaHostAlloc((void **)&htt, 3 * sizeof(cu_real),
                                 cudaHostAllocDefault));
    htt[0] = t0;
    htt[1] = tn;
    htt[2] = dt0;

    cu1DArray<cu_real> dtt(3, htt);

    cu1DArray<cu_real> dTci(nsp + 1, Tci);

    unsigned int dimX = 0;
    if ((nsp + 1) < 64) {
        dimX = 32;
    } else if ((nsp + 1) < 128) {
        dimX = 64;
    } else {
        dimX = 128;
    }
    dim3 blockSize(dimX, 1, 1);
    dim3 gridSize(1, 1, 1);
    const auto sharedMemSize =
        4 * (nsp + 1) * sizeof(cu_real) + nsp * sizeof(unsigned short);
    solveChemODEsMTSTest<<<gridSize, blockSize, sharedMemSize, 0>>>(
        odes, solver, dtt, dTci);

    dtt.copyToHost(htt);
    dTci.copyToHost(Tci);

    tn = htt[1];
    dt0 = htt[2];

    checkCUDAError(cudaFreeHost(htt));
}

__global__ void yiToCi(cu2DArray<cu_real> TyiRho,
                       const cuChem::reactionTable reactions,
                       const cu_integer nCells) {
    const auto nsp = reactions.nsp();

    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;
    // The dimensional density of this cell [kg/m^3]
    cu_real rhod = 0;
    // To get the dimensional density.
    if (cellId < nCells) {
        rhod = TyiRho(cellId, nsp + 1);
    }
    auto i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            TyiRho(cellId, i + 1) =
                reactions.species().yiToci(rhod, TyiRho(cellId, i + 1), i);
            i += blockDim.x;
        }
    }
}

__global__ void ciToYi(cu2DArray<cu_real> TyiRho,
                       const cuChem::reactionTable reactions,
                       const cu_integer nCells) {
    extern __shared__ cu_real myshare[];
    cu_real *ciS = myshare;
    const auto nsp = reactions.nsp();

    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

    // To get the pointer of ci (the molar concentrations) of this cell from the allocated shared memory.
    cu_real *ci = (cu_real *)(ciS + threadIdx.y * nsp);

    auto i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            ci[i] = TyiRho(cellId, i + 1) * reactions.species().Wi(i);
            i += blockDim.x;
        }
    }
    __syncthreads();

    auto dimX = blockDim.x;
    i = threadIdx.x;
    cudaReduction::reduce(ci, i, nsp, dimX);

    i = threadIdx.x;
    if (cellId < nCells) {
        if (i == 0) {
            TyiRho(cellId, nsp + 1) = ci[0];
        }
    }

    i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            TyiRho(cellId, i + 1) *= (reactions.species().Wi(i) / ci[0]);
            i += blockDim.x;
        }
    }
    __syncthreads();
}

__global__ void solveChemODEs(cu2DArray<cu_real> TciRho,
                              CUDAChemistrySourceODEs odes, EulerCUDA solver,
                              cu1DArray<cu_real> dt,
                              cu1DArray<cu_real> subdt,
                              const cu_integer nCells) {
    extern __shared__ cu_real myshare[];
    auto dimX = blockDim.x;
    cu_real *Tcit = myshare;
    cu_real *tmpt = (cu_real *)&Tcit[blockDim.y * odes.nEqns()];
    cu_real *dcidtt = (cu_real *)&tmpt[blockDim.y * odes.nEqns()];
    cu_real *y2tmpt = (cu_real *)&dcidtt[blockDim.y * odes.nEqns()];

    cu_real *Tci = (cu_real *)(Tcit + threadIdx.y * odes.nEqns());
    cu_real *tmp = (cu_real *)(tmpt + threadIdx.y * odes.nEqns());
    cu_real *dcidt = (cu_real *)(dcidtt + threadIdx.y * odes.nEqns());
    cu_real *y2tmp = (cu_real *)(y2tmpt + threadIdx.y * odes.nEqns());

    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

    auto tid = threadIdx.x;

    cu_real t0 = 0;
    auto tn = dt(cellId);
    auto dt0 = subdt(cellId);
    auto leftTime = tn;
    while (leftTime > 0.0) {
        auto j = threadIdx.x;
        while (j < odes.nEqns()) {
            /*if (cellId == 0)
            {
                    printf("before Tci[%d] = %.15e tn = %e dt0 = %e\n", j, TciRho(cellId, j),tn, dt0);
            }*/
            j += dimX;
        }
        auto dtn = leftTime;
        solver.solve(odes, t0, dtn, dt0, TciRho.row(cellId), Tci, y2tmp, dcidt,
                     tmp, tid, dimX);
        leftTime -= dtn;

        auto i = threadIdx.x;
        while (i < odes.nEqns()) {
            /*if (cellId == 0)
            {
                    printf("Tci[%d] = %.15e tn = %e dt0 = %e leftTime = %e dcidt = %e\n", i, TciRho(cellId, i), tn, dt0, leftTime, dcidt[i]);
            }*/
            //TciRho(cellId, i) = Tci[i];
            i += dimX;
        }
    };
    if (tid == 0) {
        subdt(cellId) = dt0;
    }
}

__global__ void solveChemODEsFactoring(
    cu2DArray<cu_real> TciRho, CUDAChemistrySourceODEs odes,
    EulerCUDA solver, cu1DArray<cu_real> dt, cu1DArray<cu_real> subdt,
    cu1DArray<cu_real> odeFactor, const cu_integer nCells) {
    extern __shared__ cu_real myshare[];
    auto dimX = blockDim.x;
    cu_real *Tci = myshare;
    cu_real *tmp = (cu_real *)&Tci[odes.nEqns()];
    cu_real *dcidt = (cu_real *)&tmp[odes.nEqns()];
    cu_real *y2tmp = (cu_real *)&dcidt[odes.nEqns()];

    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

    auto tid = threadIdx.x;
    if (cellId < nCells) {
        cu_real t0 = 0;
        auto tn = dt(cellId);
        auto dt0 = subdt(cellId);
        auto leftTime = tn;
        while (leftTime > 0.0) {
            auto dtn = leftTime;
            solver.solveFactoring(odes, t0, dtn, dt0, odeFactor(cellId),
                                  TciRho.row(cellId), Tci, y2tmp, dcidt, tmp,
                                  tid, dimX);
            leftTime -= dtn;
        };
        if (tid == 0) {
            subdt(cellId) = dt0;
        }
    }

    tid = threadIdx.x;
    if (cellId < nCells) {
        while (tid < odes.nEqns()) {
            TciRho(cellId, tid) = Tci[tid];
            tid += dimX;
        }
    }
}

__global__ void solveChemODEsMTS(cu2DArray<cu_real> TciRho,
                                 CUDAChemistrySourceODEs odes, EulerCUDA solver,
                                 cu1DArray<cu_real> dt,
                                 cu1DArray<cu_real> subdt,
                                 const cu_integer nCells) {
    extern __shared__ cu_real myshare[];

    cu_real *tcit = myshare;
    cu_real *tmpt = (cu_real *)&tcit[blockDim.y * odes.nEqns()];
    cu_real *dcidtt = (cu_real *)&tmpt[blockDim.y * odes.nEqns()];
    cu_real *y2tmpt = (cu_real *)&dcidtt[blockDim.y * odes.nEqns()];
    unsigned short *whichGroupt =
        (unsigned short *)&y2tmpt[blockDim.y * odes.nEqns()];

    cu_real *tci = (cu_real *)(tcit + threadIdx.y * odes.nEqns());
    cu_real *tmp = (cu_real *)(tmpt + threadIdx.y * odes.nEqns());
    cu_real *dcidt = (cu_real *)(dcidtt + threadIdx.y * odes.nEqns());
    cu_real *y2tmp = (cu_real *)(y2tmpt + threadIdx.y * odes.nEqns());
    unsigned short *whichGroup =
        (unsigned short *)(whichGroupt + threadIdx.y * odes.nEqns());

    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

    //cu_real t0 = 0;
    auto tn = dt(cellId); //auto dt0 = subdt(cellId);
    auto tid = threadIdx.x;
    auto dimX = blockDim.x;
    auto i = threadIdx.x;
    while (i < odes.nEqns()) {
        y2tmp[i] = TciRho.row(cellId)[i];
        i += dimX;
    }

    odes.timescale(TciRho.row(cellId)[0], TciRho.row(cellId) + 1, tmp, tci, tid,
                   dimX);
    odes.MTS_whichGroup(tn, tci, whichGroup, tid, dimX);
    const auto Nm = odes.MTS_GroupDT(tn, tci, whichGroup, tid, dimX);
    __syncthreads();

    cu_real dtDid = 0;
    const cu_real beta = 0.8;
    for (unsigned short im = 0; im < Nm; ++im) {
        while (dtDid < tn) {
            cu_real dt = beta * tci[im];
            if (dtDid + dt > tn) {
                dt = tn - dtDid;
            }
            odes.DyDt(0, y2tmp, tmp, dcidt, tid, dimX);
            i = threadIdx.x;
            while (i < odes.nEqns() - 1) {
                if (Nm - whichGroup[i] < im) {
                    dcidt[i + 1] = 0;
                }
                i += dimX;
            }
            dtDid += dt;
            auto i = threadIdx.x;
            while (i < odes.nEqns()) {
                auto err = dt * dcidt[i];
                auto ay0 = fabs(y2tmp[i]);
                y2tmp[i] += err;
                auto tol = solver.ATOL() +
                           solver.RTOL() * cu_max(fabs(y2tmp[i]), ay0);
                tmp[i] = i == 0 ? 0 : fabs(err) / tol;
                i += dimX;
            }
            if (im != Nm - 1) {
                i = threadIdx.x;
                while (i < odes.nEqns() - 1) {
                    if (Nm - whichGroup[i] != im) {
                        tmp[i + 1] = 0;
                    }
                    i += dimX;
                }
                __syncthreads();
                cudaReduction::reduceMax(tmp, tid, odes.nEqns(), dimX);
                if (tmp[0] < 1) {
                    break;
                }
            }
        }
    }

    /*if (tid == 0)
    {
            subdt(cellId) = dt0;
    }*/

    tid = threadIdx.x;
    while (tid < odes.nEqns()) {
        TciRho(cellId, tid) = y2tmp[tid];
        tid += dimX;
    }
}

void getBlockAndGridSizeODEs(const cu_integer nsp, const cu_integer nCells,
                             dim3 &blockSize, dim3 &gridSize) {
    unsigned int blocksize_x, blocksize_y;
    unsigned int gridsize_x;
    if (nsp < 8) {
        blocksize_x = 4;
        blocksize_y = 32;
        //blocksize_y = 64;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 16) {
        blocksize_x = 8;
        blocksize_y = 16;
        //blocksize_y = 32;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 32) {
        blocksize_x = 16;
        blocksize_y = 8;
        //blocksize_y = 16;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 64) {
        blocksize_x = 32;
        blocksize_y = 4;
        //blocksize_y = 8;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else if (nsp < 128) {
        blocksize_x = 64;
        blocksize_y = 2;
        //blocksize_y = 4;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    } else {
        blocksize_x = 128;
        blocksize_y = 1;
        //blocksize_y = 2;
        gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
    }

    dim3 blockSizeTmp(blocksize_x, blocksize_y, 1);
    dim3 gridSizeTmp(gridsize_x, 1, 1);

    blockSize = blockSizeTmp;
    gridSize = gridSizeTmp;
}

void calcChemODEsArray(cu_real *__restrict__ hostTYiRho,
                       const cu_integer nsp, const cu_integer nCells,
                       CUDAChemistrySourceODEs &odes, EulerCUDA &solver,
                       const nThreadsAndBlocks &nTBa,
                       const cu_real *__restrict__ hostDt,
                       cu_real *__restrict__ hostSubDt, const bool returnCi) {
    unsigned int nn = nTBa.gridSize().x * nTBa.blockSize().y;
    cu2DArray<cu_real> dTYiRho(nsp + 2, nn, hostTYiRho);
    cu1DArray<cu_real> dDt(nn, hostDt);
    cu1DArray<cu_real> dSubDt(nn, hostSubDt);

    yiToCi<<<nTBa.gridSize(), nTBa.blockSize(), 0, 0>>>(
        dTYiRho, odes.reactions(), nCells);

    const auto sharedMemSize1 =
        nTBa.blockSize().y * 4 * (nsp + 1) * sizeof(cu_real);
    solveChemODEs<<<nTBa.gridSize(), nTBa.blockSize(), sharedMemSize1, 0>>>(
        dTYiRho, odes, solver, dDt, dSubDt, nCells);

    if (!returnCi) {
        const auto sharedMemSize2 =
            nTBa.blockSize().y * nsp * sizeof(cu_real);
        ciToYi<<<nTBa.gridSize(), nTBa.blockSize(), sharedMemSize2, 0>>>(
            dTYiRho, odes.reactions(), nCells);
    }

    dTYiRho.copyToHost(hostTYiRho);
    dSubDt.copyToHost(hostSubDt);
    dTYiRho.clear();
    dDt.clear();
    dSubDt.clear();
}

void calcChemODEsArrayMTS(cu_real *__restrict__ hostTYiRho,
                          const cu_integer nsp, const cu_integer nCells,
                          CUDAChemistrySourceODEs &odes, EulerCUDA &solver,
                          const nThreadsAndBlocks &nTBa,
                          const cu_real *__restrict__ hostDt,
                          cu_real *__restrict__ hostSubDt,
                          const bool returnCi) {
    unsigned int nn = nTBa.gridSize().x * nTBa.blockSize().y;
    cu2DArray<cu_real> dTYiRho(nsp + 2, nn, hostTYiRho);
    cu1DArray<cu_real> dDt(nn, hostDt);
    cu1DArray<cu_real> dSubDt(nn, hostSubDt);

    yiToCi<<<nTBa.gridSize(), nTBa.blockSize(), 0, 0>>>(
        dTYiRho, odes.reactions(), nCells);

    const auto sharedMemSize1 =
        nTBa.blockSize().y *
        (4 * (nsp + 1) * sizeof(cu_real) + nsp * sizeof(unsigned short));
    solveChemODEsMTS<<<nTBa.gridSize(), nTBa.blockSize(), sharedMemSize1, 0>>>(
        dTYiRho, odes, solver, dDt, dSubDt, nCells);

    if (!returnCi) {
        const auto sharedMemSize2 =
            nTBa.blockSize().y * nsp * sizeof(cu_real);
        ciToYi<<<nTBa.gridSize(), nTBa.blockSize(), sharedMemSize2, 0>>>(
            dTYiRho, odes.reactions(), nCells);
    }

    dTYiRho.copyToHost(hostTYiRho);
    dSubDt.copyToHost(hostSubDt);
    dTYiRho.clear();
    dDt.clear();
    dSubDt.clear();
}

void calcChemODEsArray(cu_real *__restrict__ hostTYiRho,
                       const cu_integer nsp, const cu_integer nCells,
                       CUDAChemistrySourceODEs &odes, EulerCUDA &solver,
                       const cu_real *__restrict__ hostDt,
                       cu_real *__restrict__ hostSubDt,
                       const cu_real *__restrict__ hostOdeFactor) {
    cu2DArray<cu_real> dTYiRho(nsp + 2, nCells, hostTYiRho);
    cu1DArray<cu_real> dDt(nCells, hostDt);
    cu1DArray<cu_real> dSubDt(nCells, hostSubDt);
    cu1DArray<cu_real> dOdeFactor(nCells, hostOdeFactor);

    dim3 blockSize;
    dim3 gridSize;
    getBlockAndGridSizeODEs(nsp, nCells, blockSize, gridSize);
    yiToCi<<<gridSize, blockSize, 0, 0>>>(dTYiRho, odes.reactions(), nCells);

    const auto sharedMemSize1 = blockSize.y * 4 * (nsp + 1) * sizeof(cu_real);
    solveChemODEsFactoring<<<gridSize, blockSize, sharedMemSize1, 0>>>(
        dTYiRho, odes, solver, dDt, dSubDt, dOdeFactor, nCells);

    const auto sharedMemSize2 = blockSize.y * nsp * sizeof(cu_real);
    ciToYi<<<gridSize, blockSize, sharedMemSize2, 0>>>(
        dTYiRho, odes.reactions(), nCells);

    dTYiRho.copyToHost(hostTYiRho);
    dSubDt.copyToHost(hostSubDt);
    dTYiRho.clear();
    dDt.clear();
    dSubDt.clear();
    dOdeFactor.clear();
}
} // namespace OpenHurricane

#endif // CUDA_PARALLEL