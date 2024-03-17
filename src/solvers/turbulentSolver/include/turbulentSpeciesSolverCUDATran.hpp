/*!
 * \file turbulentSpeciesSolverCUDATran.hpp
 * \brief Headers of turbulent species transport Solver using CUDA accelerating
 *        calculation of chemical source terms and transport properties.
 *        The subroutines and functions are in the <i>turbulentSpeciesSolverCUDATran.cpp</i> file.
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
#include "combustionModel.hpp"
#include "finiteRate.hpp"
#include "solver.hpp"
#include "thermoListCUDA.hpp"
#include "transportListCUDA.hpp"

namespace OpenHurricane {
    /*!
     * \brief The class of turbulent species solver using CUDA accelerating
     *        calculation of chemical source terms and transport properties.
     */
    class turbulentSpeciesSolverCUDATran : public solver {
    protected:
        /*!\brief Delete the pointer.*/
        template <class Type> inline void deletePointer(geometryArray<Type, cellMesh> *_gPtr) const;

        uniquePtr<combustionModel> chemtryPtr_;

        integer rhoId_;
        integer rhouId_;
        integer rhoEId_;
        integer rhoYi0Id_;
        integer rhoTurb0Id_;

    public:
        /**
         * \brief The convective flux Jacobian.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[in] Vt - The contravariant velocity of the face of the control volume is set to zero for stationary grids
         * \param[out] AC - The convective flux Jacobian
         */
        virtual void Ac(const integer celli, const vector &normal, const real Vt,
                        realSquareMatrix &AC) const;

        /**
         * \brief The derivative of pressure with respect to the conservative variables.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[out] pq - The derivative of pressure with respect to the conservative variables
         */
        virtual void dpdq(const integer celli, const vector &normal, realArray &pq) const;

        /**
         * \brief The convective flux Jacobian.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[in] dq - The increment of the conservative variables
         * \param[out] adq - The convective flux Jacobian multiplied by the increment of the conservative variables
         */
        virtual void Acdq(const integer celli, const vector &normal, const realArray &dq,
                          realArray &adq) const;

    public:
        declareClassName(turbulentSpeciesSolverCUDATran);

        turbulentSpeciesSolverCUDATran(iteration &iter, const runtimeMesh &mesh);

        /*!\brief Destructor.*/
        virtual ~turbulentSpeciesSolverCUDATran() noexcept;

        virtual void solving();

        virtual void BDFSolve();

        /*!\brief Clear the solver.*/
        virtual void clear() noexcept;

        /**
         * \brief Update boundary conditions.
         */
        virtual void bc();

        void timeStep(realArray &dt);

        /**
         * \brief Update mixture physical and chemical properties.
         */
        virtual void updateProperties();

        /**
         * \brief Initialize the calculation.
         */
        virtual void initialize();

        /**
         * \brief Compute convective flux.
         */
        virtual void calculateFc();

        /**
         * \brief Compute viscous flux.
         */
        virtual void calculateFv();

        /**
         * \brief Compute viscous flux.
         */
        virtual void actualCalFv();

        /**
         * \brief Compute source terms.
         */
        virtual void calculateSource();

        virtual void updatePrimitives(const bool shouldUpdateTemp = false);

        virtual void calculateOtherTerms();

        virtual void updateFlowOld();

        /**
         * \brief Write the residuals and results.
         */
        virtual void write();

    protected:
        nThreadsAndBlocks nTBForTran_;
        nThreadsAndBlocks nTBForHaiCp_;

        void getYiTP(real *hur_restrict hostYiTPPtr_) const;
        void getYiTPSlice(real *hur_restrict hostYiTPPtr_, const integer pos,
                          const integer offset) const;

        mutable uniquePtr<integerList> fzBoundStartPtr_;

        hur_nodiscard inline const integerList &fzBoundStart() const;
        void makeFzBoundStart() const;

        void muKappaDiffBoundary();

        void getYiTPBoundary(real *hur_restrict hostYiTPPtr) const;

        void setExtrapolate(const real *hur_restrict hostDimmMuKap,
                            const real *hur_restrict hostHai);

        void transferTranP();

        /**
         * \brief Update gas properties by slice.
         */
        void gasPropertiesBySlice(const integer pos, const integer offset, const bool firstCall);

        void updatePropertiesBySlice(const integer nSlice);

        void chemicalSourceBySlice(const integer pos, const integer offset, const bool firstCall);
        void chemicalSourceImpBySlice(const integer pos, const integer offset, const bool firstCall,
                                      const bool isHybridPrecision);

        void updateChemSourceBySlice(const integer nSlice);
        void updateChemSourceImpBySlice(const integer nSlice, const bool isHybridPrecision);
#ifdef TEST_PROCESS_TIME
                
        class gasProProcessTime {
        public:
            // Unit: [s]

            real CPU_bc_gamma_E_;
            real CPU_sync_;
            real CPU_convert_data_;
            real CPU_bc_HaiMuKDim_;
            real GPU_HaiMuKDim_;
            real Total_calc_Propeeer_;
            real Total_before_calc_;
            real Total_GPUCal_to_sync_;
            real Total_begin_to_sync_;
            real Total_Fvturb_to_sync_;

            gasProProcessTime() = default;

            ~gasProProcessTime() noexcept {}

            inline void writeTecplotHead(fileOsstream &fos) const {
                if (HurMPI::master()) {
                    fos.os() << "variables = "
                                "\"step\",\t\"CPU_bc_gamma_E[s]\",\t\"CPU_sync[s]"
                                "\",\t\"CPU_convert_data[s]\""
                             << ",\t\"CPU_bc_HaiMuKDim[s]\",\t\"GPU_HaiMuKDim[s]\","
                             << "\t\"Total_calc_Propeeer[s]\","
                                "\t\"Total_before_calc[s]\""
                             << ",\t\"Total_GPUCal_to_sync[s]\",\t\"Total_begin_to_sync[s]\","
                                "\t\"Total_Fvturb_to_sync[s]\""
                             << std::endl;
                }
            }

            inline void writeTecplot(const integer cStep, fileOsstream &fos) const {
                real CPU_bc_gamma_E = CPU_bc_gamma_E_;
                real CPU_sync = CPU_sync_;
                real CPU_convert_data = CPU_convert_data_;
                real CPU_bc_HaiMuKDim = CPU_bc_HaiMuKDim_;
                real GPU_HaiMuKDim = GPU_HaiMuKDim_;
                real Total_calc_Propeeer = Total_calc_Propeeer_;
                real Total_before_calc = Total_before_calc_;
                real Total_GPUCal_to_sync = Total_GPUCal_to_sync_;
                real Total_begin_to_sync = Total_begin_to_sync_;
                real Total_Fvturb_to_sync = Total_Fvturb_to_sync_;
                HurMPI::allReduce(CPU_bc_gamma_E, MPI_MAX);
                HurMPI::allReduce(CPU_sync, MPI_MAX);
                HurMPI::allReduce(CPU_convert_data, MPI_MAX);
                HurMPI::allReduce(CPU_bc_HaiMuKDim, MPI_MAX);

                HurMPI::allReduce(GPU_HaiMuKDim, MPI_MAX);
                HurMPI::allReduce(Total_before_calc, MPI_MAX);
                HurMPI::allReduce(Total_GPUCal_to_sync, MPI_MAX);
                HurMPI::allReduce(Total_begin_to_sync, MPI_MAX);
                HurMPI::allReduce(Total_Fvturb_to_sync, MPI_MAX);

                if (HurMPI::master()) {
                    fos.setRealPrecision();
                    fos.os() << cStep << "\t" << CPU_bc_gamma_E << "\t" << CPU_sync << "\t"
                             << CPU_convert_data << "\t" << CPU_bc_HaiMuKDim << "\t"
                             << GPU_HaiMuKDim << "\t" << Total_calc_Propeeer
                             << "\t" << Total_before_calc << "\t"
                             << Total_GPUCal_to_sync << "\t" << Total_begin_to_sync << "\t"
                             << Total_Fvturb_to_sync << std::endl;
                    fos.unsetRealPrecision();
                }
            }

            void clear() noexcept {
                CPU_bc_gamma_E_ = 0;
                CPU_sync_ = 0;
                CPU_convert_data_ = 0;
                CPU_bc_HaiMuKDim_ = 0;
                GPU_HaiMuKDim_ = 0;
                Total_calc_Propeeer_ = 0;
                Total_before_calc_ = 0;
                Total_GPUCal_to_sync_ = 0;
                Total_begin_to_sync_ = 0;
                Total_Fvturb_to_sync_ = 0;
            }
        };
        fileOsstream gasProFOS_;
        gasProProcessTime gasProProcessTime_;

        class chemSorProcessTime {
        public:
            // Unit: [s]

            real CPU_FV_turbSource_;
            real CPU_sync_;
            real CPU_convert_data_;
            real GPU_ChemicalSource_;
            real Total_calc_source_;
            real Total_before_calc_;
            real Total_GPUCal_to_sync_;
            real Total_begin_to_sync_;
            real Total_Fvturb_to_sync_;

            chemSorProcessTime() = default;

            ~chemSorProcessTime() noexcept {}

            inline void writeTecplotHead(fileOsstream &fos) const {
                if (HurMPI::master()) {
                    fos.os() << "variables = "
                                "\"step\",\t\"CPU_FV_turbSource[s]\",\t\"CPU_sync[s]"
                                "\",\t\"CPU_convert_data[s]\""
                             << ",\t\"GPU_ChemicalSource[s]\",\t\"Total_calc_source[s]\","
                                "\t\"Total_before_calc[s]\""
                             << ",\t\"Total_GPUCal_to_sync[s]\",\t\"Total_begin_to_sync[s]\","
                                "\t\"Total_Fvturb_to_sync[s]\""
                             << std::endl;
                }
            }

            inline void writeTecplot(const integer cStep, fileOsstream &fos) const {
                real CPU_FV_turbSource = CPU_FV_turbSource_;
                real CPU_sync = CPU_sync_;
                real CPU_convert_data = CPU_convert_data_;
                real GPU_ChemicalSource = GPU_ChemicalSource_;
                real Total_calc_source = Total_calc_source_;
                real Total_before_calc = Total_before_calc_;
                real Total_GPUCal_to_sync = Total_GPUCal_to_sync_;
                real Total_begin_to_sync = Total_begin_to_sync_;
                real Total_Fvturb_to_sync = Total_Fvturb_to_sync_;
                HurMPI::allReduce(CPU_FV_turbSource, MPI_MAX);
                HurMPI::allReduce(CPU_sync, MPI_MAX);
                HurMPI::allReduce(CPU_convert_data, MPI_MAX);
                HurMPI::allReduce(GPU_ChemicalSource, MPI_MAX);

                HurMPI::allReduce(Total_calc_source, MPI_MAX);
                HurMPI::allReduce(Total_before_calc, MPI_MAX);
                HurMPI::allReduce(Total_GPUCal_to_sync, MPI_MAX);
                HurMPI::allReduce(Total_begin_to_sync, MPI_MAX);
                HurMPI::allReduce(Total_Fvturb_to_sync, MPI_MAX);

                if (HurMPI::master()) {
                    fos.setRealPrecision();
                    fos.os() << cStep << "\t" << CPU_FV_turbSource << "\t" << CPU_sync << "\t"
                             << CPU_convert_data << "\t" << GPU_ChemicalSource << "\t"
                             << Total_calc_source << "\t" << Total_before_calc << "\t"
                             << Total_GPUCal_to_sync << "\t" << Total_begin_to_sync << "\t"
                             << Total_Fvturb_to_sync << std::endl;
                    fos.unsetRealPrecision();
                }
            }

            void clear() noexcept {
                CPU_FV_turbSource_ = 0;
                CPU_sync_ = 0;
                CPU_convert_data_ = 0;
                GPU_ChemicalSource_ = 0;
                Total_calc_source_ = 0;
                Total_before_calc_ = 0;
                Total_GPUCal_to_sync_ = 0;
                Total_begin_to_sync_ = 0;
                Total_Fvturb_to_sync_ = 0;
            }
        };
        fileOsstream chemSorFOS_;
        chemSorProcessTime chemSorProcessTime_;

#endif // TEST_PROCESS_TIME
    };
} // namespace OpenHurricane

#include "turbulentSpeciesSolverCUDATran.inl"
#endif // CUDA_PARALLEL