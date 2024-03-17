/*!
 * \file iteration.hpp
 * \brief Header of CFD time advance iteration
 *       The subroutines and functions are in the <i>iteration.cpp</i> file.
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
#include "argParse.hpp"
#include "clock.hpp"
#include "controller.hpp"
#include "dataStructure.hpp"
#include "fileOsstream.hpp"
#include "physicalTimeStep.hpp"
#include "registerTable.hpp"
#include "subIteration.hpp"
#include "version.hpp"

#include "ArrayInclude.hpp"

#include "referenceValues.hpp"

namespace OpenHurricane {
    class mixture;
    class flowModel;
    class geometryMesh;
    class faceZone;
    class monitors;
    class writeFaceZone;

    class solutionWrite;

    /*!\brief The class of CFD time advance iteration.*/
    class iteration : public registerTable {

    public:
        enum class flowStateType : short { steady, unsteady };

        enum class timestepType : short { localTimeStep, globalTimeStep };

    private:
        // Private data

        /**\brief The computation controller.*/
        controller cont_;

        /*!\brief Name of the iteration.*/
        string name_;

    public:
        /*!\brief Name of the config file.*/
        hur_nodiscard inline const fileName &configName() const noexcept {
            return argParse::caseFileName();
        }

        hur_nodiscard inline const fileName &caseFileName() const noexcept {
            return argParse::caseFileName();
        }

        hur_nodiscard inline const fileName &dataName() const noexcept {
            return argParse::dataFileName();
        }

        /*!\brief Name of the output file.*/
        hur_nodiscard inline const fileName &outputName() const noexcept {
            return argParse::outputFileName();
        }

    private:
        uniquePtr<monitors> myMonitorPtr_;

        /**
         * \brief The pointer of solution IO.
         */
        uniquePtr<solutionWrite> solWritePtr_;

        /**\brief write field variables.*/
        uniquePtr<writeFaceZone> writeFaceZonePtr_;

    private:
        mutable bool isResConvergence_;

    public:
        hur_nodiscard inline const fileName &executePath() const noexcept {
            return argParse::executePath();
        }

        hur_nodiscard inline const fileName &inputPath() const noexcept {
            return argParse::inputPath();
        }
        hur_nodiscard inline const fileName &outputPath() const noexcept {
            return argParse::outputPath();
        }

        /**\brief If the program is under GUI mode.*/
        hur_nodiscard inline bool isGUIMode() const noexcept { return argParse::isGUIMode(); }

        /*!\brief Set isResConvergence_ true.*/
        inline void setResConveg() const;

    private:
        /**\brief Current iteration step.*/
        integer cStep_;

        /**\brief Tha max steps that the time advance should be iterated.*/
        integer maxStep_;

        /**\brief Total steps that the time advance has been iterated.*/
        integer totalStep_;

        integer writeOutputStep_;

        bool restart_;
        bool restartFromUnsteady_;
        fileName restartFrom_;

    public:
        hur_nodiscard inline bool restart() const noexcept;
        hur_nodiscard inline bool restartFromUnsteady() const noexcept;

    private:
        /*!\brief Physical time step pointer.*/
        mutable uniquePtr<physicalTimeStep> pTStepPtr_;

        /*!\brief Sub-iteration pointer.*/
        mutable uniquePtr<subIteration> subIterPtr_;

        /** \brief Reference values. */
        referenceValues refValues_;

        flowStateType flowState_;

        timestepType timestepType_;

        /**\brief Should read last time field from relay file.*/
        bool readLastFromRelay_;

        /**\brief Should write last time field to relay file.*/
        bool writeLastToRelay_;

    public:
        hur_nodiscard static const iteration &nullObject();

        iteration(const iteration &) = delete;
        iteration &operator=(const iteration &) = delete;

        iteration(const char *_c, const argParse &arg);

        iteration(const char *_c, const argParse &arg, const std::string &contStr);

        /*!\brief Destructor.*/
        virtual ~iteration() noexcept;

        /**\brief The computation controller.*/
        hur_nodiscard inline const controller &cont() const noexcept { return cont_; }

        /**\brief The computation controller.*/
        hur_nodiscard inline controller &cont() noexcept { return cont_; }

        /*!\brief Name of the iteration.*/
        hur_nodiscard inline const string &name() const noexcept;

        /*!\brief Name of the program.*/
        hur_nodiscard inline const string &pName() const noexcept { return programName::name; }

        hur_nodiscard fileName meshFileName() const;
        hur_nodiscard fileName getRelayFileName() const;
        hur_nodiscard fileName getCaseConfigName() const;

        /*!\brief Changing name of the iteration.*/
        inline void changeName(const string &nN);

        /*!\brief Current iteration step.*/
        hur_nodiscard inline integer cStep() const noexcept;

        /*!\brief Tha max steps that the time advance should be iterated.*/
        hur_nodiscard inline integer maxStep() const noexcept;

        /*!\brief Total steps that the time advance has been iterated.*/
        hur_nodiscard inline integer totalStep() const noexcept;

        /*!\brief Total steps that the time advance has been iterated.*/
        inline integer setTotalStep(const integer newTS) noexcept;

        hur_nodiscard inline integer writeOutputStep() const noexcept;

        /*!\brief If the iteration has physical time step.*/
        hur_nodiscard inline bool hasPhysicalTimeStep() const noexcept;

        /*!\brief If the iteration has sub-iteration.*/
        hur_nodiscard inline bool hasSubIteration() const noexcept;

        /*!\brief Physical time step.*/
        hur_nodiscard const physicalTimeStep &pTStep() const;

        /*!\brief Physical time step.*/
        hur_nodiscard physicalTimeStep &pTStep();

        /*!\brief Sub-iteration.*/
        hur_nodiscard const subIteration &subIter() const;

        hur_nodiscard inline const referenceValues &refValues() const noexcept;

        hur_nodiscard inline referenceValues &refValues() noexcept;

        hur_nodiscard inline flowStateType flowState() const noexcept { return flowState_; }

        hur_nodiscard inline bool isSteadyFlow() const noexcept;
        hur_nodiscard inline bool isUnsteadyFlow() const noexcept;

        hur_nodiscard inline timestepType timeStepType() const noexcept { return timestepType_; }

        hur_nodiscard inline bool isLoacalTimeStep() const noexcept {
            return timestepType_ == timestepType::localTimeStep;
        }

        inline void setLoacalTimeStep() noexcept { timestepType_ = timestepType::localTimeStep; }

        hur_nodiscard inline bool isGlobalTimeStep() const noexcept {
            return timestepType_ == timestepType::globalTimeStep;
        }

        inline void setGlobalTimeStep() noexcept { timestepType_ = timestepType::globalTimeStep; }

        /*!\brief If the iteration is ended.*/
        hur_nodiscard bool end() const noexcept;

        /*!\brief Begin sub-iteration loop.*/
        hur_nodiscard inline bool subLoop() noexcept;

        /*!\brief reset sub-iteration.*/
        inline void resetSubiter();

        /*!\brief The steps left.*/
        hur_nodiscard inline integer leftStep() const noexcept;

        /*\!brief Reset the max physical time steps.*/
        void setMaxTime(const real mT) noexcept;

        /*!\brief Reset the physical time step.*/
        void setTimeStep(const real tS) noexcept;

        inline void setWriteLastToRelay(const bool isw) noexcept;
        inline void setReadLastFromRelay(const bool isw) noexcept;

        /**\brief Should read last time field from relay file.*/
        hur_nodiscard inline bool isReadLastFromRelay() const noexcept;

        /**\brief Should write last time field to relay file.*/
        hur_nodiscard inline bool isWriteLastToRelay() const noexcept;

        /*!\brief Clear the iteration.*/
        void clear() noexcept;

        /*!\brief Iterating.*/
        hur_nodiscard bool iterating() noexcept;

        /**
         * \brief Read from controller.
         */
        void readCont();

        /**
         * \brief Read controller from string.
         */
        void readCont(std::string contStr);

        /**
         * \brief Read controller from file: contN.
         */
        void readAndAddCont(const fileName &contN);

    private:
        void getIterSteps(const controller &interCont);
        void getFlowState(const controller &interCont);
        void getTimeStepType(const controller &interCont);
        void getPhysicalTimeStep(const controller &interCont);
        void getSubIterationStep(const controller &interCont);
        void getRestartFrom(const controller &interCont);
        void getProfile(controller &interCont);

        void writeDualTimePhyInfo() const;

    public:
        // Read
        void readRelay();

        /*!\brief Write relay file and output file.*/
        void write() const;

        /*!\brief Write mesh file.*/
        void writeMesh(const fileName &meshFileName) const;

        /*!\brief Write mesh and case controller to file.*/
        void writeCaseConfig(const fileName &meshFileName) const;

    public:
        void writeRelay(const fileName &relayFileName, const fileName &meshFileName) const;

        void writeResult() const;

        void writeFieldVarResult(const fileName &outN, const stringList &varList) const;

        /*!\brief set the name of output file.*/
        hur_nodiscard fileName setFileName(const string &type) const;

    public:
        /*!\brief Write residuals.*/
        void writeSubIterResiduals() const;

        /*!\brief Write mesh file to tecplot.*/
        void writeMeshToTecplot() const;

        hur_nodiscard monitors &myMonitor() noexcept;
        hur_nodiscard solutionWrite &solWrite() noexcept;
        hur_nodiscard const solutionWrite &solWrite() const noexcept;

    private:
        void setMyMonitorPtr(monitors *newPtr);
        void setSolWrite(const flowModel &flows);
        void setWriteFaceZonePtr(writeFaceZone *newPtr);

    public:
        void initializing(const flowModel &flows);

        void monitoring() const;
        void subMonitoring() const;

    public:
        hur_nodiscard inline bool isBeginLimit() const;

        iteration &operator++() noexcept;
        inline iteration &operator++(int) noexcept { return operator++(); }
    };
} // namespace OpenHurricane

#include "iteration.inl"