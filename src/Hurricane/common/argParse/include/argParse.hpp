/*!
 * \file argParse.hpp
 * \brief Header of argv parsing
 *       The subroutines and functions are in the <i>argParse.cpp</i> file.
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

#include "stdMaps.hpp"
#include "fileName.hpp"
#include "logFile.hpp"
#include "Lists.hpp"
#include "version.hpp"

namespace OpenHurricane {
    /**
     * \brief The namespace of argument options.
     */
    namespace argOptions {
        /**
         * \brief The enum of argument options.
         */
        enum options : short {
            NON = 0,     //!< None argument option.
            HELP,        //!< Help option.
            DATAFORMAT,  //!< The format of data to be output.
            CASE,        //!< To specify the case file.
            DATA,        //!< To specify the data file
            CD,          //!< To specify the case and data file which have the same name.
            LOGFILE,     //!< To create log file option.
            VERSION,     //!< To print the version.
            INPUT,       //!< To specify the input file.
            INPUTPATH,   //!< To specify the inputpath.
            OUT,         //!< To specify the result output file name.
            OUTPUTPATH,  //!< To specify the result output path.
            RESIDUAL,    //!< To specify the residual results output file name.
            GUI,         //!< To specify the GUI mode, that is the program is running at GUI mode.
            nGPU,        //!< To specify the GPU number per machine.
            GPUZEROCOPY, //!< To specify mapping CPU host memory to CUDA GPU.
            TECPLOTFILEWRITEONMASTER, //!< To specify if tecplot file written by master process.
            noWRITEAtEND, //!< To specify do not write the output and relay files at the end of the computation.
            SKEWNESSCORRECTION, //!< To specify the implementation of a correction for the interpolated value at the face centre when the face is skewed.
            CHECKWRITESOLUTION //!< Should check the solution for writting out is NaN or INF?. Default is false.
        };
    } // End namespace argOptions

    /**
     * \brief The class of argument parsing.
     */
    class argParse {
    public:
        using mapType = std::map<argOptions::options, std::string>;

    private:
        stdStringList args_;
        mapType options_;

        static string executable_;
        std::string argListStr_;

        static fileName caseFileName_;
        static fileName dataFileName_;
        static fileName outputFileName_;
        static fileName residualFileName_;

        static fileName executePath_;
        static fileName inputPath_;
        static fileName outputPath_;

        static bool isASCII_;

        /**\brief If the program is under GUI mode.*/
        static bool isGUIMode_;

        /**\brief If the program is using GPU.*/
        static bool hasGPU_;

        /**\brief The number of GPUs per machine.*/
        static int nGPUPerMachine_;

        void defaultName() {
            caseFileName_ = "OpenHurricane.cont";
            dataFileName_ = "OpenHurricane.h5";
            outputFileName_ = "HurricaneOut.dat";
            residualFileName_ = "HurricaneResidual.dat";
        }

        void defaultName(const string &prefix) {
            caseFileName_ = prefix + ".cont";
            dataFileName_ = prefix + ".h5";
            outputFileName_ = prefix + "Out.dat";
            residualFileName_ = prefix + "Residual.dat";
        }

        void printHelp() const;
        void printVersion() const;

        // Static parameters

        static bool gpuZeroCopy_;
        static bool tecplotFileWriteOnMaster_;

        static bool noWriteAtEnd_;
        static bool isSkewCorrect_;

        /**
         * \brief Should check the solution for writting out is NaN or INF?. Default is false.
         */
        static bool checkWriteSolution_;

    public:
        static std::map<std::string, argOptions::options> validOptions;

        static void addValidOptions(const std::string &opt, const argOptions::options param);

        //! \cond internalClass
        class initValidTables {
        public:
            initValidTables();
        };

        /*!\brief Construct from argc and argv.
         */
        argParse(int &argc, char **&argv);

        /*!\brief Construct from argc and argv.
         */
        argParse(int &argc, char **&argv, const char *prefix);

        /*!\brief Construct from argc and argv.
         */
        argParse(int &argc, char **&argv, const string &prefix);

        /**\brief Destructor*/
        virtual ~argParse() noexcept;

        // Member Functions

        void parse();

        void clear() noexcept;

        // Access

        static hur_nodiscard inline const string &executable() noexcept { return executable_; }

        hur_nodiscard inline const std::string &argListStr() noexcept { return argListStr_; }

        static hur_nodiscard inline const fileName &caseFileName() noexcept {
            return caseFileName_;
        }
        static hur_nodiscard inline const fileName &dataFileName() noexcept {
            return dataFileName_;
        }
        static hur_nodiscard inline const fileName &outputFileName() noexcept {
            return outputFileName_;
        }
        static hur_nodiscard inline const fileName &residualFileName() noexcept {
            return residualFileName_;
        }

        static hur_nodiscard inline const fileName &executePath() noexcept { return executePath_; }
        static hur_nodiscard inline const fileName &inputPath() noexcept { return inputPath_; }
        static hur_nodiscard inline const fileName &outputPath() noexcept { return outputPath_; }

        static hur_nodiscard inline bool isASCII() noexcept { return isASCII_; }

        /**\brief If the program is under GUI mode.*/
        static hur_nodiscard inline bool isGUIMode() noexcept { return isGUIMode_; }

        /**\brief If the program is using GPU.*/
        static hur_nodiscard inline bool hasGPU() noexcept { return hasGPU_; }

        /**\brief The number of GPUs per machine.*/
        static hur_nodiscard inline int nGPUPerMachine() noexcept { return nGPUPerMachine_; }

        // Static functions

        static hur_nodiscard inline bool gpuZeroCopy() noexcept { return gpuZeroCopy_; }

        static hur_nodiscard inline bool tecplotFileWriteOnMaster() noexcept {
            return tecplotFileWriteOnMaster_;
        }

        static hur_nodiscard inline bool noWriteAtEnd() noexcept { return noWriteAtEnd_; }

        static hur_nodiscard inline bool isSkewCorrect() noexcept { return isSkewCorrect_; }

        /**
         * \brief Should check the solution for writting out is NaN or INF?. Default is false.
         */
        static hur_nodiscard inline bool checkWriteSolution() noexcept {
            return checkWriteSolution_;
        }
    };

} // namespace OpenHurricane
