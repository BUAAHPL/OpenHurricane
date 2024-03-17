/*!
 * \file solverSelector.cpp
 * \brief creator for solver.
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

#include "solver.hpp"

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::solver>
OpenHurricane::solver::creator(iteration &iter, const runtimeMesh &mesh) {
    string solverType;
    const auto &cont = iter.cont();
    const auto &flowCont = cont.subController("flow");
    string flowModel = flowCont.findWord("flowModel");
    if (flowModel == "EulerFlow") {
        const auto &mixtureCont = flowCont.subController("mixture");
        string cheicalModel = mixtureCont.findWord("chemical");

        if (cheicalModel == "singleSp") {
            solverType = "EulerSolver";
        } else if (cheicalModel == "mixing" || cheicalModel == "reaction") {
            solverType = "EulerSpeciesSolver";
        }
    } else {
        const auto &mixtureCont = flowCont.subController("mixture");
        string cheicalModel = mixtureCont.findWord("chemical");

        if (cheicalModel == "singleSp") {
            if (flowModel == "laminarFlow") {
                solverType = "laminarPerfectGasSolver";
            } else if (flowModel == "eddyViscosity") {
                solverType = "turbulentPerfectGasSolver";
            }
        } else if (cheicalModel == "mixing" || cheicalModel == "reaction") {
            string chemistryOption = "coupled";
            if (mixtureCont.found("chemistryOption")) {
                chemistryOption = mixtureCont.findWord("chemistryOption");
            }
            if (flowModel == "laminarFlow") {
                solverType = "laminarSpeciesSolver";
#ifdef CUDA_PARALLEL

                if (mixtureCont.found("reactions")) {
                    const auto &reactionCont = mixtureCont.subController("reactions");
                    if (reactionCont.found("chemistrySource")) {
                        if (reactionCont.isController("chemistrySource")) {
                            const auto &chemCont = reactionCont.subController("chemistrySource");
                            if (chemCont.found("chemistrySource")) {
                                const auto cw = chemCont.findWord("chemistrySource");
                                if (cw == "CUDAChemistrySourceNoSolver") {
                                    if (HurGPU::useGPU()) {
                                        solverType = "laminarSpeciesSolverCUDATran";
                                    } else {
                                        const_cast<controller &>(chemCont).remove(
                                            "chemistrySource");
                                        const_cast<controller &>(chemCont).addWord(
                                            "chemistrySource", "chemistryNoSolver");
                                    }
                                }
                            }
                        } else {
                            const auto cw = reactionCont.findWord("chemistrySource");
                            if (cw == "CUDAChemistrySourceNoSolver") {
                                if (HurGPU::useGPU()) {
                                    solverType = "laminarSpeciesSolverCUDATran";
                                } else {
                                    const_cast<controller &>(reactionCont)
                                        .remove("chemistrySource");
                                    const_cast<controller &>(reactionCont)
                                        .addWord("chemistrySource", "chemistryNoSolver");
                                }
                            }
                        }
                    }
                }

#endif // CUDA_PARALLEL
            } else if (flowModel == "eddyViscosity") {
                if (chemistryOption == "coupled") {
                    solverType = "turbulentSpeciesSolver";
                } else if (chemistryOption == "splitted") {
                    solverType = "turbulentSpeciesSplitSolver";
                } else {
                    errorAbortStr(("Unknown chemistry option: " + chemistryOption));
                }
#ifdef CUDA_PARALLEL

                if (mixtureCont.found("reactions")) {
                    const auto &reactionCont = mixtureCont.subController("reactions");
                    if (reactionCont.found("chemistrySource")) {
                        if (reactionCont.isController("chemistrySource")) {
                            const auto &chemCont = reactionCont.subController("chemistrySource");
                            if (chemCont.found("chemistrySource")) {
                                const auto cw = chemCont.findWord("chemistrySource");
                                if (cw == "CUDAChemistrySourceNoSolver") {
                                    if (HurGPU::useGPU()) {
                                        solverType = "turbulentSpeciesSolverCUDATran";
                                    } else {
                                        const_cast<controller &>(chemCont).remove(
                                            "chemistrySource");
                                        const_cast<controller &>(chemCont).addWord(
                                            "chemistrySource", "chemistryNoSolver");
                                    }
                                }
                            }
                        } else {
                            const auto cw = reactionCont.findWord("chemistrySource");
                            if (cw == "CUDAChemistrySourceNoSolver") {
                                if (HurGPU::useGPU()) {
                                    solverType = "turbulentSpeciesSolverCUDATran";
                                } else {
                                    const_cast<controller &>(reactionCont)
                                        .remove("chemistrySource");
                                    const_cast<controller &>(reactionCont)
                                        .addWord("chemistrySource", "chemistryNoSolver");
                                }
                            }
                        }
                    }
                }

#endif // CUDA_PARALLEL
            }
        }
    }
    Pout << "    Info: setting solver: " << solverType << std::endl;
    defineInObjCreator(solver, solverType, controller, (iter, mesh));
}