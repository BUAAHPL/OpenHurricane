/*!
 * \file chemSolvers.cpp
 * \brief Main subroutines for chemical solvers main program.
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

#include "commonInclude.hpp"

#include "dataStructure.hpp"
#include <iostream>

#include "chemSolver.hpp"
#include "chemSolvers.hpp"
#include <iostream>
#include <string>

using namespace OpenHurricane;

int main(int argc, char *argv[]) {
    // Initializing the program
    initializing(&argc, &argv);
    argParse myArg(argc, argv, "chemSolvers");
    programName::printVersion();
    std::string timeName = "     Starting time: ";
    timeName += clock::dateTime();

    //testCUDA();
    // Create the iteration of the whole computation
    Pout << "\n    Info: Creating iteration...\n" << std::endl;
    iteration CFDIteration("CFDIter", myArg, contStr);
    myArg.clear();

    // Create mesh for computation
    runtimeMesh gm(object(CFDIteration.name(), CFDIteration, object::WRITE_RELAY_OUTPUT),
                   CFDIteration.cont(), meshStr);
    CFDIteration.readAndAddCont(CFDIteration.configName());

    uniquePtr<solver> mySolver(new chemSolver(CFDIteration, gm));

    // Solve the equations
    Pout << "    Begin iteration..." << std::endl;
    mySolver->solve();

    Pout << "    Info: calculation finishes at " << clock::dateTime() << std::endl;

    // Finish the calculation
    finish();
    return 0;
}
