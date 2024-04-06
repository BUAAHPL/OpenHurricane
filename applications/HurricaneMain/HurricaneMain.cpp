/*!
 * \file HurricaneMain.cpp
 * \brief Main subroutines for OpenHurricane main program.
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

#include "HurricaneMain.hpp"
#include "clock.hpp"

using namespace OpenHurricane;

int main(int argc, char *argv[]) {
    // Initializing the program
    initializing(&argc, &argv);

    argParse myArg(argc, argv);
    greeting(myArg);

    // Create the iteration of the whole computation
    PLInfo("\n    Info: Creating iteration...\n\n");
    iteration CFDIteration("CFDIter", myArg);

    myArg.clear();

    // Create mesh for computation
    runtimeMesh gm(object(CFDIteration.name(), CFDIteration, object::WRITE_RELAY_OUTPUT),
                   CFDIteration.cont());

    uniquePtr<solver> mySolver = solver::creator(CFDIteration, gm);

    // Solve the equations
    mySolver->solve();
    PLInfo("    Info: calculation finishes at %s\n", clock::dateTime().c_str());

    // Finish the calculation
    finish();
    return 0;
}
