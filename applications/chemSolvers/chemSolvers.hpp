/*!
 * \file chemSolvers.hpp
 * \brief Header of chemistry solvers
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

#include <string>
namespace OpenHurricane {
std::string meshStr("(0 \"Created by : Fluent Interface 20.2.0\")\n"
                    "(0 \"Units: Meters\")\n"
                    "(2 3)\n"
                    "(0 \"Node Section\")\n"
                    "(10 (0 1 8 0 3))\n"
                    "(10 (5 1 8 1 3)\n"
                    "	(			\n"
                    "		0 0 0	\n"
                    "		1 0 0	\n"
                    "		0 1 0	\n"
                    "		1 1 0	\n"
                    "		0 0 1	\n"
                    "		1 0 1	\n"
                    "		0 1 1	\n"
                    "		1 1 1	\n"
                    "		))		\n"
                    "(12 (0 1 1 0 0))\n"
                    "(12 (6 1 1 1 4))\n"
                    "(13 (0 1 6 0 0))\n"
                    "(13 (0 1 6 0 0))\n"
                    "(0 \"Faces of zone WALL\")\n"
                    "(13 (7 1 6 3 4)(\n"
                    "	1 3 7 5 1 0	 \n"
                    "	1 5 6 2 1 0	 \n"
                    "	1 2 4 3 1 0	 \n"
                    "	5 7 8 6 1 0	 \n"
                    "	3 4 8 7 1 0	 \n"
                    "	2 6 8 4 1 0	 \n"
                    "	)			 \n"
                    "	)			 \n"
                    "(0 \"Zone Sections\")\n"
                    "(39 (6 fluid FLUID)())\n"
                    "(39 (7 wall WALL)())");

std::string contStr(
    "<Config>"
    "<meshSet>"
    "	<meshFile>fluent.msh</meshFile>"
    "    <meshFormat>fluent</meshFormat>"
    "	<unit>mm</unit>"
    "	<reordering>RCMReordering</reordering>"
    "</meshSet>"
    "<iteration>"
    "	<maxStep>200000</maxStep>"
    "	<totalStep>0</totalStep>"
    "	<writeOutputStep>1</writeOutputStep>"
    "	<flowState>steady</flowState>"
    "	<pseudoTime>"
    "		<timeStepMethod>cellBasedMethod</timeStepMethod>"
    "		<CFLSetting>"
    "			<cflMin>1.1</cflMin>"
    "			<cflMax>10.0</cflMax>"
    "			<stepForPrintCFL>20</stepForPrintCFL>"
    "			<CFLType>linearCFL</CFLType>"
    "			<linearCFL>"
    "				<cflConst>200</cflConst>"
    "				<cflLinear>400</cflLinear>"
    "			</linearCFL>"
    "		</CFLSetting>"
    "		<timeStepType>localTimeStep</timeStepType>"
    "		<CForTimeStep>2.0</CForTimeStep>"
    "	</pseudoTime>"
    "	<timeMethod>"
    "		<timeMarching>LUSGS</timeMarching>"
    "		<LUSGS>"
    "			<omegaForLUSGS>1.050</omegaForLUSGS>"
    "			<betaT>0.2</betaT>"
    "		</LUSGS>"
    "	</timeMethod>"
    "</iteration>"
    "<spatialScheme>"
    "	<spatialScheme>upwindSchemes</spatialScheme>"
    "	<upwindSchemes>"
    "		<upwindScheme>HLLC</upwindScheme>"
    "		<AUSMPlusUP>"
    "			<Kp>0.25</Kp>"
    "			<Ku>0.75</Ku>"
    "			<sigma>1.0</sigma>"
    "			<beta>0.125</beta>"
    "		</AUSMPlusUP>"
    "	</upwindSchemes>"
    "	<isLowMachCorrection>off</isLowMachCorrection>"
    "	<reconstruction>firstOrder</reconstruction>"
    "	<gradient>cellGaussGrad</gradient>"
    "</spatialScheme>"
    "<flow>"
    "	<flowModel>EulerFlow</flowModel>"
    "   <Prl>0.72</Prl>"
    "   <limits>EulerFlow"
    "		<limitMethod>directly</limitMethod>EulerFlow"
    "        <THigh>4000.0</THigh>"
    "        <TLow>1</TLow>"
    "        <pHigh>6930000</pHigh>"
    "        <pLow>0.10</pLow>"
    "    </limits>"
    "    <mixture>"
    "		<equationOfState>perfectGas</equationOfState>"
    "		<chemical>reaction</chemical>"
    "		<thermo>"
    "			<type>JANAF</type>"
    "		</thermo>"
    "		<transport>"
    "			<type>kinetic</type>"
    "		</transport>"
    "		<reactions>"
    "			<chemistrySource>"
    "				<solver>"
    "					<ODEsSolver>SIEuler</ODEsSolver>"
    "					<maxstep>10000</maxstep>"
    "					<atol>1e-11</atol>"
    "					<rtol>1e-4</rtol>"
    "				</solver>"
    "				"
    "<chemistrySource>chemistryODEs</chemistrySource>"
    "				<chemicalMaxStep>1e10</chemicalMaxStep>"
    "				<chemistryReduction>"
    "					"
    "<chemistryReduction>DRG</chemistryReduction>"
    "				</chemistryReduction>"
    "			</chemistrySource>"
    "		</reactions>"
    "	</mixture>"
    "</flow>"
    "<boundaryCondition>"
    "	<WALL>"
    "        <bcType>wall</bcType>"
    "        <momentum>"
    "            <shearCondition>invSlip</shearCondition>"
    "        </momentum>"
    "        <thermal>    "
    "            <thermalCondition>adiabatic</thermalCondition>"
    "        </thermal>"
    "    </WALL>"
    "</boundaryCondition>"
    "<ref>"
    "	<length>1.0</length> "
    "	<area>1</area> "
    "	<density>1.225</density> "
    "	<enthalpy>0.0</enthalpy> "
    "	<pressure>0.0</pressure>"
    "	<temperature>288.16</temperature>"
    "	<velocity>1</velocity>"
    "	<viscosity>1.7894e-5</viscosity>"
    "	<specificHeatRatio>1.4</specificHeatRatio>"
    "	<gasConstantNumber>287.06</gasConstantNumber>"
    "</ref>"
    "<freeStream>"
    "	<specifiedType>auto</specifiedType>"
    "</freeStream>"
    "<initialization>"
    "    <initFromBoundary>WALL</initFromBoundary>"
    "</initialization>"
    "<writeControl>"
    "	<chemSolverOut>"
    "	<writeSpecies>massFraction</writeSpecies>"
    "	</chemSolverOut>"
    "</writeControl>"
    "</Config>");

} // namespace OpenHurricane