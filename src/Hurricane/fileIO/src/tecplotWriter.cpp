/*!
 * \file tecplotWriter.cpp
 * \brief Main subroutines of tecplot writer.
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
#include "tecplotWriter.hpp"

const int OpenHurricane::tecplotWriter::intergerValue_ = 1;
const char OpenHurricane::tecplotWriter::magicNumber_[] = "#!TDV112";
const float OpenHurricane::tecplotWriter::EOHMARKER_ = 357.0f;
const float OpenHurricane::tecplotWriter::ZONEMARKER_ = 299.0f;

void OpenHurricane::tecplotWriter::checkNumVar(const integer nVar1, const stdStringList &namelist) {
    if (nVar1 != namelist.size()) {
        LFatal("The number of variables is not euqal to the size of name list");
    }
}

void OpenHurricane::tecplotWriter::writeASCIIToBinary(const char *str, fileOsstream &fos) {
    int value = 0;
    while ((*str) != '\0') {
        value = (int)*str;
        fos.write(reinterpret_cast<const char *>(&value), sizeof(int));
        str++;
    }

    char null_char[] = "";
    value = (int)*null_char;
    fos.write(reinterpret_cast<const char *>(&value), sizeof(int));
}

void OpenHurricane::tecplotWriter::writeInit(fileOsstream &fos) {
    //------i. Magic number, Version number
    fos.write(magicNumber_, 8);
    //------ii. Integer value of 1.
    fos.write(reinterpret_cast<const char *>(&intergerValue_), sizeof(intergerValue_));
}

void OpenHurricane::tecplotWriter::writeFileHeader(fileOsstream &fos, const std::string &title,
                                               const integer nVar, const stdStringList &varName,
                                               const tecplotFormat::FILETYPE fileType) {
    //==============Header Secontion =================//
    //------i. Magic number, Version number
    //------ii. Integer value of 1.
    writeInit(fos);

    /*
    iii. Title and variable names.
    +-----------+
    | INT32		|		FileType: 0 = FULL,
    +-----------+				  1 = GRID,
                                                              2 = SOLUTION
    +-----------+
    | INT32*N	|		The TITLE.
    +-----------+
    +-----------+
    | INT32		|		Number of variables (NumVar) in the datafile.
    +-----------+
    +-----------+
    | INT32*N	|		Variable names.
    +-----------+				N = L[1] + L[2] + .... L[NumVar]
                                            where:
                                            L[i] = length of the ith variable name + 1
                                            (for the terminating 0 value).
    */
    int FileType = 0;
    switch (fileType) {
    case tecplotFormat::FULL:
        FileType = 0;
        break;
    case tecplotFormat::GRID:
        FileType = 1;
        break;
    case tecplotFormat::SOLUTION:
        FileType = 2;
        break;
    default:
        break;
    }
    fos.write(reinterpret_cast<const char *>(&FileType), sizeof(FileType));

    writeASCIIToBinary(title.c_str(), fos);

    int NumVar = (int)nVar;
    fos.write(reinterpret_cast<const char *>(&NumVar), sizeof(NumVar));
    checkNumVar(nVar, varName);

    for (integer i = 0; i < nVar; ++i) {
        writeASCIIToBinary(varName[i].c_str(), fos);
    }
}

void OpenHurricane::tecplotWriter::writeZoneHeader(fileOsstream &fos, const geometryMesh &mesh,
                                               const real solutionTime, const integer nVar) {
    /*
            +-----------+
            | FLOAT32	| Zone marker. Value = 299.0
            +-----------+
    */
    fos.write(reinterpret_cast<const char *>(&ZONEMARKER_), sizeof(ZONEMARKER_));

    /*
            +-----------+
            | INT32*N	|	Zone name. (See note 1.)
            +-----------+	        N = (length of zone name) + 1.
    */
    std::string zoneName = "process ";
    zoneName += toString(HurMPI::getProcRank());
    writeASCIIToBinary(zoneName.c_str(), fos);

    /*
            +-----------+
            | INT32		|	 ParentZone: Zero-based zone number within this
            +-----------+			     datafile to which this zone is
                                                                     a child.
    */
    int parantZone = -1;
    fos.write(reinterpret_cast<const char *>(&parantZone), sizeof(parantZone));

    /*
            +-----------+
            | INT32		|	StrandID: -2 = pending strand ID for assignment
            +-----------+			  by Tecplot
                                                             -1 = static strand ID
                                                          0 <= N < 32700 valid strand ID
    */
    int StrandID = -1;
    fos.write(reinterpret_cast<const char *>(&StrandID), sizeof(StrandID));

    /*	+-----------+
            | FLOAT64	| Solution time.
            +-----------+
    */
    double sT = (double)solutionTime;
    fos.write(reinterpret_cast<const char *>(&sT), sizeof(sT));

    /*
            +-----------+
            | INT32		| Not used. Set to -1.
            +-----------+
    */
    int notUsed = -1;
    fos.write(reinterpret_cast<const char *>(&notUsed), sizeof(notUsed));

    /*
    *	+-----------+
            | INT32		|	 ZoneType 0=ORDERED,		1=FELINESEG,
            +-----------+			  2=FETRIANGLE,		3=FEQUADRILATERAL,
                                                              4=FETETRAHEDRON,	5=FEBRICK,
                                                              6=FEPOLYGON,		7=FEPOLYHEDRON
    */
    int ZoneType = 5;
    fos.write(reinterpret_cast<const char *>(&ZoneType), sizeof(ZoneType));

    /*
            +-----------+
            | INT32		|	Data packing.
            +-----------+	     0 = Block
                                                     1 = Point
    */
    /*int DataPacking = 0;
    fos.write
    (
            reinterpret_cast<const char*>(&DataPacking),
            sizeof(DataPacking)
    );*/

    /*
            +-----------+
            | INT32		|	Specify Var Location.
            +-----------+			0 = Don¡¯t specify, all data is located
                                                            at the nodes.
                                                            1 = Specify
    */
    int SpecifyVarLocation = 1;
    fos.write(reinterpret_cast<const char *>(&SpecifyVarLocation), sizeof(SpecifyVarLocation));

    /*
            if ¡°specify var location¡± == 1
                    +-----------+
                    | INT32*NV	|	Variable Location (only specify if above is 1).
                    +-----------+	0 = Node, 1 = Cell Centered (See note 5.)
    */
    int VariableLocation = 1;
    for (integer i = 0; i < nVar; ++i) {
        if (i < 3) {
            VariableLocation = 0;
        } else {
            VariableLocation = 1;
        }
        fos.write(reinterpret_cast<const char *>(&VariableLocation), sizeof(VariableLocation));
    }
    /*
            +-----------+
            | INT32		|	Are raw local 1-to-1 face neighbors supplied?
            +-----------+	(0=FALSE 1=TRUE). These raw values are a
                                            compact form of the local 1-to-1 face neighbors.
                                            If supplied, Tecplot assumes that the face
                                            neighbors are fully specified. As such, it
                                            will not perform auto face neighbor assignment.
                                            This improves Tecplot¡¯s time to first plot.
                                            See the data section below for format details.
                                            ORDERED and FELINESEG zones must specify 0 for
                                            this value because raw face neighbors are not
                                            defined for these zone types. FEPOLYGON and
                                            FEPOLYHEDRON zones must specify 0 for this value
                                            since face neighbors are defined in the face map
                                            for these zone types.
    */
    int rawLocal = 0;
    fos.write(reinterpret_cast<const char *>(&rawLocal), sizeof(rawLocal));

    /*
            +-----------+
            | INT32		|	Number of miscellaneous user-defined face
            +-----------+	neighbor connections (value >= 0). This value
                                            is in addition to the face neighbors supplied
                                            in the raw section. FEPOLYGON and FEPOLYHEDRON
                                            zones must specify 0.
    */
    int nFaceNeighbor = 0;
    if (HurMPI::parRun()) {
        nFaceNeighbor = mesh.globalMeshInfo().globalFaceIndeces().sharedFaceSize();
    }
    fos.write(reinterpret_cast<const char *>(&nFaceNeighbor), sizeof(nFaceNeighbor));

    if (HurMPI::parRun()) {
        /*
        if ¡°number of miscellaneous user-defined
                face neighbor connections¡± != 0
                +-----------+
                | INT32		|	User defined face neighbor mode
                +-----------+	(0=Local 1-to-1, 1=Local 1-to-many,
                                                2=Global 1-to-1, 3=Global 1-to-many)
                if FE Zone:
                +-----------+
                | INT32		|	Indicates if the finite element face neighbors
                +-----------+	are completely specified by the miscellaneous
                                                face neighbors given: (0=NO, 1=YES). If yes,
                                                then Tecplot will not perform auto assignment
                                                of face neighbors otherwise all faces not
                                                specified are considered boundaries. If no,
                                                then Tecplot will perform auto-assignment of
                                                the face neighbors unless the raw face neighbor
                                                array was supplied. This option is not valid
                                                for ORDERED zones.
        */
        int neighborMode = 2;
        fos.write(reinterpret_cast<const char *>(&neighborMode), sizeof(neighborMode));
        int isGiven = 0;
        fos.write(reinterpret_cast<const char *>(&isGiven), sizeof(isGiven));
    }

    /*
    if FE Zone:
            +-----------+
            | INT32		|	 NumPts
            +-----------+
    */
    int nPs = mesh.nPoints();
    fos.write(reinterpret_cast<const char *>(&nPs), sizeof(nPs));

    /*
    +-----------+
    | INT32		|	NumElements
    +-----------+
    +-----------+
    | INT32*3	|	ICellDim,JCellDim,
    +-----------+	KCellDim (for future use; set to zero)
    */
    int NumElements = mesh.nCells();
    fos.write(reinterpret_cast<const char *>(&NumElements), sizeof(NumElements));
    int ICellDim = 0;
    fos.write(reinterpret_cast<const char *>(&ICellDim), sizeof(ICellDim));
    int JCellDim = 0;
    fos.write(reinterpret_cast<const char *>(&JCellDim), sizeof(JCellDim));
    int KCellDim = 0;
    fos.write(reinterpret_cast<const char *>(&KCellDim), sizeof(KCellDim));

    /*
    For all zone types (repeat for each Auxiliary data name/value pair):
    +-----------+
    | INT32 | 1=Auxiliary name/value pair to follow
    +-----------+ 0=No more Auxiliary name/value pairs.
    */
    int AuxiliaryName = 0;
    fos.write(reinterpret_cast<const char *>(&AuxiliaryName), sizeof(AuxiliaryName));
}

void OpenHurricane::tecplotWriter::writeZoneHeader(fileOsstream &fos, const string &zoneName,
                                               const real solutionTime, const integer nVar,
                                               const int ip, const integer nPts, const integer nCls,
                                               const integer nFNgh, const int ZoneType,
                                               const int StrandID) {
    /*
            +-----------+
            | FLOAT32	| Zone marker. Value = 299.0
            +-----------+
    */
    fos.write(reinterpret_cast<const char *>(&ZONEMARKER_), sizeof(ZONEMARKER_));
    /*
            +-----------+
            | INT32*N	|	Zone name. (See note 1.)
            +-----------+	        N = (length of zone name) + 1.
    */
    /*std::string zoneName = "process ";
    zoneName += toString(ip);*/
    writeASCIIToBinary(zoneName.c_str(), fos);

    /*
            +-----------+
            | INT32		|	 ParentZone: Zero-based zone number within this
            +-----------+			     datafile to which this zone is
                                                                     a child.
    */
    int parantZone = -1;
    fos.write(reinterpret_cast<const char *>(&parantZone), sizeof(parantZone));

    /*
            +-----------+
            | INT32		|	StrandID: -2 = pending strand ID for assignment
            +-----------+			  by Tecplot
                                                             -1 = static strand ID
                                                              0 <= N < 32700 valid strand ID
    */
    //int StrandID = -1;
    fos.write(reinterpret_cast<const char *>(&StrandID), sizeof(StrandID));

    /*	+-----------+
            | FLOAT64	| Solution time.
            +-----------+
    */
    double sT = (double)solutionTime;
    fos.write(reinterpret_cast<const char *>(&sT), sizeof(sT));

    /*
            +-----------+
            | INT32		| Not used. Set to -1.
            +-----------+
    */
    int notUsed = -1;
    fos.write(reinterpret_cast<const char *>(&notUsed), sizeof(notUsed));

    /*
    *	+-----------+
            | INT32		|	 ZoneType 0=ORDERED,		1=FELINESEG,
            +-----------+			  2=FETRIANGLE,		3=FEQUADRILATERAL,
                                                              4=FETETRAHEDRON,	5=FEBRICK,
                                                              6=FEPOLYGON,		7=FEPOLYHEDRON
    */
    //int ZoneType = 5;
    fos.write(reinterpret_cast<const char *>(&ZoneType), sizeof(ZoneType));

    /*
            +-----------+
            | INT32		|	Data packing.
            +-----------+	     0 = Block
                                                     1 = Point
    */
    /*int DataPacking = 0;
    fos.write
    (
            reinterpret_cast<const char*>(&DataPacking),
            sizeof(DataPacking)
    );*/

    /*
            +-----------+
            | INT32		|	Specify Var Location.
            +-----------+			0 = Don¡¯t specify, all data is located
                                                            at the nodes.
                                                            1 = Specify
    */
    int SpecifyVarLocation = 1;
    fos.write(reinterpret_cast<const char *>(&SpecifyVarLocation), sizeof(SpecifyVarLocation));

    /*
            if ¡°specify var location¡± == 1
                    +-----------+
                    | INT32*NV	|	Variable Location (only specify if above is 1).
                    +-----------+	0 = Node, 1 = Cell Centered (See note 5.)
    */
    int VariableLocation = 1;
    for (integer i = 0; i < nVar; ++i) {
        if (i < 3) {
            VariableLocation = 0;
        } else {
            VariableLocation = 1;
        }
        fos.write(reinterpret_cast<const char *>(&VariableLocation), sizeof(VariableLocation));
    }

    /*
            +-----------+
            | INT32		|	Are raw local 1-to-1 face neighbors supplied?
            +-----------+	(0=FALSE 1=TRUE). These raw values are a
                                            compact form of the local 1-to-1 face neighbors.
                                            If supplied, Tecplot assumes that the face
                                            neighbors are fully specified. As such, it
                                            will not perform auto face neighbor assignment.
                                            This improves Tecplot¡¯s time to first plot.
                                            See the data section below for format details.
                                            ORDERED and FELINESEG zones must specify 0 for
                                            this value because raw face neighbors are not
                                            defined for these zone types. FEPOLYGON and
                                            FEPOLYHEDRON zones must specify 0 for this value
                                            since face neighbors are defined in the face map
                                            for these zone types.
    */
    int rawLocal = 0;
    fos.write(reinterpret_cast<const char *>(&rawLocal), sizeof(rawLocal));

    /*
            +-----------+
            | INT32		|	Number of miscellaneous user-defined face
            +-----------+	neighbor connections (value >= 0). This value
                                            is in addition to the face neighbors supplied
                                            in the raw section. FEPOLYGON and FEPOLYHEDRON
                                            zones must specify 0.
    */
    int nFaceNeighbor = nFNgh;
    fos.write(reinterpret_cast<const char *>(&nFaceNeighbor), sizeof(nFaceNeighbor));

    if (HurMPI::parRun() && nFNgh != 0) {
        /*
        if ¡°number of miscellaneous user-defined
                face neighbor connections¡± != 0
                +-----------+
                | INT32		|	User defined face neighbor mode
                +-----------+	(0=Local 1-to-1, 1=Local 1-to-many,
                                                2=Global 1-to-1, 3=Global 1-to-many)
                if FE Zone:
                +-----------+
                | INT32		|	Indicates if the finite element face neighbors
                +-----------+	are completely specified by the miscellaneous
                                                face neighbors given: (0=NO, 1=YES). If yes,
                                                then Tecplot will not perform auto assignment
                                                of face neighbors otherwise all faces not
                                                specified are considered boundaries. If no,
                                                then Tecplot will perform auto-assignment of
                                                the face neighbors unless the raw face neighbor
                                                array was supplied. This option is not valid
                                                for ORDERED zones.
        */
        int neighborMode = 2;
        fos.write(reinterpret_cast<const char *>(&neighborMode), sizeof(neighborMode));
        int isGiven = 0;
        fos.write(reinterpret_cast<const char *>(&isGiven), sizeof(isGiven));
    }

    /*
    if FE Zone:
            +-----------+
            | INT32		|	 NumPts
            +-----------+
    */
    int nPs = nPts;
    fos.write(reinterpret_cast<const char *>(&nPs), sizeof(nPs));

    /*
    +-----------+
    | INT32		|	NumElements
    +-----------+
    +-----------+
    | INT32*3	|	ICellDim,JCellDim,
    +-----------+	KCellDim (for future use; set to zero)
    */
    int NumElements = nCls;
    fos.write(reinterpret_cast<const char *>(&NumElements), sizeof(NumElements));
    int ICellDim = 0;
    fos.write(reinterpret_cast<const char *>(&ICellDim), sizeof(ICellDim));
    int JCellDim = 0;
    fos.write(reinterpret_cast<const char *>(&JCellDim), sizeof(JCellDim));
    int KCellDim = 0;
    fos.write(reinterpret_cast<const char *>(&KCellDim), sizeof(KCellDim));

    /*
    For all zone types (repeat for each Auxiliary data name/value pair):
    +-----------+
    | INT32 | 1=Auxiliary name/value pair to follow
    +-----------+ 0=No more Auxiliary name/value pairs.
    */
    int AuxiliaryName = 0;
    fos.write(reinterpret_cast<const char *>(&AuxiliaryName), sizeof(AuxiliaryName));
}

void OpenHurricane::tecplotWriter::writeDataSection(fileOsstream &fos, const geometryMesh &mesh,
                                                const integer nVar) {
    /*DATA SECTION (don¡¯t forget to separate the header from the data
    with an EOHMARKER). The data section contains all of the data
    associated with the zone definitions in the header.
    */

    /*
    +-----------+
    | FLOAT32	| Zone marker Value = 299.0
    +-----------+
    */
    fos.write(reinterpret_cast<const char *>(&ZONEMARKER_), sizeof(ZONEMARKER_));
    /*
    +-----------+
    | INT32*N	|	Variable data format, N=Total number of vars
    +-----------+		1=Float,		2=Double,	3=LongInt,
                                            4=ShortInt,		5=Byte,		6=Bit
    */
    for (integer i = 0; i < nVar; ++i) {
        fos.write(reinterpret_cast<const char *>(&feature<real>::dataFormat), sizeof(int));
    }
    /*
    +-----------+
    | INT32		|	Has passive variables: 0 = no, 1 = yes.
    +-----------+
    */
    int hasPassiveVar = 0;
    fos.write(reinterpret_cast<const char *>(&hasPassiveVar), sizeof(hasPassiveVar));

    /*
    +-----------+
    | INT32		|	Has variable sharing 0 = no, 1 = yes.
    +-----------+
    */
    int hasVarSharing = 0;
    fos.write(reinterpret_cast<const char *>(&hasVarSharing), sizeof(hasVarSharing));

    /*
    +-----------+
    | INT32		|	Zero based zone number to share connectivity
    +-----------+	list with (-1 = no sharing).
    */
    int zeroBasedZN = -1;
    fos.write(reinterpret_cast<const char *>(&zeroBasedZN), sizeof(zeroBasedZN));

    /*
    Compressed list of min/max pairs for each non-shared and non-passive
    variable. For each non-shared and non-passive variable (as specified
    above):
            +-----------+
            | FLOAT64	|	 Min value
            +-----------+
            +-----------+
            | FLOAT64	|	 Max value
            +-----------+
    */
    real minX = veryLarge;
    real maxX = -veryLarge;
    real minY = veryLarge;
    real maxY = -veryLarge;
    real minZ = veryLarge;
    real maxZ = -veryLarge;

    for (integer i = 0; i < mesh.nPoints(); ++i) {
        minX = min(minX, mesh.points()[i].x());
        maxX = max(maxX, mesh.points()[i].x());
        minY = min(minY, mesh.points()[i].y());
        maxY = max(maxY, mesh.points()[i].y());
        minZ = min(minZ, mesh.points()[i].z());
        maxZ = max(maxZ, mesh.points()[i].z());
    }
    fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&minY), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&maxY), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&minZ), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&maxZ), sizeof(double));
}

void OpenHurricane::tecplotWriter::writeDataSectionByMaster(fileOsstream &fos, const geometryMesh &mesh,
                                                        const integer nVar) {
    /*DATA SECTION (don¡¯t forget to separate the header from the data
    with an EOHMARKER). The data section contains all of the data
    associated with the zone definitions in the header.
    */
    if (HurMPI::master()) {
        /*
         +-----------+
         | FLOAT32	| Zone marker Value = 299.0
         +-----------+
         */
        fos.write(reinterpret_cast<const char *>(&ZONEMARKER_), sizeof(ZONEMARKER_));

        /*
        +-----------+
        | INT32*N	|	Variable data format, N=Total number of vars
        +-----------+		1=Float,		2=Double,	3=LongInt,
                                                4=ShortInt,		5=Byte,		6=Bit
        */
        for (integer i = 0; i < nVar; ++i) {
            fos.write(reinterpret_cast<const char *>(&feature<real>::dataFormat), sizeof(int));
        }
        /*
        +-----------+
        | INT32		|	Has passive variables: 0 = no, 1 = yes.
        +-----------+
        */
        int hasPassiveVar = 0;
        fos.write(reinterpret_cast<const char *>(&hasPassiveVar), sizeof(hasPassiveVar));

        /*
        +-----------+
        | INT32		|	Has variable sharing 0 = no, 1 = yes.
        +-----------+
        */
        int hasVarSharing = 0;
        fos.write(reinterpret_cast<const char *>(&hasVarSharing), sizeof(hasVarSharing));

        /*
        +-----------+
        | INT32		|	Zero based zone number to share connectivity
        +-----------+	list with (-1 = no sharing).
        */
        int zeroBasedZN = -1;
        fos.write(reinterpret_cast<const char *>(&zeroBasedZN), sizeof(zeroBasedZN));
    }
    /*
    Compressed list of min/max pairs for each non-shared and non-passive
    variable. For each non-shared and non-passive variable (as specified
    above):
            +-----------+
            | FLOAT64	|	 Min value
            +-----------+
            +-----------+
            | FLOAT64	|	 Max value
            +-----------+
    */
    real minX = veryLarge;
    real maxX = -veryLarge;
    real minY = veryLarge;
    real maxY = -veryLarge;
    real minZ = veryLarge;
    real maxZ = -veryLarge;

    for (integer i = 0; i < mesh.nPoints(); ++i) {
        minX = min(minX, mesh.points()[i].x());
        maxX = max(maxX, mesh.points()[i].x());
        minY = min(minY, mesh.points()[i].y());
        maxY = max(maxY, mesh.points()[i].y());
        minZ = min(minZ, mesh.points()[i].z());
        maxZ = max(maxZ, mesh.points()[i].z());
    }
    HurMPI::reduce(minX, MPI_MIN);
    HurMPI::reduce(maxX, MPI_MAX);
    HurMPI::reduce(minY, MPI_MIN);
    HurMPI::reduce(maxY, MPI_MAX);
    HurMPI::reduce(minZ, MPI_MIN);
    HurMPI::reduce(maxZ, MPI_MAX);
    if (HurMPI::master()) {
        fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&minY), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&maxY), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&minZ), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&maxZ), sizeof(double));
    }
}

void OpenHurricane::tecplotWriter::writeDataSection(fileOsstream &fos, const geometryMesh &mesh,
                                                const integer nVar, const integer fzId) {
    real minX = veryLarge;
    real maxX = -veryLarge;
    real minY = veryLarge;
    real maxY = -veryLarge;
    real minZ = veryLarge;
    real maxZ = -veryLarge;

    const auto &nodel = mesh.globalFaceZoneInfo(fzId).facePoints();

    for (integer i = 0; i < nodel.size(); ++i) {
        minX = min(minX, nodel[i].x());
        maxX = max(maxX, nodel[i].x());
        minY = min(minY, nodel[i].y());
        maxY = max(maxY, nodel[i].y());
        minZ = min(minZ, nodel[i].z());
        maxZ = max(maxZ, nodel[i].z());
    }

    if (HurMPI::parRun()) {
        HurMPI::reduce(minX, MPI_MIN);
        HurMPI::reduce(maxX, MPI_MAX);
        HurMPI::reduce(minY, MPI_MIN);
        HurMPI::reduce(maxY, MPI_MAX);
        HurMPI::reduce(minZ, MPI_MIN);
        HurMPI::reduce(maxZ, MPI_MAX);
    }
    if (!HurMPI::master()) {
        return;
    }
    /*DATA SECTION (don¡¯t forget to separate the header from the data
    with an EOHMARKER). The data section contains all of the data
    associated with the zone definitions in the header.
    */

    /*
    +-----------+
    | FLOAT32	| Zone marker Value = 299.0
    +-----------+
    */
    fos.write(reinterpret_cast<const char *>(&ZONEMARKER_), sizeof(ZONEMARKER_));
    /*
    +-----------+
    | INT32*N	|	Variable data format, N=Total number of vars
    +-----------+		1=Float,		2=Double,	3=LongInt,
                                            4=ShortInt,		5=Byte,		6=Bit
    */
    for (integer i = 0; i < nVar; ++i) {
        fos.write(reinterpret_cast<const char *>(&feature<real>::dataFormat), sizeof(int));
    }
    /*
    +-----------+
    | INT32		|	Has passive variables: 0 = no, 1 = yes.
    +-----------+
    */
    int hasPassiveVar = 0;
    fos.write(reinterpret_cast<const char *>(&hasPassiveVar), sizeof(hasPassiveVar));

    /*
    +-----------+
    | INT32		|	Has variable sharing 0 = no, 1 = yes.
    +-----------+
    */
    int hasVarSharing = 0;
    fos.write(reinterpret_cast<const char *>(&hasVarSharing), sizeof(hasVarSharing));

    /*
    +-----------+
    | INT32		|	Zero based zone number to share connectivity
    +-----------+	list with (-1 = no sharing).
    */
    int zeroBasedZN = -1;
    fos.write(reinterpret_cast<const char *>(&zeroBasedZN), sizeof(zeroBasedZN));

    /*
    Compressed list of min/max pairs for each non-shared and non-passive
    variable. For each non-shared and non-passive variable (as specified
    above):
            +-----------+
            | FLOAT64	|	 Min value
            +-----------+
            +-----------+
            | FLOAT64	|	 Max value
            +-----------+
    */

    fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&minY), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&maxY), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&minZ), sizeof(double));
    fos.write(reinterpret_cast<const char *>(&maxZ), sizeof(double));
}

void OpenHurricane::tecplotWriter::writeCell(fileOsstream &fos, const geometryMesh &mesh) {
    const integer nCells = mesh.nCells();
    const auto &cells_ = mesh.cells();
    int N[8];
    for (integer n = 0; n < nCells; ++n) {
        const integerList &nl = cells_[n].nodesList();
        if (cells_[n].isTetrahedral()) {
            // For tetrahedral cell the cell connectivity is: N1,N2,N3,N3,N4,N4,N4,N4
            /*fos.os() << nl[0] + 1 << " " << nl[1] + 1 << " " << nl[2] + 1 << " " << nl[2] + 1 << " "
                    << nl[3] + 1 << " " << nl[3] + 1 << " " << nl[3] + 1 << " " << nl[3] + 1 << std::endl;*/
            /*N[0] = nl[0] + 1;
            N[1] = nl[1] + 1;
            N[2] = nl[2] + 1;
            N[3] = nl[2] + 1;
            N[4] = nl[3] + 1;
            N[5] = nl[3] + 1;
            N[6] = nl[3] + 1;
            N[7] = nl[3] + 1;*/
            N[0] = nl[0];
            N[1] = nl[1];
            N[2] = nl[2];
            N[3] = nl[2];
            N[4] = nl[3];
            N[5] = nl[3];
            N[6] = nl[3];
            N[7] = nl[3];
        } else if (cells_[n].isPyramid()) {
            //// For pyramid cell: N1,N2,N3,N4,N5,N5,N5,N5
            //fos.os() << nl[0] + 1 << " " << nl[1] + 1 << " " << nl[2] + 1 << " " << nl[3] + 1 << " "
            //	<< nl[4] + 1 << " " << nl[4] + 1 << " " << nl[4] + 1 << " " << nl[4] + 1 << std::endl;
            /*N[0] = nl[0] + 1;
            N[1] = nl[1] + 1;
            N[2] = nl[2] + 1;
            N[3] = nl[3] + 1;
            N[4] = nl[4] + 1;
            N[5] = nl[4] + 1;
            N[6] = nl[4] + 1;
            N[7] = nl[4] + 1;*/
            N[0] = nl[0];
            N[1] = nl[1];
            N[2] = nl[2];
            N[3] = nl[3];
            N[4] = nl[4];
            N[5] = nl[4];
            N[6] = nl[4];
            N[7] = nl[4];
        } else if (cells_[n].isWedge()) {
            // For wedge cell: N1,N2,N3,N3,N4,N5,N6,N6
            /*fos.os() << nl[0] + 1 << " " << nl[1] + 1 << " " << nl[2] + 1 << " " << nl[2] + 1 << " "
                    << nl[3] + 1 << " " << nl[4] + 1 << " " << nl[5] + 1 << " " << nl[5] + 1 << std::endl;*/
            /*N[0] = nl[0] + 1;
            N[1] = nl[1] + 1;
            N[2] = nl[2] + 1;
            N[3] = nl[2] + 1;
            N[4] = nl[3] + 1;
            N[5] = nl[4] + 1;
            N[6] = nl[5] + 1;
            N[7] = nl[5] + 1;*/
            N[0] = nl[0];
            N[1] = nl[1];
            N[2] = nl[2];
            N[3] = nl[2];
            N[4] = nl[3];
            N[5] = nl[4];
            N[6] = nl[5];
            N[7] = nl[5];
        } else if (cells_[n].isHexahedral()) {
            //// For hexahedral cell: N1,N2,N3,N4,N5,N6,N7,N8
            //fos.os() << nl[0] + 1 << " " << nl[1] + 1 << " " << nl[2] + 1 << " " << nl[3] + 1 << " "
            //	<< nl[4] + 1 << " " << nl[5] + 1 << " " << nl[6] + 1 << " " << nl[7] + 1 << std::endl;
            /*N[0] = nl[0] + 1;
            N[1] = nl[1] + 1;
            N[2] = nl[2] + 1;
            N[3] = nl[3] + 1;
            N[4] = nl[4] + 1;
            N[5] = nl[5] + 1;
            N[6] = nl[6] + 1;
            N[7] = nl[7] + 1;*/
            N[0] = nl[0];
            N[1] = nl[1];
            N[2] = nl[2];
            N[3] = nl[3];
            N[4] = nl[4];
            N[5] = nl[5];
            N[6] = nl[6];
            N[7] = nl[7];
        } else {
            LFatal("Not support cell type : %d", integer(cells_[n].shapeType()));
        }
        fos.write(reinterpret_cast<const char *>(N), 8 * sizeof(int));
    }
}

void OpenHurricane::tecplotWriter::writeCellByMaster(fileOsstream &fos, const geometryMesh &mesh) {
    if (!HurMPI::parRun()) {
        writeCell(fos, mesh);
        return;
    }
    const integer nCells = mesh.nCells();
    const auto &cells = mesh.cells();
    integerList nSizeL;
    integerList displs;
    integer allSize = 0;
    nSizeL.resize(HurMPI::getProcSize(), Zero);
    nSizeL[HurMPI::getProcRank()] = 8 * nCells;
    if (HurMPI::master()) {
        displs.resize(HurMPI::getProcSize(), Zero);
    }
    HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());
    if (HurMPI::master()) {
        for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
            displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
        }
        for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
            allSize += nSizeL[ip];
        }
    }
    HurMPI::barrier(HurMPI::getComm());
    integerArray nnl(8 * nCells);
    for (integer n = 0; n < nCells; ++n) {
        const integerList &nl = cells[n].nodesList();
        if (cells[n].isTetrahedral()) {
            // For tetrahedral cell the cell connectivity is: N1,N2,N3,N3,N4,N4,N4,N4
            const auto N1 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[0]);
            const auto N2 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[1]);
            const auto N3 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[2]);
            const auto N4 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[3]);
            nnl[n * 8 + 0] = N1;
            nnl[n * 8 + 1] = N2;
            nnl[n * 8 + 2] = N3;
            nnl[n * 8 + 3] = N3;
            nnl[n * 8 + 4] = N4;
            nnl[n * 8 + 5] = N4;
            nnl[n * 8 + 6] = N4;
            nnl[n * 8 + 7] = N4;
        } else if (cells[n].isPyramid()) {
            //// For pyramid cell: N1,N2,N3,N4,N5,N5,N5,N5
            const auto N1 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[0]);
            const auto N2 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[1]);
            const auto N3 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[2]);
            const auto N4 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[3]);
            const auto N5 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[4]);
            nnl[n * 8 + 0] = N1;
            nnl[n * 8 + 1] = N2;
            nnl[n * 8 + 2] = N3;
            nnl[n * 8 + 3] = N4;
            nnl[n * 8 + 4] = N5;
            nnl[n * 8 + 5] = N5;
            nnl[n * 8 + 6] = N5;
            nnl[n * 8 + 7] = N5;
        } else if (cells[n].isWedge()) {
            // For wedge cell: N1,N2,N3,N3,N4,N5,N6,N6
            const auto N1 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[0]);
            const auto N2 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[1]);
            const auto N3 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[2]);
            const auto N4 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[3]);
            const auto N5 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[4]);
            const auto N6 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[5]);
            nnl[n * 8 + 0] = N1;
            nnl[n * 8 + 1] = N2;
            nnl[n * 8 + 2] = N3;
            nnl[n * 8 + 3] = N3;
            nnl[n * 8 + 4] = N4;
            nnl[n * 8 + 5] = N5;
            nnl[n * 8 + 6] = N6;
            nnl[n * 8 + 7] = N6;
        } else if (cells[n].isHexahedral()) {
            // For hexahedral cell: N1,N2,N3,N4,N5,N6,N7,N8
            const auto N1 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[0]);
            const auto N2 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[1]);
            const auto N3 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[2]);
            const auto N4 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[3]);
            const auto N5 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[4]);
            const auto N6 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[5]);
            const auto N7 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[6]);
            const auto N8 = mesh.globalMeshInfo().globalPointIndeces().toGlobalIndex(nl[7]);
            nnl[n * 8 + 0] = N1;
            nnl[n * 8 + 1] = N2;
            nnl[n * 8 + 2] = N3;
            nnl[n * 8 + 3] = N4;
            nnl[n * 8 + 4] = N5;
            nnl[n * 8 + 5] = N6;
            nnl[n * 8 + 6] = N7;
            nnl[n * 8 + 7] = N8;
        } else {
            LFatal("Not support cell type : %d", integer(cells[n].shapeType()));
        }
    }
    integerArray allNNL;
    if (HurMPI::master()) {
        allNNL.resize(allSize);
    }
    /*HurMPI::gatherv
    (
            nnl.data(),
            8 * nCells,
            feature<integer>::MPIType,
            allNNL.data(),
            nSizeL.data(),
            displs.data(),
            feature<integer>::MPIType,
            HurMPI::masterNo(),
            HurMPI::getComm()
    );*/

    HurMPI::Request request;
    HurMPI::igatherv(nnl.data(), 8 * nCells, feature<integer>::MPIType, allNNL.data(),
                     nSizeL.data(), displs.data(), feature<integer>::MPIType, HurMPI::masterNo(),
                     HurMPI::getComm(), &request);
    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

    if (HurMPI::master()) {
        fos.write(reinterpret_cast<const char *>(&allNNL[0]), allNNL.byteSize());
    }
}

void OpenHurricane::tecplotWriter::writeFaceConnectTecplot(fileOsstream &fos,
                                                       const geometryMesh &mesh) {
    if (HurMPI::parRun()) {
        const auto &cells_ = mesh.cells();
        const auto &faces_ = mesh.faces();
        const auto &globalMeshInfo = mesh.globalMeshInfo();
        // Face connectivity information(for global-one-to-one):
        //       cz, fz, zr, cr
        //       cz - the cell number in the current zone
        //       fz - the number of the cell face in the current zone
        //       zr - the remote zone number
        //       cr - the cell number of the neighboring cell in the remote zone

        const globalFaceIndex &gfi = globalMeshInfo.globalFaceIndeces();
        const typename globalFaceIndex::facePairMapType &gfiMap = gfi.facePairMap();
        for (typename globalFaceIndex::facePairMapType::const_iterator iter = gfiMap.begin();
             iter != gfiMap.end(); ++iter) {
            integer fi = iter->first;
            integer cz = faces_[fi].leftCell();
            integer fz = -1;
            for (integer i = 0; i < cells_[cz].faceSize(); i++) {
                if (fi == cells_[cz].facei(i)) {
                    fz = i;
                    break;
                }
            }

            if (fz == -1) {
                LFatal("Shared face map is not right!");
            }

            fz++; // fz+1

            if (cells_[cz].isHexahedral()) {
                // nothing to be done
            } else if (cells_[cz].isTetrahedral()) {
                if (fz == 1) // f1 -> f5
                {
                    fz = 5;
                } else if (fz == 2) // f2 -> f3
                {
                    fz = 3;
                } else if (fz == 3) // f3 -> f1
                {
                    fz = 1;
                } else if (fz == 4) // f4 -> f2
                {
                    fz = 2;
                }
            } else if (cells_[cz].isPyramid()) {
                // nothing to be done
                // Because the face number order of cell is the same as that of hexahedral cell.
            } else if (cells_[cz].isWedge()) {
                if (fz > 3) {
                    fz++; // f4 -> f5; f5 -> f6
                }
            }

            // Because the start of array in C/C++ is zero.
            integer zr = iter->second[0] + 1;
            integer cr = iter->second[3] + 1;

            // Writting cz, fz, zr, cr.
            //fos.os() << cz + 1 << " " << fz << " " << zr << " " << cr << std::endl;
            int CZ = cz;
            fos.write(reinterpret_cast<const char *>(&CZ), sizeof(CZ));
            int FZ = fz - 1;
            fos.write(reinterpret_cast<const char *>(&FZ), sizeof(FZ));
            int ZR = zr - 1;
            fos.write(reinterpret_cast<const char *>(&ZR), sizeof(ZR));
            int CR = cr - 1;
            fos.write(reinterpret_cast<const char *>(&CR), sizeof(CR));
        }
    }
}

void OpenHurricane::tecplotWriter::writeMesh(fileOsstream &fos, const geometryMesh &mesh) {
    if (HurMPI::master()) {
        string title = mesh.Iteration().configName().name(true);
        title += "-Grid";
        stdStringList varName(3);
        varName[0] = "x";
        varName[1] = "y";
        varName[2] = "z";
        // Tecplot file header
        writeFileHeader(fos, title, 3, varName, tecplotFormat::GRID);
    }
    if (!HurMPI::parRun()) {
        writeZoneHeader(fos, mesh, 0.0, 3);
    } else {
        integerList nPts(HurMPI::getProcSize(), Zero);
        nPts[HurMPI::getProcRank()] = mesh.nPoints();
        HurMPI::gatherList(nPts, HurMPI::masterNo(), HurMPI::getComm());

        integerList nCls(HurMPI::getProcSize(), Zero);
        nCls[HurMPI::getProcRank()] = mesh.nCells();
        HurMPI::gatherList(nCls, HurMPI::masterNo(), HurMPI::getComm());

        integerList nFNgh(HurMPI::getProcSize(), Zero);
        nFNgh[HurMPI::getProcRank()] = mesh.globalMeshInfo().globalFaceIndeces().sharedFaceSize();
        HurMPI::gatherList(nFNgh, HurMPI::masterNo(), HurMPI::getComm());
        /*
        *	+-----------+
                | INT32		|	 ZoneType 0=ORDERED,		1=FELINESEG,
                +-----------+			  2=FETRIANGLE,		3=FEQUADRILATERAL,
                                                                  4=FETETRAHEDRON,	5=FEBRICK,
                                                                  6=FEPOLYGON,		7=FEPOLYHEDRON
        */
        int ZoneType = 5;
        if (HurMPI::master()) {
            for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                string zoneName = "P";
                zoneName += toString(ip);
                // Zone header.
                writeZoneHeader(fos, zoneName, 0.0, 3, ip, nPts[ip], nCls[ip], nFNgh[ip], ZoneType,
                                -1);
            }
        }
    }
    // Wait
    HurMPI::barrier(HurMPI::getComm());

    if (HurMPI::master()) {
        fos.write(reinterpret_cast<const char *>(&EOHMARKER_), sizeof(float));
    }
    HurMPI::barrier(HurMPI::getComm());

    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (HurMPI::isThisProc(ip)) {
            writeDataSection(fos, mesh, 3);

            // Points
            writeData(fos, mesh.points(), mesh.nPoints());

            // Cells
            writeCell(fos, mesh);

            // Face connectivity if necessarily
            writeFaceConnectTecplot(fos, mesh);
            fos.close();
        }
        // Wait
        HurMPI::barrier(HurMPI::getComm());
    }
}

void OpenHurricane::tecplotWriter::faceConnectList(integerList &faceMap, const geometryMesh &mesh,
                                               const integer fzId) {
    const auto &gfz = mesh.globalFaceZoneInfo(fzId);
    const auto &pm = gfz.antiOrderMap();
    const auto &fz = gfz.fZ();
    const auto &f = mesh.faces();

    integer count = 0;
    faceMap.resize(fz.size() * 4);
    for (integer fi = fz.firstIndex(); fi <= fz.lastIndex(); ++fi) {
        if (f[fi].isTriangular()) {
            /*fos.os() << pm.at(f[fi][0]) + 1 << " " << pm.at(f[fi][1]) + 1 << " "
                    << pm.at(f[fi][2]) + 1 << " " << pm.at(f[fi][2]) + 1 << std::endl;*/
            integer i = count * 4 + 0;
#ifdef HUR_DEBUG
            for (integer ij = 0; ij < 3; ++ij) {
                if (pm.find(f[fi][ij]) == pm.end()) {
                    errorAbortStr(("Cannot find point: " + toString(f[fi][ij]) + " in point map"));
                }
            }
#endif // HUR_DEBUG
            faceMap[i++] = pm.at(f[fi][0]);
            faceMap[i++] = pm.at(f[fi][1]);
            faceMap[i++] = pm.at(f[fi][2]);
            faceMap[i] = pm.at(f[fi][2]);
        } else if (f[fi].isQuadrilateral()) {
            /*fos.os() << pm.at(f[fi][0]) + 1 << " " << pm.at(f[fi][1]) + 1 << " "
                    << pm.at(f[fi][2]) + 1 << " " << pm.at(f[fi][3]) + 1 << std::endl;*/
            integer i = count * 4 + 0;
#ifdef HUR_DEBUG
            for (integer ij = 0; ij < 4; ++ij) {
                if (pm.find(f[fi][ij]) == pm.end()) {
                    errorAbortStr(("Cannot find point: " + toString(f[fi][ij]) + " in point map"));
                }
            }
#endif // HUR_DEBUG

            faceMap[i++] = pm.at(f[fi][0]);
            faceMap[i++] = pm.at(f[fi][1]);
            faceMap[i++] = pm.at(f[fi][2]);
            faceMap[i] = pm.at(f[fi][3]);
        } else {
            LFatal("Other type of face not supported yet");
        }
        count++;
    }
    if (!HurMPI::parRun()) {
        return;
    }
    List<integer> nSizeL(HurMPI::getProcSize(), Zero);
    List<integer> displs;
    if (HurMPI::master()) {
        displs.resize(HurMPI::getProcSize(), Zero);
    }
    nSizeL[HurMPI::getProcRank()] = faceMap.size();
    HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());
    integer totalSize = Zero;
    List<integer> rootT;

    if (HurMPI::master()) {
        for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
            displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
        }
        for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
            totalSize += nSizeL[ip];
        }
        rootT.resize(totalSize);
    }
    HurMPI::barrier(HurMPI::getComm());

    /*HurMPI::gatherv
    (
            faceMap.data(),
            faceMap.size(),
            feature<integer>::MPIType,
            rootT.data(),
            nSizeL.data(),
            displs.data(),
            feature<integer>::MPIType,
            HurMPI::masterNo(),
            HurMPI::getComm()
    );*/

    HurMPI::Request request;
    HurMPI::igatherv(faceMap.data(), faceMap.size(), feature<integer>::MPIType, rootT.data(),
                     nSizeL.data(), displs.data(), feature<integer>::MPIType, HurMPI::masterNo(),
                     HurMPI::getComm(), &request);
    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

    if (HurMPI::master()) {
        faceMap.clear();
        faceMap.transfer(rootT);
    }
}

void OpenHurricane::tecplotWriter::writeZoneMesh(fileOsstream &fos, const geometryMesh &mesh,
                                             const integer fzId) {
    if (HurMPI::master()) {
        const auto &name = mesh.globalFaceZoneInfo(fzId).fZ().name();

        string title = mesh.Iteration().configName().name(true);
        title += "-";
        title += name;
        title += "-Grid";
        stdStringList varName(3);
        varName[0] = "x";
        varName[1] = "y";
        varName[2] = "z";
        // Tecplot file header
        writeFileHeader(fos, title, 3, varName, tecplotFormat::GRID);
        const auto nPts = mesh.globalFaceZoneInfo(fzId).totalNodes();
        const auto nEles = mesh.globalFaceZoneInfo(fzId).totalFaces();
        /*
        *	+-----------+
                | INT32		|	 ZoneType 0=ORDERED,		1=FELINESEG,
                +-----------+			  2=FETRIANGLE,		3=FEQUADRILATERAL,
                                                                  4=FETETRAHEDRON,	5=FEBRICK,
                                                                  6=FEPOLYGON,		7=FEPOLYHEDRON
        */
        int ZoneType = 3;
        writeZoneHeader(fos, name, 0.0, 3, HurMPI::masterNo(), nPts, nEles, 0, ZoneType, -1);

        fos.write(reinterpret_cast<const char *>(&EOHMARKER_), sizeof(float));
    }

    writeDataSection(fos, mesh, 3, fzId);
    if (HurMPI::master()) {
        const auto &nodel = mesh.globalFaceZoneInfo(fzId).facePoints();
        // Points
        writeData(fos, nodel, nodel.size());
    }
    integerList faceMap;
    faceConnectList(faceMap, mesh, fzId);

    if (HurMPI::master()) {
        fos.write(reinterpret_cast<const char *>(faceMap.data()), faceMap.byteSize());
    }
}

void OpenHurricane::tecplotWriter::writeEOHMARKER(fileOsstream &fos) {
    fos.write(reinterpret_cast<const char *>(&EOHMARKER_), sizeof(float));
}

void OpenHurricane::tecplotWriter::writeResultToTecplot(fileOsstream &fos, const geometryMesh &mesh) {
    int nVar = mesh.outputTitleSize();
    nVar += 3;
    if (HurMPI::master()) {
        string title = mesh.Iteration().configName().name(true);
        stdStringList varName(nVar);
        varName[0] = "x";
        varName[1] = "y";
        varName[2] = "z";
        auto nlist = mesh.outputTitleNameDocList();
        for (integer i = 0; i < nlist.size(); ++i) {
            integer j = i + 3;
            varName[j] = nlist[i];
        }

        // Tecplot file header
        writeFileHeader(fos, title, nVar, varName, tecplotFormat::FULL);
    }
    real st = 0.0;
    if (mesh.Iteration().hasPhysicalTimeStep()) {
        st = mesh.Iteration().pTStep().totalTime();
    }
    if (!HurMPI::parRun()) {
        int ZoneType = 5;
        string zoneName = "P0";
        // Zone header.
        int strandid = -1;
        if (mesh.Iteration().hasPhysicalTimeStep()) {
            strandid = 0;
        }
        writeZoneHeader(fos, zoneName, st, nVar, 0, mesh.nPoints(), mesh.nCells(), 0, ZoneType,
                        strandid);
    } else {
        integerList nPts(HurMPI::getProcSize(), Zero);
        nPts[HurMPI::getProcRank()] = mesh.nPoints();
        HurMPI::gatherList(nPts, HurMPI::masterNo(), HurMPI::getComm());

        integerList nCls(HurMPI::getProcSize(), Zero);
        nCls[HurMPI::getProcRank()] = mesh.nCells();
        HurMPI::gatherList(nCls, HurMPI::masterNo(), HurMPI::getComm());

        integerList nFNgh(HurMPI::getProcSize(), Zero);
        nFNgh[HurMPI::getProcRank()] = mesh.globalMeshInfo().globalFaceIndeces().sharedFaceSize();
        HurMPI::gatherList(nFNgh, HurMPI::masterNo(), HurMPI::getComm());
        /*
        *	+-----------+
                | INT32		|	 ZoneType 0=ORDERED,		1=FELINESEG,
                +-----------+			  2=FETRIANGLE,		3=FEQUADRILATERAL,
                                                                  4=FETETRAHEDRON,	5=FEBRICK,
                                                                  6=FEPOLYGON,		7=FEPOLYHEDRON
        */
        int ZoneType = 5;
        if (HurMPI::master()) {
            for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                string zoneName = "P";
                zoneName += toString(ip);
                int strandid = -1;
                if (mesh.Iteration().hasPhysicalTimeStep()) {
                    strandid = 0;
                }
                // Zone header.
                writeZoneHeader(fos, zoneName, st, nVar, ip, nPts[ip], nCls[ip], nFNgh[ip],
                                ZoneType, strandid);
            }
        }
    }
    HurMPI::barrier(HurMPI::getComm());
    if (HurMPI::master()) {
        writeEOHMARKER(fos);
    }
    HurMPI::barrier(HurMPI::getComm());
    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (HurMPI::isThisProc(ip)) {
            writeDataSection(fos, mesh, nVar);

            mesh.writeMinMaxOutput(fos);

            // Points
            writeData(fos, mesh.points(), mesh.nPoints());

            // Results
            mesh.writeOutput(fos);

            // Cells
            writeCell(fos, mesh);

            // Face connectivity if necessarily
            writeFaceConnectTecplot(fos, mesh);
            fos.close();
        }
        // Wait
        HurMPI::barrier(HurMPI::getComm());
    }
}

void OpenHurricane::tecplotWriter::writeResultToTecplot(fileOsstream &fos, const geometryMesh &mesh,
                                                    const stringList &outVarName) {
    int nVar = mesh.outputTitleSize(outVarName);
    nVar += 3;
    if (HurMPI::master()) {
        string title = mesh.Iteration().configName().name(true);
        stdStringList varName(nVar);
        varName[0] = "x";
        varName[1] = "y";
        varName[2] = "z";
        auto nlist = mesh.outputTitleNameDocList(outVarName);
        for (integer i = 0; i < nlist.size(); ++i) {
            integer j = i + 3;
            varName[j] = nlist[i];
        }

        // Tecplot file header
        writeFileHeader(fos, title, nVar, varName, tecplotFormat::FULL);
    }
    real st = 0.0;
    if (mesh.Iteration().hasPhysicalTimeStep()) {
        st = mesh.Iteration().pTStep().totalTime();
    }
    if (!HurMPI::parRun()) {
        int ZoneType = 5;
        string zoneName = "P0-";
        zoneName += toString(mesh.Iteration().cStep());
        // Zone header.
        int strandid = -1;
        if (mesh.Iteration().hasPhysicalTimeStep()) {
            strandid = 0;
        }
        writeZoneHeader(fos, zoneName, st, nVar, 0, mesh.nPoints(), mesh.nCells(), 0, ZoneType,
                        strandid);
    } else {
        integerList nPts(HurMPI::getProcSize(), Zero);
        nPts[HurMPI::getProcRank()] = mesh.nPoints();
        HurMPI::gatherList(nPts, HurMPI::masterNo(), HurMPI::getComm());

        integerList nCls(HurMPI::getProcSize(), Zero);
        nCls[HurMPI::getProcRank()] = mesh.nCells();
        HurMPI::gatherList(nCls, HurMPI::masterNo(), HurMPI::getComm());

        integerList nFNgh(HurMPI::getProcSize(), Zero);
        nFNgh[HurMPI::getProcRank()] = mesh.globalMeshInfo().globalFaceIndeces().sharedFaceSize();
        HurMPI::gatherList(nFNgh, HurMPI::masterNo(), HurMPI::getComm());
        /*
        *	+-----------+
                | INT32		|	 ZoneType 0=ORDERED,		1=FELINESEG,
                +-----------+			  2=FETRIANGLE,		3=FEQUADRILATERAL,
                                                                  4=FETETRAHEDRON,	5=FEBRICK,
                                                                  6=FEPOLYGON,		7=FEPOLYHEDRON
        */
        int ZoneType = 5;
        if (HurMPI::master()) {
            for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                string zoneName = "P";
                zoneName += toString(ip);
                zoneName += "-";
                zoneName += toString(mesh.Iteration().cStep());
                int strandid = -1;
                if (mesh.Iteration().hasPhysicalTimeStep()) {
                    strandid = 0;
                }
                // Zone header.
                writeZoneHeader(fos, zoneName, st, nVar, ip, nPts[ip], nCls[ip], nFNgh[ip],
                                ZoneType, strandid);
            }
        }
    }
    HurMPI::barrier(HurMPI::getComm());
    if (HurMPI::master()) {
        writeEOHMARKER(fos);
    }
    HurMPI::barrier(HurMPI::getComm());
    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (HurMPI::isThisProc(ip)) {
            writeDataSection(fos, mesh, nVar);

            mesh.writeMinMaxOutput(fos, outVarName);

            // Points
            writeData(fos, mesh.points(), mesh.nPoints());

            // Results
            mesh.writeOutput(fos, outVarName);

            // Cells
            writeCell(fos, mesh);

            // Face connectivity if necessarily
            writeFaceConnectTecplot(fos, mesh);
            fos.close();
        }
        // Wait
        HurMPI::barrier(HurMPI::getComm());
    }
}

void OpenHurricane::tecplotWriter::writeResultToTecplotByMaster(fileOsstream &fos,
                                                            const geometryMesh &mesh) {
    if (!HurMPI::parRun()) {
        writeResultToTecplot(fos, mesh);
        return;
    }

    int nVar = mesh.outputTitleSize();
    nVar += 3;
    if (HurMPI::master()) {
        string title = mesh.Iteration().configName().name(true);
        stdStringList varName(nVar);
        varName[0] = "x";
        varName[1] = "y";
        varName[2] = "z";
        auto nlist = mesh.outputTitleNameDocList();
        for (integer i = 0; i < nlist.size(); ++i) {
            integer j = i + 3;
            varName[j] = nlist[i];
        }

        // Tecplot file header
        writeFileHeader(fos, title, nVar, varName, tecplotFormat::FULL);
    }
    if (HurMPI::master()) {
        real st = 0.0;
        int strandid = -1;
        if (mesh.Iteration().hasPhysicalTimeStep()) {
            st = mesh.Iteration().pTStep().totalTime();
            strandid = 0;
        }
        int ZoneType = 5;
        string zoneName = "P0-";
        zoneName += toString(mesh.Iteration().cStep());
        // Zone header.
        tecplotWriter::writeZoneHeader(fos, zoneName, st, nVar, 0, mesh.allMeshNodeNumber(),
                                       mesh.allMeshCellNumber(), 0, ZoneType, strandid);

        tecplotWriter::writeEOHMARKER(fos);
    }
    writeDataSectionByMaster(fos, mesh, nVar);
    mesh.writeMinMaxOutputByMaster(fos);
    writeMeshPointsToTecplotByMaster(fos, mesh);

    mesh.writeOutputByMaster(fos);
    writeCellByMaster(fos, mesh);
}

void OpenHurricane::tecplotWriter::writeResultToTecplotByMaster(fileOsstream &fos,
                                                            const geometryMesh &mesh,
                                                            const stringList &outVarName) {
    if (!HurMPI::parRun()) {
        writeResultToTecplot(fos, mesh, outVarName);
        return;
    }

    int nVar = mesh.outputTitleSize(outVarName);
    nVar += 3;
    if (HurMPI::master()) {
        string title = mesh.Iteration().configName().name(true);
        stdStringList varName(nVar);
        varName[0] = "x";
        varName[1] = "y";
        varName[2] = "z";
        auto nlist = mesh.outputTitleNameDocList(outVarName);
        for (integer i = 0; i < nlist.size(); ++i) {
            integer j = i + 3;
            varName[j] = nlist[i];
        }

        // Tecplot file header
        writeFileHeader(fos, title, nVar, varName, tecplotFormat::FULL);
    }
    if (HurMPI::master()) {
        real st = 0.0;
        int strandid = -1;
        if (mesh.Iteration().hasPhysicalTimeStep()) {
            st = mesh.Iteration().pTStep().totalTime();
            strandid = 0;
        }
        int ZoneType = 5;
        string zoneName = "P0-";
        zoneName += toString(mesh.Iteration().cStep());
        // Zone header.
        tecplotWriter::writeZoneHeader(fos, zoneName, st, nVar, 0, mesh.allMeshNodeNumber(),
                                       mesh.allMeshCellNumber(), 0, ZoneType, strandid);

        tecplotWriter::writeEOHMARKER(fos);
    }
    writeDataSectionByMaster(fos, mesh, nVar);
    mesh.writeMinMaxOutputByMaster(fos, outVarName);
    writeMeshPointsToTecplotByMaster(fos, mesh);

    mesh.writeOutputByMaster(fos, outVarName);

    writeCellByMaster(fos, mesh);
}

void OpenHurricane::tecplotWriter::writeMeshPointsToTecplotByMaster(fileOsstream &fos,
                                                                const geometryMesh &mesh) {
    List<pointField> tpp(mesh.pointZones().size());
    for (integer pzid = 0; pzid < mesh.pointZones().size(); ++pzid) {
        mesh.globalMeshInfo().globalPointIndeces().getNodeList(pzid, tpp[pzid]);
    }
    if (HurMPI::master()) {
        for (int j = 0; j < feature<point>::nElements_; ++j) {
            for (integer i = 0; i < tpp.size(); ++i) {
                Array<real> components = tpp[i].component(j);

                fos.write(reinterpret_cast<const char *>(&components[0]), components.byteSize());
            }
        }
    }
}
