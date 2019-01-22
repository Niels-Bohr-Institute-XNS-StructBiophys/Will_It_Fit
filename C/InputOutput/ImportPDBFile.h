/*
 1- 6   - Chars Record name "ATOM  " or "HETATM ".
 7-11   - Integer Atom serial number.
12      - ? /empty
13-16   - Chars Atom name.
17      - Char altLoc Alternate location indicator.
18-20   - Chars Residue name resName Residue name.
21      - ? /empty
22      - Char chainID Chain identifier.
23-26   - Integer resSeq Residue sequence number.
27      - Char iCode Code for insertion of residues.
28-30   - ? /empty
31-38   - Real(8.3) x Orthogonal coordinates for X in Angstroms.
39-46   - Real(8.3) y Orthogonal coordinates for Y in Angstroms.
47-54   - Real(8.3) z Orthogonal coordinates for Z in Angstroms.
55-60   - Real(6.2) Occupancy.
61-66   - Real(6.2) Temperature factor.
67-76   - ? /empty
77-78   - Chars Element symbol, right-justified.
79-80   - Chars Charge on the atom.

.....|....|?...||..|?|...||???.......|.......|.......|.....|.....|??????????.|.|
         11111111112222222222333333333344444444445555555555666666666677777777778
12345678901234567890123456789012345678901234567890123456789012345678901234567890
.....|....|?...||..|?|...||???.......|.......|.......|.....|.....|??????????.|.|
ATOM      1  N   GLY A   1      -0.959  20.540 -13.302  1.00 34.83           N  
ATOM      2  CA  GLY A   1      -0.916  20.008 -11.944  1.00 32.41           C  

HETATM 4894  N   UZ9 L1029      15.220 -10.635 -48.148  1.00 25.35           N  
HETATM 4895  CA  UZ9 L1029      15.613 -11.666 -49.109  1.00 27.85           C  

*/



/* function delallspc deletes all spaces in the input string str */
/* http://www.cs.bu.edu/teaching/cpp/string/array-vs-ptr/ */
void delallspc(char *str)
{
	char dummy[strlen(str)+1] ;

	for (; *str; ++str)
	{
		if (*str == ' ')
		{
			strcpy( dummy, str+1) ;
			strcpy( str, dummy) ;
			--str ;
		}
	}
}

// writes all chars between the last '/' and '.' from str to new_str
// OS dependent directory separator ?!
// str = file => new_str = file
// str = file.pdb => new_str = file
// str = /usr/lib/file.pdb => new_str = file
// str = /usr/lib/file => new_str = file
// str = /usr/lib.sdd/file.pdb => new_str = file
// str = /usr/lib.sdd/file => new_str => error
// str = /usr/lib/file.pdb/ss => new_str => error
void GetBaseFileNameFromPath(char *str, char *new_str, int WriteLog, FILE* logfile)
{
	int len ;
	char *dir_delim_pos = strrchr(str, '/') ; // last occurence of "/", can be NULL if not found
	char *file_delim_pos = strrchr(str, '.') ; // last occurence of ".", can be NULL if not found

	if ( dir_delim_pos == 0 ) { dir_delim_pos = str ; /* printf("dir_delim_pos is NULL\n"); fflush(NULL) ; */ }
	else { ++dir_delim_pos ; }
	if ( file_delim_pos == 0 ) { file_delim_pos = str + strlen(str) - 1 ; /* printf("file_delim_pos is NULL\n"); fflush(NULL) ; */ }
	else { --file_delim_pos ; }
	// printf("*dir_delim_pos = %c, *file_delim_pos = %c\n", *dir_delim_pos, *file_delim_pos) ;

	len = file_delim_pos - dir_delim_pos ;
	// printf("len = %d\n", len) ; fflush(NULL) ;

	if ( len < 1 )
	{
		if ( abs(WriteLog) > 0 )
		{	fprintf( logfile, "An error occured when finding the basefilename in path %s. Exit.\n", str) ;
			fflush(logfile) ;
		}
		exit(0) ;
	}

	for ( int i=0; i<=len; ++i ) { new_str[i] = *(dir_delim_pos+i) ; }
	new_str[len+1] = 0 ;

	// printf("%s\n", new_str) ; fflush(NULL) ;
	// printf("%s\n", str) ; fflush(NULL) ;
}



int CheckNumberOfResiduesInPDBFile(char Filename[256])
{
	// Declarations
	FILE *PointerToPDBFile ;
	char Linebuffer[82] ; // length of PDB file incl newline (\n) and terminating NULL (\0)
	char dummy[82] ; // dummy char array
	int NumberOfResidues = 0 ;
	int PreviousResidueID = 0 ;
	int ResidueID = 0 ;

	// I/O
	PointerToPDBFile = fopen( Filename, "r") ;

	if( PointerToPDBFile == 0) { return -1 ; }

	while( fgets( Linebuffer, sizeof(Linebuffer), PointerToPDBFile) != NULL )
	{
		// extract lines starting with ATOM or HETATM pattern
		if ( strncmp( Linebuffer, "ATOM", 4) == 0 || strncmp( Linebuffer, "HETATM", 6) == 0 )
		{
			// read ResidueID from line
			memcpy( dummy, Linebuffer + 22, 5 * sizeof Linebuffer[0]) ; // read ResidueID (pos 23-26)
			dummy[5] = 0 ;
			delallspc( dummy ) ;
			ResidueID = strtol( dummy, NULL, 10) ;

			// count ResidueIDs
			if( ResidueID != PreviousResidueID && ResidueID != 0)
			{
				memcpy( dummy, Linebuffer + 16, 1 * sizeof Linebuffer[0]) ; // read altLoc from pos 17
				dummy[1] = 0 ;
				// in case of multiple occupancies (alternative locations) count only once
				if ( strcmp( " ", dummy) == 0 || strcmp( "A", dummy) == 0) { ++NumberOfResidues; }
			}

			PreviousResidueID = ResidueID;
		}
	}

	fclose(PointerToPDBFile) ;

	return NumberOfResidues;
}



int CheckNumberOfAtomsInPDBFile(char Filename[256])
{
	// Declarations
	FILE *PointerToPDBFile ;
	char Linebuffer[82] ; // length of PDB file incl newline (\n) and terminating NULL (\0)
	char dummy[82] ; // dummy char array
	int NumberOfAtoms = 0 ;

	// I/O
	PointerToPDBFile = fopen( Filename, "r") ;

	if( PointerToPDBFile == 0) { return -1 ; }

	while( fgets( Linebuffer, sizeof(Linebuffer), PointerToPDBFile) != NULL )
	{
		// count lines starting with ATOM or HETATM pattern
		if ( strncmp( Linebuffer, "ATOM", 4) == 0 || strncmp( Linebuffer, "HETATM", 6) == 0 )
		{
			memcpy( dummy, Linebuffer + 16, 1 * sizeof Linebuffer[0]) ; // read altLoc from pos 17
			dummy[1] = 0 ;
			// in case of multiple occupancies (alternative locations) count atoms only once
			if ( strcmp( " ", dummy) == 0 || strcmp( "A", dummy) == 0) { ++NumberOfAtoms; }
		}
	}

	fclose(PointerToPDBFile) ;

	return NumberOfAtoms;
}


void AssignAtom( struct Protein * ProteinStructure, int IndexAtomID, int WriteLog, FILE* logfile)
{
	// function for assigning the scattering length to the correct atoms
	int AtomRecg = 0;

	double r0 = 2.818e-13 ; // classical electron radius [cm]

	// Displaced volumes [Å^3]
	const double HVolume = 5.15;
	const double DVolume = 5.15;
	const double CVolume = 16.44;
	const double NVolume = 2.49;
	const double OVolume = 9.13;
	const double PVolume = 5.73;
	const double SVolume = 19.86;

	const double CLVolume = 28.81;
	const double ZNVolume = 9.85;
	const double CAVolume = 31.89;
	const double FEVolume = 7.99;
	const double NAVolume = 49.0 ; // r = 2.27 Å

	const double QVolume = 4.133 * 30.0;


	// Isotope/Compound weights [u]
	const double HWeight = 1.008;
	const double DWeight = 2.014;
	const double CWeight = 12.011;
	const double NWeight = 14.007;
	const double OWeight = 15.999;
	const double PWeight = 30.974;
	const double SWeight = 32.06 ;

	const double CLWeight = 35.45;
	const double ZNWeight = 65.38;
	const double CAWeight = 40.08;
	const double FEWeight = 55.85;
	const double NAWeight = 126.9;

	const double QWeight = 4.133 * ( 2.0 * HWeight + OWeight ) ;


	// X-ray scattering lengths [cm]
	const double HXRayScatteringLength =  1.0 * r0;
	const double DXRayScatteringLength =  1.0 * r0;
	const double CXRayScatteringLength =  6.0 * r0;
	const double NXRayScatteringLength =  7.0 * r0;
	const double OXRayScatteringLength =  8.0 * r0;
	const double PXRayScatteringLength = 15.0 * r0;
	const double SXRayScatteringLength = 16.0 * r0;

	const double CLXRayScatteringLength = 17.0 * r0;
	const double ZNXRayScatteringLength = 30.0 * r0;
	const double FEXRayScatteringLength = 26.0 * r0;
	const double CAXRayScatteringLength = 20.0 * r0;
	const double NAXRayScatteringLength = 11.0* r0;

	const double QXRayScatteringLength =  4.133 * 10.0 * r0; // 4.133 * 2.818e-12


	// neutron scattering lengths [cm]
	const double HNeutronScatteringLength = -3.742e-13;
	const double DNeutronScatteringLength = 6.674e-13;
	const double CNeutronScatteringLength = 6.6484e-13;
	const double NNeutronScatteringLength = 9.36e-13;
	const double ONeutronScatteringLength = 5.803e-13;
	const double PNeutronScatteringLength = 5.13E-13;
	const double SNeutronScatteringLength = 2.847e-13;

	const double CLNeutronScatteringLength = 9.577e-13;
	const double ZNNeutronScatteringLength = 5.680e-13;
	const double FENeutronScatteringLength = 9.450e-13;
	const double CANeutronScatteringLength = 4.700e-13;
	const double NANeutronScatteringLength = 3.630e-13;

	const double QNeutronScatteringLength = QXRayScatteringLength ; // Note that this is a dummy value. For neutron contrasts the scattering lengths of the dummy waters are set in Model.h


	// the following if statements are to read the atoms. Note this was done in a switch statement before but that does not work for multi character atom names such as "Zn"

	// single char atom names
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] ==  ' ' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'H' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = HXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = HNeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = HVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = HWeight ;
		ProteinStructure->NumberOfHAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] ==  ' ' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'D' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = DXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = DNeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = DVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = DWeight ;
		ProteinStructure->NumberOfDAtoms += 1 ;
		AtomRecg = 1;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == ' ' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'C' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = CXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = CNeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = CVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = CWeight ;
		ProteinStructure->NumberOfCAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == ' ' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'N' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = NXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = NNeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = NVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = NWeight ;
		ProteinStructure->NumberOfNAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == ' ' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'O' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = OXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = ONeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = OVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = OWeight ;
		ProteinStructure->NumberOfOAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == ' ' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'P' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = PXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = PNeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = PVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = PWeight ;
		ProteinStructure->NumberOfPAtoms += 1 ;
		AtomRecg = 1;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == ' ' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'S' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = SXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = SNeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = SVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = SWeight ;
		ProteinStructure->NumberOfSAtoms += 1 ;
		AtomRecg = 1 ;
	}
	/*
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == ' ' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'I' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = IXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = INeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = IVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = IWeight ;
		ProteinStructure->NumberOfIAtoms += 1 ;
		AtomRecg = 1 ;
	}
	*/

	// double char atom names
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == 'Z' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'N' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = ZNXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = ZNNeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = ZNVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = ZNWeight ;
		ProteinStructure->NumberOfZNAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == 'C' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'L' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = CLXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = CLNeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = CLVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = CLWeight ;
		ProteinStructure->NumberOfCLAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == 'N' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'A' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = NAXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = NANeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = NAVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = NAWeight ;
		ProteinStructure->NumberOfNAAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == 'C' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'A' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = CAXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = CANeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = CAVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = CAWeight ;
		ProteinStructure->NumberOfCAAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == 'F' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'E' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = FEXRayScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = FENeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = FEVolume ;
		ProteinStructure->Atoms[IndexAtomID].Weight = FEWeight ;
		ProteinStructure->NumberOfFEAtoms += 1 ;
		AtomRecg = 1 ;
	}

	// special "atoms" like dummy waters (residues WAT with atom names Q)
	// note that crystal waters (residues HOH) consist of and are treated like normal atoms, i.e. 2xH and 1xO atoms
	if ( ProteinStructure->Atoms[IndexAtomID].Type[0] == ' ' && ProteinStructure->Atoms[IndexAtomID].Type[1] == 'Q' )
	{
		ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength = QXRayScatteringLength ; // 1.165E-11
		ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength = QNeutronScatteringLength ;
		ProteinStructure->Atoms[IndexAtomID].Volume = QVolume ; // 123.990
		ProteinStructure->Atoms[IndexAtomID].Weight = QWeight ; // Note that this is the value for dummy-H2O, for neutron contrasts H2O/D2O mixtures the correct weights are not considered in Model.h
		ProteinStructure->NumberOfQAtoms += 1 ;
		AtomRecg = 1 ;
	}


	// if not found error and exit
	if ( AtomRecg == 0 )
	{
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "AssignAtom(): Atom name %s not found in database. Exit.\n", ProteinStructure->Atoms[IndexAtomID].Type) ; }
		fflush(logfile) ;
		exit(0) ;
	}

	return ;
}




void CenterPDBAtom( struct Protein * ProteinStructure, int WriteLog, FILE* logfile)
{
	// Function for finding the geometric center of PDB Structure and translating it so that the center becomes (0,0,0)
	// Find geometric center for protein structure
	// Shift all coordinates (independent of weighting scheme) by the coordinates of the geometric center (none of their centers will match the geometric center)

	double xDum = 0.0 ;
	double yDum = 0.0 ;
	double zDum = 0.0 ;
	for ( int j = 0; j < ProteinStructure->NumberOfAtoms; j++)
	{
		xDum = xDum + ProteinStructure->Atoms[j].x ;
		yDum = yDum + ProteinStructure->Atoms[j].y ;
		zDum = zDum + ProteinStructure->Atoms[j].z ;
	}
	xDum = xDum / (double)ProteinStructure->NumberOfAtoms ;
	yDum = yDum / (double)ProteinStructure->NumberOfAtoms ;
	zDum = zDum / (double)ProteinStructure->NumberOfAtoms ;

	if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\t\tFound center at (x,y,z) = (%f , %f , %f)\n", xDum, yDum, zDum) ; }

	for ( int j = 0; j < ProteinStructure->NumberOfAtoms; j++)
	{
		ProteinStructure->Atoms[j].x = ProteinStructure->Atoms[j].x - xDum ;
		ProteinStructure->Atoms[j].y = ProteinStructure->Atoms[j].y - yDum ;
		ProteinStructure->Atoms[j].z = ProteinStructure->Atoms[j].z - zDum ;
	}

	// Reset dum and recalculate new center
	xDum = 0.0 ;
	yDum = 0.0 ;
	zDum = 0.0 ;
	for ( int j = 0; j < ProteinStructure->NumberOfAtoms; j++)
	{
		xDum = xDum + ProteinStructure->Atoms[j].x ;
		yDum = yDum + ProteinStructure->Atoms[j].y ;
		zDum = zDum + ProteinStructure->Atoms[j].z ;
	}

	if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\t\tNew center is at (x,y,z) = (%f , %f , %f)\n\n", xDum, yDum, zDum) ; }
}


void CenterPDBResidue( struct Protein * ProteinStructure, int WriteLog, FILE* logfile)
{
	// Function for finding the geometric center of PDB Structure and translating it so that the center becomes (0,0,0)
	// Find geometric center for protein structure
	// Shift all coordinates (independent of weighting scheme) by the coordinates of the geometric center (none of their centers will match the geometric center)

	double xDum = 0.0 ;
	double yDum = 0.0 ;
	double zDum = 0.0 ;
	for ( int j = 0; j < ProteinStructure->NumberOfResidues; j++)
	{
		xDum = xDum + ProteinStructure->Residues[j].xVolume ;
		yDum = yDum + ProteinStructure->Residues[j].yVolume ;
		zDum = zDum + ProteinStructure->Residues[j].zVolume ;
	}
	xDum = xDum / (double)ProteinStructure->NumberOfResidues ;
	yDum = yDum / (double)ProteinStructure->NumberOfResidues ;
	zDum = zDum / (double)ProteinStructure->NumberOfResidues ;

	if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\t\tFound center at (x,y,z) = (%f , %f , %f)\n", xDum, yDum, zDum) ; }

	for ( int j = 0; j < ProteinStructure->NumberOfResidues; j++)
	{
		ProteinStructure->Residues[j].xVolume = ProteinStructure->Residues[j].xVolume - xDum ;
		ProteinStructure->Residues[j].yVolume = ProteinStructure->Residues[j].yVolume - yDum ;
		ProteinStructure->Residues[j].zVolume = ProteinStructure->Residues[j].zVolume - zDum ;

		ProteinStructure->Residues[j].xXRayScattering = ProteinStructure->Residues[j].xXRayScattering - xDum ;
		ProteinStructure->Residues[j].yXRayScattering = ProteinStructure->Residues[j].yXRayScattering - yDum ;
		ProteinStructure->Residues[j].zXRayScattering = ProteinStructure->Residues[j].zXRayScattering - zDum ;

		ProteinStructure->Residues[j].xNeutronScattering = ProteinStructure->Residues[j].xNeutronScattering - xDum ;
		ProteinStructure->Residues[j].yNeutronScattering = ProteinStructure->Residues[j].yNeutronScattering - yDum ;
		ProteinStructure->Residues[j].zNeutronScattering = ProteinStructure->Residues[j].zNeutronScattering - zDum ;
	}

	// Reset dum and recalculate new center
	xDum = 0.0 ;
	yDum = 0.0 ;
	zDum = 0.0 ;
	for ( int j = 0; j < ProteinStructure->NumberOfResidues; j++)
	{
		xDum = xDum + ProteinStructure->Residues[j].xVolume ;
		yDum = yDum + ProteinStructure->Residues[j].yVolume ;
		zDum = zDum + ProteinStructure->Residues[j].zVolume ;
	}

	if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\t\tNew center is at (x,y,z) = (%f , %f , %f)\n\n", xDum, yDum, zDum) ; }
}





// Fill fields from struct ProteinStructure->Residues[] and ProteinStructure->Atoms[]
// ProteinStructure has been allocated priorly in function AllocateProteinStructure(...) ../Auxillary/Allocation.h
void ImportAtomsAndResiduesFromPDBFile( char Filename[256], struct Protein * ProteinStructure, char *ResultsDirectory, int WriteLog, FILE* logfile)
{
	// Declarations
	FILE *PointerToPDBFile ;
	FILE *PointerToSkippedLinesPDBFile ;

	char BaseFilename[256] ;
	char SkippedLinesPDBFileLocation[256] ;

	char Linebuffer[82] ; // length of PDB file incl newline (\n) and terminating NULL (\0)
	char dummy[82] ; // dummy char array

	int PreviousResidueID = 0 ;
	int ResidueID = 0 ;
	int IndexResidueID = 0 ;
	char PreviousResidueName[4] ; // includes a terminating NULL (\0)
	char CurrentResidueName[4] ; // includes a terminating NULL (\0)

	int NumberOfModificationAtoms = 0 ;

	int AtomID = 0 ;
	int IndexAtomID = 0 ;
//	char AtomName[5] ; // includes a terminating NULL (\0)
//	char AtomType[3] ; // includes a terminating NULL (\0)


	double xDummy ;
	double yDummy ;
	double zDummy ;

	double VolumeOfResidue = 0.0 ;
	double WeightOfResidue = 0.0 ;
	double XRayScatteringLengthOfResidue = 0.0 ;
	double NeutronScatteringLengthOfResidue = 0.0 ;

	double xCenterOfXRayScattering = 0.0 ;
	double yCenterOfXRayScattering = 0.0 ;
	double zCenterOfXRayScattering = 0.0 ;

	double xCenterOfNeutronScattering = 0.0 ;
	double yCenterOfNeutronScattering = 0.0 ;
	double zCenterOfNeutronScattering = 0.0 ;

	double xCenterOfVolume = 0.0 ;
	double yCenterOfVolume = 0.0 ;
	double zCenterOfVolume = 0.0 ;


	// I/O

	// skip file-existence check for PDB file
	PointerToPDBFile = fopen( Filename, "r") ;

	// get base filename of PDB file and open PDB file to write skipped lines
	GetBaseFileNameFromPath( Filename, BaseFilename, WriteLog, logfile) ;
	sprintf( SkippedLinesPDBFileLocation, "%s%s%s", ResultsDirectory, BaseFilename, "_LinesNotRead.pdb") ; // sprintf appends automatically 0

	if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\t\tWrite skipped lines from PDB file to %s\n", SkippedLinesPDBFileLocation) ; }

	PointerToSkippedLinesPDBFile = fopen( SkippedLinesPDBFileLocation, "w") ;
	if ( PointerToSkippedLinesPDBFile == NULL )
	{
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "Could not write to file %s. Exit.\n", SkippedLinesPDBFileLocation) ; fflush(logfile) ; }
		exit(0) ;
	}

	// Reading of Amino acids:
	// We read each line that starts with either "ATOM" or "HETATM"
	// For the "ATOM" lines residues are assigned based on ResidueID.
	// Each residue has a weight, a volume and a scattering length density (in X-rays and neutrons)
	// As the lines are read these properties are updated according the atom on the current line until the ResidueID changes
	while( fgets( Linebuffer, sizeof(Linebuffer), PointerToPDBFile) != NULL )
	{
		// extract lines starting with ATOM or HETATM pattern
		if ( strncmp( Linebuffer, "ATOM", 4) == 0 || strncmp( Linebuffer, "HETATM", 6) == 0 )
		{
			// current line processing to read current AtomID, AtomName, altLoc, ResidueName, ResidueID, coordinates, AtomType


			// read AtomID from line
			memcpy( dummy, Linebuffer + 6, 5 * sizeof Linebuffer[0]) ; // read AtomID (pos 7-11) -> 5
			dummy[5] = 0 ;
			delallspc( dummy ) ;
			AtomID = strtol( dummy, NULL, 10) ;


			// read AtomName (here called (Current)AtomName)
//			memcpy( AtomName, Linebuffer + 12, 4 * sizeof Linebuffer[0]) ; // read AtomName (pos 13-16) -> 4
//			AtomName[4] = 0 ;
			memcpy( ProteinStructure->Atoms[IndexAtomID].Name, Linebuffer + 12, 4 * sizeof Linebuffer[0]) ; // read AtomName (pos 13-16) -> 4
			ProteinStructure->Atoms[IndexAtomID].Name[4] = 0 ;


			// read altLoc from line
			memcpy( dummy, Linebuffer + 16, 1 * sizeof Linebuffer[0]) ; // read altLoc (pos 17) -> 1
			dummy[1] = 0 ;
			// in case of multiple occupancies (alternative locations) print line (Linebuffer incl newline char) and continue
			// NOTE THAT IN NICHOLAS VERSION MULTIPLE OCCUPANCIES WERE NOT COUNTED IN CheckNumberOfResiduesInPDBFile BUT THEN READ AND TREATED HERE -> ERROR
			if ( strcmp( " ", dummy) != 0 && strcmp( "A", dummy) != 0)
			{
				fprintf( PointerToSkippedLinesPDBFile, "%s", Linebuffer) ;
				continue ;
			}


			// read ResidueName (here called (Current)ResidueName)
//			memcpy( CurrentResidueName, Linebuffer + 17, 3 * sizeof Linebuffer[0]) ; // read ResidueName (pos 18-20) -> 3
//			CurrentResidueName[3] = 0 ;
			memcpy( ProteinStructure->Atoms[IndexAtomID].ResidueName, Linebuffer + 17, 3 * sizeof Linebuffer[0]) ; // read ResidueName (pos 18-20) -> 3
			ProteinStructure->Atoms[IndexAtomID].ResidueName[3] = 0 ;

			// read ResidueID from line
			memcpy( dummy, Linebuffer + 22, 4 * sizeof Linebuffer[0]) ; // read ResidueID (pos 23-26) -> 4
			dummy[4] = 0 ;
			delallspc( dummy ) ;
			ResidueID = strtol( dummy, NULL, 10) ;


			// read x, y, z coordinates and element symbol (here called (Current)AtomName)
			memcpy( dummy, Linebuffer + 30, 8 * sizeof Linebuffer[0]) ; // read x coo (pos 31-38) -> 8
			dummy[8] = 0 ;
			delallspc( dummy ) ;
			xDummy = strtod( dummy, NULL) ;

			memcpy( dummy, Linebuffer + 38, 8 * sizeof Linebuffer[0]) ; // read y coo (pos 39-46) -> 8
			dummy[8] = 0 ;
			delallspc( dummy ) ;
			yDummy = strtod( dummy, NULL) ;

			memcpy( dummy, Linebuffer + 46, 8 * sizeof Linebuffer[0]) ; // read y coo (pos 47-54) -> 8
			dummy[8] = 0 ;
			delallspc( dummy ) ;
			zDummy = strtod( dummy, NULL) ;


			// read element symbol (here called (Current)AtomType)
//			memcpy( AtomType, Linebuffer + 76, 2 * sizeof Linebuffer[0]) ; // read element symbol, keep spaces (pos 77-78) -> 2
//			AtomType[2] = 0 ;
			memcpy( ProteinStructure->Atoms[IndexAtomID].Type, Linebuffer + 76, 2 * sizeof Linebuffer[0]) ; // read element symbol, keep spaces (pos 77-78) -> 2
			ProteinStructure->Atoms[IndexAtomID].Type[2] = 0 ;

			// assign already available information to Atoms
			//strcpy( ProteinStructure->Atoms[IndexAtomID].Name, AtomName) ;
			//strcpy( ProteinStructure->Atoms[IndexAtomID].Type, AtomType) ;
			ProteinStructure->Atoms[IndexAtomID].ID = AtomID ;

			//strcpy( ProteinStructure->Atoms[IndexAtomID].ResidueName, CurrentResidueName) ;
			ProteinStructure->Atoms[IndexAtomID].ResidueID = ResidueID ;

			ProteinStructure->Atoms[IndexAtomID].x = xDummy ;
			ProteinStructure->Atoms[IndexAtomID].y = yDummy ;
			ProteinStructure->Atoms[IndexAtomID].z = zDummy ;


			strcpy( CurrentResidueName, ProteinStructure->Atoms[IndexAtomID].ResidueName) ;



			// if the ResidueID is different from the previous one, do final assignment of residue stuff (ProteinStructure->Residues[IndexResidueID].*), before continuing with the next residue
			// exception is if first Residue ID is read, which differs PreviousResidueID (0)
			// after that "local" residue variables are reset to 0.0 and IndexResidueID is incremented for the new residue, and PreviousResidueID is reset
			if ( ResidueID != PreviousResidueID )
			{
				if (PreviousResidueID != 0)
				{
					// assign
					ProteinStructure->Residues[IndexResidueID].Volume = VolumeOfResidue ;
					ProteinStructure->Residues[IndexResidueID].Weight = WeightOfResidue ;

					ProteinStructure->Residues[IndexResidueID].XRayScatteringLength = XRayScatteringLengthOfResidue ;
					ProteinStructure->Residues[IndexResidueID].NeutronScatteringLength = NeutronScatteringLengthOfResidue ;

					ProteinStructure->Residues[IndexResidueID].xVolume = xCenterOfVolume / VolumeOfResidue ;
					ProteinStructure->Residues[IndexResidueID].yVolume = yCenterOfVolume / VolumeOfResidue ;
					ProteinStructure->Residues[IndexResidueID].zVolume = zCenterOfVolume / VolumeOfResidue ;

					ProteinStructure->Residues[IndexResidueID].xXRayScattering = xCenterOfXRayScattering / XRayScatteringLengthOfResidue ;
					ProteinStructure->Residues[IndexResidueID].yXRayScattering = yCenterOfXRayScattering / XRayScatteringLengthOfResidue ;
					ProteinStructure->Residues[IndexResidueID].zXRayScattering = zCenterOfXRayScattering / XRayScatteringLengthOfResidue ;

					ProteinStructure->Residues[IndexResidueID].xNeutronScattering = xCenterOfNeutronScattering / NeutronScatteringLengthOfResidue ;
					ProteinStructure->Residues[IndexResidueID].yNeutronScattering = yCenterOfNeutronScattering / NeutronScatteringLengthOfResidue ;
					ProteinStructure->Residues[IndexResidueID].zNeutronScattering = zCenterOfNeutronScattering / NeutronScatteringLengthOfResidue ;

					strcpy( ProteinStructure->Residues[IndexResidueID].Name, PreviousResidueName) ; //
					ProteinStructure->Residues[IndexResidueID].ID = PreviousResidueID ;

					// reset 
					VolumeOfResidue = 0.0 ;
					WeightOfResidue = 0.0 ;

					XRayScatteringLengthOfResidue = 0.0 ;
					NeutronScatteringLengthOfResidue = 0.0 ;

					xCenterOfVolume = 0.0 ;
					yCenterOfVolume = 0.0 ;
					zCenterOfVolume = 0.0 ;

					xCenterOfXRayScattering = 0.0 ;
					yCenterOfXRayScattering = 0.0 ;
					zCenterOfXRayScattering = 0.0 ;

					xCenterOfNeutronScattering = 0.0 ;
					yCenterOfNeutronScattering = 0.0 ;
					zCenterOfNeutronScattering = 0.0 ;

					// update index for a new Residue except for the first
					++IndexResidueID ;
				}

				// backup the ID and name of the new residue for the case when a new residue is found and the residue stuff needs to be assigned for the previous one
				strcpy( PreviousResidueName, CurrentResidueName ) ;
				PreviousResidueID = ResidueID ;
			}

			// in case current atom belongs to modification, count them
			if ( strcmp( CurrentResidueName, ProteinStructure->ModificationName) == 0 ) { ++NumberOfModificationAtoms ; }


			// assign information such as scattering lengths, volume and weight for current atom
			AssignAtom( ProteinStructure, IndexAtomID, WriteLog, logfile) ;


			// update values for residue with the values of the current (het)atom
			VolumeOfResidue += ProteinStructure->Atoms[IndexAtomID].Volume ;
			WeightOfResidue += ProteinStructure->Atoms[IndexAtomID].Weight ;

			XRayScatteringLengthOfResidue += ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength ;
			NeutronScatteringLengthOfResidue += ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength ;

			xCenterOfVolume += ProteinStructure->Atoms[IndexAtomID].Volume * xDummy ;
			yCenterOfVolume += ProteinStructure->Atoms[IndexAtomID].Volume * yDummy ;
			zCenterOfVolume += ProteinStructure->Atoms[IndexAtomID].Volume * zDummy ;

			xCenterOfXRayScattering += ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength * xDummy ;
			yCenterOfXRayScattering += ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength * yDummy ;
			zCenterOfXRayScattering += ProteinStructure->Atoms[IndexAtomID].XRayScatteringLength * zDummy ;

			xCenterOfNeutronScattering += ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength * xDummy ;
			yCenterOfNeutronScattering += ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength * yDummy ;
			zCenterOfNeutronScattering += ProteinStructure->Atoms[IndexAtomID].NeutronScatteringLength * zDummy ;


			// in case the last residue is processed the information for ProteinStructure->Residues[IndexResidueID].* must be assigned/updated each time
			if ( IndexResidueID == ProteinStructure->NumberOfResidues - 1 )
			{
				ProteinStructure->Residues[IndexResidueID].Volume = VolumeOfResidue ;
				ProteinStructure->Residues[IndexResidueID].Weight = WeightOfResidue ;

				ProteinStructure->Residues[IndexResidueID].XRayScatteringLength = XRayScatteringLengthOfResidue ;
				ProteinStructure->Residues[IndexResidueID].NeutronScatteringLength = NeutronScatteringLengthOfResidue ;

				ProteinStructure->Residues[IndexResidueID].xVolume = xCenterOfVolume / VolumeOfResidue ;
				ProteinStructure->Residues[IndexResidueID].yVolume = yCenterOfVolume / VolumeOfResidue ;
				ProteinStructure->Residues[IndexResidueID].zVolume = zCenterOfVolume / VolumeOfResidue ;

				ProteinStructure->Residues[IndexResidueID].xXRayScattering = xCenterOfXRayScattering / XRayScatteringLengthOfResidue ;
				ProteinStructure->Residues[IndexResidueID].yXRayScattering = yCenterOfXRayScattering / XRayScatteringLengthOfResidue ;
				ProteinStructure->Residues[IndexResidueID].zXRayScattering = zCenterOfXRayScattering / XRayScatteringLengthOfResidue ;

				ProteinStructure->Residues[IndexResidueID].xNeutronScattering = xCenterOfNeutronScattering / NeutronScatteringLengthOfResidue ;
				ProteinStructure->Residues[IndexResidueID].yNeutronScattering = yCenterOfNeutronScattering / NeutronScatteringLengthOfResidue ;
				ProteinStructure->Residues[IndexResidueID].zNeutronScattering = zCenterOfNeutronScattering / NeutronScatteringLengthOfResidue ;

				strcpy( ProteinStructure->Residues[IndexResidueID].Name, PreviousResidueName) ; //
				ProteinStructure->Residues[IndexResidueID].ID = PreviousResidueID ;
			}


			// finally update index for Atom
			++IndexAtomID ;

		}
		else
		{
			// print lines not starting with ATOM or HETATM (Linebuffer incl newline char)
			fprintf( PointerToSkippedLinesPDBFile, "%s", Linebuffer) ;
		}
	}

	ProteinStructure->NumberOfModificationAtoms = NumberOfModificationAtoms ;

	if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\t\tFound %d atoms belonging to modification: %s\n", NumberOfModificationAtoms, ProteinStructure->ModificationName) ; }


	fclose( PointerToSkippedLinesPDBFile ) ;
	fclose( PointerToPDBFile ) ;


	CenterPDBAtom( ProteinStructure, WriteLog, logfile ) ;
}
