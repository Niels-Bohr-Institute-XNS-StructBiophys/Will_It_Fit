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
			memcpy( dummy, Linebuffer + 6, 5 * sizeof Linebuffer[0]) ; // here read atom ID (pos 7-11) instead of residue ID (pos 23-26)
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


/*
	Currently supported isotopes:
	" H"," D"," C"," N"," O"," P"," N"," S"," Q" (water),"Cl","Ca","Zn","Fe","Na"

	atom radii for displaced solvent volume from
	https://www.saxier.org/forum/viewtopic.php?t=93

		The radii is a consensus from several sources:
		Radii from Fraser et al.,
		http://www.sicyon.com/
		http://www.webelements.com/
		C CH CH2 CH3 N NH NH2 NH3 O OH
		1.577,1.73,1.85,1.97,0.84,1.22,1.45,1.62,1.30,1.50,

		S SH P FE CU CA MG MN ZN NHH NHM*
		1.68,1.81,1.11,1.24,1.28,1.97,1.60,1.30,1.33,1.22,1.22,

		MO BR RB AR K XE NI CO KR NA HG
		2.01, 1.85, 2.98, 1.88, 2.75, 2.16, 1.63, 1.67, 2.02, 2.27, 1.75,

		Ag Au U SE CL F Al Si Cr Sr Pd
		1.75, 1.79, 1.86, 1.90, 1.75, 1.47, 1.82, 2.10, 1.85, 2.45, 1.63,

		I Cs Ba Pt Ga Ge Cd H
		1.98, 3.34, 2.78, 1.75, 1.87, 1.52, 1.58, 1.07

		*) NHH is NH Histidine, NHM is NH main chain


	H, C, O, N from Fraser et al paper / Crysol paper -> usually atoms bound in proteins
	P, S, Ca, Mn, Fe, Cu, Zn from Crysol paper / IUCr Tables (not found there) -> usually atoms bound / coordinated in proteins
	Na, K, Ni, Pd, Pt, Cd, Ga, Si, Se, F, Cl, Br, I, Ar, Kr, Xe, U are exactly the VdW radii from Wikipedia or webelements.com -> usually atoms not bound in proteins
	Rb, Cs, Sr, Ba, Ag, Al (very) similar to VdW radii from Wikipedia -> usually atoms not bound in proteins
	Mg, Cr, Co, Au, Hg, Ge no good reference found yet

	see also
	http://lorentz.dynstr.pasteur.fr/suny/index.php?id0=aquasaxs
	http://scripts.iucr.org/cgi-bin/paper?S0021889895007047 (Crysol paper 1995)
	http://scripts.iucr.org/cgi-bin/paper?S0021889878014296 (Fraser paper 1978)

	Note the radii are not listed in the IUCr Tables as stated in the Crysol paper.
*/
/*
void AssignAtom(char AtomName[3], double *XRayScatteringLengthOfCurrentAtom, double *NeutronScatteringLengthOfCurrentAtom, double *VolumeOfCurrentAtom, double *WeightOfCurrentAtom, int WriteLog, FILE* logfile)
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

	const double H2OVolume =  30.0;
	// const double D2OVolume ???


	// Isotope/Compound weights [u]
	const double HWeight = 1.008;
	// const double DWeight ???
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
	// const double H2OWeight ???
	// const double D2OWeight ???


	// X-ray scattering lengths [cm]]
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

	const double H2OXRayScatteringLength =  10.0 * r0; // 2.818e-12
	// const double D2OXRayScatteringLength ???


	// neutron scattering lengths [cm]]
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

	// const double H2OXNeutronScatteringLength ???
	// const double D2OXNeutronScatteringLength ???



	// the following if statements are to read the atoms. Note this was done in a switch statement before
	// but that does not work for multi character atom names such as "Zn"

	// single char atom names
	if (AtomName[0] ==  ' ' && AtomName[1] == 'H' ){
		*XRayScatteringLengthOfCurrentAtom = HXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = HNeutronScatteringLength;
		*VolumeOfCurrentAtom = HVolume;
		*WeightOfCurrentAtom = HWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] ==  ' ' && AtomName[1] == 'D' ){
		*XRayScatteringLengthOfCurrentAtom = DXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = DNeutronScatteringLength;
		*VolumeOfCurrentAtom = DVolume;
		*WeightOfCurrentAtom = HWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'C' ){
		*XRayScatteringLengthOfCurrentAtom = CXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = CNeutronScatteringLength;
		*VolumeOfCurrentAtom = CVolume;
		*WeightOfCurrentAtom = CWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'N' ){
		*XRayScatteringLengthOfCurrentAtom = NXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = NNeutronScatteringLength;
		*VolumeOfCurrentAtom = NVolume;
		*WeightOfCurrentAtom = NWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'O' ){
		*XRayScatteringLengthOfCurrentAtom = OXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = ONeutronScatteringLength;
		*VolumeOfCurrentAtom = OVolume;
		*WeightOfCurrentAtom = OWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'P' ){
		*XRayScatteringLengthOfCurrentAtom = PXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = PNeutronScatteringLength;
		*VolumeOfCurrentAtom = PVolume;
		*WeightOfCurrentAtom = PWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'S' ){
		*XRayScatteringLengthOfCurrentAtom = SXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = SNeutronScatteringLength;
		*VolumeOfCurrentAtom = SVolume;
		*WeightOfCurrentAtom = SWeight;
		AtomRecg = 1;
	}


	// double char atom names
	if (AtomName[0] == 'Z' && AtomName[1] == 'N' ){
		*XRayScatteringLengthOfCurrentAtom = ZNXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = ZNNeutronScatteringLength;
		*VolumeOfCurrentAtom = ZNVolume;
		*WeightOfCurrentAtom = ZNWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] == 'C' && AtomName[1] == 'L' ){
		*XRayScatteringLengthOfCurrentAtom = CLXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = CLNeutronScatteringLength;
		*VolumeOfCurrentAtom = CLVolume;
		*WeightOfCurrentAtom = CLWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] == 'C' && AtomName[1] == 'A' ){
		*XRayScatteringLengthOfCurrentAtom = CAXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = CANeutronScatteringLength;
		*VolumeOfCurrentAtom = CAVolume;
		*WeightOfCurrentAtom = CAWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] == 'F' && AtomName[1] == 'E' ){
		*XRayScatteringLengthOfCurrentAtom = FEXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = FENeutronScatteringLength;
		*VolumeOfCurrentAtom = FEVolume;
		*WeightOfCurrentAtom = FEWeight;
		AtomRecg = 1;
	}
	if (AtomName[0] == 'N' && AtomName[1] == 'A' ){
		*XRayScatteringLengthOfCurrentAtom = NAXRayScatteringLength;
		*NeutronScatteringLengthOfCurrentAtom = NANeutronScatteringLength;
		*VolumeOfCurrentAtom = NAVolume;
		*WeightOfCurrentAtom = NAWeight;
		AtomRecg = 1;
	}

	// special "atoms" like waters
	if (AtomName[0] == ' ' && AtomName[1] == 'Q' ){
		*XRayScatteringLengthOfCurrentAtom = 4.133*H2OXRayScatteringLength; // 1.165E-11
		*NeutronScatteringLengthOfCurrentAtom = 4.133*H2OXRayScatteringLength; // Note that this is a dummy value. for neutron contrasts the scattering lengths of the dummy waters are set in Model.h
		*VolumeOfCurrentAtom = 4.133*(H2OVolume); // 123.990
		*WeightOfCurrentAtom = 0.0; //4.133*(2*HWeight+OWeight);
		AtomRecg = 1;
	}


	// if not found error and exit
	if (AtomRecg == 0)
	{
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "AssignAtom(): Atom name %s not found in database. Exit.\n", AtomName) ; }
		fflush(logfile) ;
		exit(0) ;
	}
	return;
}
*/



void AssignAtom( char AtomName[3], struct Protein * ProteinStructure, int CurrentAtomID, int WriteLog, FILE* logfile)
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

	const double QVolume =  30.0;


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

	const double QWeight = 2*HWeight + OWeight ;


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

	const double QXRayScatteringLength =  10.0 * r0; // 2.818e-12


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

	// const double QNeutronScatteringLength ??? // assigned later by contrast, since it depends on sample




	// the following if statements are to read the atoms. Note this was done in a switch statement before
	// but that does not work for multi character atom names such as "Zn"

	// single char atom names
	if (AtomName[0] ==  ' ' && AtomName[1] == 'H' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = HXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = HNeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = HVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = HWeight ;
		ProteinStructure->NumberOfHAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if (AtomName[0] ==  ' ' && AtomName[1] == 'D' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = DXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = DNeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = DVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = DWeight ;
		ProteinStructure->NumberOfDAtoms += 1 ;
		AtomRecg = 1;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'C' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = CXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = CNeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = CVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = CWeight ;
		ProteinStructure->NumberOfCAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'N' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = NXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = NNeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = NVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = NWeight ;
		ProteinStructure->NumberOfNAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'O' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = OXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = ONeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = OVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = OWeight ;
		ProteinStructure->NumberOfOAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'P' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = PXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = PNeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = PVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = PWeight ;
		ProteinStructure->NumberOfPAtoms += 1 ;
		AtomRecg = 1;
	}
	if (AtomName[0] == ' ' && AtomName[1] == 'S' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = SXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = SNeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = SVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = SWeight ;
		ProteinStructure->NumberOfSAtoms += 1 ;
		AtomRecg = 1 ;
	}
	/*
	if (AtomName[0] == ' ' && AtomName[1] == 'I' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = IXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = INeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = IVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = IWeight ;
		ProteinStructure->NumberOfIAtoms += 1 ;
		AtomRecg = 1 ;
	}
	*/

	// double char atom names
	if (AtomName[0] == 'Z' && AtomName[1] == 'N' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = ZNXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = ZNNeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = ZNVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = ZNWeight ;
		ProteinStructure->NumberOfZNAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if (AtomName[0] == 'C' && AtomName[1] == 'L' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = CLXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = CLNeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = CLVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = CLWeight ;
		ProteinStructure->NumberOfCLAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if (AtomName[0] == 'N' && AtomName[1] == 'A' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = NAXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = NANeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = NAVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = NAWeight ;
		ProteinStructure->NumberOfNAAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if (AtomName[0] == 'C' && AtomName[1] == 'A' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = CAXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = CANeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = CAVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = CAWeight ;
		ProteinStructure->NumberOfCAAtoms += 1 ;
		AtomRecg = 1 ;
	}
	if (AtomName[0] == 'F' && AtomName[1] == 'E' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = FEXRayScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = FENeutronScatteringLength ;
		ProteinStructure->Atoms[CurrentAtomID].Volume = FEVolume ;
		ProteinStructure->Atoms[CurrentAtomID].Weight = FEWeight ;
		ProteinStructure->NumberOfFEAtoms += 1 ;
		AtomRecg = 1 ;
	}

	// special "atoms" like dummy waters (residues WAT with atom names Q)
	// note that crystal waters (residues HOH) consist of and are treated like normal atoms, i.e. O and 2*H atoms
	if (AtomName[0] == ' ' && AtomName[1] == 'Q' )
	{
		ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength = 4.133*QXRayScatteringLength ; // 1.165E-11
		ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength = 4.133*QXRayScatteringLength ; // Note that this is a dummy value. For neutron contrasts the scattering lengths of the dummy waters are set in Model.h
		ProteinStructure->Atoms[CurrentAtomID].Volume = 4.133 * QVolume ; // 123.990
		ProteinStructure->Atoms[CurrentAtomID].Weight = 4.133 * QWeight ; // Note that this is the value for dummy-H2O, For neutron contrasts H2O/D2O mixtures the correct weights are not considered in Model.h
		ProteinStructure->NumberOfQAtoms += 1 ;
		AtomRecg = 1 ;
	}


	// if not found error and exit
	if (AtomRecg == 0)
	{
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "AssignAtom(): Atom name %s not found in database. Exit.\n", AtomName) ; }
		fflush(logfile) ;
		exit(0) ;
	}
	return;
}











// Fill fields from struct ProteinStructure->Atoms[] from each line starting with ATOM / HETATM label in PDB file
// ProteinStructure has been allocated priorly in function AllocateProteinStructure(...) ../Auxillary/Allocation.h
// .x
// .y
// .z
// .Name
// .XRayScatteringLength
// .NeutronScatteringLength
// .Volume
// .Weight
//
// The number of each atom type will be counted, too
void ImportAtomsFromPDBFile(char Filename[256], struct Protein * ProteinStructure, int WriteLog, FILE *logfile)
{
	// Declarations
	FILE *PointerToPDBFile ;
	char Linebuffer[82] ; // length of PDB file incl newline (\n) and terminating NULL (\0)
	char dummy[82] ; // dummy char array
	char AtomName[3] ; // includes a terminating NULL
	int CurrentAtomID = 0 ;

	// I/O
	// skip file-existence check, done before
	PointerToPDBFile = fopen( Filename, "r") ;

	while( fgets( Linebuffer, sizeof(Linebuffer), PointerToPDBFile) != NULL )
	{
		// extract lines starting with ATOM or HETATM pattern
		if ( strncmp( Linebuffer, "ATOM", 4) == 0 || strncmp( Linebuffer, "HETATM", 6) == 0 )
		{
			memcpy( dummy, Linebuffer + 16, 1 * sizeof Linebuffer[0]) ; // read altLoc from pos 17
			dummy[1] = 0 ;
			// in case of multiple occupancies (alternative locations) count only once
			if ( strcmp( " ", dummy) == 0 || strcmp( "A", dummy) == 0)
			{
				memcpy( dummy, Linebuffer + 30, 8 * sizeof Linebuffer[0]) ; // read x coo (pos 31-38)
				dummy[8] = 0 ;
				delallspc( dummy ) ;
				ProteinStructure->Atoms[CurrentAtomID].x = strtod( dummy, NULL) ;

				memcpy( dummy, Linebuffer + 38, 8 * sizeof Linebuffer[0]) ; // read y coo (pos 39-46)
				dummy[8] = 0 ;
				delallspc( dummy ) ;
				ProteinStructure->Atoms[CurrentAtomID].y = strtod( dummy, NULL) ;

				memcpy( dummy, Linebuffer + 46, 8 * sizeof Linebuffer[0]) ; // read y coo (pos 47-54)
				dummy[8] = 0 ;
				delallspc( dummy ) ;
				ProteinStructure->Atoms[CurrentAtomID].z = strtod( dummy, NULL) ;

				memcpy( AtomName, Linebuffer + 76, 2 * sizeof Linebuffer[0]) ; // read element symbol, keep spaces (pos 77-78)
				AtomName[2] = 0 ;
				// printf( "'%s'\n", AtomName) ; // fflush( NULL) ;

//				AssignAtom( AtomName, &ProteinStructure->Atoms[CurrentAtomID].XRayScatteringLength, &ProteinStructure->Atoms[CurrentAtomID].NeutronScatteringLength, &ProteinStructure->Atoms[CurrentAtomID].Volume, &ProteinStructure->Atoms[CurrentAtomID].Weight, WriteLog, logfile) ;
				AssignAtom( AtomName, ProteinStructure, CurrentAtomID, WriteLog, logfile) ;


				++CurrentAtomID ;
			}
		}
	}

	fclose( PointerToPDBFile );
}



void CenterPDB( struct Protein * ProteinStructure, int WriteLog, FILE* logfile)
{
	// Function for finding the geometric center of PDB Structure and translating it so that the center becomes (0,0,0)
	// Find geometric center for protein structure
	//

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



// Fill fields from struct ProteinStructure->Residues[] for each residue in PDB, calculate from all ATOM / HETATM entries belonging to a residue
// ProteinStructure has been allocated priorly in function AllocateProteinStructure(...) ../Auxillary/Allocation.h
void ImportResiduesFromPDBFile( char Filename[256], struct Protein * ProteinStructure, char *ResultsDirectory, int WriteLog, FILE* logfile)
{
	// Declarations
	FILE *PointerToPDBFile ;
	FILE *PointerToSkippedLinesPDBFile ;

	char BaseFilename[256] ;
	char SkippedLinesPDBFileLocation[256] ;

	char Linebuffer[82] ; // length of PDB file incl newline (\n) and terminating NULL (\0)
	char dummy[82] ; // dummy char array

	char AtomName[3] ; // includes a terminating NULL
	char CurrentAtomName[3]; // includes a terminating NULL

	int PreviousResidueID = 0 ;
	int ResidueID = 0 ;
	int IndexResidueID = 0 ;

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

	struct Protein DummyProteinStructure ; // DummyProteinStructure used for AssignAtom
	DummyProteinStructure.NumberOfAtoms = 1 ;
	DummyProteinStructure.NumberOfResidues = 0 ;
	AllocateProteinStructure( &DummyProteinStructure, DummyProteinStructure.NumberOfResidues, DummyProteinStructure.NumberOfAtoms) ;


	int NumberOfModificationAtoms = 0 ;


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
		// fprintf( logfile, "%s", Linebuffer) ; fflush(logfile) ;

		if ( strncmp( Linebuffer, "ATOM", 4) == 0 || strncmp( Linebuffer, "HETATM", 6) == 0 )
		{
			// current line processing to read current ResidueID, altLoc, coordinates, atom name

			// extract Residue ID from line
			memcpy( dummy, Linebuffer + 6, 5 * sizeof Linebuffer[0]) ; // here read atom ID (pos 7-11) instead of residue ID (pos 23-26)
			dummy[5] = 0 ;
			delallspc( dummy ) ;
			ResidueID = strtol( dummy, NULL, 10) ;

			// read altLoc from line
			memcpy( dummy, Linebuffer + 16, 1 * sizeof Linebuffer[0]) ; // read altLoc from pos 17
			dummy[1] = 0 ;
			// printf("'%s' ", dummy) ; fflush( NULL) ; //

			// in case of multiple occupancies (alternative locations) print line (Linebuffer incl newline char) and continue
			// NOTE THAT IN NICHOLAS VERSION MULTIPLE OCCUPANCIES WERE NOT COUNTED IN CheckNumberOfResiduesInPDBFile BUT THEN READ AND TREATED HERE -> ERROR
			if ( strcmp( " ", dummy) != 0 && strcmp( "A", dummy) != 0)
			{
				fprintf( PointerToSkippedLinesPDBFile, "%s", Linebuffer) ;
				continue ;
			}

			// read x, y, z coordinates and element symbol (here called (Current)AtomName)
			memcpy( dummy, Linebuffer + 30, 8 * sizeof Linebuffer[0]) ; // read x coo (pos 31-38)
			dummy[8] = 0 ;
			delallspc( dummy ) ;
			xDummy = strtod( dummy, NULL) ;

			memcpy( dummy, Linebuffer + 38, 8 * sizeof Linebuffer[0]) ; // read y coo (pos 39-46)
			dummy[8] = 0 ;
			delallspc( dummy ) ;
			yDummy = strtod( dummy, NULL) ;

			memcpy( dummy, Linebuffer + 46, 8 * sizeof Linebuffer[0]) ; // read y coo (pos 47-54)
			dummy[8] = 0 ;
			delallspc( dummy ) ;
			zDummy = strtod( dummy, NULL) ;

			memcpy( CurrentAtomName, Linebuffer + 76, 2 * sizeof Linebuffer[0]) ; // read element symbol, keep spaces (pos 77-78)
			CurrentAtomName[2] = 0 ;
			// printf( "'%s'\n", CurrentAtomName) ; // fflush( NULL) ;

			// in case Residue ID is different from previous one, calculate residue level stuff (ProteinStructure->Residues[IndexResidueID].*) for previous residue
			// exception is if first Residue ID is read, which differs PreviousResidueID (0)
			// after that "local" residue variables are reset to 0.0 and IndexResidueID is incremented for the new residue, and PreviousResidueID is reset
			if ( ResidueID != PreviousResidueID ) {
				if (PreviousResidueID != 0) {

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

					// ProteinStructure->Residues[].AtomName[3] and AtomName[3] -> okay for strcpy with incl terminating NULL
					strcpy( ProteinStructure->Residues[IndexResidueID].AtomName, AtomName) ;

					xCenterOfVolume = 0.0 ;
					yCenterOfVolume = 0.0 ;
					zCenterOfVolume = 0.0 ;
					VolumeOfResidue = 0.0 ;
					WeightOfResidue = 0.0 ;

					xCenterOfXRayScattering = 0.0 ;
					yCenterOfXRayScattering = 0.0 ;
					zCenterOfXRayScattering = 0.0 ;
					XRayScatteringLengthOfResidue = 0.0 ;

					xCenterOfNeutronScattering = 0.0 ;
					yCenterOfNeutronScattering = 0.0 ;
					zCenterOfNeutronScattering = 0.0 ;
					NeutronScatteringLengthOfResidue = 0.0 ;

					++IndexResidueID ;
				}
				// grab here residue name
				// NOTE THIS IS AS IN NICHOLAS VERSION BUT DIFFERENT AS IN "ORIGINAL" WIF VERSION
				// IT IS ACTUALLY BETTER TO UPDATE RESIDUE NAME HERE (I.E. IF ResidueID CHANGES) AND NOT FOR EVERY LINE AS IN "ORIGINAL" WIF
				memcpy( ProteinStructure->Residues[IndexResidueID].Name, Linebuffer + 17, 3 * sizeof Linebuffer[0]) ; // keep spaces (pos 18-20)
				ProteinStructure->Residues[IndexResidueID].Name[3] = 0 ;

				// printf("'%s' ", ProteinStructure->Residues[IndexResidueID].Name) ; fflush( NULL) ;


				// does residue name match modification name? if yes rename current residue name to "  X"
				if ( strcmp( ProteinStructure->Residues[IndexResidueID].Name, ProteinStructure->ModificationName) == 0 )
				{
					strcpy( ProteinStructure->Residues[IndexResidueID].Name, "  X") ;
				}

				// printf("'%s' ", ProteinStructure->Residues[IndexResidueID].Name) ; fflush( NULL) ;

				PreviousResidueID = ResidueID ;
			}


			// further processing of information from current line

			// BEGIN SECTION

			// NOTE THAT THIS SECTION HAS BEEN MOVED HERE (IS NOW AS IN THE "ORIGINAL" WIF VERSION), INSTEAD OF KEEPING IT INSIDE THE PREVIOUS IF-CLAUSE AS IN NICHOLAS VERSION
			// THIS SHOULD HAVE NO EFFECT ON THE COMPUTATIONS SINCE HERE RESIDUES == ATOMS, BUT MAKES MORE SENSE AND PROVIDES BETTER COMPATIBILITY WITH THE "ORIGINAL" VERSION WHERE RESIDUES != ATOMS

			// grab information such as scattering lengths, volume and weight from current atom name
//			AssignAtom( CurrentAtomName, &XRayScatteringLengthOfCurrentAtom, &NeutronScatteringLengthOfCurrentAtom, &VolumeOfCurrentAtom, &WeightOfCurrentAtom, WriteLog, logfile) ;
			AssignAtom( CurrentAtomName, &DummyProteinStructure, 0, WriteLog, logfile) ;

			// in case current atom belongs to modification, count them
			if ( strcmp( ProteinStructure->Residues[IndexResidueID].Name, "  X") == 0 ) { ++NumberOfModificationAtoms ; }
			// printf("'%s' ", ProteinStructure->Residues[IndexResidueID].Name) ; fflush( NULL) ; //

			// END SECTION


			// printf("'%s'\n", CurrentAtomName) ; // fflush( NULL) ;
			strcpy(AtomName, CurrentAtomName);
			// printf("'%s'\n", AtomName) ; // fflush( NULL) ;

			VolumeOfResidue = DummyProteinStructure.Atoms[0].Volume ;
			WeightOfResidue = DummyProteinStructure.Atoms[0].Weight ;
			XRayScatteringLengthOfResidue = DummyProteinStructure.Atoms[0].XRayScatteringLength ;
			NeutronScatteringLengthOfResidue = DummyProteinStructure.Atoms[0].NeutronScatteringLength ;

			xCenterOfVolume = DummyProteinStructure.Atoms[0].Volume * xDummy ;
			yCenterOfVolume = DummyProteinStructure.Atoms[0].Volume * yDummy ;
			zCenterOfVolume = DummyProteinStructure.Atoms[0].Volume * zDummy ;

			xCenterOfXRayScattering = DummyProteinStructure.Atoms[0].XRayScatteringLength * xDummy ;
			yCenterOfXRayScattering = DummyProteinStructure.Atoms[0].XRayScatteringLength * yDummy ;
			zCenterOfXRayScattering = DummyProteinStructure.Atoms[0].XRayScatteringLength * zDummy ;

			xCenterOfNeutronScattering = DummyProteinStructure.Atoms[0].NeutronScatteringLength * xDummy ;
			yCenterOfNeutronScattering = DummyProteinStructure.Atoms[0].NeutronScatteringLength * yDummy ;
			zCenterOfNeutronScattering = DummyProteinStructure.Atoms[0].NeutronScatteringLength * zDummy ;


			// in case the last residue is processed the information for ProteinStructure->Residues[IndexResidueID].* must be updated for each read (het)atom
			// compared to Nicholas WIF code, this section has been moved here, since the information of the local variables (VolumeOfResidue, ...) must be known
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

				// ProteinStructure->Residues[]
				//	.Name[4]
				//	.ModificationName[4]
				//	.AtomName[3]
				// ProteinStructure->Residues[].AtomName[3] and AtomName[3] -> okay for strcpy with incl terminating NULL
				strcpy( ProteinStructure->Residues[IndexResidueID].AtomName, AtomName) ;
			}
		}
		else
		{
			// print lines not starting with ATOM or HETATM (Linebuffer incl newline char)
			fprintf( PointerToSkippedLinesPDBFile, "%s", Linebuffer) ;
		}
	}

	ProteinStructure->NumberOfModificationAtoms = NumberOfModificationAtoms ;

	if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\t\tFound %d atoms belonging to modification: %s\n", NumberOfModificationAtoms, ProteinStructure->ModificationName) ; }


	// free DummyProteinStructure
	free(DummyProteinStructure.Residues) ;
	free(DummyProteinStructure.Atoms) ;



	fclose( PointerToSkippedLinesPDBFile ) ;
	fclose( PointerToPDBFile ) ;

	CenterPDB( ProteinStructure, WriteLog, logfile ) ;
}
