/// Declaration of Dataset-structure
struct Dataset {
	double * QValues;
	double * IValues;
	double * SigmaValues;
	double * SigmaQValues;
	double ** ResolutionWeights;
	double * FitValues;
	double Concentration;
	double Contrast;
	double * Constraints;
	double * ScatteringLengths;
	bool IncludeResolutionEffects;
	int NMin;
	int NMax;
	int NumberOfDatapoints;
	char Filename[256];
};

/// Declaration of parameter-structure
struct Parameter
{
	double Value;
	bool iParameter;
	double MinValue;
	double MaxValue;
	double Error;
	char Name[64];
};

/// Declaration of structures used for PDB-files
struct Atom
{
	// Coordinates
	double x;
	double y;
	double z;

	// Physical properties
	double XRayScatteringLength;
	double NeutronScatteringLength;
	double Volume;
	double Weight;

	// Book-keeping
//	int ID ;
	char Name; // not assigned / used in ImportAtomsFromPDBFile, needs to be [5]
//	char Type[3] ;
//	int ResidueID ;
//	char ResidueName[4] ;
};

struct Residue
{
	// Coordinates of center of volume
	double xVolume;
	double yVolume;
	double zVolume;

	double xXRayScattering;
	double yXRayScattering;
	double zXRayScattering;

	double xNeutronScattering;
	double yNeutronScattering;
	double zNeutronScattering;

	// Physical properties
	double XRayScatteringLength;
	double NeutronScatteringLength;
	double Volume;
	double Weight;

	// Book-keeping
	int ResidueID;
	char Name[4]; // rename to ResidueName or rename ResidueID to ID
	char AtomName[3]; // Fake name for atom. works when each atom is considered a residue // remove
};

struct Protein
{
	struct Atom * Atoms ;
	int NumberOfAtoms ;

	int NumberOfHAtoms ;
	int NumberOfDAtoms ;
	int NumberOfCAtoms ;
	int NumberOfNAtoms ;
	int NumberOfOAtoms ;
	int NumberOfPAtoms ;
	int NumberOfSAtoms ;
	// int NumberOfIAtoms ;

	int NumberOfZNAtoms ;
	int NumberOfCLAtoms ;
	int NumberOfNAAtoms ;
	int NumberOfCAAtoms ;
	int NumberOfFEAtoms ;

	int NumberOfQAtoms ;

	struct Residue * Residues ;
	int NumberOfResidues ;

	char ModificationName[4] ;
	int NumberOfModificationAtoms ;

	char PDBFileLocation[256] ;
	double Weight ;
} ;
