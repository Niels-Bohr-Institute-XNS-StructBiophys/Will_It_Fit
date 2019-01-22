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
	// Book-keeping
	char Name[5] ; // atom name incl terminating NULL (\0)
	char Type[3] ; // element symbol incl terminating NULL (\0)
	int ID ; // normally unqiue or can be easily made unique in PDBs

	char ResidueName[4] ; // associated residue name incl terminating NULL (\0)
	int ResidueID ; // associated residue ID, not necessarily unique when using flattened PDB structures

	// Physical properties
	double Volume ;
	double Weight ;
	double XRayScatteringLength ;
	double NeutronScatteringLength ;

	// Coordinates
	double x ;
	double y ;
	double z ;
} ;

struct Residue
{
	// Book-keeping
	char Name[4] ; // residue name incl terminating NULL (\0)
	int ID ; // residue ID

	// Physical properties
	double Volume ;
	double Weight ;
	double XRayScatteringLength ;
	double NeutronScatteringLength ;

	// Coordinates of center of volume
	double xVolume ;
	double yVolume ;
	double zVolume ;

	double xXRayScattering ;
	double yXRayScattering ;
	double zXRayScattering ;

	double xNeutronScattering ;
	double yNeutronScattering ;
	double zNeutronScattering ;
} ;

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

	char ModificationName[4] ; // residue name of modification incl terminating NULL (\0)
	int NumberOfModificationAtoms ;

	char PDBFileLocation[256] ;
	double Weight ;
} ;

/// Function used to assign fitting ranges
void CopyProteinStructure( struct Protein * ProteinStructureCopy, struct Protein * ProteinStructure)
{

	ProteinStructureCopy->NumberOfHAtoms = ProteinStructure->NumberOfHAtoms ;
	ProteinStructureCopy->NumberOfDAtoms = ProteinStructure->NumberOfDAtoms ;
	ProteinStructureCopy->NumberOfCAtoms = ProteinStructure->NumberOfCAtoms ;
	ProteinStructureCopy->NumberOfNAtoms = ProteinStructure->NumberOfNAtoms ;
	ProteinStructureCopy->NumberOfOAtoms = ProteinStructure->NumberOfOAtoms ;
	ProteinStructureCopy->NumberOfPAtoms = ProteinStructure->NumberOfPAtoms ;
	ProteinStructureCopy->NumberOfSAtoms = ProteinStructure->NumberOfSAtoms ;
	// ProteinStructureCopy->NumberOfIAtoms = ProteinStructure->NumberOfIAtoms ;

	ProteinStructureCopy->NumberOfZNAtoms = ProteinStructure->NumberOfZNAtoms ;
	ProteinStructureCopy->NumberOfCLAtoms = ProteinStructure->NumberOfCLAtoms ;
	ProteinStructureCopy->NumberOfNAAtoms = ProteinStructure->NumberOfNAAtoms ;
	ProteinStructureCopy->NumberOfCAAtoms = ProteinStructure->NumberOfCAAtoms ;
	ProteinStructureCopy->NumberOfFEAtoms = ProteinStructure->NumberOfFEAtoms ;

	ProteinStructureCopy->NumberOfQAtoms = ProteinStructure->NumberOfQAtoms ;


	strcpy( ProteinStructureCopy->ModificationName, ProteinStructure->ModificationName) ;
	ProteinStructureCopy->NumberOfModificationAtoms = ProteinStructure->NumberOfModificationAtoms ;


	strcpy( ProteinStructureCopy->PDBFileLocation, ProteinStructure->PDBFileLocation) ;
	ProteinStructureCopy->Weight = ProteinStructure->Weight ;


	// copy Residues stuff
	for ( int k = 0; k < ProteinStructureCopy->NumberOfResidues; ++k)
	{
		strcpy( ProteinStructureCopy->Residues[k].Name, ProteinStructure->Residues[k].Name) ;
		ProteinStructureCopy->Residues[k].ID = ProteinStructure->Residues[k].ID ;

		ProteinStructureCopy->Residues[k].Volume                  = ProteinStructure->Residues[k].Volume ;
		ProteinStructureCopy->Residues[k].Weight                  = ProteinStructure->Residues[k].Weight ;
		ProteinStructureCopy->Residues[k].XRayScatteringLength    = ProteinStructure->Residues[k].XRayScatteringLength ;
		ProteinStructureCopy->Residues[k].NeutronScatteringLength = ProteinStructure->Residues[k].NeutronScatteringLength ;


		ProteinStructureCopy->Residues[k].xVolume = ProteinStructure->Residues[k].xVolume ;
		ProteinStructureCopy->Residues[k].yVolume = ProteinStructure->Residues[k].yVolume ;
		ProteinStructureCopy->Residues[k].zVolume = ProteinStructure->Residues[k].zVolume ;

		ProteinStructureCopy->Residues[k].xXRayScattering = ProteinStructure->Residues[k].xXRayScattering ;
		ProteinStructureCopy->Residues[k].yXRayScattering = ProteinStructure->Residues[k].yXRayScattering ;
		ProteinStructureCopy->Residues[k].zXRayScattering = ProteinStructure->Residues[k].zXRayScattering ;

		ProteinStructureCopy->Residues[k].xNeutronScattering = ProteinStructure->Residues[k].xNeutronScattering ;
		ProteinStructureCopy->Residues[k].yNeutronScattering = ProteinStructure->Residues[k].yNeutronScattering ;
		ProteinStructureCopy->Residues[k].zNeutronScattering = ProteinStructure->Residues[k].zNeutronScattering ;
	}

	// copy Atoms stuff
	for ( int k = 0; k < ProteinStructureCopy->NumberOfAtoms; ++k)
	{
		strcpy( ProteinStructureCopy->Atoms[k].Name, ProteinStructure->Atoms[k].Name) ;
		strcpy( ProteinStructureCopy->Atoms[k].Type, ProteinStructure->Atoms[k].Type) ;
		ProteinStructureCopy->Atoms[k].ID = ProteinStructure->Atoms[k].ID ;

		strcpy( ProteinStructureCopy->Atoms[k].ResidueName, ProteinStructure->Atoms[k].ResidueName) ;
		ProteinStructureCopy->Atoms[k].ResidueID = ProteinStructure->Atoms[k].ResidueID ;

		ProteinStructureCopy->Atoms[k].Volume                  = ProteinStructure->Atoms[k].Volume ;
		ProteinStructureCopy->Atoms[k].Weight                  = ProteinStructure->Atoms[k].Weight ;
		ProteinStructureCopy->Atoms[k].XRayScatteringLength    = ProteinStructure->Atoms[k].XRayScatteringLength ;
		ProteinStructureCopy->Atoms[k].NeutronScatteringLength = ProteinStructure->Atoms[k].NeutronScatteringLength ;

		ProteinStructureCopy->Atoms[k].x = ProteinStructure->Atoms[k].x ;
		ProteinStructureCopy->Atoms[k].y = ProteinStructure->Atoms[k].y ;
		ProteinStructureCopy->Atoms[k].z = ProteinStructure->Atoms[k].z ;
	}
}
