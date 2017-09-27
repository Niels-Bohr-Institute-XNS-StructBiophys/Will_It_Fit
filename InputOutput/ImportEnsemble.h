int CheckSizeOfEnsemble(char EnsembleLocation[256], int ** ProteinResidues, int ** ProteinAtoms, bool CMD){

	int i;
	int NumberOfPDBs;
	int AtomsInCurrentPDB;
	int MostAtomsInAPDB = 0;
	int ResiduesInCurrentPDB;
	int MostResiduesInAPDB = 0;
	int linenum;
    int * AtomDummy;
    int * ResidueDummy;

	//Is the .pdb files located in a subdirectory
	bool SubDir = false;

	char Buffer[128];
	char Directory[256];
	memset(Directory, 0, sizeof(Directory)); //Avoid random values in Directory path


	FILE *fp;
	char Path[sizeof(Buffer)];
	char Dummy[1+10*2][sizeof(Buffer)]; //Room for 10 datasets

	fp = fopen(EnsembleLocation, "r");

	if (fp == NULL) {
        printf("An error occured when attempting to open the .ens-file. \n");
        return -1;
    }

    printf(".ens-file found!\n");



	SubDir = DirectoryFinder(EnsembleLocation, Directory);

    linenum = 0;

    while(fgets(Buffer, sizeof(Buffer), fp) != 0){
    	sscanf(Buffer, "%s", Dummy[linenum]);
    	++linenum;
    }

    fclose(fp);

    NumberOfPDBs = atoi(Dummy[0]);

    //Initialize Atom and Residue arrays
    Initialize1DIntegerArray(&AtomDummy, NumberOfPDBs);
    Initialize1DIntegerArray(&ResidueDummy, NumberOfPDBs);
    for(i = 0; i <= NumberOfPDBs; ++i){
    	if(SubDir){
    		sprintf(Path, "%s%s", Directory, Dummy[1+i*2]);
    	}
    	else{
    		sprintf(Path, "%s", Dummy[1+i*2]);
    	}
    	ResiduesInCurrentPDB = CheckNumberOfResiduesInPDBFile(Path);
    	AtomsInCurrentPDB = CheckNumberOfAtomsInPDBFile(Path);
        ResidueDummy[i] = ResiduesInCurrentPDB;
    	AtomDummy[i] = AtomsInCurrentPDB;
    }
    *ProteinAtoms = AtomDummy;
    *ProteinResidues = ResidueDummy;
    return NumberOfPDBs;
}

void ImportEnsemble(struct Protein * Ensemble, double ** Weights, char EnsembleLocation[256], int NumberOfPDBs, int * ProteinResidues, int * ProteinAtoms, bool CMD){
	int i;
    int j;
	//Is the .pdb files located in a subdirectory
    bool SubDir;

    char Buffer[128];
    char Directory[256];
    memset(Directory, 0, sizeof(Directory)); //Avoid random values in Directory path
    double * WeightDummy;

    int linenum;

    FILE *fp;

    char Path[sizeof(Buffer)];
    char Dummy[1+10*2][sizeof(Buffer)]; //Room for 10 datasets 

    fp  = fopen(EnsembleLocation, "r");

    if (fp == NULL) {
        printf("An error occured when attempting to open the .ens-file. \n");
    }

    printf(".ens-file found!\n");

    SubDir = DirectoryFinder(EnsembleLocation, Directory);

    while(fgets(Buffer, sizeof(Buffer), fp) != 0){
    	sscanf(Buffer, "%s", Dummy[linenum]);
    	++linenum;
    }

    fclose(fp);
    Initialize1DArray(&WeightDummy, NumberOfPDBs);
    for (i = 0; i < NumberOfPDBs; ++i)
    {
    	if (SubDir)
    	{
    		sprintf(Path, "%s%s", Directory, Dummy[1+i*2]);
    	}
    	else
    	{
    		sprintf(Path, "%s", Dummy[1+i*2]);
    	}
        WeightDummy[i] = atof(Dummy[2+i*2]);
        
        ImportResiduesFromPDBFile(Path, Ensemble[i], ProteinResidues[i]);
        ImportAtomsFromPDBFile(Path, Ensemble[i], ProteinAtoms[i]);
        Ensemble[i].NumberOfResidues = ProteinResidues[i];
        Ensemble[i].NumberOfAtoms = ProteinAtoms[i];
        *Weights = WeightDummy;
    }
}
//Simple hack for cleaning -nan values from the first element of the ensemble
void CleanEnsemble(struct Protein * Ensemble, int NumberOfProteins){
    int i;
    //Clean first residue in each protein
    for(i = 0; i < NumberOfProteins; i++){
        Ensemble[i].Residues[0].xVolume = 0.0;
        Ensemble[i].Residues[0].yVolume = 0.0;
        Ensemble[i].Residues[0].zVolume = 0.0;
        Ensemble[i].Residues[0].xXRayScattering = 0.0;
        Ensemble[i].Residues[0].yXRayScattering = 0.0;
        Ensemble[i].Residues[0].zXRayScattering = 0.0;
        Ensemble[i].Residues[0].xNeutronScattering = 0.0;
        Ensemble[i].Residues[0].yNeutronScattering = 0.0;
        Ensemble[i].Residues[0].zNeutronScattering = 0.0;

        // Physical properties
        Ensemble[i].Residues[0].XRayScatteringLength = 0.0;
        Ensemble[i].Residues[0].NeutronScatteringLength = 0.0;
        Ensemble[i].Residues[0].Volume = 1.0;
        Ensemble[i].Residues[0].Weight = 0.0;
    }


}