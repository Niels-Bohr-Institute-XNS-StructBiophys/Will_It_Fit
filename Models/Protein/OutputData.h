void OutputData(double ChiSquare, double QMin, double QMax, struct Parameter * Parameters, int NumberOfParameters,
                struct Dataset * Data, int NumberOfSpectra, char cardfilename[128], struct Protein ProteinStructure, 
                struct UserDefined UserDefinedStructure, char SampleFilename[256], char *ResultsDirectory)
{
    /// Declarations
    // Declare dummy variables used in function
    int i;

    // Variables describing the output file
    FILE *fp;
    //const char *filename;
    char Filename[256];
    struct stat status = { 0 };
    sprintf(Filename, "%s/Results.wif", ResultsDirectory);
    //filename = "Results.wif";

    /// I/O
    // Create the output file
    printf("%s\n", Filename);
    fp = fopen(Filename, "w+");
    //Check if results directory exists, and create it if not
    //Create directory for output
    //This might be platform specific
    if(stat (ResultsDirectory, &status) ==-1){
        mkdir(ResultsDirectory, 0700);
        fp = fopen(Filename, "w+");
    }
    // Print filfilenames to file
    
    fprintf(fp, ".card-file: \n");
    fprintf(fp, "%s \n", cardfilename);

    fprintf(fp, "\n");
    fprintf(fp, "Datafiles: \n");
    for(i = 0; i < NumberOfSpectra; ++i){
        fprintf(fp, "%s \n", Data[i].Filename);
    }

    fprintf(fp, "\n");

    // Print location of sample file
    fprintf(fp, "%s \n", SampleFilename);

    fprintf(fp, "\n");

    // Print fit quality to file
    fprintf(fp, "Final chisq = %g \n", ChiSquare);

    fprintf(fp, "\n");

    // Print range of q to file
    fprintf(fp, "Lower limit on q = %g \n", QMin);
    fprintf(fp, "Upper limit on q = %g \n", QMax);

    fprintf(fp, "\n");

    // Print parameters and properties of parameters to file
    fprintf(fp, "Parameters: \n");
    fprintf(fp, "Name:             Value:            Error: \n");
    for(i = 0; i < NumberOfParameters; ++i){
        if (Parameters[i].iParameter == true) {
            fprintf(fp, "%-25s   %-15f   %-15g \n", Parameters[i].Name, Parameters[i].Value, Parameters[i].Error);
        } else {
            fprintf(fp, "%-25s   %-15f   Fixed \n", Parameters[i].Name, Parameters[i].Value);
        }
    }

    fprintf(fp, "\n");
    // Print info on protein file
    //
//    int AtomsToPrint = 10;
//    fprintf(fp,"Residue name, Residue ID. Note: when doing all atom calculations each atom is considered a residue and there will be multiple \"residues\" with the same name but different IDs\n");
//    for(i = 0; i < AtomsToPrint; ++i){
//    fprintf(fp, "%s %d\n",ProteinStructure.Residues[i].Name,i);
//    }
    // Calculate Molecular Weigth
    double Weight = 0.0;
    double WeightModification = 0.0;
    double VolumeModification = 0.0; 
    for (i = 0; i < ProteinStructure.NumberOfResidues; ++i){
	    if(strcmp(ProteinStructure.Residues[i].Name, "W") != 0 || strcmp(ProteinStructure.Residues[i].Name, "HOH") != 0 ){ // Skip waters in calculation of molecular weight
	   
	Weight += ProteinStructure.Residues[i].Weight;
	    }
	 if(strcmp(ProteinStructure.Residues[i].Name, "  X") == 0 ){ // Calculate Molecular Weight of Fatty Acids
	VolumeModification +=  ProteinStructure.Residues[i].Volume;
	WeightModification += ProteinStructure.Residues[i].Weight;
	    }    
    }
    fprintf(fp, "Molecular Weight: %f\n",Weight);
    fprintf(fp, "Molecular Weight of Modification (Total): %f\n",WeightModification);
    fprintf(fp, "Volume of Modification (Total): %f\n",VolumeModification*(1./(1.+Parameters[8].Value)));
    fprintf(fp, "Specific Volume of Modification (Total): %f\n",VolumeModification * 1e-24 *6.022e23* (1./(1.+Parameters[8].Value))/(WeightModification ));


    // Close file and end program
    fclose(fp);
}
