void OutputData(double ChiSquare, double QMin, double QMax, struct Parameter * Parameters, int NumberOfParameters,
                struct Dataset * Data, int NumberOfSpectra, char cardfilename[128], struct Protein ProteinStructure, 
                struct UserDefined UserDefinedStructure, char SampleFilename[256])
{
    /// Declarations
    // Declare dummy variables used in function
    int i;

    // Variables describing the output file
    FILE *fp;
    const char *filename;
    filename = "Results.wif";

    /// I/O
    // Create the output file
    fp = fopen(filename, "w+");

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
    fprintf(fp, "Sample info-file: \n");
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
    int AtomsToPrint = 10;
    fprintf(fp,"Residue name, Residue ID. Note: when doing all atom calculations each atom is considered a residue and there will be multiple \"residues\" with the same name but different IDs\n")
    for(i = 0; i < AtomsToPrint; ++i){
    fprintf(fp, "%c%c%c%d\n",ProteinStructure.Residues[i].Name[0],ProteinStructure.Residues[i].Name[1],ProteinStructure.Residues[i].Name[2],i);
    }
    // Close file and end program
    fclose(fp);
}
