void OutputData(double ChiSquare, double QMin, double QMax, struct Parameter * Parameters, int NumberOfParameters,
                struct Dataset * Data, int NumberOfSpectra, char cardfilename[128], struct Protein ProteinStructure, struct UserDefined UserDefinedStructure,
                char SampleFilename[256])
{
    // Declare dummy variables used in function
    int i;

    // Variables describing the output file
    FILE *Outputfile;
    const char *Filename;
    Filename = "Results.wif";

    // Create the output file
    Outputfile = fopen(Filename, "w+");

    // Print filfilenames to file
    fprintf(Outputfile, ".card-file: \n");
    fprintf(Outputfile, "%s \n", cardfilename);

    fprintf(Outputfile, "\n");
    fprintf(Outputfile, "Datafiles: \n");

    for(i = 0; i < NumberOfSpectra; ++i){
        fprintf(Outputfile, "%s \n", Data[i].Filename);
    }

    fprintf(Outputfile, "\n");

    // Print location of sample file
    fprintf(Outputfile, "Sample info-file: \n");
    fprintf(Outputfile, "%s \n", SampleFilename);
    fprintf(Outputfile, "\n");

    // Print fit quality to file
    fprintf(Outputfile, "Final chisquare = %g \n", ChiSquare);
    fprintf(Outputfile, "\n");

    // Print range of q to file
    fprintf(Outputfile, "Lower limit on q = %g \n", QMin);
    fprintf(Outputfile, "Upper limit on q = %g \n", QMax);

    fprintf(Outputfile, "\n");

    // Print parameters and properties of parameters to file
    fprintf(Outputfile, "Parameters: \n");
    fprintf(Outputfile, "Name:             Value:            Error: \n");

    for(i = 0; i < NumberOfParameters; ++i){

        if (Parameters[i].iParameter == true) {
            fprintf(Outputfile, "%-15s   %-15f   %-15g \n", Parameters[i].Name, Parameters[i].Value, Parameters[i].Error);
        } else {
            fprintf(Outputfile, "%-15s   %-15f   Fixed \n", Parameters[i].Name, Parameters[i].Value);
        }
    }

    fprintf(Outputfile, "\n");

    fprintf(Outputfile, "Location of inital .pdb-file: %s \n", ProteinStructure.PDBFileLocation);

    // Close file 
    fclose(Outputfile);
}
