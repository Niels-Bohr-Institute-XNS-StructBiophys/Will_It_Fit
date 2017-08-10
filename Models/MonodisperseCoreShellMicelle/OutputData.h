void OutputData(double ChiSquare, double ChiSquareRed, double QMin, double QMax, struct Parameter * Parameters, int NumberOfParameters,
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
    fprintf(fp, "Final chisqred = %g \n", ChiSquareRed);

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

    // Print constraints to file
    for (i = 0; i < NumberOfSpectra; ++i) {
        fprintf(fp, "For spectrum %d: \n", i);
        fprintf(fp, "Scattering length density of heads = %g \n", Data[i].Constraints[0]);
        fprintf(fp, "Scattering length density of tails = %g \n", Data[i].Constraints[1]);
        fprintf(fp, "\n");
    }

    // Print global parameters to file
    fprintf(fp, "Global parameteres: \n");
    fprintf(fp, "Total radius                      = %g \n", Data[0].Constraints[2]);
    fprintf(fp, "Radius of core                    = %g \n", Data[0].Constraints[3]);
    fprintf(fp, "Thickness of shell                = %g \n", Data[0].Constraints[2] - Data[0].Constraints[3]);
    fprintf(fp, "Number of molecules per aggregate = %g \n", Data[0].Constraints[5]);
    fprintf(fp, "Volume of head                    = %g \n", Data[0].Constraints[6]);
    fprintf(fp, "Volume of tail                    = %g \n", Data[0].Constraints[7]);
    fprintf(fp, "\n");

    // Close file and end program
    fclose(fp);
}
