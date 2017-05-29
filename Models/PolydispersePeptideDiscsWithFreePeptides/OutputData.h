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

    // Print global parameters to file
    for (i = 0; i < NumberOfSpectra; ++i) {
        fprintf(fp, "For spectrum %d: \n", i);
        fprintf(fp, "Excess scattering length density of headgroups = %g \n", Data[i].Constraints[0]);
        fprintf(fp, "Excess scattering length density of core       = %g \n", Data[i].Constraints[1]);
        fprintf(fp, "Excess scattering length density of methyl     = %g \n", Data[i].Constraints[2]);
        fprintf(fp, "Excess scattering length density of belt       = %g \n", Data[i].Constraints[3]);
        fprintf(fp, "\n");
    }

    // Print global parameters to file
    fprintf(fp, "Global parameteres: \n");
    fprintf(fp, "Height of bilayer             = %g \n", Data[0].Constraints[6]);
    fprintf(fp, "Height of hydrophobic bilayer = %g \n", Data[0].Constraints[7]);
    fprintf(fp, "Height of methyl layer        = %g \n", Data[0].Constraints[8]);
    fprintf(fp, "Mol. volume of dry lipid head = %g \n", Data[0].Constraints[13]);
    fprintf(fp, "Mol. volume of lipid tail     = %g \n", Data[0].Constraints[14]);
    fprintf(fp, "Mol. volume of methyl         = %g \n", Data[0].Constraints[15]);
    fprintf(fp, "Mol. volume of belt           = %g \n", Data[0].Constraints[16]);
    fprintf(fp, "\n");
    fprintf(fp, "Number density                   = %g \n", Data[0].Constraints[21]);
    fprintf(fp, "Mean number of peptides per disc = %g \n", Data[0].Constraints[17]);
    fprintf(fp, "Ratio of peptides in discs       = %g \n", Data[0].Constraints[18]);
    fprintf(fp, "Volume of core                   = %g \n", Data[0].Constraints[25]);
    fprintf(fp, "Volume of headgroups             = %g \n", Data[0].Constraints[26]);
    fprintf(fp, "Volume of methyl                 = %g \n", Data[0].Constraints[28]);
    fprintf(fp, "Height of peptide trimer         = %g \n", Data[0].Constraints[23]);

    // Close file and end program
    fclose(fp);
}
