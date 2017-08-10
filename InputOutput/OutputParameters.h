void OutputParameters(struct Parameter * Parameters, int NumberOfParameters, int ChooseFittingRoutine, double ChiSquare)
{
    /// Declarations
    int i;
    int j;

    FILE *fp;

    /// Print to file
    //Create the output file
    fp = fopen("ParametersAfter.mcp", "w+");

    // Print parameters and properties of parameters to file
    for(i = 0; i < NumberOfParameters; ++i){

        if (Parameters[i].iParameter == true) {
            j = 0;
        } else {
            j = 1;
        }

        fprintf(fp, "%-15g   %-15g   %-15g   %d   %-15s \n", Parameters[i].MinValue, Parameters[i].Value, Parameters[i].MaxValue, j, Parameters[i].Name);
    }

    // Close file
    fclose(fp);

    /// Print to screen
    printf("Parameters: \n");
    printf("Final Chisquare = %g \n", ChiSquare);
    printf("\n");

    printf("                Value          Error          Name \n");

    for(i = 0; i < NumberOfParameters; ++i){

        if (Parameters[i].iParameter == true) {
            printf("Parameter %2d:   %-12g   %-12g   %s \n", i, Parameters[i].Value, Parameters[i].Error, Parameters[i].Name);
        } else {
            printf("Parameter %2d:   %-12g   Fixed          %s \n", i, Parameters[i].Value, Parameters[i].Name);
        }
    }

    printf("\n");
}

