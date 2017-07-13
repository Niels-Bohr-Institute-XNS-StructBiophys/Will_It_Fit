void OutputSpectra(struct Dataset * Data, int NumberOfSpectra, char *ResultsDirectory)
{
    // Dummy variables used in function
    int i;
    int j;
    bool PrintResolution = false;

    // Variables describing the files
    FILE *Outputfile;
    char Filename[256];


    // Begin outputting data
    for(i = 0; i < NumberOfSpectra; ++i){

        // Output dummy data
        sprintf(Filename, "%s/dat%d.dat",ResultsDirectory, i + 1);

        printf("%s\n", Filename);
        Outputfile = fopen(Filename, "w+");

        fprintf(Outputfile, "Datafile   = %s \n", Data[i].Filename);
        fprintf(Outputfile, "Datapoints = %d \n", Data[i].NumberOfDatapoints);
        fprintf(Outputfile, "Contrast   = %g \n", Data[i].Contrast);

        for(j = 0; j < Data[i].NumberOfDatapoints; ++j){
            fprintf(Outputfile, " %15g   %15g   %15g   1.0 \n", Data[i].QValues[j], Data[i].IValues[j], 1.0 / sqrt(Data[i].SigmaValues[j]));
        }

        fclose(Outputfile);

        // Output dummy fitting data
        sprintf(Filename, "%s/fit%ds.dat",ResultsDirectory, i + 1);

        Outputfile = fopen(Filename, "w+");

        fprintf(Outputfile, "Fit of data from file   = %s \n", Data[i].Filename);
        fprintf(Outputfile, "Number of points in fit = %d \n", Data[i].NMax - Data[i].NMin);
        fprintf(Outputfile, "Contrast                = %g \n", Data[i].Contrast);

        for(j = Data[i].NMin; j < Data[i].NMax; ++j){
            fprintf(Outputfile, " %15g   %15g          1.0          1.0 \n", Data[i].QValues[j], Data[i].FitValues[j]);
        }

        fclose(Outputfile);

        if (PrintResolution == true) {

            if (Data[i].IncludeResolutionEffects == true) {
                sprintf(Filename, "Resolution%d.mcp", i + 1);

                Outputfile = fopen(Filename, "w+");

                for(j = Data[i].NMin; j < Data[i].NMax; ++j){
                    fprintf(Outputfile, " %15g   %15g \n", Data[i].QValues[j], Data[i].SigmaQValues[j] / Data[i].QValues[j]);
                }

                fclose(Outputfile);
            }
        }
    }
}
