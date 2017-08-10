int CheckSizeOfData(char CardFileLocation[256], int *NumberOfSpectra)
{
    /// Declarations
    double Dummy1;
    double Dummy2;
    double Dummy3;
    double Dummy4;

    int i;
    int NumberOfDatasetsDummy;
    int PointsInCurrentSpectrum;
    int MostPointsInASpectrum = 0;

    char Buffer[128];
    int linenum;

    // Variables describing the files
    FILE *fp;
    char Path[sizeof(Buffer)];
    char Dummy[2 + 10 * 6][sizeof(Buffer)]; // Room for 10 datasets

    /// Open from .card-file
    fp = fopen(CardFileLocation, "r");

    if (fp == NULL) {
        printf("An error occured when attempting to open the .card-file. \n");
        return -1;
    }

    printf(".card-file found!\n");

    // Write the data and names into dummy variable
    linenum = 0;

    while (fgets(Buffer, sizeof(Buffer), fp) != NULL) {
        sscanf(Buffer, "%s", Dummy[linenum]);
        ++linenum;
    }

    // Close .card-file
    fclose(fp);

    // Set number of datasets
    NumberOfDatasetsDummy = atoi(Dummy[0]);

    /// Begin reading from .rad-files
    for (i = 0; i < NumberOfDatasetsDummy; ++i) {
        sprintf(Path, "%s", Dummy[2 + i * 6]);
        fp = fopen(Path,"r");

        if(fp == NULL){
            printf("An error occured when attempting to open the .rad-files. \n");
            return -1;
        }

        printf("Datafile found! - %s \n", Path);

        // Read lines into dummy variables
        linenum = 0;

        while (fgets(Buffer, sizeof(Buffer), fp) != NULL) {

            if (sscanf(Buffer, "%lf  %lf  %lf  %lf", &Dummy1, &Dummy2, &Dummy3, &Dummy4) == 4) {
                ++linenum;
            } else if (sscanf(Buffer, "%lf  %lf  %lf", &Dummy1, &Dummy2, &Dummy3) == 3) {
                ++linenum;
            }
        }

        PointsInCurrentSpectrum = linenum;

        if (PointsInCurrentSpectrum > MostPointsInASpectrum) {
            MostPointsInASpectrum = PointsInCurrentSpectrum;
        }

        // Close .rad-file
        fclose(fp);
    }

    // Return highest number of datapoints in file and number of spectra
    *NumberOfSpectra = NumberOfDatasetsDummy;
    return MostPointsInASpectrum;
}

void ImportSpectra(struct Dataset * Data, char CardFileLocation[256], int NumberOfSpectra)
{
    // Variables used in iterations
    int i;

    double Dummy1;
    double Dummy2;
    double Dummy3;
    double Dummy4;

    char Buffer[128];
    int linenum;

    // Variables describing the files
    FILE *fp;
    char Path[sizeof(Buffer)];
    char Dummy[MaxSizeOfCardfile][sizeof(Buffer)];

    // Variables used to set the scaling, concentration, contrast, and background for each spectrum
    double Contrast;
    double Concentration;
    double Scaling;
    double Background;

    /// Open from .card-file
    fp = fopen(CardFileLocation, "r");

    if (fp == NULL) {
        printf("An error occured when attempting to open the .card-file. \n");
        return;
    }

    // Write the data and names into dummy variable
    linenum = 0;

    while (fgets(Buffer, sizeof(Buffer), fp) != NULL) {
        sscanf(Buffer, "%s", Dummy[linenum]);
        ++linenum;
    }

    // Close .card-file
    fclose(fp);

    /// Begin reading from .rad-files
    for (i = 0; i < NumberOfSpectra; ++i) {
        // Prepare .rad-file
        sprintf(Path, "%s", Dummy[2 + i * 6]);
        fp = fopen(Path,"r");

        Contrast      = atof(Dummy[3 + i * 6]);
        Concentration = atof(Dummy[4 + i * 6]);
        Scaling       = atof(Dummy[5 + i * 6]);
        Background    = atof(Dummy[6 + i * 6]);

        if(fp == NULL){
            printf("An error occured when attempting to open the .rad-files. \n");
            return;
        }

        // Read lines into dummy variables
        linenum = 0;

        while (fgets(Buffer, sizeof(Buffer), fp) != NULL) {

            if (sscanf(Buffer, "%lf  %lf  %lf  %lf", &Dummy1, &Dummy2, &Dummy3, &Dummy4) == 4) {
                Data[i].QValues[linenum]     = Dummy1;
                Data[i].IValues[linenum]     = Scaling * Dummy2 - Background;
                Data[i].SigmaValues[linenum] = 1.0 / pow(Scaling * Dummy3, 2);

                ++linenum;
            } else if (sscanf(Buffer, "%lf  %lf  %lf", &Dummy1, &Dummy2, &Dummy3) == 3) {
                Data[i].QValues[linenum]     = Dummy1;
                Data[i].IValues[linenum]     = Scaling * Dummy2 - Background;
                Data[i].SigmaValues[linenum] = 1.0 / pow(Scaling * Dummy3, 2);

                ++linenum;
            }
        }

        // Close .rad-file
        fclose(fp);

        Data[i].Concentration      = Concentration * 6.02 * 1e17;
        Data[i].Contrast           = Contrast;
        Data[i].NumberOfDatapoints = linenum;

        sprintf(Data[i].Filename, "%s", Path);
    }
}
