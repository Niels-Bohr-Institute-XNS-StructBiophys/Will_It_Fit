/// Check size of sample file
int CheckSizeOfSampleInformation(char SamplesFileLocation[256]) {
    // Declarations
    FILE *fp;
    char Line[1024];

    // Open file
    fp = fopen(SamplesFileLocation, "r");

    if (fp == NULL){
        return -1;
    }

    // Fetch values from file
    if(fgets(Line, sizeof(Line), fp));

    fclose(fp);

    return atoi(Line);
}

/// Import sample information
void ImportSampleInformation(struct Dataset * Data, double * VolumesOfMolecules, char SamplesFileLocation[256], int NumberOfSampleInformations, int NumberOfSpectra)
{
    /// Declarations
    // Dummy variables used in function
    int i;
    int j;
    char Line[1024];
    int LineIterator = 0;

    // Dummy variables used to store inputs
    double ScatteringLengthsNeutrons100[NumberOfSampleInformations];
    double ScatteringLengthsNeutrons0[NumberOfSampleInformations];
    double ScatteringLengthsXrays[NumberOfSampleInformations];
    char Dummy[2 + 4 * (NumberOfSampleInformations + 2)][128];

    /// I/O
    // Variables describing the file
    FILE *fp;

    // Open file
    printf("Loading samples file. \n");

    fp = fopen(SamplesFileLocation, "r");

    if (fp == NULL){
        printf("An error occured when attempting to open samples file.");
        return;
    }

    printf("File opened. \n");

    // Fetch values from file
    while (fgets(Line, sizeof(Line), fp) != NULL) {
        sscanf(Line, "%s ", Dummy[LineIterator]);
        ++LineIterator;
    }

    for (i = 0; i < NumberOfSampleInformations; ++i) {
        VolumesOfMolecules[i]           = atof(Dummy[3 + i]);
        ScatteringLengthsNeutrons100[i] = atof(Dummy[5 + NumberOfSampleInformations + i]);
        ScatteringLengthsNeutrons0[i]   = atof(Dummy[7 + 2 * NumberOfSampleInformations + i]);
        ScatteringLengthsXrays[i]       = atof(Dummy[9 + 3 * NumberOfSampleInformations + i]);
    }

    // Close file and dump values into ScatteringLengths and VolumesOfMolecules
    fclose(fp);

    /// Compute scattering lengths for the various contrasts
    for (i = 0; i < NumberOfSpectra; ++i) {

        if (Data[i].Contrast >= 0.0 && Data[i].Contrast <= 100.0) {

            for (j = 0; j < NumberOfSampleInformations; ++j) {
                Data[i].ScatteringLengths[j] = ScatteringLengthsNeutrons100[j] * Data[i].Contrast / 100.0 +
                                               ScatteringLengthsNeutrons0[j]   * (100.0 - Data[i].Contrast) / 100.0;
            }

        } else {

            for (j = 0; j < NumberOfSampleInformations; ++j) {
                Data[i].ScatteringLengths[j] = ScatteringLengthsXrays[j];
            }
        }
    }
}
