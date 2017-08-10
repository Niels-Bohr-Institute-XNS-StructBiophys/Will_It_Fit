int CheckNumberOfParameters(char ParameterFileLocation[256])
{
    // Declare variables used in function
    char Line[128];
    int LineIterator;

    // Dummy variables used to read file
    double Dummy1;
    double Dummy2;
    double Dummy3;
    int Dummy4;
    char Dummy5[sizeof(Line)];

    // Specify filename with parameters:
    FILE *ParameterFile;

    // Open file and read parameters
    printf("Loading parameter file. \n");
    printf("%s\n", ParameterFileLocation);
    ParameterFile = fopen(ParameterFileLocation, "r");

    if (ParameterFile == NULL){
        printf("An error occured when attempting to open parameter file.");
        return -1;
    }

    printf("File opened. \n");

    // Fetch parameters from file (as strings)
    LineIterator = 0;

    while (fgets(Line, sizeof(Line), ParameterFile) != NULL) {

        if (sscanf(Line, "%lf %lf %lf %d %s", &Dummy1, &Dummy2, &Dummy3, &Dummy4, Dummy5) == 5) {
            ++LineIterator;
        }
    }

    if (LineIterator == 0) {
        return -1;
    }

    // Close and return
    fclose(ParameterFile);

    return LineIterator;
}

void ImportParameters(struct Parameter * Parameters, char ParameterFileLocation[256])
{
    // Declare variables used in function
    char Line[128];
    int LineIterator;

    // Dummy variables used to read file
    double Dummy1;
    double Dummy2;
    double Dummy3;
    int Dummy4;
    char Dummy5[sizeof(Line)];

    // Specify filename
    FILE *ParameterFile;

    // Fetch parameters from file (as strings)
    ParameterFile = fopen(ParameterFileLocation, "r");
    LineIterator = 0;

    while (fgets(Line, sizeof(Line), ParameterFile) != NULL) {

        if (sscanf(Line, "%lf %lf %lf %d %s", &Dummy1, &Dummy2, &Dummy3, &Dummy4, Dummy5) == 5) {
            Parameters[LineIterator].MinValue = Dummy1;
            Parameters[LineIterator].Value    = Dummy2;
            Parameters[LineIterator].MaxValue = Dummy3;
            Parameters[LineIterator].Error    = 0.0;

            if (Dummy4 == 0) {
                Parameters[LineIterator].iParameter = true;
            } else {
                Parameters[LineIterator].iParameter = false;
            }

            sprintf(Parameters[LineIterator].Name, "%s", Dummy5);

            ++LineIterator;
        }
    }

    // Close and return
    fclose(ParameterFile);
}
