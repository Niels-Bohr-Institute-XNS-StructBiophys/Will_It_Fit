/// This function is used to clear the terminal window in a non-OS-dependent way
void ClearScreen()
{
    printf("******************************************************************************** \n");
}

/// Function used to communicate a return message to interface
void ReturnMessage(char Message[])
{
    FILE *fp;

    fp = fopen(".data/ReturnMessage.mcp", "w");

    fprintf(fp, "%s", Message);

    fclose(fp);
}

/// Errorcatching function
void Errorcheck(int Errorcode, char Message[]) {
    char ErrorMessage[256];

    if (Errorcode == -1) {
        sprintf(ErrorMessage, "An error occured when %s...", Message);
        ReturnMessage(ErrorMessage);

        exit(0);
    }
}

/// Function used to assign fitting ranges
void AssignFittingRanges(struct Dataset * Data, double QMin, double QMax, int NumberOfSpectra)
{
    // Declarations
    int i;
    int j;

    int NMinDummy;
    int NMaxDummy;

    int TotalNumberOfDatapoints = 0;

    // Computation
    for (i = 0; i < NumberOfSpectra; ++i) {
        NMinDummy = 0;
        NMaxDummy = Data[i].NumberOfDatapoints;

        for (j = 0; j < Data[i].NumberOfDatapoints; ++j) {

            if (Data[i].QValues[j] < QMin) {
                NMinDummy = j + 1;
            }
        }

        for (j = Data[i].NumberOfDatapoints - 1; j >= 0; --j) {

            if (Data[i].QValues[j] > QMax){
                NMaxDummy = j;
            }
        }

        TotalNumberOfDatapoints += NMaxDummy - NMinDummy;

        printf("Dataset nr. %d contains %d points. \n", i + 1, Data[i].NumberOfDatapoints);
        printf("Fit will be from point %d to point %d. \n", NMinDummy + 1, NMaxDummy);
        printf("\n");

        Data[i].NMin = NMinDummy;
        Data[i].NMax = NMaxDummy;
    }

    printf("The total number of fitted points is %d. \n", TotalNumberOfDatapoints);

    return;
}

/// Function used to assign arguments from the console to the appropriate variables
void AssignArgumentsBash(int NumberOfArguments, char *Arguments[], char CardFileLocation[256], char SamplesFileLocation[256], char ParameterFileLocation[256],
                     double *QMin, double *QMax,
                     bool *IncludeResolutionEffects, int *NumberOfSmearingFolds, char ResolutionFileLocation[256],
                     char PDBFileLocation[256], int *ComputeModel, bool *CMD)
{
    // Declarations
    int i;

    char Category;
    char Argument[256];

    // Main loop
    for (i = 1; i < NumberOfArguments; ++i) {

        // Scan and assign
        sscanf(Arguments[i], "-%c=%s", &Category, Argument);

        switch (Category) {

            // .card-file location
            case 'c':
                sprintf(CardFileLocation, "%s", Argument);
            break;

            // Samples file location
            case 's':
                sprintf(SamplesFileLocation, "%s", Argument);
            break;

            // Parameter file location
            case 'p':
                sprintf(ParameterFileLocation, "%s", Argument);
            break;

            // Fitting range
            case 'n':
                *QMin = atof(Argument);
            break;

            case 'x':
                *QMax = atof(Argument);
            break;

            // Obtain resolution information
            case 't':
                if (atoi(Argument) != 0) {
                    *IncludeResolutionEffects = true;
                    *NumberOfSmearingFolds = atoi(Argument);
                }
            break;

            case 'e':
                sprintf(ResolutionFileLocation, "%s", Argument);
            break;

            // .pdb-file location
            case 'd':
                sprintf(PDBFileLocation, "%s", Argument);
            break;
	       // Will we fit or compute model? 
	        case 'f':
               *ComputeModel = atoi(Argument);
            break;
            case 'z':
                *CMD = atoi;
                printf("CMD is %d \n", *CMD);

            // Unknown argument
            default:
                printf("Unknown arguments passed to algorithm! \n");
            break;
        }
    }
}
void AssignArguments(int NumberOfArguments, char *Arguments[], char CardFileLocation[256], char SamplesFileLocation[256], char ParameterFileLocation[256],
                     double *QMin, double *QMax, int *ChooseFittingRoutine, int *FittingRoutineArgument2,
                     bool *IncludeResolutionEffects, int *NumberOfSmearingFolds, char ResolutionFileLocation[256], bool *PrintCorrelationMatrix,
                     char PDBFileLocation[256], double *ChiSquareFractile, int *FittingRoutineArgument3, bool * CMD)
{
    // Declarations
    int i;

    char Category;
    char Argument[256];

    // Main loop
    for (i = 1; i < NumberOfArguments; ++i) {

        // Scan and assign
        sscanf(Arguments[i], "-%c=%s", &Category, Argument);

        switch (Category) {

            // .card-file location
            case 'c':
                sprintf(CardFileLocation, "%s", Argument);
            break;

            // Samples file location
            case 's':
                sprintf(SamplesFileLocation, "%s", Argument);
            break;

            // Parameter file location
            case 'p':
                sprintf(ParameterFileLocation, "%s", Argument);
            break;

            // Fitting range
            case 'n':
                *QMin = atof(Argument);
            break;

            case 'x':
                *QMax = atof(Argument);
            break;

            // Choose fitting routine
            case 'r':
                *ChooseFittingRoutine = atoi(Argument);
            break;

            case 'y':
                *FittingRoutineArgument2 = atoi(Argument);
            break;

            // Obtain resolution information
            case 't':
                if (atoi(Argument) != 0) {
                    *IncludeResolutionEffects = true;
                    *NumberOfSmearingFolds = atoi(Argument);
                }
            break;

            case 'e':
                sprintf(ResolutionFileLocation, "%s", Argument);
            break;

            // Print correlation matrix
            case 'o':
                if (atoi(Argument) != 0) {
                    *PrintCorrelationMatrix = true;
                }
            break;

            // .pdb-file location
            case 'd':
                sprintf(PDBFileLocation, "%s", Argument);
            break;

            // Profile likelihood fractile
            case 'h':
                *ChiSquareFractile = atof(Argument);
            break;

            // Number of agents in swarm
            case 'a':
                *FittingRoutineArgument3 = atoi(Argument);
            break;

            case 'z':
                *CMD = atoi(Argument);
                printf("CMD is %d\n", *CMD);
            break;
            // Unknown argument
            default:
                printf("Unknown arguments passed to algorithm! \n");
            break;
        }
    }
}

void DirectoryFinder(char FileLocation[256], char Directory[256]){
    int length;
    int i;
    int j;
    char curr;
    char check = '/';
    //loop backwards through FileLocation string to find the end of the directory
    length = strlen(FileLocation);
    for (i= length -1; i >= 0; --i)
    {
        curr = FileLocation[i];
        if(curr == check)
        {
            break;
        }
    }
    //Second loop assigns characters in the directory string
    for(j = 0; j <= i; ++j)
    {
        Directory[j] = FileLocation[j];
    }


}
