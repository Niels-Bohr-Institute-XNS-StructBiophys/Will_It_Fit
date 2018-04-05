inline int max ( int a, int b ) { return a > b ? a : b ; }
inline int min ( int a, int b ) { return a < b ? a : b ; }

void ClearScreen(FILE *logfile)
{
	fprintf( logfile, "********************************************************************************\n") ;
}


/// Write return message to ReturnMessage.mcp in the .card-file folder
void ReturnMessage(char Message[], char *ResultsDirectory)
{
	FILE *fp ;
	char ReturnMessageFileLocation[256] ;

	sprintf( ReturnMessageFileLocation, "%s%s", ResultsDirectory, "ReturnMessage.mcp") ;
	fp = fopen( ReturnMessageFileLocation, "w") ;

	fprintf(fp, "%s", Message) ;

	fclose(fp) ;
}

/// Errorcatching function
// writes ErrorMessage to ReturnMessage.mcp in ResultsDirectory and if applicable to logfile / stdout
void ErrorCheck(int Errorcode, char const * Message, char *ResultsDirectory, int WriteLog, FILE* logfile)
{
	char ErrorMessage[256] ;

	if (Errorcode == -1)
	{
		sprintf( ErrorMessage, "An error occured %s...\n", Message) ;

		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "%s", ErrorMessage) ; }

		if (ResultsDirectory != NULL ) { ReturnMessage( ErrorMessage, ResultsDirectory) ; }

		exit(0) ;
	}
}

/// Function used to assign fitting ranges
void AssignFittingRanges(struct Dataset * Data, double QMin, double QMax, int NumberOfSpectra, int WriteLog, FILE* logfile)
{
	// Declarations
	int i ;
	int j ;

	int NMinDummy ;
	int NMaxDummy ;

	int TotalNumberOfDatapoints = 0 ;

	// Computation
	for (i = 0; i < NumberOfSpectra; ++i)
	{
		NMinDummy = 0 ;
		NMaxDummy = Data[i].NumberOfDatapoints ;

		for (j = 0; j < Data[i].NumberOfDatapoints; ++j)
		{
			if (Data[i].QValues[j] < QMin) { NMinDummy = j + 1 ; }
		}

		for (j = Data[i].NumberOfDatapoints - 1; j >= 0; --j)
		{
			if (Data[i].QValues[j] > QMax) { NMaxDummy = j ; }
		}

		TotalNumberOfDatapoints += NMaxDummy - NMinDummy ;

		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\tDataset no. %d contains %d points.\n", i + 1, Data[i].NumberOfDatapoints) ;
			fprintf( logfile, "\t\tFit will be from point %d to point %d.\n", NMinDummy + 1, NMaxDummy) ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}

		Data[i].NMin = NMinDummy ;
		Data[i].NMax = NMaxDummy ;
	}

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\t\tThe total number of fitted points is %d.\n", TotalNumberOfDatapoints) ;
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;
	}

	return ;
}


int AssignArguments(int NumberOfArguments, char *Arguments[], char CardFileLocation[256], char SamplesFileLocation[256], char ParameterFileLocation[256], char PDBFileLocation[256], double *QMin, double *QMax, int *ChooseFittingRoutine, int *FittingRoutineArgument1, int *FittingRoutineArgument2, int *FittingRoutineArgument3, bool *IncludeResolutionEffects, int *NumberOfSmearingFolds, char ResolutionFileLocation[256], bool *PrintCovarianceMatrix, double *ChiSquareFractile, bool *CMD, int *WriteLog)
{
	// Declarations
	char Category ;
	char Argument[256] ;

	// Main loop
	for ( int i = 1; i < NumberOfArguments; ++i)
	{
		// Scan and assign
		sscanf( Arguments[i], "-%c=%s", &Category, Argument) ;

		switch ( Category ) {

			// .card-file location
			case 'c':
				sprintf(CardFileLocation, "%s", Argument) ;
			break ;

			// .dat samples file location
			case 's':
				sprintf(SamplesFileLocation, "%s", Argument) ;
			break ;

			// .par parameter file location
			case 'p':
				sprintf(ParameterFileLocation, "%s", Argument) ;
			break ;

			// .pdb file location
			case 'd':
				sprintf(PDBFileLocation, "%s", Argument) ;
			break ;


			// Fitting range
			case 'n':
				*QMin = atof(Argument) ;
			break ;

			case 'x':
				*QMax = atof(Argument) ;
			break ;


			// Choose fitting routine and its optional arguments
			case 'r':
				*ChooseFittingRoutine = atoi(Argument) ;
			break ;

			case 'i':
				*FittingRoutineArgument1 = atoi(Argument) ;
			break ;

			case 'j':
				*FittingRoutineArgument2 = atoi(Argument) ;
			break ;

			case 'k':
				*FittingRoutineArgument3 = atoi(Argument) ;
			break ;


			// Obtain resolution information
			case 't':
				if (atoi(Argument) != 0)
				{
					*IncludeResolutionEffects = true ;
					*NumberOfSmearingFolds = atoi(Argument) ;
				}
			break ;

			case 'e':
				sprintf(ResolutionFileLocation, "%s", Argument) ;
			break ;


			// Print covariance matrix
			case 'o':
				if (atoi(Argument) != 0) { *PrintCovarianceMatrix = true ; }
			break ;

			// Profile likelihood fractile
			case 'h':
				*ChiSquareFractile = atof(Argument) ;
			break ;


			// CMD mode
			case 'z':
				*CMD = atoi(Argument) ;
			break ;


			// write logfile
			case 'l':
				*WriteLog = atoi(Argument) ;
			break ;


			// Unknown argument
			default:
				printf("\n") ;
				printf("Unknown argument %s passed to algorithm\n", Arguments[i]) ;
				printf("\n") ;
				return -1 ;
			break ;
		}
	}
	return 0 ;
}

bool DirectoryFinder(char* FileLocation, char* Directory)
{
	int length ;
	int i ;
	int j ;
	char curr ;

	// potential OS problem with folder binders / (*nix, Mac/BSD) vs \ (Windows)
	char check = '/' ;

	bool found = false ;
	//loop backwards through FileLocation string to find the end of the directory
	length = strlen(FileLocation) ;
	for (i= length-2; i >= 0; --i)
	{
		curr = FileLocation[i] ;
		if(curr == check)
		{
			found = true ;
			break ;
		}
	}
	if ( found )
	{
		//Second loop assigns characters in the directory string
		for(j = 0; j <= i; ++j)
		{
			Directory[j] = FileLocation[j] ;
		}
		Directory[i+1] = 0 ;
	}

	return found ;
}
