int OutputParameters(struct Parameter * Parameters, int NumberOfParameters, int ChooseFittingRoutine, double ChiSquare, char * ResultsDirectory, int WriteLog, FILE* logfile)
{
	/// Declarations
	int j ;

	char Filename[256] ;
	FILE *Outputfile ;

	/// Print to file
	//Create the output file
	sprintf( Filename, "%sParametersAfter.mcp", ResultsDirectory) ;

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\t\tWrite fitted parameters to parameter file at %s\n", Filename) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}
	Outputfile = fopen(Filename, "w") ;
	if ( Outputfile == NULL)
	{
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\tAn error occured when attempting to open the parameter file\n") ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}
		return -1 ;
	}

	// Print parameters and properties of parameters to file
	for( int i = 0; i < NumberOfParameters; ++i)
	{
		if (Parameters[i].iParameter == true) { j = 0 ; }
		else { j = 1 ; }

		fprintf( Outputfile, "%-15g   %-15g   %-15g   %d   %-15s \n", Parameters[i].MinValue, Parameters[i].Value, Parameters[i].MaxValue, j, Parameters[i].Name) ;
	}

	// Close file
	fclose( Outputfile ) ;


	/// Print to logfile
	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\t\tSummarize fitted parameters\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\t\t\tFinal Chisquare = %g\n", ChiSquare) ;
		printf("\n") ;
		fprintf( logfile, "\t\t\tParameters:\n") ;
		printf("\n") ;
		fprintf( logfile, "\t\t\t     Value             Error             Name\n") ;

		for( int i = 0; i < NumberOfParameters; ++i)
		{
			if ( Parameters[i].iParameter == true ) { fprintf( logfile, "\t\t\t%2d   %-15g   %-15g   %s\n", i, Parameters[i].Value, Parameters[i].Error, Parameters[i].Name) ; }
			else { fprintf( logfile, "\t\t\t%2d   %-15g   Fixed             %s\n", i, Parameters[i].Value, Parameters[i].Name) ; }
		}

		fflush( logfile ) ;
	}

	return 0 ;
}

