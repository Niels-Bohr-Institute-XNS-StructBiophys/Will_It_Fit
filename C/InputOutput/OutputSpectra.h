int OutputSpectra(struct Dataset * Data, int NumberOfSpectra, char *ResultsDirectory, int WriteLog, FILE* logfile)
{
	// Dummy variables used in function
	bool PrintResolution = false ;

	// Variables describing the files
	FILE *Outputfile ;
	char Filename[256] ;


	// Begin outputting data
	for( int i = 0; i < NumberOfSpectra; ++i)
	{
		// Output dummy data
		sprintf( Filename, "%sdat%d.dat", ResultsDirectory, i+1) ;
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\tWrite spectra no. %d to file %s\n", i+1, Filename) ;

			fflush( logfile ) ;
		}
		Outputfile = fopen( Filename, "w") ;
		if ( Outputfile == NULL)
		{
			if ( abs(WriteLog) > 0 )
			{
				fprintf( logfile, "\t\tAn error occured when attempting to open the file\n") ;
				fprintf( logfile, "\n") ;

				fflush( logfile ) ;
			}
			return -1 ;
		}

		fprintf( Outputfile, "# Data file             = %s\n", Data[i].Filename) ;
		fprintf( Outputfile, "# Number of data points = %d\n", Data[i].NumberOfDatapoints) ;
		fprintf( Outputfile, "# Contrast              = %g\n", Data[i].Contrast) ;
		for( int j = 0; j < Data[i].NumberOfDatapoints; ++j)
		{
			fprintf( Outputfile, " %15g   %15g   %15g   1.0 \n", Data[i].QValues[j], Data[i].IValues[j], 1.0 / sqrt(Data[i].SigmaValues[j])) ;
		}

		fclose( Outputfile) ;


		// Output dummy fitting data
		sprintf( Filename, "%sfit%d.dat", ResultsDirectory, i+1) ;
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\tWrite fitted spectra no. %d to file %s\n", i+1, Filename) ;

			fflush( logfile ) ;
		}
		Outputfile = fopen( Filename, "w") ;
		if ( Outputfile == NULL)
		{
			if ( abs(WriteLog) > 0 )
			{
				fprintf( logfile, "\t\tAn error occured when attempting to open the file\n") ;
				fprintf( logfile, "\n") ;

				fflush( logfile ) ;
			}
			return -1 ;
		}

		fprintf( Outputfile, "# Fit of data from file   = %s\n", Data[i].Filename) ;
		fprintf( Outputfile, "# Number of points in fit = %d\n", Data[i].NMax - Data[i].NMin) ;
		fprintf( Outputfile, "# Contrast                = %g\n", Data[i].Contrast) ;

		for( int j = Data[i].NMin; j < Data[i].NMax; ++j)
		{
			fprintf( Outputfile, "%15g   %15g          1.0          1.0\n", Data[i].QValues[j], Data[i].FitValues[j]) ;
		}

		fclose( Outputfile) ;


		// Output dummy resolution files with Q vs calculated dQ/Q
		if ( PrintResolution == true )
		{
			if ( Data[i].IncludeResolutionEffects == true )
			{
				sprintf( Filename, "%sres%d.dat", ResultsDirectory, i+1) ;
				if ( abs(WriteLog) > 0 )
				{
					fprintf( logfile, "\t\tWrite Q-resolution for spectra no. %d to file %s\n", i+1, Filename) ;

					fflush( logfile ) ;
				}
				Outputfile = fopen( Filename, "w") ;
				if ( Outputfile == NULL)
				{
						if ( abs(WriteLog) > 0 )
						{
							fprintf( logfile, "\t\tAn error occured when attempting to open the file\n") ;
							fprintf( logfile, "\n") ;

							fflush( logfile ) ;
						}
						return -1 ;
				}

				for( int j = Data[i].NMin; j < Data[i].NMax; ++j)
				{
					fprintf( Outputfile, "%15g   %15g\n", Data[i].QValues[j], Data[i].SigmaQValues[j] / Data[i].QValues[j]) ;
				}

				fclose( Outputfile) ;
			}
		}

		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\n") ; }
	}

	return 0 ;
}
