int CheckSizeOfData(char* CardFileLocation, int *NumberOfSpectra, bool CMD, int WriteLog, FILE* logfile)
{
	/// Declarations
	double Dummy1 ;
	double Dummy2 ;
	double Dummy3 ;
	double Dummy4 ;

	int i ;
	int NumberOfDatasetsDummy ;
	int PointsInCurrentSpectrum ;
	int MostPointsInASpectrum = 0 ;

	//Is the .rad files located in a subdirectory
	bool SubDir = false ; 

	char Buffer[256] ;
	char Directory[256] ;
	int linenum ;

	// Variables describing the files
	FILE *fp ;

	char Path[sizeof(Buffer)] ;

	char Dummy[2 + 10 * 6][sizeof(Buffer)] ; // Room for 10 datasets

	/// Open from .card-file
	fp = fopen(CardFileLocation, "r") ;
	if (fp == NULL)
	{
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\tAn error occured when attempting to open the .card-file\n") ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}
		return -1 ;
	}

	if (CMD) { SubDir = DirectoryFinder(CardFileLocation, Directory)  ; }


	// Write the data and names into dummy variable
	linenum = 0 ;

	while (fgets(Buffer, sizeof(Buffer), fp) != NULL)
	{
		sscanf(Buffer, "%s", Dummy[linenum]) ;
		++linenum ;
	}

	// Close .card-file
	fclose(fp) ;

	// Set number of datasets
	NumberOfDatasetsDummy = atoi(Dummy[0]) ;
	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\t\tAttempt to read %d spectra\n", NumberOfDatasetsDummy) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}


	/// Begin reading from .rad-files
	for (i = 0; i < NumberOfDatasetsDummy; ++i)
	{
		if (CMD && SubDir ) { sprintf(Path, "%s%s", Directory, Dummy[2 + i * 6]) ; }
		else{ sprintf(Path, "%s", Dummy[2 + i * 6]); }

		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\t\tRead file no. %d from %s\n", i+1, Path) ;

			fflush( logfile ) ;
		}

		fp = fopen(Path,"r") ;
		if( fp == NULL )
		{
			if ( abs(WriteLog) > 0 )
			{
				fprintf( logfile, "\t\t\tAn error occured when attempting to open the file\n") ;
				fprintf( logfile, "\n") ;

				fflush( logfile ) ;
			}
			return -1  ;
		}

		// Read lines into dummy variables
		linenum = 0 ;
		while (fgets(Buffer, sizeof(Buffer), fp) != NULL)
		{
			if (sscanf(Buffer, "%lf  %lf  %lf  %lf", &Dummy1, &Dummy2, &Dummy3, &Dummy4) == 4) { ++linenum ; }
			else if (sscanf(Buffer, "%lf  %lf  %lf", &Dummy1, &Dummy2, &Dummy3) == 3) { ++linenum ; }
		}

		// Close .rad-file
		fclose(fp) ;

		PointsInCurrentSpectrum = linenum ;
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\t\tFound %d data points\n", PointsInCurrentSpectrum) ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}

		if (PointsInCurrentSpectrum > MostPointsInASpectrum) { MostPointsInASpectrum = PointsInCurrentSpectrum ; }
	}

	// Return highest number of datapoints in file and number of spectra
	*NumberOfSpectra = NumberOfDatasetsDummy ;

	return MostPointsInASpectrum ;
}

int ImportSpectra(struct Dataset* Data, char* CardFileLocation, int NumberOfSpectra, bool CMD, int WriteLog, FILE* logfile)
{
	// Variables used in iterations
	int i ;

	double Dummy1 ;
	double Dummy2 ;
	double Dummy3 ;
	double Dummy4 ;
	//Is the .rad files located in a subdirectory
	bool SubDir = false ; 

	char Buffer[256] ;
	char Directory[256] ;
	int linenum ;

	// Variables describing the files
	FILE *fp ;
	char Path[sizeof(Buffer)] ;
	char Dummy[MaxSizeOfCardfile][sizeof(Buffer)] ;

	// Variables used to set the scaling, concentration, contrast, and background for each spectrum
	double Contrast ;
	double Concentration ;
	double Scaling ;
	double Background ;

	/// Open from .card-file
	fp = fopen(CardFileLocation, "r") ;
	if (fp == NULL)
	{
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\tAn error occured when attempting to open the .card-file\n") ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}
		return -1 ;
	}

	if (CMD) { SubDir = DirectoryFinder(CardFileLocation, Directory) ; }


	// Write the data and names into dummy variable
	linenum = 0 ;

	while (fgets(Buffer, sizeof(Buffer), fp) != NULL)
	{
		sscanf(Buffer, "%s", Dummy[linenum]) ;
		++linenum ;
	}

	// Close .card-file
	fclose(fp) ;

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\t\tAttempt to import %d spectra\n", NumberOfSpectra) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}


	/// Begin reading from .rad-files
	for (i = 0; i < NumberOfSpectra; ++i)
	{
		// Prepare .rad-file
		if(CMD && SubDir){ sprintf(Path, "%s%s",Directory, Dummy[2 + i * 6]) ; }
		else { sprintf(Path, "%s", Dummy[2 + i * 6]) ; }

		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\t\tImport data from file no. %d from %s\n", i+1, Path) ;

			fflush( logfile ) ;
		}

		fp = fopen(Path,"r") ;
		if( fp == NULL )
		{
			if ( abs(WriteLog) > 0 )
			{
				fprintf( logfile, "\t\t\tAn error occured when attempting to open the file\n") ;
				fprintf( logfile, "\n") ;

				fflush( logfile ) ;
			}
			return -1  ;
		}

		Contrast      = atof(Dummy[3 + i * 6]) ;
		Concentration = atof(Dummy[4 + i * 6]) ;
		Scaling       = atof(Dummy[5 + i * 6]) ;
		Background    = atof(Dummy[6 + i * 6]) ;


		// Read lines into dummy variables
		linenum = 0 ;
		while (fgets(Buffer, sizeof(Buffer), fp) != NULL)
		{
			if (sscanf(Buffer, "%lf  %lf  %lf  %lf", &Dummy1, &Dummy2, &Dummy3, &Dummy4) == 4)
			{
				Data[i].QValues[linenum]     = Dummy1 ;
				Data[i].IValues[linenum]     = Scaling * Dummy2 - Background ;
				Data[i].SigmaValues[linenum] = 1.0 / pow(Scaling * Dummy3, 2) ;

				++linenum ;
			}
			else if (sscanf(Buffer, "%lf  %lf  %lf", &Dummy1, &Dummy2, &Dummy3) == 3)
			{
				Data[i].QValues[linenum]     = Dummy1 ;
				Data[i].IValues[linenum]     = Scaling * Dummy2 - Background ;
				Data[i].SigmaValues[linenum] = 1.0 / pow(Scaling * Dummy3, 2) ;

				++linenum ;
			}
		}

		// Close .rad-file
		fclose(fp) ;

		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\t\tProcessed %d data points\n", linenum) ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}

		Data[i].Concentration      = Concentration * 6.022 * 1e17 ; // [1/cm^3] conversion mM -> M solution & L->cm^3, i.e. mmol/L to mol/cm^3, thus 6.022 e23 * 1e-3 * 1e-3 = 6.022 e17
		Data[i].Contrast           = Contrast ;
		Data[i].NumberOfDatapoints = linenum ;

		sprintf( Data[i].Filename, "%s", Path) ;
	}

	return 1  ;
}
