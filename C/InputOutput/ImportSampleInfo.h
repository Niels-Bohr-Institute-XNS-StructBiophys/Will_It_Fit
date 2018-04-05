/// Check size of sample file
int CheckSizeOfSampleInformation(char SamplesFileLocation[256])
{
	// Declarations
	FILE *fp ;
	char Line[1024] ;

	// Open file
	fp = fopen(SamplesFileLocation, "r") ;

	if (fp == NULL) { return -1 ; }

	// Fetch values from file
	if(fgets(Line, sizeof(Line), fp)) ;

	fclose(fp) ;

	return atoi(Line) ;
}


/// Import sample information
int ImportSampleInformation(struct Dataset * Data, double * VolumesOfMolecules, char SamplesFileLocation[256], int NumberOfSampleInformations, int NumberOfSpectra, char *ModificationName)
{
	/// Declarations
	char Line[1024] ;
	int LineIterator = 0 ;

	// Dummy variables used to store inputs
	double ScatteringLengthsNeutrons100[NumberOfSampleInformations] ;
	double ScatteringLengthsNeutrons0[NumberOfSampleInformations] ;
	double ScatteringLengthsXrays[NumberOfSampleInformations] ;
	char Dummy[2 + 6 * (NumberOfSampleInformations + 2)][228] ;
	

	/// I/O
	// Open file
	FILE *fp = fopen( SamplesFileLocation, "r") ;

	if (fp == NULL) { return -1 ; }

	// Fetch values from file
	while (fgets(Line, sizeof(Line), fp) != NULL)
	{
		sscanf(Line, "%s ", Dummy[LineIterator]) ;
		++LineIterator ;
	}

	for ( int i = 0; i < NumberOfSampleInformations; ++i)
	{
		VolumesOfMolecules[i]           = atof(Dummy[3 + i]) ;
		ScatteringLengthsNeutrons100[i] = atof(Dummy[5 + NumberOfSampleInformations + i]) ;
		ScatteringLengthsNeutrons0[i]   = atof(Dummy[7 + 2 * NumberOfSampleInformations + i]) ;
		ScatteringLengthsXrays[i]       = atof(Dummy[9 + 3 * NumberOfSampleInformations + i]) ;
		ModificationName[0]       = Dummy[11 + 4 * NumberOfSampleInformations + i][0] ;
		ModificationName[1]       = Dummy[11 + 4 * NumberOfSampleInformations + i][1] ;
		ModificationName[2]       = Dummy[11 + 4 * NumberOfSampleInformations + i][2] ;
		ModificationName[3]       = 0 ;

		//printf("%s\n",Dummy[11 + 4 * NumberOfSampleInformations + i]) ;
	}

	// Close file and dump values into ScatteringLengths and VolumesOfMolecules
	fclose(fp) ;


	/// Compute scattering lengths for the various contrasts
	for ( int i = 0; i < NumberOfSpectra; ++i)
	{

		if (Data[i].Contrast >= 0.0 && Data[i].Contrast <= 100.0) {

			for ( int j = 0; j < NumberOfSampleInformations; ++j) {
				Data[i].ScatteringLengths[j] = ScatteringLengthsNeutrons100[j] * Data[i].Contrast / 100.0 +
				                               ScatteringLengthsNeutrons0[j]   * (100.0 - Data[i].Contrast) / 100.0 ;
			}

		}
		else
		{
			for ( int j = 0; j < NumberOfSampleInformations; ++j)
			{
				Data[i].ScatteringLengths[j] = ScatteringLengthsXrays[j] ;
			}
		}
	}

	return 0 ;
}
