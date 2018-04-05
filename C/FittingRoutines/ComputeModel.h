int ComputeModel(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, double *ChiXX, int NumberOfSmearingFolds, double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, int TotalNumberOfDatapoints, int NumberOfFreeParameters, bool WriteLog, FILE* logfile)
{
	/// Declarations
	double ChiSquare ;

	/// Comptutation
	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tComputing model function for given parameters\n") ;
		fprintf( logfile, "\t\n") ;

		fflush( logfile ) ;
	}

	ChiSquare = ComputeChiSquare( Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters) ;

	*ChiXX = ChiSquare ;

	return 0 ;
}
