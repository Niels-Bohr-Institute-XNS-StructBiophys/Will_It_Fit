// int LevenbergMarquardt(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int MaxIterations, double *ChiXX, int NumberOfSmearingFolds, double * VolumesOfMolecules, bool UseGridsearchRoutine, bool PrintCovarianceMatrix, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, double DeltaForDifferentiations, int NumberOfSampleInformations, int TotalNumberOfDatapoints, int NumberOfFreeParameters, int HighestNumberOfDatapoints)
// keep UseGridsearchRoutine variable (is used with false for LevenbergMarquardt to print iteration steps and used with true for GridsearchLM, ProfileLikelihood & ProfileLikelihoodSingleParameter to suppress output of iteration steps)
int LevenbergMarquardt(struct Dataset * Data, int NumberOfSpectra, int TotalNumberOfDatapoints, int HighestNumberOfDatapoints, struct Parameter * Parameters, int NumberOfParameters, int NumberOfFreeParameters, double * VolumesOfMolecules, int NumberOfSampleInformations, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, int MaxIterations, int NumberOfSmearingFolds, bool UseGridsearchRoutine, double DeltaForDifferentiations, double *ChiXX, bool PrintCovarianceMatrix, char *ResultsDirectory, int WriteLog, FILE* logfile)
{
	// Declarations
	double Chisquare ;
	double Lambda = 1.0 ;
	double ** AlphaMatrix ;
	double ** CovarianceMatrix ;
	double * Beta ;
	double DummyLambda = 0.0 ;
	double PreviousChiSquare = 1e100 ;

	Initialize2DArray( &AlphaMatrix, NumberOfParameters, NumberOfParameters) ;
	Initialize2DArray( &CovarianceMatrix, NumberOfParameters, NumberOfParameters) ;
	Initialize1DArray( &Beta, NumberOfParameters) ;

	/* option to write input arguments to logfile */
	if (  UseGridsearchRoutine == false && abs(WriteLog) > 1 )
	{
		fprintf( logfile, "\tSummarizing initial arguments for LevenbergMarquardt algorithm\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\t\tMaxIterations = %d\n", MaxIterations) ;
		fprintf( logfile, "\t\tNumberOfSpectra = %d\n", NumberOfSpectra) ;
		fprintf( logfile, "\t\tNumberOfParameters = %d\n", NumberOfParameters) ;
		fprintf( logfile, "\t\tNumberOfFreeParameters = %d\n", NumberOfFreeParameters) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\t\tNumberOfSampleInformations = %d\n", NumberOfSampleInformations) ;
		fprintf( logfile, "\t\tHighestNumberOfDatapoints = %d\n", HighestNumberOfDatapoints) ;
		fprintf( logfile, "\t\tTotalNumberOfDatapoints = %d\n", TotalNumberOfDatapoints) ;
		fprintf( logfile, "\t\tNumberOfSmearingFolds = %d\n", NumberOfSmearingFolds) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\t\tUseGridsearchRoutine = %d\n", (int)UseGridsearchRoutine) ;
		fprintf( logfile, "\t\tPrintCovarianceMatrix = %d\n", (int)PrintCovarianceMatrix) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\t\tDeltaForDifferentiations = %lf\n", DeltaForDifferentiations) ;
		fprintf( logfile, "\t\t*ChiXX = %lf\n", *ChiXX) ;
		fprintf( logfile, "\n\n") ;

		fprintf( logfile, "\t\tInitial Fit Parameters[i]\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\t\t\ti\n") ;
		fprintf( logfile, "\t\t\t.Name\n") ;
		fprintf( logfile, "\t\t\t.MinValue\n") ;
		fprintf( logfile, "\t\t\t.Value\n") ;
		fprintf( logfile, "\t\t\t.MaxValue\n") ;
		fprintf( logfile, "\t\t\t.(int)iParameter\n") ;
		fprintf( logfile, "\t\t\t.Error\n") ;
		fprintf( logfile, "\n") ;
		for (int i = 0; i < NumberOfParameters; ++i)
		{
			fprintf( logfile, "\t\t\t%2d %10s %-15g %-15g %-15g %d %-15g\n", i, Parameters[i].Name, Parameters[i].MinValue, Parameters[i].Value, Parameters[i].MaxValue, (int)Parameters[i].iParameter, Parameters[i].Error) ;
		}
		fprintf( logfile, "\n\n") ;

		for ( int i = 0; i < NumberOfSpectra; ++i)
		{
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\t\tSpectra no. %d\n", i+1) ;
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\t\t\tData[%d]\n", i) ;
			fprintf( logfile, "\t\t\t\t.Filename                 = %s\n", Data[i].Filename) ;
			fprintf( logfile, "\t\t\t\t.NumberOfDatapoints       = %d\n", Data[i].NumberOfDatapoints) ;
			fprintf( logfile, "\t\t\t\t.NMin                     = %d\n", Data[i].NMin) ;
			fprintf( logfile, "\t\t\t\t.NMax                     = %d\n", Data[i].NMax) ;
			fprintf( logfile, "\t\t\t\t.IncludeResolutionEffects = %d\n", (int)Data[i].IncludeResolutionEffects) ;
			fprintf( logfile, "\t\t\t\t.Concentration            = %lf\n", Data[i].Concentration) ;
			fprintf( logfile, "\t\t\t\t.Contrast                 = %lf\n", Data[i].Contrast) ;
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\t\t\t\tj\n") ;
			fprintf( logfile, "\t\t\t\t.QValues[j]\n") ;
			fprintf( logfile, "\t\t\t\t.IValues[j]\n") ;
			fprintf( logfile, "\t\t\t\t.SigmaValues[j]\n") ;
			fprintf( logfile, "\t\t\t\t.SigmaQValues[j]\n") ;

			for ( int j = 0; j < Data[i].NumberOfDatapoints; ++j)
			{
				fprintf( logfile, "\t\t\t\t%3d %-15f %-15f %-15f %-15f\n", j, Data[i].QValues[j], Data[i].IValues[j], Data[i].SigmaValues[j], Data[i].SigmaQValues[j]) ;
			}
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\t\t\t\tj\n") ;
			fprintf( logfile, "\t\t\t\t.ScatteringLengths[j]\n") ;
			fprintf( logfile, "\n") ;
			for ( int j = 0; j < NumberOfSampleInformations; ++j)
			{
				fprintf( logfile, "\t\t\t\t%2d %-15g\n", j, Data[i].ScatteringLengths[j]) ;
			}
		}
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;
	}



	int it = 0 ;
	// Run main algorithm
	while ( it < MaxIterations)
	{
		if ( UseGridsearchRoutine == false && abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t----------------------------------------------------------------\n") ;
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\tIteration = %d\n", it) ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}

		// initiate Lambda to -1 ; by this ComputeValueAndCovarianceMatrix() will be called twice in RunLevenbergMarquardt()
		if ( it == 0 )
		{
			if ( abs(WriteLog) > 0 )
			{
				fprintf( logfile, "\t\tInitiating the Levenberg-Marquardt-algorithm\n") ;
				fprintf( logfile, "\n") ;

				fflush( logfile ) ;
			}

			Lambda = -1.0 ;
		}

		// RunLevenbergMarquardt(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure, DeltaForDifferentiations, CovarianceMatrix, AlphaMatrix, &Lambda, &Chisquare, Beta, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints) ;
		RunLevenbergMarquardt( Data, NumberOfSpectra, TotalNumberOfDatapoints, HighestNumberOfDatapoints, Parameters, NumberOfParameters, NumberOfFreeParameters, VolumesOfMolecules, NumberOfSampleInformations, ProteinStructure, &*UserDefinedStructure, NumberOfSmearingFolds, DeltaForDifferentiations, &Chisquare, CovarianceMatrix, &Lambda, AlphaMatrix, Beta, WriteLog, logfile) ;

		if ( UseGridsearchRoutine == false && abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t\tChiSquare = %g\n", Chisquare) ;
			fprintf( logfile, "\t\tLambda    = %g\n", Lambda) ;
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\t\tParameters[i]\n") ;
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\t\t\ti\n") ;
			fprintf( logfile, "\t\t\t.Name\n") ;
			fprintf( logfile, "\t\t\t.MinValue\n") ;
			fprintf( logfile, "\t\t\t.Value\n") ;
			fprintf( logfile, "\t\t\t.MaxValue\n") ;
			fprintf( logfile, "\t\t\t.(int)iParameter\n") ;
			// fprintf( logfile, "\t\t\t.Error\n") ;
			fprintf( logfile, "\n") ;
			for ( int i = 0; i < NumberOfParameters; ++i)
			{
				// fprintf( logfile, "\t\t%2d %10s %-15g %-15g %-15g %d %-15g\n", i, Parameters[i].Name, Parameters[i].MinValue, Parameters[i].Value, Parameters[i].MaxValue, (int)Parameters[i].iParameter , Parameters[i].Error) ;
				fprintf( logfile, "\t\t\t%2d %10s %-15g %-15g %-15g %d\n", i, Parameters[i].Name, Parameters[i].MinValue, Parameters[i].Value, Parameters[i].MaxValue, (int)Parameters[i].iParameter) ;
			}
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}

		if (Lambda > 10 && PreviousChiSquare - Chisquare < 1e-3) { break ; }

		PreviousChiSquare = Chisquare ;

		++it ;
	}

	if (abs(WriteLog) > 0 ) 
	{
		fprintf( logfile, "\t----------------------------------------------------------------\n") ;
		fprintf( logfile, "\n") ;
		if ( it == MaxIterations ) { fprintf( logfile, "\tAlgorithm terminated at maximal iteration step it = %d, maybe increase number of iteration steps\n", it) ; }
		else { fprintf( logfile, "\tAlgorithm terminated at iteration step it = %d, since ChiSquare does not change anymore and Lambda > 10\n", it) ; }
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}


	// Compute errors
	if ( UseGridsearchRoutine == false && abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tEstimating errors on refined parameters with CovarianceMatrix\n") ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}

	// RunLevenbergMarquardt(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure, DeltaForDifferentiations, CovarianceMatrix, AlphaMatrix, &DummyLambda, &Chisquare, Beta, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints) ;
	RunLevenbergMarquardt( Data, NumberOfSpectra, TotalNumberOfDatapoints, HighestNumberOfDatapoints, Parameters, NumberOfParameters, NumberOfFreeParameters, VolumesOfMolecules, NumberOfSampleInformations, ProteinStructure, &*UserDefinedStructure, NumberOfSmearingFolds, DeltaForDifferentiations, &Chisquare, CovarianceMatrix, &DummyLambda, AlphaMatrix, Beta, WriteLog, logfile) ;

	for ( int i = 0; i < NumberOfParameters; ++i) { Parameters[i].Error = sqrt(CovarianceMatrix[i][i]) ; }

	// Conclude algorithm
	Chisquare = ComputeChiSquare( Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters) ;

	*ChiXX = Chisquare ;

	if ( PrintCovarianceMatrix == true )
	{
		PrintCovarianceMatrixFnc( NumberOfParameters, AlphaMatrix, ResultsDirectory, WriteLog, logfile) ;
	}

	if ( UseGridsearchRoutine == false && abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tConcluding algorithm\n") ;
		fprintf( logfile, "\n") ;
	}

	free(Beta) ;
	Free2DArray( CovarianceMatrix, NumberOfParameters) ;
	Free2DArray( AlphaMatrix, NumberOfParameters) ;

	return 0 ;
}
