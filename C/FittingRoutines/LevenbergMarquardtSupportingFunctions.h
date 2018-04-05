double ComputeValue(double q, int DatapointID, int DatasetID, struct Dataset * Data, struct Parameter * Parameters, int NumberOfParameters, int NumberOfSmearingFolds, double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
	// Declarations
	int i ;
	double DummyQ ;
	double DummyParameters[NumberOfParameters] ;
	double Intensity ;
	double Stepsize ;
	double SigmaOfQ ;

	// Computation
	for (i = 0; i < NumberOfParameters; ++i) { DummyParameters[i] = Parameters[i].Value ; }

	// with Data[DatasetID] (Data[i]) data it will calculate and return Data[DatasetID].Constraints (e.g. SLD of solvent based on Data[DatasetID].ScatteringLengths and VolumesOfMolecules for the current dataset)
	ComputeConstraints( DummyParameters, VolumesOfMolecules, Data[DatasetID].ScatteringLengths, Data[DatasetID].Contrast, Data[DatasetID].Concentration, Data[DatasetID].Constraints, ProteinStructure, &*UserDefinedStructure) ;

	// fprintf( stdout, "ComputeValue: %d %d %-.4G %-.4G\n", DatasetID, DatapointID, Data[DatasetID].Constraints[0], Data[DatasetID].Constraints[1]) ;

	if (Data[DatasetID].IncludeResolutionEffects == true)
	{
		Intensity = 0.0 ;
		SigmaOfQ  = Data[DatasetID].SigmaQValues[DatapointID] ;
		Stepsize  = 6.0 * SigmaOfQ / (1.0 * NumberOfSmearingFolds) ;

		for (i = 0; i < NumberOfSmearingFolds; ++i)
		{
			DummyQ = q + (i + 0.5 - NumberOfSmearingFolds / 2.0) * Stepsize ;

			if (DummyQ < 1e-5) { DummyQ = 1e-5 ; }

			Intensity += Model(DummyQ, DummyParameters, Data[DatasetID].Constraints, Data[DatasetID].Contrast, ProteinStructure, &*UserDefinedStructure) * Data[DatasetID].ResolutionWeights[DatapointID][i] ;
		}
	}
	else
	{
		Intensity = Model(q, DummyParameters, Data[DatasetID].Constraints, Data[DatasetID].Contrast, ProteinStructure, &*UserDefinedStructure) ;
	}

	return Intensity ;
}


void ComputeTheGradient(double q, int DatapointID, int DatasetID, double Intensity, struct Dataset * Data, struct Parameter * Parameters, int NumberOfParameters, int NumberOfSmearingFolds, double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, double DeltaForDifferentiations, double * Gradient)
{
	// Declarations
	int i;
	int j;
	double DummyQ;
	double dParameter;
	double DummyIntensity;
	double Stepsize;
	double OriginalParameterValue;
	double SigmaOfQ;
	double DummyParameters[NumberOfParameters];

	// Computations
	for (i = 0; i < NumberOfParameters; ++i) { DummyParameters[i] = Parameters[i].Value ; }

	for (i = 0; i < NumberOfParameters; ++i)
	{ 
		if (Parameters[i].iParameter == true)
		{
			dParameter = DeltaForDifferentiations * DummyParameters[i] ;
			OriginalParameterValue = DummyParameters[i] ;

			if (dParameter < 1e-5) { dParameter = 1e-5 ; }

			DummyParameters[i] += dParameter ;

			ComputeConstraints(DummyParameters, VolumesOfMolecules, Data[DatasetID].ScatteringLengths, Data[DatasetID].Contrast, Data[DatasetID].Concentration, Data[DatasetID].Constraints, ProteinStructure, &*UserDefinedStructure) ;

			if (Data[DatasetID].IncludeResolutionEffects == true)
			{
				DummyIntensity = 0.0 ;
				SigmaOfQ       = Data[DatasetID].SigmaQValues[DatapointID] ;
				Stepsize       = 6.0 * SigmaOfQ / (1.0 * NumberOfSmearingFolds) ;

				for (j = 0; j < NumberOfSmearingFolds; ++j)
				{
					DummyQ = q + (j + 0.5 - NumberOfSmearingFolds / 2.0) * Stepsize ;

					DummyIntensity += Model(DummyQ, DummyParameters, Data[DatasetID].Constraints, Data[DatasetID].Contrast, ProteinStructure, &*UserDefinedStructure) * Data[DatasetID].ResolutionWeights[DatapointID][j] ;
				}
			}
			else
			{
				DummyIntensity = Model(q, DummyParameters, Data[DatasetID].Constraints, Data[DatasetID].Contrast, ProteinStructure, &*UserDefinedStructure);
			}

			Gradient[i] = (Intensity - DummyIntensity) / dParameter ;
			DummyParameters[i] = OriginalParameterValue ;
		}
		else
		{
			Gradient[i] = 0.0 ;
		}
	}
}


void ComputeValueAndCovarianceMatrix(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int NumberOfSmearingFolds, double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, double ** AlphaMatrix, double * Beta, double *ChisquareFinal, double DeltaForDifferentiations, int NumberOfSampleInformations, int TotalNumberOfDatapoints, int NumberOfFreeParameters, int HighestNumberOfDatapoints, int WriteLog, FILE* logfile)
{
	// Declarations
	int i ;
	int j ;
	int k ;
	int l ;
	char sdummy[256] ;
	double q ;
	double Intensity ;
	double DifferenceInIntensity ;
	double Chisquare = 0.0 ;
	double DeltaChisquare ;
	double Weight ;
	double * Gradient ;
	struct Parameter * DummyParameters ;
	struct Dataset * DummyData ;
// /*
	struct UserDefined UserDefinedCopy ;
	struct Protein ProteinStructureCopy ;
// */

	// Initialisation
	for (i = 0; i < NumberOfParameters; ++i)
	{
		for (j = 0; j < NumberOfParameters; ++j)
		{
			if (i == j && Parameters[i].iParameter == false) { AlphaMatrix[i][j] = 1.0 ; }
			else { AlphaMatrix[i][j] = 0.0 ; }
		}

		Beta[i] = 0.0 ;
	}

	if ( abs(WriteLog) > 1 )
	{
		fprintf( logfile, "\t\tPrint modeled (difference) intensities and incremental ChiSquares in ComputeValueAndCovarianceMatrix()\n") ;
		fprintf( logfile, "\n") ;
		sprintf( sdummy, "\t\t\t%3s %3s  %-15s  %-15s  %-15s  %-15s  %-17s", "i", "j", "Q_ij", "I_exp_ij", "I_mod_ij", "DeltaI_ij", "DeltaChiSquare_ij") ;
		fprintf( logfile, "%s\n", sdummy) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}

	// Computation
	for (i = 0; i < NumberOfSpectra; ++i)
	{
// /*
		#pragma omp parallel for schedule(dynamic) default(shared) private(q, j, k, l, DeltaChisquare, Intensity, DifferenceInIntensity, Weight, Gradient, DummyData, DummyParameters, ProteinStructureCopy, UserDefinedCopy)
// */
//		#pragma omp parallel for schedule(dynamic) default(none) shared(i, Data, NumberOfSpectra, HighestNumberOfDatapoints, TotalNumberOfDatapoints, NumberOfParameters, NumberOfFreeParameters, NumberOfSampleInformations, VolumesOfMolecules, NumberOfSmearingFolds, Parameters, ProteinStructure, UserDefinedStructure, DeltaForDifferentiations, AlphaMatrix, Chisquare, Beta, WriteLog, logfile) private(q, j, k, l, DeltaChisquare, Intensity, DifferenceInIntensity, Weight, Gradient, DummyData, DummyParameters)
		for (j = Data[i].NMin; j < Data[i].NMax; ++j)
		{
			// Initialization of copies for Data (DummyData), Parameters (DummyParameters) arrays
			AllocateParameters(&DummyParameters, NumberOfParameters) ;
			AllocateData(&DummyData, NumberOfSpectra) ;

// /*
			InitializeUserDefinedStructure(&UserDefinedCopy) ;
			CopyUserDefined(&*UserDefinedStructure, &UserDefinedCopy) ;
// */

			Initialize1DArray(&DummyData[i].Constraints,       MaxNumberOfConstraints) ;
			Initialize1DArray(&DummyData[i].ScatteringLengths, NumberOfSampleInformations) ;

			Initialize1DArray(&DummyData[i].SigmaQValues,      HighestNumberOfDatapoints) ;
			Initialize2DArray(&DummyData[i].ResolutionWeights, HighestNumberOfDatapoints, NumberOfSmearingFolds) ;

			// copy resolution information, q-range and .Concentration & .Contrast only
			DummyData[i].SigmaQValues[j] = Data[i].SigmaQValues[j] ;

			for (k = 0; k < NumberOfSmearingFolds; ++k)
			{
				DummyData[i].ResolutionWeights[j][k] = Data[i].ResolutionWeights[j][k] ;
			}

			DummyData[i].IncludeResolutionEffects = Data[i].IncludeResolutionEffects ;
			DummyData[i].Concentration            = Data[i].Concentration ;
			DummyData[i].Contrast                 = Data[i].Contrast ;
			DummyData[i].NMax                     = Data[i].NMax ;
			DummyData[i].NMin                     = Data[i].NMin ;

			for (k = 0; k < NumberOfSampleInformations; ++k)
			{
				DummyData[i].ScatteringLengths[k] = Data[i].ScatteringLengths[k] ;
			}

			// copy .Value, .iParameter, .MinValue, .MaxValue only
			for (k = 0; k < NumberOfParameters; ++k)
			{
				DummyParameters[k].Value      = Parameters[k].Value ;
				DummyParameters[k].iParameter = Parameters[k].iParameter ;
				DummyParameters[k].MinValue   = Parameters[k].MinValue ;
				DummyParameters[k].MaxValue   = Parameters[k].MaxValue ;
			}

/// *
			// Initialization and full copy of ProteinStructure (ProteinStructureCopy)
			if (ProteinStructure.NumberOfAtoms != 0)
			{
				ProteinStructureCopy.NumberOfAtoms    = ProteinStructure.NumberOfAtoms ;
				ProteinStructureCopy.NumberOfResidues = ProteinStructure.NumberOfResidues ;
				sprintf(ProteinStructureCopy.PDBFileLocation, "%s", ProteinStructure.PDBFileLocation) ;

				AllocateProteinStructure(&ProteinStructureCopy, ProteinStructureCopy.NumberOfResidues, ProteinStructureCopy.NumberOfAtoms) ;

				for (k = 0; k < ProteinStructureCopy.NumberOfResidues; ++k)
				{
					ProteinStructureCopy.Residues[k].xVolume = ProteinStructure.Residues[k].xVolume ;
					ProteinStructureCopy.Residues[k].yVolume = ProteinStructure.Residues[k].yVolume ;
					ProteinStructureCopy.Residues[k].zVolume = ProteinStructure.Residues[k].zVolume ;

					ProteinStructureCopy.Residues[k].xXRayScattering = ProteinStructure.Residues[k].xXRayScattering ;
					ProteinStructureCopy.Residues[k].yXRayScattering = ProteinStructure.Residues[k].yXRayScattering ;
					ProteinStructureCopy.Residues[k].zXRayScattering = ProteinStructure.Residues[k].zXRayScattering ;

					ProteinStructureCopy.Residues[k].xNeutronScattering = ProteinStructure.Residues[k].xNeutronScattering ;
					ProteinStructureCopy.Residues[k].yNeutronScattering = ProteinStructure.Residues[k].yNeutronScattering ;
					ProteinStructureCopy.Residues[k].zNeutronScattering = ProteinStructure.Residues[k].zNeutronScattering ;

					ProteinStructureCopy.Residues[k].XRayScatteringLength    = ProteinStructure.Residues[k].XRayScatteringLength ;
					ProteinStructureCopy.Residues[k].NeutronScatteringLength = ProteinStructure.Residues[k].NeutronScatteringLength ;
					ProteinStructureCopy.Residues[k].Volume                  = ProteinStructure.Residues[k].Volume ;
					ProteinStructureCopy.Residues[k].ResidueID               = ProteinStructure.Residues[k].ResidueID ;

					strcpy(ProteinStructureCopy.Residues[k].Name, ProteinStructure.Residues[k].Name) ;
				}

				for (k = 0; k < ProteinStructureCopy.NumberOfAtoms; ++k)
				{
					ProteinStructureCopy.Atoms[k].x = ProteinStructure.Atoms[k].x ;
					ProteinStructureCopy.Atoms[k].y = ProteinStructure.Atoms[k].y ;
					ProteinStructureCopy.Atoms[k].z = ProteinStructure.Atoms[k].z ;

					ProteinStructureCopy.Atoms[k].XRayScatteringLength    = ProteinStructure.Atoms[k].XRayScatteringLength ;
					ProteinStructureCopy.Atoms[k].NeutronScatteringLength = ProteinStructure.Atoms[k].NeutronScatteringLength ;
					ProteinStructureCopy.Atoms[k].Volume                  = ProteinStructure.Atoms[k].Volume ;
					ProteinStructureCopy.Atoms[k].Name                    = ProteinStructure.Atoms[k].Name ;
				}
			}
//*/


			// Computation
			q = Data[i].QValues[j] ;
			Initialize1DArray(&Gradient, NumberOfParameters) ;


//			fprintf( logfile, "using WIF with shared structs\n") ;
//			Intensity = ComputeValue(q, j, i, DummyData, DummyParameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure) ;


///*
			Intensity = ComputeValue(q, j, i, DummyData, DummyParameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &UserDefinedCopy) ;
//*/

			Data[i].FitValues[j] = Intensity ;
			DifferenceInIntensity = Data[i].IValues[j] - Intensity ;

			DeltaChisquare = DifferenceInIntensity * DifferenceInIntensity * Data[i].SigmaValues[j] / (1.0 * (TotalNumberOfDatapoints - NumberOfFreeParameters)) ;

			#pragma omp atomic
			Chisquare += DeltaChisquare ;

			if ( abs(WriteLog) > 1 )
			{
				// sheared: i, Data[i].IValues[j] ; all other are thread-private variables
				fprintf( logfile, "\t\t\t% 3d %3d  %-.4G  %-.4G  %-.4G  %-.4G  %-.4G\n", i, j, q, Data[i].IValues[j], Intensity, DifferenceInIntensity, DeltaChisquare) ;

				fflush( logfile ) ;
			}


			// Construct alphamatrix
//			ComputeTheGradient(q, j, i, Intensity, DummyData, DummyParameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure, DeltaForDifferentiations, Gradient) ;

// /*
			ComputeTheGradient(q, j, i, Intensity, DummyData, DummyParameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &UserDefinedCopy, DeltaForDifferentiations, Gradient) ;
// */

			// Free variables
			free(DummyData[i].Constraints) ;
			free(DummyData[i].ScatteringLengths) ;
			free(DummyData[i].SigmaQValues) ;

			for (k = 0; k < HighestNumberOfDatapoints; ++k) { free(DummyData[i].ResolutionWeights[k]) ; }

			free(DummyData[i].ResolutionWeights) ;

			free(DummyData) ;
			free(DummyParameters) ;
// /*
			FreeUserDefined(&UserDefinedCopy) ;

			if (ProteinStructure.NumberOfAtoms != 0)
			{
				free(ProteinStructureCopy.Residues) ;
				free(ProteinStructureCopy.Atoms) ;
			}
// */


			// Generate alpha and beta
			for (k = 0; k < NumberOfParameters; ++k)
			{
				if (Parameters[k].iParameter == true)
				{
					Weight = Gradient[k] * Data[i].SigmaValues[j] / (1.0 * (TotalNumberOfDatapoints - NumberOfFreeParameters)) ;

					for (l = 0; l <= k; ++l)
					{

						if (Parameters[l].iParameter == true)
						{
							#pragma omp atomic
							AlphaMatrix[k][l] += Weight * Gradient[l] ;
					 	}
					}

					#pragma omp atomic
					Beta[k] += DifferenceInIntensity * Weight ;
				}
			}

			// Free more memory
			free(Gradient) ;
		}
		if ( abs(WriteLog) > 1 )
		{
			fprintf( logfile, "\n\n") ;

			fflush( logfile ) ;
		}
	}

	// Fill out the other half of the matrix
	for (i = 0; i < NumberOfParameters; ++i)
	{
		for (j = 0; j < i; ++j) { AlphaMatrix[j][i] = AlphaMatrix[i][j] ; }
	}

	*ChisquareFinal = Chisquare ;
}


int GaussJordanElimination(double ** Matrix, int SizeOfMatrix, double * Vector)
{
	// Define macro used to swap to elements
	double Temp ;
	#define SWAP(a, b) { Temp = (a) ; (a) = (b) ; (b) = Temp ; }

	// Integers used for iteration
	int i ;
	int j ;
	int k ;

	// Dummy variables
	double BigDummy ;
	double Dummy ;

	double InversePivot ;
	int Pivot[SizeOfMatrix] ;

	// Locations in the matrix
	int Column = 0 ;
	int Row = 0 ;

	int ColumnIndex[SizeOfMatrix] ;
	int RowIndex[SizeOfMatrix] ;

	// Function Body
	for (i = 0; i < SizeOfMatrix; ++i) { Pivot[i] = 0 ; }

	for (i = 0; i < SizeOfMatrix; ++i)
	{
		BigDummy = 0.0 ;

		// Search for largest element in matrix
		for (j = 0; j < SizeOfMatrix; ++j)
		{
			if (Pivot[j] != 1)
			{
				for (k = 0; k < SizeOfMatrix; ++k)
				{
					if (Pivot[k] == 0)
					{
						if (fabs(Matrix[j][k]) >= BigDummy)
						{
							BigDummy = fabs(Matrix[j][k]) ;
							Row = j ;
							Column = k ;
						}
					// Test for singularity
					}
					else if (Pivot[k] > 1) { return -1 ; }
				}
			}
		}

		++Pivot[Column] ;

		// Put largest element in diagonal
		if (Row != Column)
		{
			for (j = 0; j < SizeOfMatrix; ++j)
			{ 
				SWAP(Matrix[Row][j], Matrix[Column][j]) ;
			}

			SWAP(Vector[Row], Vector[Column]) ;
		}

		RowIndex[i] = Row ;
		ColumnIndex[i] = Column ;

		// Test for singularity
		if (Matrix[Column][Column] == 0.0) { return -1 ; }

		// Invert the pivot and divide the column and the vector element by it
		InversePivot = 1.0 / Matrix[Column][Column] ;
		Matrix[Column][Column] = 1.0 ;

		for (j = 0; j < SizeOfMatrix; ++j) { Matrix[Column][j] *= InversePivot ; }

		Vector[Column] *= InversePivot ;

		// Reduce the rest of the matrix
		for (j = 0; j < SizeOfMatrix; ++j)
		{
			if (j != Column)
			{
				Dummy = Matrix[j][Column] ;
				Matrix[j][Column] = 0.0 ;

				for (k = 0; k < SizeOfMatrix; ++k)
				{
					Matrix[j][k] -= Matrix[Column][k] * Dummy ;
				}

				Vector[j] -= Vector[Column] * Dummy ;
			}
		}
	}

	// Reconstruct matrix to original order
	for (i = SizeOfMatrix - 1; i >= 0; --i)
	{
		if (RowIndex[i] != ColumnIndex[i])
		{
			for (j = 0; j < SizeOfMatrix; ++j)
			{
				SWAP(Matrix[j][RowIndex[i]], Matrix[j][ColumnIndex[i]]) ;
			}
		}
	}

	return 0 ;
}


// void RunLevenbergMarquardt(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int NumberOfSmearingFolds, double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure,  double DeltaForDifferentiations, double ** CovarianceMatrix, double ** AlphaMatrix, double * Lambda, double *Chisquare, double * Beta, int NumberOfSampleInformations, int TotalNumberOfDatapoints, int NumberOfFreeParameters, int HighestNumberOfDatapoints)
void RunLevenbergMarquardt(struct Dataset * Data, int NumberOfSpectra, int TotalNumberOfDatapoints, int HighestNumberOfDatapoints, struct Parameter * Parameters, int NumberOfParameters, int NumberOfFreeParameters, double * VolumesOfMolecules, int NumberOfSampleInformations, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, int NumberOfSmearingFolds, double DeltaForDifferentiations, double *Chisquare, double ** CovarianceMatrix, double * Lambda, double ** AlphaMatrix, double * Beta, int WriteLog, FILE* logfile)
{
	// Declarations
	int i;
	int j;

	int GaussError;
	struct Parameter * DummyParameters;
	static double * ParameterSteps;
	static double * ParameterStepsDummy;
	static double PreviousChisquare;

	// Initialize
	AllocateParameters(&DummyParameters, NumberOfParameters);
	Initialize1DArray(&ParameterSteps, NumberOfParameters);
	Initialize1DArray(&ParameterStepsDummy, NumberOfParameters);

	for (i = 0; i < NumberOfParameters; ++i)
	{
		DummyParameters[i].iParameter = Parameters[i].iParameter ;
	}

	if (*Lambda < 0.0)
	{
		*Lambda = 0.001 ;

		ComputeValueAndCovarianceMatrix(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure, AlphaMatrix, Beta, &*Chisquare, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints, WriteLog, logfile) ;

		PreviousChisquare = *Chisquare ;
	}

	// Prepare CovarianceMatrix matrix
	for (i = 0; i < NumberOfParameters; ++i)
	{
		for (j = 0; j < NumberOfParameters; ++j)
		{
			CovarianceMatrix[i][j] = AlphaMatrix[i][j] ;
		}

		if (Parameters[i].iParameter == true)
		{
			CovarianceMatrix[i][i] = AlphaMatrix[i][i] * (1.0 + (*Lambda)) ;
		}

		ParameterStepsDummy[i] = Beta[i] ;
	}

	// Solve the matrix problem
	GaussError = GaussJordanElimination(CovarianceMatrix, NumberOfParameters, ParameterStepsDummy) ;

	if (GaussError < 0)
	{
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n") ;
			fprintf( logfile, "\t! Singularity in covariance matrix - unable to invert          !\n") ;
			fprintf( logfile, "\t! Perhaps algorithm tries to refine a parameter without impact !\n") ;
			fprintf( logfile, "\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n") ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}

		free(ParameterStepsDummy) ;
		free(ParameterSteps) ;
		free(DummyParameters) ;

		return ;
	}

	// Use the derived stepsizes
	if (*Lambda == 0.0)
	{
		free(ParameterStepsDummy) ;
		free(ParameterSteps) ;
		free(DummyParameters) ;

		return ;
	}

	for (i = 0; i < NumberOfParameters; ++i)
	{
		ParameterSteps[i] = ParameterStepsDummy[i] ;
	}

	// Compute new parameter values
	for (i = 0; i < NumberOfParameters; ++i)
	{
		DummyParameters[i].Value = Parameters[i].Value - ParameterSteps[i] ;

		if (Parameters[i].iParameter == true)
		{
			if (DummyParameters[i].Value > Parameters[i].MaxValue)
			{
				DummyParameters[i].Value = Parameters[i].MaxValue ;
			}

			if (DummyParameters[i].Value < Parameters[i].MinValue)
			{
				DummyParameters[i].Value = Parameters[i].MinValue ;
			}
		}
	}

	ComputeValueAndCovarianceMatrix(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure, CovarianceMatrix, ParameterSteps, &*Chisquare, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints, WriteLog, logfile) ;

	if (*Chisquare < PreviousChisquare)
	{
		*Lambda *= 0.1 ; 
		PreviousChisquare = *Chisquare ;

		for (i = 0; i < NumberOfParameters; ++i)
		{
			for (j = 0; j < NumberOfParameters; ++j)
			{
				AlphaMatrix[i][j] = CovarianceMatrix[i][j] ;
			}

			Beta[i] = ParameterSteps[i] ;
		}

		for (i = 0; i < NumberOfParameters; ++i)
		{
			Parameters[i].Value = DummyParameters[i].Value ;
		}
	}
	else
	{
		*Lambda *= 10.0 ;
		*Chisquare = PreviousChisquare ;
	}

	free(ParameterStepsDummy) ;
	free(ParameterSteps) ;
	free(DummyParameters) ;

	return ;
}


void PrintCovarianceMatrixFnc(int NumberOfParameters, double **CovarianceMatrix, char *ResultsDirectory, int WriteLog, FILE* logfile)
{
	// Declarations
	int i ;
	int j ;
	char CovarianceMatrixFileLocation[256] ;
	FILE *Outputfile ;

	// Print-out
	sprintf( CovarianceMatrixFileLocation, "%s%s", ResultsDirectory, "CovarianceMatrix.dat") ;
	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tWriting covariance matrix to file at %s\n", CovarianceMatrixFileLocation) ;
		fprintf( logfile, "\n") ;
	}
	Outputfile = fopen( CovarianceMatrixFileLocation, "w") ;

	fprintf( Outputfile, "Covariance matrix:\n") ;
	fprintf( Outputfile, "      ") ;

	for (i = 0; i < NumberOfParameters; ++i)
	{
		fprintf( Outputfile, "%10d:  ", i) ;
	}

	fprintf( Outputfile, "\n") ;

	for (i = 0; i < NumberOfParameters; ++i)
	{
		fprintf( Outputfile, "%2d:  ", i) ;

		for (j = 0; j < NumberOfParameters; ++j)
		{
			fprintf( Outputfile, "%12g ", CovarianceMatrix[j][i]) ;
		}

		fprintf( Outputfile, "\n") ;
	}

	fclose( Outputfile ) ;
}
