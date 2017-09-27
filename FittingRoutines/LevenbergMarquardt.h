int LevenbergMarquardt(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int MaxIterations, double *ChiXX,
                       int NumberOfSmearingFolds, double * VolumesOfMolecules, bool UseGridsearchRoutine, bool PrintCovarianceMatrix, struct Protein * Ensemble, int NumberOfProteins, 
                       double * ProteinWeights, struct UserDefined * UserDefinedStructure, double DeltaForDifferentiations, int NumberOfSampleInformations, int TotalNumberOfDatapoints,
                       int NumberOfFreeParameters, int HighestNumberOfDatapoints)
{
    // Declarations
    int i = 0;
    int j;
    double Chisquare;
    double Lambda = 1.0;
    double ** AlphaMatrix;
    double ** CovarianceMatrix;
    double * Beta;
    double DummyLambda = 0.0;
    double PreviousChiSquare = 1e100;

    Initialize2DArray(&AlphaMatrix, NumberOfParameters, NumberOfParameters);
    Initialize2DArray(&CovarianceMatrix, NumberOfParameters, NumberOfParameters);
    Initialize1DArray(&Beta, NumberOfParameters);

    // Run main algorithm
    while (i < MaxIterations) {

        if (i == 0) {
            Lambda = -1.0;

            if (UseGridsearchRoutine == false) {
                printf("Initiating the Levenberg-Marquardt-algorithm...\n\n");
            }
        }

        RunLevenbergMarquardt(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, Ensemble, NumberOfProteins, ProteinWeights,
                            &*UserDefinedStructure,DeltaForDifferentiations, CovarianceMatrix, AlphaMatrix, &Lambda, &Chisquare, Beta, NumberOfSampleInformations,
                            TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints);

        ++i;

        if (UseGridsearchRoutine == false) {
            ClearScreen();
            printf("Iteration = %d \n", i);
            printf("Chisquare = %g \n", Chisquare);
            printf("Lambda    = %g \n", Lambda);
            printf("\n");

            printf("Parameters: \n");

            for (j = 0; j < NumberOfParameters; ++j) {
                printf("Parameter %2d = %15g          %s \n", j, Parameters[j].Value, Parameters[j].Name);
            }
            printf("\n");
        }

        if (Lambda > 10 && PreviousChiSquare - Chisquare < 1e-3) {
            break;
        }

        PreviousChiSquare = Chisquare;
    }

    // Compute errors
    if (UseGridsearchRoutine == false) {
        ClearScreen();
        printf("Estimating errors on refined parameters...\n\n");
    }

    RunLevenbergMarquardt(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, Ensemble, NumberOfProteins, ProteinWeights,
                            &*UserDefinedStructure,DeltaForDifferentiations, CovarianceMatrix, AlphaMatrix, &Lambda, &Chisquare, Beta, NumberOfSampleInformations,
                            TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints);

    for (i = 0; i < NumberOfParameters; ++i) {
        Parameters[i].Error = sqrt(CovarianceMatrix[i][i]);
    }

    // Conclude algorithm
    Chisquare = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, Ensemble, NumberOfProteins,
                                ProteinWeights, &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters);

    *ChiXX = Chisquare;

    if (PrintCovarianceMatrix == true) {
        printf("Printing covariance matrix...\n");
        PrintCovarianceMatrixFnc(NumberOfParameters, AlphaMatrix);
    }

    if (UseGridsearchRoutine == false) {
        printf("Concluding algorithm...\n");
    }

    free(Beta);
    Free2DArray(CovarianceMatrix, NumberOfParameters);
    Free2DArray(AlphaMatrix, NumberOfParameters);

    return 0;
}
