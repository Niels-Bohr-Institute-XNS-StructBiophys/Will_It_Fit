int ProfileLikelihood(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int MaxIterations, double *ChiXX, int NumberOfSmearingFolds,
                    double * VolumesOfMolecules, struct Protein * Ensemble, int NumberOfProteins, double * ProteinWeights, struct UserDefined * UserDefinedStructure,
                    double DeltaForDifferentiations, double ChiSquareFractile, int NumberOfCycles, char CardFileLocation[256], int NumberOfSampleInformations,
                    int TotalNumberOfDatapoints, int NumberOfFreeParameters, int HighestNumberOfDatapoints)
{
    // Variables used in iterations
    int i;
    int j;
    int k;
    int NumberOfParametersToSearch;
    const int NumberOfSteps = EstimatedNumberOfStepsInLikelihoodProfile;

    double ChiXXDummy;
    double ChiXXInitial;
    double dParameter;
    double ArbitraryStepSize;
    double ArbitraryChiSquare;

    struct Parameter * DummyParameters;
    AllocateParameters(&DummyParameters, NumberOfParameters);

    int * ParametersToSearch;
    Initialize1DIntegerArray(&ParametersToSearch, NumberOfParameters);

    bool Continue;
    bool LimitReached;

    char Outputfile[256];
    FILE *fp;

    // Parameters used in linear interpolation
    double * LowerLikelihoodBoundary;
    Initialize1DArray(&LowerLikelihoodBoundary, NumberOfParameters);
    double * UpperLikelihoodBoundary;
    Initialize1DArray(&UpperLikelihoodBoundary, NumberOfParameters);

    double Slope;
    double Intersection;

    double PreviousChiSquare = 0.0;
    double PreviousParameterValue = 0.0;

    double CurrentChiSquare = 0.0;
    double CurrentParameterValue = 0.0;

    // Use the unreduced chisquare for these computations
    NumberOfFreeParameters = TotalNumberOfDatapoints - 1;

    // Initial chisquare computation
    printf("Initializing Profile Likelihood-computation... \n");

    ChiXXInitial = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, Ensemble, NumberOfProteins,
                                ProteinWeights, &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters);

    printf("Initial chisquare is %g...\n", ChiXXInitial);
    printf("Likelihood boundary is set to %g...\n", ChiXXInitial + ChiSquareFractile);

    // Identify unfixed parameters
    NumberOfParametersToSearch = 0;

    for(i = 0; i < NumberOfParameters; ++i) {

    	if (Parameters[i].iParameter == true) {
     		ParametersToSearch[NumberOfParametersToSearch] = i;
     		++NumberOfParametersToSearch;
     	}

     	DummyParameters[i].iParameter = Parameters[i].iParameter;
        DummyParameters[i].Value = Parameters[i].Value;
        DummyParameters[i].MinValue = Parameters[i].MinValue;
        DummyParameters[i].MaxValue = Parameters[i].MaxValue;
    }

    printf("Investigating likelihood of %d parameters...\n\n", NumberOfParametersToSearch);

    // Main loop
    for (i = 0; i < NumberOfParametersToSearch; ++i) {
        /// Initialize search
        DummyParameters[ParametersToSearch[i]].iParameter = false;

        // Reset parameters
        for (j = 0; j < NumberOfParameters; ++j) {
            DummyParameters[j].Value = Parameters[j].Value;
        }

        ClearScreen();
        printf("Initializing search of parameter %d...\n", ParametersToSearch[i]);

        // Prepare output file and print first datapoint
        sprintf(Outputfile, "ProfileLikelihood/%s.plh", Parameters[ParametersToSearch[i]].Name);

        fp = fopen(Outputfile, "w+");

        fprintf(fp, "%10g %10g \n", Parameters[ParametersToSearch[i]].Value, ChiXXInitial);

        /// Begin search in positive direction
        printf("Beginning search in positive direction...\n");

        ArbitraryStepSize = 0.05 * Parameters[ParametersToSearch[i]].Value;

        if (ArbitraryStepSize <= 0.001) {
            ArbitraryStepSize = 0.001;
        }

        if (DummyParameters[ParametersToSearch[i]].Value + ArbitraryStepSize >= Parameters[ParametersToSearch[i]].MaxValue) {
            ArbitraryStepSize = fabs(Parameters[ParametersToSearch[i]].MaxValue - DummyParameters[ParametersToSearch[i]].Value);
        }

        DummyParameters[ParametersToSearch[i]].Value += ArbitraryStepSize;

        LevenbergMarquardt(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, MaxIterations, &ChiXXDummy, NumberOfSmearingFolds, VolumesOfMolecules, true, false,
                           Ensemble, NumberOfProteins, ProteinWeights, &*UserDefinedStructure, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints,
                           NumberOfFreeParameters, HighestNumberOfDatapoints);

        ArbitraryChiSquare = ChiXXDummy;
        dParameter = sqrt(ChiSquareFractile) * ArbitraryStepSize / (NumberOfSteps * sqrt(fabs(ArbitraryChiSquare - ChiXXInitial)));

        for (j = 0; j < NumberOfParameters; ++j) {
            DummyParameters[j].Value = Parameters[j].Value;
        }

        printf("Stepsize is %g...\n\n", dParameter);

        ChiXXDummy = ChiXXInitial;
        Continue = true;
        j = 0;

        while (Continue) {
            ++j;
            PreviousChiSquare = ChiXXDummy;
            PreviousParameterValue = DummyParameters[ParametersToSearch[i]].Value;

            // Reset parameters
            for (k = 0; k < NumberOfParameters; ++k) {
                DummyParameters[k].Value = Parameters[k].Value;
            }

            DummyParameters[ParametersToSearch[i]].Value = PreviousParameterValue;

            // Modify parameter of interest
            DummyParameters[ParametersToSearch[i]].Value += dParameter;
            LimitReached = false;

            if (DummyParameters[ParametersToSearch[i]].Value >= Parameters[ParametersToSearch[i]].MaxValue) {
                DummyParameters[ParametersToSearch[i]].Value = Parameters[ParametersToSearch[i]].MaxValue;
                LimitReached = true;
            }

            // Search
            LevenbergMarquardt(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, MaxIterations, &ChiXXDummy, NumberOfSmearingFolds, VolumesOfMolecules, true, false,
                            Ensemble, NumberOfProteins, ProteinWeights, &*UserDefinedStructure, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints,
                            NumberOfFreeParameters, HighestNumberOfDatapoints);

            fprintf(fp, "%10g %10g \n", DummyParameters[ParametersToSearch[i]].Value, ChiXXDummy);

            printf("Profile minimum found at %g with a chisquare of %g...\n", DummyParameters[ParametersToSearch[i]].Value, ChiXXDummy);

            CurrentChiSquare = ChiXXDummy;
            CurrentParameterValue = DummyParameters[ParametersToSearch[i]].Value;

            // Check for reasonably smooth convergence
            if (LimitReached == true) {
                Continue = false;
            } else if (CurrentChiSquare > ChiXXInitial + ChiSquareFractile) {
                Continue = false;
            }
        }

        // Compute positive likelihood boundary
        if (LimitReached == false || CurrentChiSquare > ChiXXInitial + ChiSquareFractile) {
            Slope = (CurrentChiSquare - PreviousChiSquare) / (CurrentParameterValue - PreviousParameterValue);
            Intersection = PreviousChiSquare - Slope * PreviousParameterValue;
            UpperLikelihoodBoundary[ParametersToSearch[i]] = ((ChiXXInitial + ChiSquareFractile) - Intersection) / Slope;
        } else {
            printf("Upper parameter boundary reached...\n");
            UpperLikelihoodBoundary[ParametersToSearch[i]] = Parameters[ParametersToSearch[i]].MaxValue;
        }

        /// Begin search in negative direction
        // Reset parameters
        for (j = 0; j < NumberOfParameters; ++j) {
            DummyParameters[j].Value = Parameters[j].Value;
        }

        printf("\n");
        printf("Beginning search in negative direction...\n");

        ArbitraryStepSize = 0.05 * Parameters[ParametersToSearch[i]].Value;

        if (ArbitraryStepSize <= 0.001) {
            ArbitraryStepSize = 0.001;
        }

        if (DummyParameters[ParametersToSearch[i]].Value - ArbitraryStepSize<= Parameters[ParametersToSearch[i]].MinValue) {
            ArbitraryStepSize = fabs(DummyParameters[ParametersToSearch[i]].Value - Parameters[ParametersToSearch[i]].MinValue);
        }

        DummyParameters[ParametersToSearch[i]].Value -= ArbitraryStepSize;

        LevenbergMarquardt(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, MaxIterations, &ChiXXDummy, NumberOfSmearingFolds, VolumesOfMolecules, true, false,
                           Ensemble, NumberOfProteins, ProteinWeights, &*UserDefinedStructure, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints,
                           NumberOfFreeParameters, HighestNumberOfDatapoints);

        ArbitraryChiSquare = ChiXXDummy;
        dParameter = sqrt(ChiSquareFractile) * ArbitraryStepSize / (NumberOfSteps * sqrt(fabs(ArbitraryChiSquare - ChiXXInitial)));

        for (j = 0; j < NumberOfParameters; ++j) {
            DummyParameters[j].Value = Parameters[j].Value;
        }

        printf("Stepsize is %g...\n\n", dParameter);

        ChiXXDummy = ChiXXInitial;
        Continue = true;
        j = 0;

        while (Continue) {
            ++j;
            PreviousChiSquare = ChiXXDummy;
            PreviousParameterValue = DummyParameters[ParametersToSearch[i]].Value;

            // Reset parameters
            for (k = 0; k < NumberOfParameters; ++k) {
                DummyParameters[k].Value = Parameters[k].Value;
            }

            DummyParameters[ParametersToSearch[i]].Value = PreviousParameterValue;

            // Modify parameter of interest
            DummyParameters[ParametersToSearch[i]].Value -= dParameter;
            LimitReached = false;

            if (DummyParameters[ParametersToSearch[i]].Value <= Parameters[ParametersToSearch[i]].MinValue) {
                DummyParameters[ParametersToSearch[i]].Value = Parameters[ParametersToSearch[i]].MinValue;
                LimitReached = true;
            }

            // Search
            LevenbergMarquardt(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, MaxIterations, &ChiXXDummy, NumberOfSmearingFolds, VolumesOfMolecules, true, false,
                            Ensemble, NumberOfProteins, ProteinWeights, &*UserDefinedStructure, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints,
                            NumberOfFreeParameters, HighestNumberOfDatapoints);

            fprintf(fp, "%10g %10g \n", DummyParameters[ParametersToSearch[i]].Value, ChiXXDummy);

            printf("Profile minimum found at %g with a chisquare of %g...\n", DummyParameters[ParametersToSearch[i]].Value, ChiXXDummy);

            CurrentChiSquare = ChiXXDummy;
            CurrentParameterValue = DummyParameters[ParametersToSearch[i]].Value;

            // Check for reasonably smooth convergence
            if (LimitReached == true) {
                Continue = false;
            } else if (CurrentChiSquare > ChiXXInitial + ChiSquareFractile) {
                Continue = false;
            }
        }

        // Compute negative likelihood boundary
        if (LimitReached == false || CurrentChiSquare > ChiXXInitial + ChiSquareFractile) {
            Slope = (CurrentChiSquare - PreviousChiSquare) / (CurrentParameterValue - PreviousParameterValue);
            Intersection = PreviousChiSquare - Slope * PreviousParameterValue;
            LowerLikelihoodBoundary[ParametersToSearch[i]] = ((ChiXXInitial + ChiSquareFractile) - Intersection) / Slope;
        } else {
            printf("Lower parameter boundary reached...\n");
            LowerLikelihoodBoundary[ParametersToSearch[i]] = Parameters[ParametersToSearch[i]].MinValue;
        }

        // Conclude search
        printf("\nConcluding search for parameter %d...\n\n", ParametersToSearch[i]);

        fclose(fp);

        DummyParameters[ParametersToSearch[i]].iParameter = true;
    }

    ClearScreen();

    // Print results to screen
    printf("Profile Likelihood-search results - chisquare fractile = %g:\n", ChiSquareFractile);
    printf("Id: Lower boundary:     Fitted value:       Upper boundary:\n");

    for (i = 0; i < NumberOfParametersToSearch; ++i) {
        printf("%2d: %-15g     %-15g     %-15g \n", ParametersToSearch[i], LowerLikelihoodBoundary[ParametersToSearch[i]],
               Parameters[ParametersToSearch[i]].Value, UpperLikelihoodBoundary[ParametersToSearch[i]]);
    }

    // Print results to file
    fp = fopen("ProfileLikelihood/ProfileLikelihood.plh", "w+");

    fprintf(fp, "Profile Likelihood-search results - chisquare fractile = %g:\n", ChiSquareFractile);
    fprintf(fp, ".card-file = %s:\n", CardFileLocation);
    fprintf(fp, "Id: Lower boundary:     Fitted value:       Upper boundary:\n");

    for (i = 0; i < NumberOfParametersToSearch; ++i) {
        fprintf(fp, "%2d: %-15g     %-15g     %-15g \n", ParametersToSearch[i], LowerLikelihoodBoundary[ParametersToSearch[i]],
                Parameters[ParametersToSearch[i]].Value, UpperLikelihoodBoundary[ParametersToSearch[i]]);
    }

    fclose(fp);

    // Generate original fits
    ChiXXInitial = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, Ensemble, NumberOfProteins,
                                ProteinWeights, &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters);

    // Conclusion
    *ChiXX = ChiXXInitial;
    printf("\n");
    ClearScreen();

    free(DummyParameters);
    free(ParametersToSearch);
    free(LowerLikelihoodBoundary);
    free(UpperLikelihoodBoundary);

    return 0;
}
