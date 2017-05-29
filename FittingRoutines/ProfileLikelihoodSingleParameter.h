int ProfileLikelihoodSingleParameter(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int MaxIterations, double *ChiXX,
                                     int NumberOfSmearingFolds, double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure,
                                     double DeltaForDifferentiations, double ChiSquareFractile, int NumberOfCycles, char CardFileLocation[256], int ParameterID,
                                     int NumberOfSampleInformations, int TotalNumberOfDatapoints, int NumberOfFreeParameters, int HighestNumberOfDatapoints)
{
    // Variables used in iterations
    int i;
    int j;
    const int NumberOfSteps = EstimatedNumberOfStepsInLikelihoodProfile;

    double ChiXXDummy;
    double ChiXXInitial;
    double dParameter;
    double ArbitraryStepSize;
    double ArbitraryChiSquare;

    struct Parameter * DummyParameters;
    AllocateParameters(&DummyParameters, NumberOfParameters);

    bool Continue;
    bool LimitReached;

    char Outputfile[256];
    FILE *fp;

    // Parameters used in linear interpolation
    double LowerLikelihoodBoundary;
    double UpperLikelihoodBoundary;

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

    ChiXXInitial = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure,
                                    TotalNumberOfDatapoints, NumberOfFreeParameters);

    printf("Initial chisquare is %g...\n", ChiXXInitial);
    printf("Likelihood boundary is set to %g...\n", ChiXXInitial + ChiSquareFractile);

    // Identify unfixed parameters
    for(i = 0; i < NumberOfParameters; ++i) {
     	DummyParameters[i].iParameter = Parameters[i].iParameter;
        DummyParameters[i].Value = Parameters[i].Value;
        DummyParameters[i].MinValue = Parameters[i].MinValue;
        DummyParameters[i].MaxValue = Parameters[i].MaxValue;
    }

    printf("Investigating likelihood of parameter number %d...\n\n", ParameterID);

    /// Initialize search
    DummyParameters[ParameterID].iParameter = false;

    // Reset parameters
    for (i = 0; i < NumberOfParameters; ++i) {
        DummyParameters[i].Value = Parameters[i].Value;
    }

    ClearScreen();

    // Prepare output file and print first datapoint
    sprintf(Outputfile, "ProfileLikelihood/%s.plh", Parameters[ParameterID].Name);

    fp = fopen(Outputfile, "w+");

    fprintf(fp, "%10g %10g \n", Parameters[ParameterID].Value, ChiXXInitial);

    /// Begin search in positive direction
    printf("Beginning search in positive direction...\n");

    ArbitraryStepSize = 0.05 * Parameters[ParameterID].Value;

    if (ArbitraryStepSize <= 0.001) {
        ArbitraryStepSize = 0.001;
    }

    if (DummyParameters[ParameterID].Value + ArbitraryStepSize >= Parameters[ParameterID].MaxValue) {
        ArbitraryStepSize = fabs(Parameters[ParameterID].MaxValue - DummyParameters[ParameterID].Value);
    }

    DummyParameters[ParameterID].Value += ArbitraryStepSize;

    LevenbergMarquardt(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, MaxIterations, &ChiXXDummy, NumberOfSmearingFolds, VolumesOfMolecules, true, false, ProteinStructure,
                       &*UserDefinedStructure, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints);

    ArbitraryChiSquare = ChiXXDummy;
    dParameter = sqrt(ChiSquareFractile) * ArbitraryStepSize / (NumberOfSteps * sqrt(fabs(ArbitraryChiSquare - ChiXXInitial)));

    for (i = 0; i < NumberOfParameters; ++i) {
        DummyParameters[i].Value = Parameters[i].Value;
    }

    printf("Stepsize is %g...\n\n", dParameter);

    ChiXXDummy = ChiXXInitial;
    Continue = true;
    i = 0;

    while (Continue) {
        ++i;
        PreviousChiSquare = ChiXXDummy;
        PreviousParameterValue = DummyParameters[ParameterID].Value;

        // Reset parameters
        for (j = 0; j < NumberOfParameters; ++j) {
            DummyParameters[j].Value = Parameters[j].Value;
        }

        DummyParameters[ParameterID].Value = PreviousParameterValue;

        // Modify parameter of interest
        DummyParameters[ParameterID].Value += dParameter;
        LimitReached = false;

        if (DummyParameters[ParameterID].Value >= Parameters[ParameterID].MaxValue) {
            DummyParameters[ParameterID].Value = Parameters[ParameterID].MaxValue;
            LimitReached = true;
        }

        // Search
        LevenbergMarquardt(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, MaxIterations, &ChiXXDummy, NumberOfSmearingFolds, VolumesOfMolecules, true, false,
                           ProteinStructure, &*UserDefinedStructure, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters,
                           HighestNumberOfDatapoints);

        fprintf(fp, "%10g %10g \n", DummyParameters[ParameterID].Value, ChiXXDummy);

        printf("Profile minimum found at %g with a chisquare of %g...\n", DummyParameters[ParameterID].Value, ChiXXDummy);

        CurrentChiSquare = ChiXXDummy;
        CurrentParameterValue = DummyParameters[ParameterID].Value;

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
        UpperLikelihoodBoundary = ((ChiXXInitial + ChiSquareFractile) - Intersection) / Slope;
    } else {
        printf("Upper parameter boundary reached...\n");
        UpperLikelihoodBoundary = Parameters[ParameterID].MaxValue;
    }

    /// Begin search in negative direction
    // Reset parameters
    for (i = 0; i < NumberOfParameters; ++i) {
        DummyParameters[i].Value = Parameters[i].Value;
    }

    printf("\n");
    printf("Beginning search in negative direction...\n");

    ArbitraryStepSize = 0.05 * Parameters[ParameterID].Value;

    if (ArbitraryStepSize <= 0.001) {
        ArbitraryStepSize = 0.001;
    }

    if (DummyParameters[ParameterID].Value - ArbitraryStepSize<= Parameters[ParameterID].MinValue) {
        ArbitraryStepSize = fabs(DummyParameters[ParameterID].Value - Parameters[ParameterID].MinValue);
    }

    DummyParameters[ParameterID].Value -= ArbitraryStepSize;

    LevenbergMarquardt(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, MaxIterations, &ChiXXDummy, NumberOfSmearingFolds, VolumesOfMolecules, true, false, ProteinStructure,
                       &*UserDefinedStructure, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints);

    ArbitraryChiSquare = ChiXXDummy;
    dParameter = sqrt(ChiSquareFractile) * ArbitraryStepSize / (NumberOfSteps * sqrt(fabs(ArbitraryChiSquare - ChiXXInitial)));

    for (i = 0; i < NumberOfParameters; ++i) {
        DummyParameters[i].Value = Parameters[i].Value;
    }

    printf("Stepsize is %g...\n\n", dParameter);

    ChiXXDummy = ChiXXInitial;
    Continue = true;
    i = 0;

    while (Continue) {
        ++i;
        PreviousChiSquare = ChiXXDummy;
        PreviousParameterValue = DummyParameters[ParameterID].Value;

        // Reset parameters
        for (j = 0; j < NumberOfParameters; ++j) {
            DummyParameters[j].Value = Parameters[j].Value;
        }

        DummyParameters[ParameterID].Value = PreviousParameterValue;

        // Modify parameter of interest
        DummyParameters[ParameterID].Value -= dParameter;
        LimitReached = false;

        if (DummyParameters[ParameterID].Value <= Parameters[ParameterID].MinValue) {
            DummyParameters[ParameterID].Value = Parameters[ParameterID].MinValue;
            LimitReached = true;
        }

        // Search
        LevenbergMarquardt(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, MaxIterations, &ChiXXDummy, NumberOfSmearingFolds, VolumesOfMolecules, true, false,
                           ProteinStructure, &*UserDefinedStructure, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters,
                           HighestNumberOfDatapoints);

        fprintf(fp, "%10g %10g \n", DummyParameters[ParameterID].Value, ChiXXDummy);

        printf("Profile minimum found at %g with a chisquare of %g...\n", DummyParameters[ParameterID].Value, ChiXXDummy);

        CurrentChiSquare = ChiXXDummy;
        CurrentParameterValue = DummyParameters[ParameterID].Value;

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
        LowerLikelihoodBoundary = ((ChiXXInitial + ChiSquareFractile) - Intersection) / Slope;
    } else {
        printf("Lower parameter boundary reached...\n");
        LowerLikelihoodBoundary = Parameters[ParameterID].MinValue;
    }

    // Conclude search
    fprintf(fp, "\n");
    fprintf(fp, "Conclusion: \n");
    fprintf(fp, "Lower limit: %g \n", LowerLikelihoodBoundary);
    fprintf(fp, "Upper limit: %g \n", UpperLikelihoodBoundary);

    fclose(fp);

    DummyParameters[ParameterID].iParameter = true;

    ClearScreen();

    // Print results to screen
    printf("Profile Likelihood-search results - chisquare fractile = %g:\n", ChiSquareFractile);
    printf("Id: Lower boundary:     Fitted value:       Upper boundary:\n");

    printf("%2d: %-15g     %-15g     %-15g \n", ParameterID, LowerLikelihoodBoundary,
           Parameters[ParameterID].Value, UpperLikelihoodBoundary);

    // Generate original fits
    ChiXXInitial = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure,
                                    TotalNumberOfDatapoints, NumberOfFreeParameters);

    // Conclusion
    *ChiXX = ChiXXInitial;
    printf("\n");
    ClearScreen();

    free(DummyParameters);

    return 0;
}
