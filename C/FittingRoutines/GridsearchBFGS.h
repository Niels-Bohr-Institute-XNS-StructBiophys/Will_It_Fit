int GridsearchBFGS(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int MaxIterations, double *ChiXX, int NumberOfSmearingFolds,
                   double * VolumesOfMolecules, int NumberOfCyclesInGridsearch, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure,
                   double DeltaForDifferentiations, bool PrintOut, int TotalNumberOfDatapoints, int NumberOfFreeParameters)
{
    // Variables used in iterations
    int i;
    int j;
    int NumberOfCycles;
    int NumberOfParametersToFit;
    double ChiXXDummy = 0.0;
    double BestChisquare = 1e40;
    double StartingValue;

    // Define variables describing fitting parameters
    struct Parameter * DummyParameters;
    AllocateParameters(&DummyParameters, NumberOfParameters);
    int * ParametersToFit;
    Initialize1DIntegerArray(&ParametersToFit, NumberOfParameters);

    // Begin computation
    for(i = 0; i < NumberOfParameters; ++i){
        DummyParameters[i].iParameter = false;
        DummyParameters[i].Value = Parameters[i].Value;
        DummyParameters[i].MinValue = Parameters[i].MinValue;
        DummyParameters[i].MaxValue = Parameters[i].MaxValue;

        if (PrintOut == true) {
            printf("Parameter %2d = %15g          %s \n", i, Parameters[i].Value, Parameters[i].Name);
        }
    }

    if (PrintOut == true) {
        printf("\n");
    }

    NumberOfCycles = NumberOfCyclesInGridsearch;
    NumberOfParametersToFit = 0;

    for (i = 0; i < NumberOfParameters; ++i) {

        if (Parameters[i].iParameter == true) {
            ParametersToFit[NumberOfParametersToFit] = i;

            ++NumberOfParametersToFit;
        }
    }

    if (PrintOut == true) {
        ClearScreen( stdout ) ;
    }

    // Begin iterating
    for (i = 0; i < NumberOfCycles; ++i) {

        for (j = 0; j < NumberOfParametersToFit; ++j) {
            DummyParameters[ParametersToFit[j]].iParameter = true;
            StartingValue = DummyParameters[ParametersToFit[j]].Value;

            BFGS(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, MaxIterations, &ChiXXDummy, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure,
                 &*UserDefinedStructure, DeltaForDifferentiations, true, TotalNumberOfDatapoints, NumberOfFreeParameters);

            DummyParameters[ParametersToFit[j]].iParameter = false;

            if (ChiXXDummy <= BestChisquare) {
                BestChisquare = ChiXXDummy;
                Parameters[ParametersToFit[j]].Error = DummyParameters[ParametersToFit[j]].Error;
                Parameters[ParametersToFit[j]].Value = DummyParameters[ParametersToFit[j]].Value;
            } else {
                DummyParameters[ParametersToFit[j]].Value = StartingValue;
            }

            if (PrintOut == true) {
                printf("Cycle number     = %d \n", i + 1);
                printf("Parameter nr. %2d = %g \n", ParametersToFit[j], DummyParameters[ParametersToFit[j]].Value);
                printf("Chisquare        = %g \n", BestChisquare);
                printf("\n");

                ClearScreen( stdout ) ;
            }
        }
    }

    // Conclusion
    *ChiXX = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure,
                              TotalNumberOfDatapoints, NumberOfFreeParameters);

    if (PrintOut == true) {
        printf("\n");
        ClearScreen( stdout ) ;
    }

    free(DummyParameters);
    free(ParametersToFit);

    return 0;
}
