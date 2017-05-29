int BFGS(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int MaxIterations, double * ChiXX, int NumberOfSmearingFolds,
         double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, double DeltaForDifferentiations, bool Gridsearching,
         int TotalNumberOfDatapoints, int NumberOfFreeParameters)
{
    // Declaration of dummy variables
    double Dummy1;
    double Dummy2;
    double Dummy3;

    int i;
    int j;
    int k;

    // Convergence declarations
    const double GradientConvergenceLimit = 1.0e-6;
    const double Epsilon = 3.0e-8;
    const double ToleranceLevel = 4.0 * Epsilon;
    const double MaxStepSize = 100.0;

    // Declarations of supporting variables
    double ChiSquare;
    double MaxStep;

    double DifferenceOfGradients[NumberOfParameters];
    double OldGradient[NumberOfParameters];
    double Gradient[NumberOfParameters];
    double HessianTimesDeltaParameters[NumberOfParameters];
    double HessianMatrix[NumberOfParameters][NumberOfParameters];
    double DeltaParameters[NumberOfParameters];

    struct Parameter * NewParameters;
    AllocateParameters(&NewParameters, NumberOfParameters);

    double ScalarProductOfDifferenceAndDelta;
    double ScalarProductOfDifferenceAndHessianTimesDelta;

    double Sum = 0.0;
    double SumOverDerivatives;
    double SumOverDeltaParameters;

    // Calculate initial value and gradient
    ChiSquare = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure,
                                 TotalNumberOfDatapoints, NumberOfFreeParameters);

    ComputeGradient(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure, Gradient,
                    DeltaForDifferentiations, TotalNumberOfDatapoints, NumberOfFreeParameters);

    // Initialize the Hessian matrix (as the identity matrix)
    for (i = 0; i < NumberOfParameters; ++i) {

        for (j = 0; j < NumberOfParameters; ++j) {
            HessianMatrix[i][j] = 0.0;
        }

        HessianMatrix[i][i] = 1.0;

        DeltaParameters[i] = - Gradient[i];
        Sum += pow(Parameters[i].Value, 2);
    }

    if (sqrt(Sum) > 1.0 * NumberOfParameters) {
        MaxStep = MaxStepSize * sqrt(Sum);
    } else {
        MaxStep = MaxStepSize * NumberOfParameters;
    }

    // Initialize main loop
    i = 0;

    while (i < MaxIterations) {
        if (Gridsearching == false) {
            ClearScreen();
            printf("Iteration nr.: %d \n", i);
            printf("ChiSquare = %g \n", ChiSquare);
            printf("\n");

            for (j = 0; j < NumberOfParameters; ++j) {
                printf("Parameter %2d = %15g          %s \n", j, Parameters[j].Value, Parameters[j].Name);
            }

            printf("\n");
        }

        ++i;

        // Search the gradient for the best step to take
        ChiSquare = PerformLinesearch(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure,
                                      &*UserDefinedStructure, ChiSquare, Gradient, DeltaParameters, NewParameters, MaxStep, TotalNumberOfDatapoints, NumberOfFreeParameters);

        for (j = 0; j < NumberOfParameters; ++j) {
            DeltaParameters[j] = NewParameters[j].Value - Parameters[j].Value;

            if (NewParameters[j].Value < Parameters[j].MinValue) {
                Parameters[j].Value = Parameters[j].MinValue;
            } else if (NewParameters[j].Value > Parameters[j].MaxValue) {
                Parameters[j].Value = Parameters[j].MaxValue;
            } else {
                Parameters[j].Value = NewParameters[j].Value;
            }
        }

        // Check, whether conclusion criteria has been reached
        Dummy1 = 0.0;

        for (j = 0; j < NumberOfParameters; ++j) {

            if (fabs(Parameters[j].Value) > 1.0) {
                Dummy2 = fabs(DeltaParameters[j]) / fabs(Parameters[j].Value);
            } else {
                Dummy2 = fabs(DeltaParameters[j]);
            }

            if (Dummy2 > Dummy1) {
                Dummy1 = Dummy2;
            }
        }

        if (Dummy1 < ToleranceLevel) {

            if (Gridsearching == false) {
                printf("Tolerance level reached - concluding algorithm. \n");
            }

            goto Conclusion;
        }

        // Store old gradient and compute gradient
        for (j = 0; j < NumberOfParameters; ++j) {
            OldGradient[j] = Gradient[j];
        }

        ComputeGradient(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure, Gradient,
                        DeltaForDifferentiations, TotalNumberOfDatapoints, NumberOfFreeParameters);

        // Compute tolerance
        Dummy1 = 0.0;

        if (ChiSquare > 1.0) {
            Dummy3 = ChiSquare;
        } else {
            Dummy3 = 1.0;
        }

        for (j = 0; j < NumberOfParameters; ++j) {

            if (fabs(Parameters[j].Value) > 1.0) {
                Dummy2 = fabs(Gradient[j]) * fabs(Parameters[j].Value) / Dummy3;
            } else {
                Dummy2 = fabs(Gradient[j]) / Dummy3;
            }

            if (Dummy2 > Dummy1) {
                Dummy1 = Dummy2;
            }
        }

        // Check, whether stopping criteria has been reached
        if (Dummy1 < GradientConvergenceLimit) {

            if (Gridsearching == false) {
                printf("Convergence limit reached - concluding algorithm.");
            }

            goto Conclusion;
        }

        // Compute difference between gradients
        for (j = 0; j < NumberOfParameters; ++j) {
            DifferenceOfGradients[j] = Gradient[j] - OldGradient[j];
        }

        // Compute DFP and BFGS terms
        for (j = 0; j < NumberOfParameters; ++j) {
            HessianTimesDeltaParameters[j] = 0.0;

            for (k = 0; k < NumberOfParameters; ++k) {
                HessianTimesDeltaParameters[j] += HessianMatrix[j][k] * DifferenceOfGradients[k];
            }
        }

        // Evaluate scalar products
        ScalarProductOfDifferenceAndDelta = 0.0;
        ScalarProductOfDifferenceAndHessianTimesDelta = 0.0;
        SumOverDerivatives = 0.0;
        SumOverDeltaParameters = 0.0;

        for (j = 0; j < NumberOfParameters; ++j) {
            ScalarProductOfDifferenceAndDelta += DifferenceOfGradients[j] * DeltaParameters[j];
            ScalarProductOfDifferenceAndHessianTimesDelta += DifferenceOfGradients[j] * HessianTimesDeltaParameters[j];
            SumOverDerivatives += pow(DifferenceOfGradients[j], 2);
            SumOverDeltaParameters += pow(DeltaParameters[j], 2);
        }

        // Evaluate update
        if (ScalarProductOfDifferenceAndDelta > sqrt(Epsilon * SumOverDerivatives * SumOverDeltaParameters)) {
            ScalarProductOfDifferenceAndDelta = 1.0 / ScalarProductOfDifferenceAndDelta;
            ScalarProductOfDifferenceAndHessianTimesDelta = 1.0 / ScalarProductOfDifferenceAndHessianTimesDelta;

            // Do BGFS update
            for (j = 0; j < NumberOfParameters; ++j) {
                DifferenceOfGradients[j] = ScalarProductOfDifferenceAndDelta * DeltaParameters[j] -
                                           ScalarProductOfDifferenceAndHessianTimesDelta * HessianTimesDeltaParameters[j];
            }

            for (j = 0; j < NumberOfParameters; ++j) {

                for (k = 0; k < NumberOfParameters; ++k) {
                    HessianMatrix[j][k] += ScalarProductOfDifferenceAndDelta * DeltaParameters[j] * DeltaParameters[k] -
                                           ScalarProductOfDifferenceAndHessianTimesDelta * HessianTimesDeltaParameters[j] *
                                           HessianTimesDeltaParameters[k] + ScalarProductOfDifferenceAndHessianTimesDelta *
                                           DifferenceOfGradients[j] * DifferenceOfGradients[k];
                }
            }
        }

        // Calculate new direction in which to go
        for (j = 0; j < NumberOfParameters; ++j) {
            DeltaParameters[j] = 0.0;

            for (k = 0; k < NumberOfParameters; ++k){
                DeltaParameters[j] -= HessianMatrix[j][k] * Gradient[k];
            }
        }
    }

    // Maximum number of iterations reached
    if (Gridsearching == false) {
        printf("Maximum number of iterations reached - concluding algorithm. \n");
    }

Conclusion:
    ChiSquare = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure,
                                 TotalNumberOfDatapoints, NumberOfFreeParameters);

    *ChiXX = ChiSquare;

    for (i = 0; i < NumberOfParameters; ++i) {
        Parameters[i].Error = 0.0;
    }

    free(NewParameters);
    return 0;
}
