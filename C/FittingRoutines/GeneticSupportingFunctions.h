/// Agent struct
struct Individual {
    double * Parameters;
    double Chisquare;
    double Fitness;
};

/// Initialize individual
void AllocateSubject(struct Individual *Subject, int NumberOfParameters) {
    // Declaration
    int i;

    // Allocation
    Subject->Parameters = (double *) calloc(NumberOfParameters, sizeof(double));

    for (i = 0; i < NumberOfParameters; ++i) {
        Subject->Parameters[i] = 0.0;
    }
}

/// Definition of fitness
double EvaluateFitness(double Chisquare) {
    return exp(-Chisquare);
}

/// Initialization of individual
void InitializeSubject(struct Individual *Subject, struct Parameter * Parameters, int NumberOfParameters)
{
    // Declaration
    int i;

    // Generation of subject
    for (i = 0; i < NumberOfParameters; ++i) {

        if (Parameters[i].iParameter == true) {
            Subject->Parameters[i] = (Parameters[i].MaxValue - Parameters[i].MinValue) * RandomFraction() + Parameters[i].MinValue;
        } else {
            Subject->Parameters[i] = Parameters[i].Value;
        }
    }
}

/// Printing best candidate so far
void PrintBestCandidate(struct Individual Subject, int NumberOfGenerations, int NumberOfParameters, struct Parameter * Parameters)
{
    // Declarations
    int i;

    // Printing
    ClearScreen( stdout ) ;
    printf("Chisquare   = %10g \n", Subject.Chisquare);
    printf("Generations = %10d \n", NumberOfGenerations);
    printf("\n");

    for (i = 0; i < NumberOfParameters; ++i) {
        printf("%-30s = %-10g \n", Parameters[i].Name, Subject.Parameters[i]);
    }

    printf("\n");
}

/// Generate an individual
void CreateIndividual(struct Individual *Subject, struct Individual ParentOne, struct Individual ParentTwo, int NumberOfParameters)
{
    // Declaration
    int i;

    // Parameter a
    for (i = 0; i < NumberOfParameters; ++i) {

        if (RandomFraction() < 0.5) {
            Subject->Parameters[i] = ParentOne.Parameters[i] ;
        } else {
            Subject->Parameters[i] = ParentTwo.Parameters[i] ;
        }
    }
}

/// Mutate an individual
void MutateSubject(struct Individual *Subject, struct Parameter * Parameters, int NumberOfParameters, double ProbabilityOfMutation)
{
    // Declaration
    int i;

    // Generation of subject
    for (i = 0; i < NumberOfParameters; ++i) {

        if (Parameters[i].iParameter == true) {

            if (RandomFraction() < ProbabilityOfMutation) {
                Subject->Parameters[i] = (Parameters[i].MaxValue - Parameters[i].MinValue) * RandomFraction() + Parameters[i].MinValue;
            }
        }
    }
}

/// Supporting functions for swarm algorithm
double ComputeChiSquareGenetic(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int NumberOfSmearingFolds,
                               double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, int TotalNumberOfDatapoints,
                               int NumberOfFreeParameters)
{
    // Declarations
    double ChiSquare = 0.0;

    int i;
    int j;
    int k;

    double q;
    double DummyQ;
    double Stepsize;
    double Intensity;
    double DummyParameters[NumberOfParameters];
    double SigmaOfQ;

    // Generate parameter array
    for (i = 0; i < NumberOfParameters; ++i) {
        DummyParameters[i] = Parameters[i].Value;
    }

    // Compute constraints and elliptic distribution for each dataset
    for (i = 0; i < NumberOfSpectra; ++i) {
        ComputeConstraints(DummyParameters, VolumesOfMolecules, Data[i].ScatteringLengths, Data[i].Contrast,
                           Data[i].Concentration, Data[i].Constraints, ProteinStructure, &*UserDefinedStructure);

        for (j = Data[i].NMin; j < Data[i].NMax; ++j) {
            q = Data[i].QValues[j];

            if (Data[i].IncludeResolutionEffects == true) {
                Intensity = 0.0;
                SigmaOfQ  = Data[i].SigmaQValues[j];
                Stepsize  = 6.0 * SigmaOfQ / (1.0 * NumberOfSmearingFolds);

                for (k = 0; k < NumberOfSmearingFolds; ++k) {
                    DummyQ = q + (k + 0.5 - NumberOfSmearingFolds / 2.0) * Stepsize;

                    if (DummyQ < 1e-5) {
                        DummyQ = 1e-5;
                    }

                    Intensity += Model(DummyQ, DummyParameters, Data[i].Constraints, Data[i].Contrast, ProteinStructure, &*UserDefinedStructure) * Data[i].ResolutionWeights[j][k];
                }
            } else {
                Intensity = Model(q, DummyParameters, Data[i].Constraints,  Data[i].Contrast, ProteinStructure, &*UserDefinedStructure);
            }

            Data[i].FitValues[j] = Intensity;

            ChiSquare += pow(Data[i].IValues[j] - Data[i].FitValues[j], 2) * Data[i].SigmaValues[j] / (1.0 * (TotalNumberOfDatapoints - NumberOfFreeParameters));
        }
    }

    return ChiSquare;
}
