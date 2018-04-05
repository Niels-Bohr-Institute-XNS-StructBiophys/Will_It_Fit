/// Agent struct
struct Agent {
    double * Parameters;
    double * Velocity;
    double Chisquare;
    double PersonalBestChisquare;
    double * PersonalBestParameters;
};

/// Initialize agent
void InitializeAgent(struct Agent *CurrentAgent, int NumberOfParameters) {
    // Declaration
    int i;

    // Allocation
    CurrentAgent->Parameters             = (double *) calloc(NumberOfParameters, sizeof(double));
    CurrentAgent->Velocity               = (double *) calloc(NumberOfParameters, sizeof(double));
    CurrentAgent->PersonalBestParameters = (double *) calloc(NumberOfParameters, sizeof(double));

    for (i = 0; i < NumberOfParameters; ++i) {
        CurrentAgent->Parameters[i]             = 0.0;
        CurrentAgent->Velocity[i]               = 0.0;
        CurrentAgent->PersonalBestParameters[i] = 0.0;
    }
}

/// RNG
double RandomFraction(void)
{
    // Declarations
    double Value;

    // Calculate pseudo-random number
    Value = rand() / (RAND_MAX * 1.0);

    return Value;
}

/// Supporting functions for swarm algorithm
double ComputeChiSquareSwarm(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int NumberOfSmearingFolds, double * VolumesOfMolecules,
                             struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, int TotalNumberOfDatapoints, int NumberOfFreeParameters)
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
    double SigmaOfQ;
    double DummyParameters[NumberOfParameters];

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
                Intensity = Model(q, DummyParameters, Data[i].Constraints, Data[i].Contrast, ProteinStructure, &*UserDefinedStructure);
            }

            Data[i].FitValues[j] = Intensity;

            ChiSquare += pow(Data[i].IValues[j] - Data[i].FitValues[j], 2) * Data[i].SigmaValues[j] / (1.0 * (TotalNumberOfDatapoints - NumberOfFreeParameters));
        }
    }

    return ChiSquare;
}
