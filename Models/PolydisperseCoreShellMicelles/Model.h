double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    int i;
    const int NumberOfStepsInPolydispersity = 51;
    double Roughness;
    double ScatteringLengthDensityOfTails;
    double ScatteringLengthDensityOfHeads;
    double TotalRadius;
    double RadiusOfCore;
    double ConcentrationValue;
    double NumberOfMoleculesPerAggregate;
    double SigmaNumberOfMoleculesPerAggregate;
    double Intensity = 0.0;
    double NumberOfMoleculesPerAggregateMin;
    double NumberOfMoleculesPerAggregateMax;
    double NumberOfMoleculesPerAggregateStep;
    double CurrentNumberOfMoleculesPerAggregate;
    double VolumeOfTail;
    double VolumeOfHead;
    double TotalVolumeOfCore;
    double TotalVolume;
    double Prefactor;
    double Weight;
    double Scaling;
    double Background;

    /// Get values from constraints
    ScatteringLengthDensityOfHeads = Constraints[0];
    ScatteringLengthDensityOfTails = Constraints[1];
    ConcentrationValue             = Constraints[2];
    VolumeOfHead                   = Constraints[3];
    VolumeOfTail                   = Constraints[4];

    NumberOfMoleculesPerAggregate      = Parameters[7];
    SigmaNumberOfMoleculesPerAggregate = Parameters[11];
    NumberOfMoleculesPerAggregateMin   = NumberOfMoleculesPerAggregate - 3.0 * NumberOfMoleculesPerAggregate * SigmaNumberOfMoleculesPerAggregate;
    NumberOfMoleculesPerAggregateMax   = NumberOfMoleculesPerAggregate + 3.0 * NumberOfMoleculesPerAggregate * SigmaNumberOfMoleculesPerAggregate;

    if (NumberOfMoleculesPerAggregateMin < 0.0) {
        NumberOfMoleculesPerAggregateMin = 0.0;
    }

    NumberOfMoleculesPerAggregateStep = (NumberOfMoleculesPerAggregateMax - NumberOfMoleculesPerAggregateMin) / (1.0 * NumberOfStepsInPolydispersity);

    /// Compute model
    Prefactor = 1.0 / (sqrt(2.0 * pi) * NumberOfMoleculesPerAggregate * SigmaNumberOfMoleculesPerAggregate);

    for (i = 0; i < NumberOfStepsInPolydispersity; ++i) {
        CurrentNumberOfMoleculesPerAggregate = NumberOfMoleculesPerAggregateMin + (i + 0.5) * NumberOfMoleculesPerAggregateStep;

        TotalVolumeOfCore = CurrentNumberOfMoleculesPerAggregate * VolumeOfTail;
        TotalVolume       = CurrentNumberOfMoleculesPerAggregate * (VolumeOfTail + VolumeOfHead);

        RadiusOfCore = cbrt(3.0 * TotalVolumeOfCore / (4.0 * pi));
        TotalRadius  = cbrt(3.0 * TotalVolume       / (4.0 * pi));

        Weight = exp(- pow(CurrentNumberOfMoleculesPerAggregate - NumberOfMoleculesPerAggregate, 2) /
                     (2.0 * pow(NumberOfMoleculesPerAggregate * SigmaNumberOfMoleculesPerAggregate, 2)));

        Intensity += Prefactor * Weight * NumberOfMoleculesPerAggregateStep * MonodisperseMicelles(q, TotalRadius, RadiusOfCore, ScatteringLengthDensityOfHeads, ScatteringLengthDensityOfTails, ConcentrationValue);
    }

    Intensity /= NumberOfMoleculesPerAggregate;

    /// Rescale result and return
    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Roughness  = exp(-pow(q * Parameters[13], 2));
        Scaling    = Parameters[14] * (Contrast / 100.0 * Parameters[0] + (100.0 - Contrast) / 100.0 * Parameters[1]);
        Background = Contrast / 100.0 * Parameters[8] + (100.0 - Contrast) / 100.0 * Parameters[9];
    } else {
        Roughness  = exp(-pow(q * Parameters[12], 2));
        Scaling    = Parameters[14] * Parameters[2];
        Background = Parameters[10];
    }

    Intensity = Intensity * Scaling * Roughness - Background;

    return Intensity;
}
