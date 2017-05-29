double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Parameters for model
    double Dummy1;
    double Dummy2;
    int i;

    // Variables describing the sample
    double Roughness;
    double ScatteringLength;
    double Thickness;
    double Radius;
    double ConcentrationValue;

    // Variables describing the polydispersity
    double RMin;
    double RMax;
    double RStep;

    const int StepsInPolydispersity = 100;
    double DistributionOfRadii[StepsInPolydispersity];
    double DistributionOfRadiiWeight[StepsInPolydispersity];
    double TotalDistribution;

    // Declarations of partial sums
    double Intensity = 0.0;
    double Scaling;
    double Background;

    /// Get values from constraints
    ScatteringLength = Constraints[0];
    Thickness = Constraints[2];
    ConcentrationValue = Constraints[4];
    Radius = Constraints[5];

    /// Computation
    // Include polydispersity of sample
    RMin = Radius - 3.0 * Radius * Parameters[7];
    RMax = Radius + 3.0 * Radius * Parameters[7];

    if (RMin < 0.0) {
        RMin = 0.0;
    }

    RStep = (RMax - RMin) / (StepsInPolydispersity - 1.0);

    TotalDistribution = 0.0;

    for (i = 0; i < StepsInPolydispersity; ++i) {
        DistributionOfRadii[i] = RMin + i * RStep;
        Dummy1 = 0.5 * pow((DistributionOfRadii[i] - Radius) / (Radius * Parameters[7]), 2);

        if (Dummy1 > 70) {
            Dummy2 = 0.0;
        } else {
            Dummy2 = 1.0 / sqrt(2.0 * pi * pow(Radius * Parameters[7], 2)) * exp(-Dummy1);
        }

        DistributionOfRadiiWeight[i] = Dummy2 * RStep;

        TotalDistribution += DistributionOfRadiiWeight[i];
    }

    for (i = 0; i < StepsInPolydispersity; ++i) {
        Intensity += DistributionOfRadiiWeight[i] / TotalDistribution * SingleShell(q, DistributionOfRadii[i], Thickness, ScatteringLength, ConcentrationValue);
    }

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
