double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Parameters for model
    double Dummy;
    double Dummy1;
    double Dummy2;
    int i;

    // Variables describing the sample
    double Roughness;
    double ScatteringLengthOfTails;
    double ScatteringLengthOfHeads;
    double ThicknessOfHeads;
    double ThicknessOfTails;
    double Radius;
    double ConcentrationValue;

    // Variables describing the polydispersity
    double RMin;
    double RMax;
    double RStep;
    const int NumberOfBins = 500;
    double DistributionOfRadii[NumberOfBins];
    double DistributionOfRadiiWeight[NumberOfBins];
    double TotalDistribution;
    double Intensity = 0.0;
    double Scaling;
    double Background;

    /// Get values from constraints
    ScatteringLengthOfHeads = Constraints[0];
    ScatteringLengthOfTails = Constraints[1];
    ThicknessOfHeads        = Constraints[2];
    ThicknessOfTails        = Constraints[3];
    ConcentrationValue      = Constraints[4];
    Radius                  = Constraints[5];

    /// Include polydispersity of sample
    // Create gaussian distribution of radii
    RMin = Radius - 3.0 * Radius * Parameters[7];
    RMax = Radius + 3.0 * Radius * Parameters[7];

    if (RMin < 0.0) {
        RMin = 0.0;
    }

    RStep = (RMax - RMin) / (NumberOfBins - 1.0);

    TotalDistribution = 0;

    for (i = 0; i < NumberOfBins; ++i) {
        DistributionOfRadii[i] = RMin + i * RStep;
        Dummy1 = pow((DistributionOfRadii[i] - Radius) / (sqrt(2) * Radius * Parameters[7]), 2);

        if (Dummy1 > 40.0) {
            Dummy2 = 0;
        } else {
            Dummy2 = 1.0 / sqrt(2.0 * pi * pow(Radius * Parameters[7], 2)) * exp(-Dummy1);
        }

        DistributionOfRadiiWeight[i] = Dummy2 * RStep;

        TotalDistribution += DistributionOfRadiiWeight[i];
    }

    for (i = 0; i < NumberOfBins; ++i) {
        Intensity += DistributionOfRadiiWeight[i] / TotalDistribution * Liposome(q, DistributionOfRadii[i], ThicknessOfTails, ThicknessOfHeads, ScatteringLengthOfHeads, ScatteringLengthOfTails, ConcentrationValue);
    }

    /// Rescale result and return
    // Scale the sum using the computed roughness
    if (Contrast == 1) {
        Intensity = Intensity * Parameters[0] * Roughness * Parameters[18] - Parameters[8];
    }

    if (Contrast == 2) {
        Intensity = Intensity * Parameters[1] * Roughness * Parameters[18] - Parameters[9];
    }

    if (Contrast == 3) {
        Intensity = Intensity * Parameters[2] * Roughness * Parameters[18] - Parameters[10];
    }




    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Roughness  = exp(-pow(q * Parameters[13], 2));
        Scaling    = Parameters[18] * (Contrast / 100.0 * Parameters[0] + (100.0 - Contrast) / 100.0 * Parameters[1]);
        Background = Contrast / 100.0 * Parameters[8] + (100.0 - Contrast) / 100.0 * Parameters[9];
    } else {
        Roughness  = exp(-pow(q * Parameters[12], 2));
        Scaling    = Parameters[18] * Parameters[2];
        Background = Parameters[10];
    }

    Intensity = Intensity * Scaling * Roughness - Background;

    return Intensity;
}
