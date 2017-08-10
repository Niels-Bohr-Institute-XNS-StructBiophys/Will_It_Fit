double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    double ScatteringLengthDensityOfTails;
    double ScatteringLengthDensityOfHeads;
    double TotalRadius;
    double RadiusOfCore;
    double ConcentrationValue;
    double NumberOfMoleculesPerAggregate;
    double Intensity;
    double Roughness;
    double Scaling;
    double Background;

    /// Get values from constraints
    ScatteringLengthDensityOfHeads = Constraints[0];
    ScatteringLengthDensityOfTails = Constraints[1];
    TotalRadius                    = Constraints[2];
    RadiusOfCore                   = Constraints[3];
    ConcentrationValue             = Constraints[4];
    NumberOfMoleculesPerAggregate  = Constraints[5];

    /// Compute model
    Intensity = MonodisperseMicelles(q, TotalRadius, RadiusOfCore, ScatteringLengthDensityOfHeads, ScatteringLengthDensityOfTails, ConcentrationValue);
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

    Intensity = Scaling * Roughness * Intensity - Background;

    return Intensity;
}
