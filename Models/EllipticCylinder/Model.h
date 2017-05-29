double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Variables describing the sample
    double ScatteringLength;
    double MinorRadius;
    double MajorRadius;
    double Height;
    double ConcentrationValue;
    double Volume;
    double Intensity;
    double Scaling;
    double Background;

    /// Computation
    ScatteringLength   = Constraints[0];
    MinorRadius        = Constraints[1];
    MajorRadius        = Constraints[2];
    Height             = Constraints[3];
    ConcentrationValue = Constraints[4];

    Volume = pi * MajorRadius * MinorRadius * Height;
    Intensity = EllipticCylinder(q, MinorRadius, MajorRadius, Height, ConcentrationValue, Volume, ScatteringLength);

    /// Rescale result and return
    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Scaling    = Parameters[9] * (Contrast / 100.0 * Parameters[0] + (100.0 - Contrast) / 100.0 * Parameters[1]);
        Background = Contrast / 100.0 * Parameters[6] + (100.0 - Contrast) / 100.0 * Parameters[7];
    } else {
        Scaling    = Parameters[9] * Parameters[2];
        Background = Parameters[8];
    }

    Intensity = Intensity * Scaling - Background;

    return Intensity;
}
