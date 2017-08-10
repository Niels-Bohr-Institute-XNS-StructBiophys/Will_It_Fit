double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    double Roughness;
    double Scaling;
    double Background;
    double ConcentrationValue;
    double Intensity = 0.0;
    double PeakValue = Parameters[3];
    double PeakWidth = Parameters[10];

    /// Computation
    Intensity = 1.0 / (pi * PeakWidth * (1.0 + pow((q - PeakValue) / PeakWidth, 2)));

    /// Rescale result and return
    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Roughness  = exp(-pow(q * Parameters[8], 2));
        Scaling    = Parameters[9] * (Contrast / 100.0 * Parameters[0] + (100.0 - Contrast) / 100.0 * Parameters[1]);
        Background = Contrast / 100.0 * Parameters[4] + (100.0 - Contrast) / 100.0 * Parameters[5];
    } else {
        Roughness  = exp(-pow(q * Parameters[7], 2));
        Scaling    = Parameters[9] * Parameters[2];
        Background = Parameters[6];
    }

    Intensity = Intensity * Scaling * Roughness - Background;

    return Intensity;
}
