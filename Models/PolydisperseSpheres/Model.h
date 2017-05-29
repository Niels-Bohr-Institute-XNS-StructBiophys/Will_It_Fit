double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    int i;

    int NumberOfSteps = 100;
	double MeanRadius;
	double MinRadius;
    double MaxRadius;
    double Radius;
    double RadiusStepsize;
    double SigmaRadius;
    double Prefactor;
    double Weight;

    double Roughness;
    double Scaling;
    double Background;
    double ScatteringLength;
    double ConcentrationValue;
    double Intensity = 0.0;

    /// Computation
    ScatteringLength   = Constraints[0];
    MeanRadius         = Constraints[1];
    ConcentrationValue = Constraints[2];
    SigmaRadius        = Parameters[10];

    MinRadius = MeanRadius - 3.0 * SigmaRadius * MeanRadius;
    MaxRadius = MeanRadius + 3.0 * SigmaRadius * MeanRadius;
   
    RadiusStepsize = (MaxRadius - MinRadius) / (1.0 * NumberOfSteps);

    Prefactor = 1.0 / (SigmaRadius * MeanRadius * sqrt(2 * pi));

    for (i = 0; i < NumberOfSteps; ++i) {
        Radius = MinRadius + i * RadiusStepsize;
        Weight = exp(- pow(MeanRadius - Radius, 2) / (2.0 * pow(SigmaRadius * MeanRadius, 2)));

        Intensity += Prefactor * Weight * Spheres(q, Radius, ConcentrationValue, ScatteringLength);
    }

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
