double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Declarations of dummy variables used in computations
    int i;
    double Dummy;

    // Parameters describing the core
    double HeightOfCore;
    double VolumeOfCore;
    double ScatteringLengthOfCore;

    // Parameters describing the lipids
    double HeightOfLipids;

    // Parameters describing the cap
    double VolumeOfCaps;
    double ScatteringLengthOfCaps;

    // Parameters describing the methyl
    double HeightOfMethyl;
    double VolumeOfMethyl;
    double ScatteringLengthOfMethyl;

    // Parameters describing the belt
    double ThicknessOfBelt;
    double HeightOfBelt;
    double VolumeOfBelt;
    double ScatteringLengthOfBelt;

    // Parameters describing the solution
    double NumberOfLipids;
    double Roughness;
    double Scaling;
    double Background;
    double ConcentrationOfSample;
    double Intensity = 0.0;

    // Variables used in integral over NLipids
    double Radius;
    double AreaPerHeadgroup;
    double Area;
    double Prefactor;
    double Weight;

    double SigmaOfNLipids;
    double NLipids;
    double NLipidsMean;
    double NLipidsMin;
    double NLipidsMax;
    double NLipidsStep;
    int NumberOfStepsInNLipids = 50;

    /// Scattering lengths
    ScatteringLengthOfCaps   = Constraints[0];
    ScatteringLengthOfCore   = Constraints[1];
    ScatteringLengthOfMethyl = Constraints[2];
    ScatteringLengthOfBelt   = Constraints[3];

    /// Get the parameters describing the core
    HeightOfLipids = Constraints[6];
    HeightOfCore   = Constraints[7];
    HeightOfMethyl = Constraints[8];

    /// Get the parameters describing the belt
    HeightOfBelt    = Parameters[2];
    VolumeOfBelt    = Constraints[16];
    ThicknessOfBelt = Parameters[0];

    /// Other parameters
    NumberOfLipids        = Constraints[20];
    ConcentrationOfSample = Constraints[21];
    VolumeOfCore          = Constraints[25];
    VolumeOfCaps          = Constraints[26];
    VolumeOfMethyl        = Constraints[28];

    /// Initialize integration over NLipids
    AreaPerHeadgroup = Parameters[1];
    SigmaOfNLipids   = Parameters[4];
    NLipidsMean      = Parameters[3];

    NLipidsMin = NLipidsMean - 3.0 * SigmaOfNLipids * NLipidsMean;
    NLipidsMax = NLipidsMean + 3.0 * SigmaOfNLipids * NLipidsMean;

    if (NLipidsMin < 0.0) {
        NLipidsMin = 0.0;
    }

    NLipidsStep = (NLipidsMax - NLipidsMin) / (1.0 * NumberOfStepsInNLipids);

    Prefactor = 1.0 / (SigmaOfNLipids * NLipidsMean * sqrt(2.0 * pi));

    /// Computation
    for (i = 0; i < NumberOfStepsInNLipids; ++i) {
        NLipids = (i + 0.5) * NLipidsStep + NLipidsMin;
        Area = NLipids * AreaPerHeadgroup / 2.0;
        Radius = sqrt(Area / pi);

        Weight = exp(- pow((NLipids - NLipidsMean) / (sqrt(2.0) * SigmaOfNLipids * NLipidsMean), 2));

        VolumeOfBelt = pi * HeightOfBelt * pow(Radius + ThicknessOfBelt, 2) - pi * HeightOfBelt * pow(Radius, 2);

        Intensity += Prefactor * Weight * NLipidsStep *
                     PeptideDiscs(q, HeightOfBelt, HeightOfLipids, HeightOfCore, HeightOfMethyl, VolumeOfBelt, ScatteringLengthOfCaps, ScatteringLengthOfCore,
                                  ScatteringLengthOfBelt, ScatteringLengthOfMethyl, NumberOfLipids, ConcentrationOfSample, Radius, VolumeOfCore, VolumeOfCaps,
                                  VolumeOfMethyl, ThicknessOfBelt);
    }

    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Roughness  = exp(-pow(q * Parameters[18], 2));
        Scaling    = Parameters[20] * (Contrast / 100.0 * Parameters[13] + (100.0 - Contrast) / 100.0 * Parameters[15]);
        Background = Contrast / 100.0 * Parameters[12] + (100.0 - Contrast) / 100.0 * Parameters[14];
    } else {
        Roughness  = exp(-pow(q * Parameters[8], 2));
        Scaling    = Parameters[20] * Parameters[17];
        Background = Parameters[16];
    }

    Intensity = Intensity * Scaling * Roughness - Background;

    return Intensity;
}
