double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Declarations of dummy variables used in computations
    int i;
    double A;
    double B;
    double C;

    // Parameters describing the core
    double HeightOfCore;
    double VolumeOfCore;
    double ScatteringLengthOfCore;

    // Parameters describing the lipids
    double HeightOfLipids;
    double AreaPerHeadgroup;

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

    // Parameters describing the disc
    double MinorAxisOfCore;
    double MajorAxisOfCore;

    // Parameters describing the solution
    double NumberOfLipids;
    double Roughness;
    double Scaling;
    double Background;
    double ConcentrationOfSample;
    double Intensity = 0.0;

    // Properties of distribution
    const int NumberInDistribution = 100;
    double AreaUnderDistribution = 0.0;
    double NLipidsMin;
    double NLipidsMax;
    double SigmaOfNLipids;
    double NLipidsMean;
    double NLipidsStep;
    bool EmptyDistribution = true;
    FILE *DistributionFile;

    /// Get scattering lengths
    ScatteringLengthOfCaps = Constraints[0];
    ScatteringLengthOfCore = Constraints[1];
    ScatteringLengthOfMethyl = Constraints[2];
    ScatteringLengthOfBelt = Constraints[3];

    /// Get the parameters describing the core
    HeightOfLipids   = Constraints[6];
    HeightOfCore     = Constraints[7];
    HeightOfMethyl   = Constraints[8];
    AreaPerHeadgroup = Constraints[9];

    /// Get the parameters describing the belt
    HeightOfBelt    = Parameters[2];
    VolumeOfBelt    = Constraints[16];
    ThicknessOfBelt = Constraints[27];

    /// More parameters
    ConcentrationOfSample = Constraints[21];
    VolumeOfCore          = Constraints[25];
    VolumeOfCaps          = Constraints[26];
    VolumeOfMethyl        = Constraints[28];
    NLipidsMean           = Parameters[3];
    SigmaOfNLipids        = Parameters[4];

    NLipidsMin = NLipidsMean - 3.0 * SigmaOfNLipids * NLipidsMean;
    NLipidsMax = NLipidsMean + 3.0 * SigmaOfNLipids * NLipidsMean;

    if (NLipidsMin < 0.0) {
        NLipidsMin = 0.0;
    }

    NLipidsStep = (NLipidsMax - NLipidsMin) / (1.0 * NumberInDistribution);

    /// Output distribution to file
    DistributionFile = fopen("Distribution.mcp", "w+");
    fprintf(DistributionFile, "Number of lipids   Frequency \n");

    /// Intensity over distribution
    for (i = 0; i < NumberInDistribution; ++i) {
        NumberOfLipids = (i + 0.5) * NLipidsStep + NLipidsMin;

        A = 1.0;
        B = (ThicknessOfBelt - VolumeOfBelt / (pi * HeightOfBelt * ThicknessOfBelt));
        C = NumberOfLipids * AreaPerHeadgroup / (2.0 * pi);

        if (pow(B, 2) - 4.0 * A * C > 0) {
            MinorAxisOfCore = (- B - sqrt(pow(B, 2) - 4.0 * A * C)) / (2.0 * A);
            MajorAxisOfCore = (- B + sqrt(pow(B, 2) - 4.0 * A * C)) / (2.0 * A);

            if (MajorAxisOfCore / MinorAxisOfCore < 1.8) {
                Intensity += Nanodisc(q, HeightOfBelt, HeightOfLipids, HeightOfCore, HeightOfMethyl, VolumeOfBelt, ScatteringLengthOfCaps, ScatteringLengthOfCore,
                                ScatteringLengthOfBelt, ScatteringLengthOfMethyl, NumberOfLipids, ConcentrationOfSample, MajorAxisOfCore, MinorAxisOfCore,
                                VolumeOfCore, VolumeOfCaps, VolumeOfMethyl, ThicknessOfBelt) *
                       exp(- pow(NumberOfLipids - NLipidsMean, 2) / (2 * pow(SigmaOfNLipids * NLipidsMean, 2)));

                AreaUnderDistribution += exp(- pow(NumberOfLipids - NLipidsMean, 2) / (2 * pow(SigmaOfNLipids * NLipidsMean, 2)));

                EmptyDistribution = false;

                fprintf(DistributionFile, "%16g   %16g \n", NumberOfLipids, exp(- pow(NumberOfLipids - NLipidsMean, 2) / (2 * pow(SigmaOfNLipids * NLipidsMean, 2))));
            }
        }
    }

    /// Close outputfile and normalize distribution
    fclose(DistributionFile);

    if (EmptyDistribution == false) {
        Intensity /= AreaUnderDistribution;
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
