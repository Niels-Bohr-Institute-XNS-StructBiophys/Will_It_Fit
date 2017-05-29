double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
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

    // Parameters describing the disc
    double MinorAxisOfCore;
    double MajorAxisOfCore;

    // Parameters describing the his-tags
    double ScatteringLengthOfTag;
    double VolumeOfTag;
    double RadiusOfGyrationForTag;

    // Parameters describing the solution
    double NumberOfLipids;
    double Roughness;
    double Scaling;
    double Background;
    double ConcentrationOfSample;
    double Intensity;

    /// Get excess scattering length densities
    ScatteringLengthOfCaps   = Constraints[0];
    ScatteringLengthOfCore   = Constraints[1];
    ScatteringLengthOfMethyl = Constraints[2];
    ScatteringLengthOfBelt   = Constraints[3];
    ScatteringLengthOfTag    = Constraints[4];

    /// Get the parameters describing the core
    HeightOfLipids = Constraints[6];
    HeightOfCore   = Constraints[7];
    HeightOfMethyl = Constraints[8];

    /// Get the parameters describing the belt
    HeightOfBelt           = Parameters[2];
    VolumeOfBelt           = Constraints[16];
    VolumeOfTag            = Constraints[17];
    RadiusOfGyrationForTag = Parameters[7];;

    /// Get more parameters
    NumberOfLipids        = Constraints[20];
    ConcentrationOfSample = Constraints[21];
    MajorAxisOfCore       = Constraints[23];
    MinorAxisOfCore       = Constraints[24];
    VolumeOfCore          = Constraints[25];
    VolumeOfCaps          = Constraints[26];
    ThicknessOfBelt       = Constraints[27];
    VolumeOfMethyl        = Constraints[28];

    /// Computation
    Intensity = NanodiscWithTags(q, HeightOfBelt, HeightOfLipids, HeightOfCore, HeightOfMethyl, VolumeOfBelt, ScatteringLengthOfCaps, ScatteringLengthOfCore,
                                 ScatteringLengthOfBelt, ScatteringLengthOfMethyl, NumberOfLipids, ConcentrationOfSample, MajorAxisOfCore, MinorAxisOfCore,
                                 VolumeOfCore, VolumeOfCaps, VolumeOfMethyl, ThicknessOfBelt, VolumeOfTag, RadiusOfGyrationForTag, ScatteringLengthOfTag);

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

    /// Return
    return Intensity;
}
