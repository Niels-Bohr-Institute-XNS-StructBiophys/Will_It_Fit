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

    // Parameters describing the solution
    double NumberOfLipids;
    double Roughness;
    double Scaling;
    double Background;
    double ConcentrationOfSample;

    // Declarations of partial sums
    double Intensity = 0.0;

    /// Get excess scattering length densities
    ScatteringLengthOfCaps   = Constraints[0];
    ScatteringLengthOfCore   = Constraints[1];
    ScatteringLengthOfMethyl = Constraints[2];
    ScatteringLengthOfBelt   = Constraints[3];

    /// Get the parameters describing the core
    HeightOfLipids = Constraints[6];
    HeightOfCore   = Constraints[7];
    HeightOfMethyl = Constraints[8];

    /// Get the parameters describing the belt
    HeightOfBelt = Parameters[2];
    VolumeOfBelt = Constraints[16];

    /// Get geometric parameters
    NumberOfLipids        = Constraints[20];
    ConcentrationOfSample = Constraints[21];
    MajorAxisOfCore       = Constraints[23];
    MinorAxisOfCore       = Constraints[24];
    VolumeOfCore          = Constraints[25];
    VolumeOfCaps          = Constraints[26];
    ThicknessOfBelt       = Constraints[27];
    VolumeOfMethyl        = Constraints[28];

    /// Model computation
    Intensity = Nanodisc(q, HeightOfBelt, HeightOfLipids, HeightOfCore, HeightOfMethyl, VolumeOfBelt, ScatteringLengthOfCaps, ScatteringLengthOfCore,
                         ScatteringLengthOfBelt, ScatteringLengthOfMethyl, NumberOfLipids, ConcentrationOfSample, MajorAxisOfCore, MinorAxisOfCore,
                         VolumeOfCore, VolumeOfCaps, VolumeOfMethyl, ThicknessOfBelt);

    Roughness  = Constraints[22];
    Scaling    = Constraints[17];
    Background = Constraints[18];

    Intensity = Intensity * Scaling * exp(- pow(q * Roughness, 2)) - Background;

    return Intensity;
}
