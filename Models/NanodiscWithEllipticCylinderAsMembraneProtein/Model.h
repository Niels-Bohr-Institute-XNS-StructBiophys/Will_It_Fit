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

    // Parameters describing the membrane protein
    double HeightOfMembraneProtein;
    double VolumeOfMembraneProtein;
    double ScatteringLengthOfMembraneProtein;
    double MinorAxisOfMembraneProtein;
    double MajorAxisOfMembraneProtein;
    double AngularOffset;
    double DisplacementOfMP;

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
    double Intensity;

    /// Get scattering length densities
    ScatteringLengthOfCaps            = Constraints[0];
    ScatteringLengthOfCore            = Constraints[1];
    ScatteringLengthOfMethyl          = Constraints[2];
    ScatteringLengthOfBelt            = Constraints[3];
    ScatteringLengthOfMembraneProtein = Constraints[5];

    /// Get the parameters describing the core
    HeightOfLipids = Constraints[6];
    HeightOfCore   = Constraints[7];
    HeightOfMethyl = Constraints[8];

    /// Get the parameters describing the belt
    HeightOfBelt = Parameters[2];
    VolumeOfBelt = Constraints[16];

    /// Parameters describing the MP
    HeightOfMembraneProtein    = Constraints[9];
    VolumeOfMembraneProtein    = Constraints[18];
    MinorAxisOfMembraneProtein = Constraints[10];
    MajorAxisOfMembraneProtein = Constraints[11];
    AngularOffset              = Constraints[19];
    DisplacementOfMP           = Constraints[20];

    /// Get geometric parameters
    NumberOfLipids        = Constraints[22];
    ConcentrationOfSample = Constraints[21];
    MajorAxisOfCore       = Constraints[23];
    MinorAxisOfCore       = Constraints[24];
    VolumeOfCore          = Constraints[25];
    VolumeOfCaps          = Constraints[26];
    ThicknessOfBelt       = Constraints[27];
    VolumeOfMethyl        = Constraints[28];

    /// Computation
    Intensity = NanodiscWithEllipticCylinderAsMembraneProtein(q, HeightOfBelt, HeightOfLipids, HeightOfCore, HeightOfMethyl, VolumeOfBelt, ScatteringLengthOfCaps,
                                                              ScatteringLengthOfCore, ScatteringLengthOfBelt, ScatteringLengthOfMethyl, NumberOfLipids, ConcentrationOfSample,
                                                              MajorAxisOfCore, MinorAxisOfCore, VolumeOfCore, VolumeOfCaps, VolumeOfMethyl, ThicknessOfBelt,
                                                              MinorAxisOfMembraneProtein, MajorAxisOfMembraneProtein, HeightOfMembraneProtein,
                                                              ScatteringLengthOfMembraneProtein, VolumeOfMembraneProtein, AngularOffset, DisplacementOfMP);

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
