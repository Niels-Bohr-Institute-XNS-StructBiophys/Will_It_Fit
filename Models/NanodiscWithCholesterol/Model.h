double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Scattering lengths
    double ExcessScatteringLengthDensityOfCaps;
    double ExcessScatteringLengthDensityOfTail;
    double ExcessScatteringLengthDensityOfMethyl;
    double ExcessScatteringLengthDensityOfBelt;
    double ExcessScatteringLengthDensityOfTag;

    // Geometry
    double HeightOfMethyl;
    double HeightOfCore;
    double HeightOfLipids;
    double HeightOfBelt;
    double HeightOfCurve;

    double ThicknessOfBelt;
    double MinorAxisOfCore;
    double MajorAxisOfCore;

    double RadiusOfGyrationForTag;

    // Volumes
    double VolumeOfMethyl;
    double VolumeOfCore;
    double VolumeOfCaps;
    double VolumeOfBelt;
    double VolumeOfTag;

    // Endcaps
    double VerticalAxisOfEllipsoid;
    double ScaleFactorOfEndcaps;
    double VerticalShiftOfEllipsoidCenter;

    // Parameters describing the solution
    double NumberOfLipids;
    double Roughness;
    double Scaling;
    double Background;
    double ConcentrationOfSample;
    double Intensity;

    /// Scattering length densities
    ExcessScatteringLengthDensityOfCaps   = Constraints[0];
    ExcessScatteringLengthDensityOfTail   = Constraints[1];
    ExcessScatteringLengthDensityOfMethyl = Constraints[2];
    ExcessScatteringLengthDensityOfBelt   = Constraints[3];
    ExcessScatteringLengthDensityOfTag    = Constraints[4];

    /// Geometry
    HeightOfMethyl = Constraints[5];
    HeightOfLipids = Constraints[6];
    HeightOfCore   = Constraints[7];
    HeightOfCurve  = Constraints[8];

    /// Endcaps
    VerticalAxisOfEllipsoid        = Constraints[9];
    ScaleFactorOfEndcaps           = Constraints[10];
    VerticalShiftOfEllipsoidCenter = Constraints[11];

    /// Belt
    HeightOfBelt = Parameters[2];
    VolumeOfBelt = Constraints[16];

    /// Tags
    VolumeOfTag = Constraints[17];
    RadiusOfGyrationForTag = fabs(Parameters[7]);

    /// More geometric parameters
    NumberOfLipids        = Constraints[20];
    ConcentrationOfSample = Constraints[21];
    MajorAxisOfCore       = Constraints[23];
    MinorAxisOfCore       = Constraints[24];
    VolumeOfCore          = Constraints[25];
    VolumeOfCaps          = Constraints[26];
    ThicknessOfBelt       = Constraints[27];
    VolumeOfMethyl        = Constraints[28];

    /// Computation
    Intensity = NanodiscWithTagsAndEllipsoidalEndcaps(q, HeightOfBelt, HeightOfLipids, HeightOfCore, HeightOfMethyl, VolumeOfBelt, VolumeOfMethyl,
                                                ExcessScatteringLengthDensityOfCaps, ExcessScatteringLengthDensityOfTail,
                                                ExcessScatteringLengthDensityOfMethyl, ExcessScatteringLengthDensityOfBelt, NumberOfLipids,
                                                ConcentrationOfSample, MajorAxisOfCore, MinorAxisOfCore, VolumeOfCore, VolumeOfCaps, ThicknessOfBelt,
                                                VolumeOfTag, RadiusOfGyrationForTag, ExcessScatteringLengthDensityOfTag, HeightOfCurve,
                                                ScaleFactorOfEndcaps, VerticalAxisOfEllipsoid, VerticalShiftOfEllipsoidCenter);

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
