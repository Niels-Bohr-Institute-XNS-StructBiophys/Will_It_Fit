void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Variables describing the lipid core
    double HeightOfCore;
    double HeightOfBelt;
    double HeightOfLipids;
    double HeightOfCurve;
    double HydrophobicMismatch;
    double HeightOfMethyl;

    // Variables describing scattering lengths
    double ScatteringLengthOfHeads;
    double ScatteringLengthOfTail;
    double ScatteringLengthOfBelt;
    double ScatteringLengthOfWater;
    double ScatteringLengthOfMethyl;

    double ScatteringLengthDensityOfWater; 

    double ExcessScatteringLengthDensityOfBelt;
    double ExcessScatteringLengthDensityOfTail;
    double ExcessScatteringLengthDensityOfCaps;
    double ExcessScatteringLengthDensityOfMethyl;

    // Variables describing the volumes
    double VolumeOfHead;
    double VolumeOfTail;
    double VolumeOfBelt;
    double VolumeOfWater;
    double VolumeOfMethyl;

    double CorrectionToVolumeOfHead;
    double CorrectionToVolumeOfTail;
    double CorrectionToVolumeOfBelt;
    double CorrectionToVolumeOfWater;

    // Variables describing water
    double NumberOfWaterAtBelt;
    double NumberOfWaterAtHead;

    // Variables describing the properties of the disc
    double RatioBetweenAxis;
    double MajorAxisOfCore;
    double MinorAxisOfCore;

    // Variables describing the lipids
    double AreaOfDisc;
    double NumberOfLipids;
    double ThicknessOfBelt;

    // Variables describing the endcaps
    double MajorAxisOfEndcaps;
    double MinorAxisOfEndcaps;
    double ScaleFactorOfEndcaps;
    double VerticalAxisOfEllipsoid;
    double VerticalShiftOfEllipsoidCenter;

    /// Get parameters from Parameters
    // Scale factors
    CorrectionToVolumeOfBelt = Parameters[9];
    CorrectionToVolumeOfTail = Parameters[10];
    CorrectionToVolumeOfHead = Parameters[10];
    CorrectionToVolumeOfWater = Parameters[19];

    // Hydration numbers
    NumberOfWaterAtHead = fabs(Parameters[5]);
    NumberOfWaterAtBelt = fabs(Parameters[6]) * VolumesOfMolecules[4] / (VolumesOfMolecules[4] + VolumesOfMolecules[5]);

    // Height of belt
    HeightOfBelt = Parameters[2];

    // Characteristics of lipids
    NumberOfLipids      = Parameters[3];
    HydrophobicMismatch = Parameters[1];
    HeightOfCurve       = fabs(Parameters[21]);
    RatioBetweenAxis    = fabs(Parameters[0]);
    HeightOfCore        = HeightOfBelt + HydrophobicMismatch;

    // Properties of endcaps
    ScaleFactorOfEndcaps = Parameters[22];
    VerticalAxisOfEllipsoid = HeightOfCurve / (1.0 - sqrt(pow(ScaleFactorOfEndcaps, 2) - 1.0) / ScaleFactorOfEndcaps);
    VerticalShiftOfEllipsoidCenter = HeightOfCurve - VerticalAxisOfEllipsoid;

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater  = VolumesOfMolecules[0];
    VolumeOfHead   = VolumesOfMolecules[1] * CorrectionToVolumeOfHead + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfTail   = VolumesOfMolecules[2] * CorrectionToVolumeOfTail;
    VolumeOfMethyl = VolumesOfMolecules[3] * CorrectionToVolumeOfTail;

    VolumeOfBelt = 2.0 * VolumesOfMolecules[4] *  CorrectionToVolumeOfBelt + 2.0 * NumberOfWaterAtBelt * VolumesOfMolecules[0];

    /// More geometry
    // Axis of lipid core
    MinorAxisOfCore = sqrt((NumberOfLipids * VolumeOfMethyl + NumberOfLipids * VolumeOfTail) /
                           (RatioBetweenAxis * pi * HeightOfCore +
                            2.0 * RatioBetweenAxis * pi * pow(ScaleFactorOfEndcaps, 2) *
                            (2.0 * VerticalAxisOfEllipsoid / 3.0 - fabs(VerticalShiftOfEllipsoidCenter) +
                             pow(fabs(VerticalShiftOfEllipsoidCenter), 3) / (3.0 * pow(VerticalAxisOfEllipsoid, 2)))));

    MajorAxisOfCore = RatioBetweenAxis * MinorAxisOfCore;

    AreaOfDisc = pi * MinorAxisOfCore * MajorAxisOfCore;

    // Compute the height of the methyl layer
    HeightOfMethyl = NumberOfLipids * VolumeOfMethyl / AreaOfDisc;

    // Endcaps
    MajorAxisOfEndcaps = ScaleFactorOfEndcaps * MinorAxisOfCore;
    MinorAxisOfEndcaps = ScaleFactorOfEndcaps * MajorAxisOfCore;

    // Compute the height of the lipids
    HeightOfLipids = NumberOfLipids * VolumeOfHead / AreaOfDisc + HeightOfCore;

    /// Obtain scattering lengths
    ScatteringLengthOfWater  = ScatteringLengths[0];
    ScatteringLengthOfHeads  = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    ScatteringLengthOfTail   = ScatteringLengths[2];
    ScatteringLengthOfMethyl = ScatteringLengths[3];
    ScatteringLengthOfBelt   = 2.0 * (ScatteringLengths[4] + NumberOfWaterAtBelt * ScatteringLengths[0]);

    /// Derive scattering length densities
    ScatteringLengthDensityOfWater  = ScatteringLengthOfWater  / VolumeOfWater;
    ExcessScatteringLengthDensityOfCaps   = ScatteringLengthOfHeads  / VolumeOfHead -   ScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfTail   = ScatteringLengthOfTail   / VolumeOfTail -   ScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfMethyl = ScatteringLengthOfMethyl / VolumeOfMethyl - ScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfBelt   = ScatteringLengthOfBelt   / VolumeOfBelt -   ScatteringLengthDensityOfWater;

    /// Assign values of the constraints
    Constraints[0]  = ExcessScatteringLengthDensityOfCaps;
    Constraints[1]  = ExcessScatteringLengthDensityOfTail;
    Constraints[2]  = ExcessScatteringLengthDensityOfMethyl;
    Constraints[3]  = ExcessScatteringLengthDensityOfBelt;
    Constraints[4]  = 0.0;
    Constraints[5]  = HeightOfMethyl;
    Constraints[6]  = HeightOfLipids;
    Constraints[7]  = HeightOfCore;
    Constraints[8]  = HeightOfCurve;
    Constraints[9]  = VerticalAxisOfEllipsoid;
    Constraints[10] = ScaleFactorOfEndcaps;
    Constraints[11] = VerticalShiftOfEllipsoidCenter;
    Constraints[12] = Concentration;
    Constraints[13] = VolumesOfMolecules[2] * CorrectionToVolumeOfHead;
    Constraints[14] = VolumeOfTail;
    Constraints[15] = VolumeOfMethyl;
    Constraints[16] = VolumeOfBelt;
    Constraints[17] = MajorAxisOfEndcaps;
    Constraints[18] = MinorAxisOfEndcaps;

    /// Geometric properties of the disc
    ThicknessOfBelt = -(MajorAxisOfCore + MinorAxisOfCore) / 2.0 +
                         sqrt(pow(MajorAxisOfCore + MinorAxisOfCore, 2) / 4.0 +
                              VolumeOfBelt / (pi * HeightOfBelt));

    Constraints[20] = Parameters[3];
    Constraints[21] = Concentration;
    Constraints[22] = 0.0;
    Constraints[23] = MajorAxisOfCore;
    Constraints[24] = MinorAxisOfCore;
    Constraints[25] = Parameters[3] * VolumeOfTail;
    Constraints[26] = Parameters[3] * VolumeOfHead;
    Constraints[27] = ThicknessOfBelt;
    Constraints[28] = Parameters[3] * VolumeOfMethyl;
}
