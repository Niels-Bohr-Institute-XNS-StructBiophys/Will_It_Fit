void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Declare dummy variables needed in function
    double Dummy1;
    double Dummy2;

    // Variables describing the lipid core
    double HeightOfCore;
    double HeightOfBelt;
    double HeightOfLipids;
    double HeightOfCurve;
    double HydrophobicMismatch;
    double HeightOfMethyl;

    // Variables describing scattering lengths
    double ExcessScatteringLengthOfHeads;
    double ExcessScatteringLengthOfTail;
    double ExcessScatteringLengthOfTag;
    double ExcessScatteringLengthOfBelt;
    double ExcessScatteringLengthOfWater;
    double ExcessScatteringLengthOfMethyl;

    double ExcessScatteringLengthDensityOfBelt;
    double ExcessScatteringLengthDensityOfTail;
    double ExcessScatteringLengthDensityOfTag;
    double ExcessScatteringLengthDensityOfWater;
    double ExcessScatteringLengthDensityOfCaps;
    double ExcessScatteringLengthDensityOfMethyl;

    // Variables describing the volumes
    double VolumeOfHead;
    double VolumeOfTail;
    double VolumeOfBelt;
    double VolumeOfWater;
    double VolumeOfTag;
    double VolumeOfMethyl;

    double CorrectionToVolumeOfHead;
    double CorrectionToVolumeOfTail;
    double CorrectionToVolumeOfBelt;
    double CorrectionToVolumeOfWater;
    double SizeOfTag;

    // Variables describing water
    double NumberOfWaterAtBelt;
    double NumberOfWaterAtHead;
    double NumberOfWaterAtTag;

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
    NumberOfWaterAtTag  = fabs(Parameters[6]) * VolumesOfMolecules[5] / (VolumesOfMolecules[4] + VolumesOfMolecules[5]);

    // Height of belt
    HeightOfBelt = Parameters[2];

    // Characteristics of lipids
    NumberOfLipids = Parameters[3];
    HydrophobicMismatch = Parameters[1];
    HeightOfCurve = fabs(Parameters[21]);
    HeightOfCore = HeightOfBelt + HydrophobicMismatch;

    // Use fabs to avoid problems with negative squareroots later in the program
    RatioBetweenAxis = fabs(Parameters[0]);

    // Properties of endcaps
    ScaleFactorOfEndcaps = Parameters[22];
    VerticalAxisOfEllipsoid = HeightOfCurve / (1.0 - sqrt(pow(ScaleFactorOfEndcaps, 2) - 1.0) / ScaleFactorOfEndcaps);
    VerticalShiftOfEllipsoidCenter = HeightOfCurve * sqrt(pow(ScaleFactorOfEndcaps, 2) - 1.0) / ScaleFactorOfEndcaps /
                                     (1.0 - sqrt(pow(ScaleFactorOfEndcaps, 2) - 1.0) / ScaleFactorOfEndcaps);

    // Assign tagsize
    SizeOfTag = Parameters[23];

    /// Get parameters from variable VolumesOfMolecules
    // Volume of water
    VolumeOfWater = VolumesOfMolecules[0];

    // Volume of the hydrophilic head groups - including hydration water
    VolumeOfHead = VolumesOfMolecules[1] * CorrectionToVolumeOfHead +
                   NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;

    // Volume of the hydrophobic tails
    VolumeOfTail   = VolumesOfMolecules[2] * CorrectionToVolumeOfTail;
    VolumeOfMethyl = VolumesOfMolecules[3] * CorrectionToVolumeOfTail;

    // Volume of two belts - i.e. the total belt volume
    VolumeOfBelt = 2.0 * (VolumesOfMolecules[4] +
                   (SizeOfTag * VolumesOfMolecules[5])) * CorrectionToVolumeOfBelt +
                   2.0 * NumberOfWaterAtBelt * VolumesOfMolecules[0];

    // Volume of one tag
    VolumeOfTag = (1.0 - SizeOfTag) * VolumesOfMolecules[5] * CorrectionToVolumeOfBelt +
                  NumberOfWaterAtTag * VolumesOfMolecules[0];

    /// More geometry
    // Compute the height of the methyl layer
    Dummy1 = sqrt(pow(ScaleFactorOfEndcaps, 2) - 1.0) / ScaleFactorOfEndcaps;
    Dummy2 = pow(ScaleFactorOfEndcaps, 2) * VerticalAxisOfEllipsoid * (2.0 / 3.0 - Dummy1 + 1.0 / 3.0 * pow(Dummy1, 3));
    HeightOfMethyl = (HeightOfCore + 2.0 * Dummy2) / (VolumeOfTail / VolumeOfMethyl + 1.0);

    // Axis of lipid core and the endcaps
    MinorAxisOfCore = sqrt(NumberOfLipids * VolumeOfMethyl / (HeightOfMethyl * pi * RatioBetweenAxis));
    MajorAxisOfCore = MinorAxisOfCore * RatioBetweenAxis;

    AreaOfDisc = pi * MinorAxisOfCore * MajorAxisOfCore;

    MajorAxisOfEndcaps = ScaleFactorOfEndcaps * MinorAxisOfCore;
    MinorAxisOfEndcaps = ScaleFactorOfEndcaps * MajorAxisOfCore;

    // Compute the height of the different parts of the lipids
    HeightOfLipids = NumberOfLipids * VolumeOfHead / AreaOfDisc + HeightOfCore;

    /// Obtain scattering lengths
    ExcessScatteringLengthOfWater = ScatteringLengths[0];
    ExcessScatteringLengthOfHeads = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    ExcessScatteringLengthOfTail = ScatteringLengths[2];
    ExcessScatteringLengthOfMethyl = ScatteringLengths[3];
    ExcessScatteringLengthOfBelt = 2.0 * (ScatteringLengths[4] + SizeOfTag * ScatteringLengths[5] + NumberOfWaterAtBelt * ScatteringLengths[0]);
    ExcessScatteringLengthOfTag = (1.0 - SizeOfTag) * ScatteringLengths[5] + NumberOfWaterAtTag * ScatteringLengths[0];


    /// Derive scattering length densities
    ExcessScatteringLengthDensityOfWater = ExcessScatteringLengthOfWater / VolumeOfWater;
    ExcessScatteringLengthDensityOfCaps = ExcessScatteringLengthOfHeads / VolumeOfHead - ExcessScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfTail = ExcessScatteringLengthOfTail / VolumeOfTail - ExcessScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfMethyl = ExcessScatteringLengthOfMethyl / VolumeOfMethyl - ExcessScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfBelt = ExcessScatteringLengthOfBelt / VolumeOfBelt - ExcessScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfTag = ExcessScatteringLengthOfTag / VolumeOfTag - ExcessScatteringLengthDensityOfWater;

    /// Assign values of the constraints
    // Values to use in computing model
    Constraints[0] = ExcessScatteringLengthDensityOfCaps;
    Constraints[1] = ExcessScatteringLengthDensityOfTail;
    Constraints[2] = ExcessScatteringLengthDensityOfMethyl;
    Constraints[3] = ExcessScatteringLengthDensityOfBelt;
    Constraints[4] = ExcessScatteringLengthDensityOfTag;
    Constraints[5] = HeightOfMethyl;
    Constraints[6] = HeightOfLipids;
    Constraints[7] = HeightOfCore;
    Constraints[8] = HeightOfCurve;
    Constraints[9] = VerticalAxisOfEllipsoid;
    Constraints[10] = ScaleFactorOfEndcaps;
    Constraints[11] = VerticalShiftOfEllipsoidCenter;
    Constraints[12] = Concentration;
    Constraints[13] = VolumesOfMolecules[2] * CorrectionToVolumeOfHead;
    Constraints[14] = VolumeOfTail;
    Constraints[15] = VolumeOfMethyl;
    Constraints[16] = VolumeOfBelt;
    Constraints[17] = VolumeOfTag;
    Constraints[18] = 0.0;

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
    Constraints[29] = MajorAxisOfEndcaps;
    Constraints[30] = MinorAxisOfEndcaps;
}
