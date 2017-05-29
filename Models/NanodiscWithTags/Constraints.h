void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Variables describing the lipid core
    double HeightOfCore;
    double HeightOfBelt;
    double HeightOfLipids;
    double HeightOfMethyl;

    // Variables describing scattering lengths
    double ScatteringLengthOfHeads;
    double ScatteringLengthOfTail;
    double ScatteringLengthOfTag;
    double ScatteringLengthOfBelt;
    double ScatteringLengthOfWater;
    double ScatteringLengthOfMethyl;

    double ExcessScatteringLengthDensityOfBelt;
    double ExcessScatteringLengthDensityOfCore;
    double ExcessScatteringLengthDensityOfTag;
    double ScatteringLengthDensityOfWater;
    double ExcessScatteringLengthDensityOfCaps;
    double ExcessScatteringLengthDensityOfMethyl;

    // Variables describing the valumes
    double VolumeOfHead;
    double VolumeOfTail;
    double VolumeOfBelt;
    double VolumeOfWater;
    double VolumeOfMethyl;
    double VolumeOfTag;

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
    double AverageAreaOfLipids;
    double NumberOfLipids;
    double AreaPerHeadgroupOfLipids;
    double ThicknessOfBelt;

    /// Get parameters from Parameters
    // Scale factor for the belt
    CorrectionToVolumeOfBelt  = Parameters[9];
    CorrectionToVolumeOfTail  = Parameters[10];
    CorrectionToVolumeOfHead  = Parameters[10];
    CorrectionToVolumeOfWater = Parameters[19];

    // Other parameters
    NumberOfWaterAtHead = Parameters[5];
    NumberOfWaterAtBelt = Parameters[6] * VolumesOfMolecules[4] / (VolumesOfMolecules[4] + VolumesOfMolecules[5]);
    NumberOfWaterAtTag  = Parameters[6] * VolumesOfMolecules[5] / (VolumesOfMolecules[4] + VolumesOfMolecules[5]);

    AreaPerHeadgroupOfLipids = Parameters[1];
    HeightOfBelt             = Parameters[2];
    NumberOfLipids           = Parameters[3];
    RatioBetweenAxis         = Parameters[0];
    SizeOfTag                = Parameters[23];

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater  = VolumesOfMolecules[0];
    VolumeOfHead   = VolumesOfMolecules[1] * CorrectionToVolumeOfHead + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfTail   = VolumesOfMolecules[2] * CorrectionToVolumeOfTail;
    VolumeOfMethyl = VolumesOfMolecules[3] * CorrectionToVolumeOfTail;
    VolumeOfBelt   = 2.0 * (VolumesOfMolecules[4] + (SizeOfTag * VolumesOfMolecules[5])) * CorrectionToVolumeOfBelt + 2.0 * NumberOfWaterAtBelt * VolumesOfMolecules[0];
    VolumeOfTag    = (1.0 - SizeOfTag) * VolumesOfMolecules[5] * CorrectionToVolumeOfBelt + NumberOfWaterAtTag * VolumesOfMolecules[0];

    /// Compute geometry of the disc
	AverageAreaOfLipids = NumberOfLipids * AreaPerHeadgroupOfLipids / 2.0;

    MinorAxisOfCore = sqrt(abs(AverageAreaOfLipids / (pi * RatioBetweenAxis)));
    MajorAxisOfCore = MinorAxisOfCore * RatioBetweenAxis;

    HeightOfLipids = 2.0 * (VolumeOfHead + VolumeOfTail + VolumeOfMethyl) / AreaPerHeadgroupOfLipids;
    HeightOfCore   = 2.0 * (VolumeOfTail + VolumeOfMethyl) / AreaPerHeadgroupOfLipids;
    HeightOfMethyl = 2.0 * VolumeOfMethyl / AreaPerHeadgroupOfLipids;

    /// Obtain scattering length densities
    ScatteringLengthOfWater  = ScatteringLengths[0];
    ScatteringLengthOfHeads  = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    ScatteringLengthOfTail   = ScatteringLengths[2];
    ScatteringLengthOfMethyl = ScatteringLengths[3];
    ScatteringLengthOfBelt   = 2.0 * (ScatteringLengths[4] + SizeOfTag * ScatteringLengths[5] + NumberOfWaterAtBelt * ScatteringLengths[0]);
    ScatteringLengthOfTag    = (1.0 - SizeOfTag) * ScatteringLengths[5] + NumberOfWaterAtTag * ScatteringLengths[0];

    /// Derive excess scattering length densitites
    ScatteringLengthDensityOfWater        = ScatteringLengthOfWater / VolumeOfWater;
    ExcessScatteringLengthDensityOfCaps   = ScatteringLengthOfHeads / VolumeOfHead - ScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfCore   = ScatteringLengthOfTail / VolumeOfTail - ScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfMethyl = ScatteringLengthOfMethyl / VolumeOfMethyl - ScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfBelt   = ScatteringLengthOfBelt / VolumeOfBelt - ScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfTag    = ScatteringLengthOfTag / VolumeOfTag - ScatteringLengthDensityOfWater;

    /// Assign values of the constraints
    Constraints[0]  = ExcessScatteringLengthDensityOfCaps;
    Constraints[1]  = ExcessScatteringLengthDensityOfCore;
    Constraints[2]  = ExcessScatteringLengthDensityOfMethyl;
    Constraints[3]  = ExcessScatteringLengthDensityOfBelt;
    Constraints[4]  = ExcessScatteringLengthDensityOfTag;
    Constraints[5]  = 0.0;
    Constraints[6]  = HeightOfLipids;
    Constraints[7]  = HeightOfCore;
    Constraints[8]  = HeightOfMethyl;
    Constraints[9]  = 0.0;
    Constraints[10] = 0.0;
    Constraints[11] = 0.0;
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
}
