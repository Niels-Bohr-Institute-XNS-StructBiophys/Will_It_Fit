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
    double ScatteringLengthDensityOfHeads;
    double ScatteringLengthDensityOfTail;
    double ScatteringLengthDensityOfBelt;
    double ScatteringLengthDensityOfWater;
    double ScatteringLengthDensityOfMethyl;

    double ScatteringLengthOfBelt;
    double ScatteringLengthOfCore;
    double ScatteringLengthOfWater;
    double ScatteringLengthOfCaps;
    double ScatteringLengthOfMethyl;

    // Variables describing the valumes
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
    double NumberOfWaterAtHead;

    // Variables describing the properties of the disc
    double Radlipidcore;
    double Npeptide;

    // Variables describing the lipids
    double AverageAreaOfLipids;
    double NumberOfLipids;
    double AreaPerHeadgroupOfLipids;
    double ThicknessOfBelt;

    /// Get parameters from Parameters
    // Scale factors
    CorrectionToVolumeOfBelt  = Parameters[9];
    CorrectionToVolumeOfTail  = Parameters[10];
    CorrectionToVolumeOfHead  = Parameters[10];
    CorrectionToVolumeOfWater = Parameters[19];

    // Other parameters
    NumberOfWaterAtHead      = fabs(Parameters[5]);
    AreaPerHeadgroupOfLipids = Parameters[1];
    HeightOfBelt             = Parameters[2];
    NumberOfLipids           = Parameters[3];
    Npeptide                 = Parameters[0];

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater = VolumesOfMolecules[0];
    VolumeOfHead = VolumesOfMolecules[1] * CorrectionToVolumeOfHead + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfTail = VolumesOfMolecules[2] * CorrectionToVolumeOfTail;
    VolumeOfMethyl = VolumesOfMolecules[3] * CorrectionToVolumeOfTail;
    VolumeOfBelt = Npeptide * VolumesOfMolecules[4] * CorrectionToVolumeOfBelt;

    /// Derive geometry of the disc
	AverageAreaOfLipids = NumberOfLipids * AreaPerHeadgroupOfLipids / 2.0;
	Radlipidcore = sqrt(AverageAreaOfLipids / pi);

    HeightOfLipids = 2.0 * (VolumeOfHead + VolumeOfTail + VolumeOfMethyl) / AreaPerHeadgroupOfLipids;
    HeightOfCore   = 2.0 * (VolumeOfTail + VolumeOfMethyl) / AreaPerHeadgroupOfLipids;
    HeightOfMethyl = 2.0 * VolumeOfMethyl / AreaPerHeadgroupOfLipids;

    /// Obtain scattering length densities
    ScatteringLengthDensityOfWater  = ScatteringLengths[0];
    ScatteringLengthDensityOfHeads  = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    ScatteringLengthDensityOfTail   = ScatteringLengths[2];
    ScatteringLengthDensityOfMethyl = ScatteringLengths[3];
    ScatteringLengthDensityOfBelt   = Npeptide * ScatteringLengths[4];

    /// Derive scattering lengths
    ScatteringLengthOfWater = ScatteringLengthDensityOfWater / VolumeOfWater;
    ScatteringLengthOfCaps = ScatteringLengthDensityOfHeads / VolumeOfHead - ScatteringLengthOfWater;
    ScatteringLengthOfCore = ScatteringLengthDensityOfTail / VolumeOfTail - ScatteringLengthOfWater;
    ScatteringLengthOfMethyl = ScatteringLengthDensityOfMethyl / VolumeOfMethyl - ScatteringLengthOfWater;
    ScatteringLengthOfBelt = ScatteringLengthDensityOfBelt / VolumeOfBelt - ScatteringLengthOfWater;

    /// Assign values of the constraints
    Constraints[0] = ScatteringLengthOfCaps;
    Constraints[1] = ScatteringLengthOfCore;
    Constraints[2] = ScatteringLengthOfMethyl;
    Constraints[3] = ScatteringLengthOfBelt;
    Constraints[4] = 0.0;
    Constraints[5] = 0.0;
    Constraints[6] = HeightOfLipids;
    Constraints[7] = HeightOfCore;
    Constraints[8] = HeightOfMethyl;
    Constraints[9] = 0.0;
    Constraints[10] = 0.0;
    Constraints[11] = 0.0;
    Constraints[12] = 0.0;
    Constraints[13] = VolumesOfMolecules[2] * CorrectionToVolumeOfHead;
    Constraints[14] = VolumeOfTail;
    Constraints[15] = VolumeOfMethyl;
    Constraints[16] = VolumeOfBelt;
    Constraints[17] = 0.0;
    Constraints[18] = 0.0;

    /// Geometric properties of the disc
    ThicknessOfBelt = sqrt(pow(Radlipidcore,2) + VolumeOfBelt / (pi*HeightOfBelt)) - Radlipidcore;

    Constraints[20] = Parameters[3];
    Constraints[21] = Concentration / Npeptide;
    Constraints[22] = 0.0;
    Constraints[23] = Radlipidcore;
    Constraints[24] = 0.0;
    Constraints[25] = Parameters[3] * VolumeOfTail;
    Constraints[26] = Parameters[3] * VolumeOfHead;
    Constraints[27] = ThicknessOfBelt;
    Constraints[28] = Parameters[3] * VolumeOfMethyl;
}
