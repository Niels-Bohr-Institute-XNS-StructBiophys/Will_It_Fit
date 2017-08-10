void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Variables describing scattering lengths
    double ScatteringLengthDensityOfHeads;
    double ScatteringLengthDensityOfTails;
    double ScatteringLengthOfHeads;
    double ScatteringLengthOfTails;

    double RadiusOfCore;
    double TotalRadius;
    double TotalVolume;
    double AreaPerHeadgroup;
    double NumberOfMoleculesPerAggregate;

    // Variables describing the volumes
    double VolumeOfHead;
    double VolumeOfTail;

    double CorrectionToVolumeOfHead;
    double CorrectionToVolumeOfTail;

    // Variables describing water
    double VolumeOfWater;
    double CorrectionToVolumeOfWater;
    double NumberOfWaterAtHead;
    double ScatteringLengthOfWater;
    double ScatteringLengthDensityOfWater;

    /// Get parameters from Parameters
    NumberOfWaterAtHead       = Parameters[3];
    CorrectionToVolumeOfHead  = Parameters[4];
    CorrectionToVolumeOfTail  = Parameters[5];
    CorrectionToVolumeOfWater = Parameters[6];
    AreaPerHeadgroup          = Parameters[7];

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater = VolumesOfMolecules[0];
    VolumeOfHead  = VolumesOfMolecules[1] * CorrectionToVolumeOfHead + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfTail  = VolumesOfMolecules[2] * CorrectionToVolumeOfTail;

    /// Derive geometric properties
    TotalRadius  = 3.0 * (VolumeOfHead + VolumeOfTail) / AreaPerHeadgroup;
    TotalVolume  = 4.0 / 3.0 * pi * pow(TotalRadius, 3);
    RadiusOfCore = cbrt(3.0 * pow(TotalRadius, 2) * VolumeOfTail / AreaPerHeadgroup);

    NumberOfMoleculesPerAggregate = TotalVolume / (VolumeOfHead + VolumeOfTail);

    /// Obtain scattering length
    ScatteringLengthOfWater = ScatteringLengths[0];
    ScatteringLengthOfHeads = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    ScatteringLengthOfTails = ScatteringLengths[2];

    /// Derive scattering length densities
    ScatteringLengthDensityOfWater = ScatteringLengthOfWater / VolumeOfWater;
    ScatteringLengthDensityOfHeads = ScatteringLengthOfHeads / VolumeOfHead - ScatteringLengthDensityOfWater;
    ScatteringLengthDensityOfTails = ScatteringLengthOfTails / VolumeOfTail - ScatteringLengthDensityOfWater;

    /// Assign values of the constraints
    Constraints[0] = ScatteringLengthDensityOfHeads;
    Constraints[1] = ScatteringLengthDensityOfTails;
    Constraints[2] = TotalRadius;
    Constraints[3] = RadiusOfCore;
    Constraints[4] = Concentration;
    Constraints[5] = NumberOfMoleculesPerAggregate;
    Constraints[6] = VolumeOfHead;
    Constraints[7] = VolumeOfTail;
}
