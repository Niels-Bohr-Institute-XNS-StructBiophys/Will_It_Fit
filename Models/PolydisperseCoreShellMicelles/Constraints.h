void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Variables describing scattering lengths
    double ScatteringLengthDensityOfHeads;
    double ScatteringLengthDensityOfTails;
    double ScatteringLengthOfHeads;
    double ScatteringLengthOfTails;

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
    CorrectionToVolumeOfHead  = Parameters[4];
    CorrectionToVolumeOfTail  = Parameters[5];
    CorrectionToVolumeOfWater = Parameters[6];
    NumberOfWaterAtHead       = fabs(Parameters[3]);

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater = VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfHead  = VolumesOfMolecules[1] * CorrectionToVolumeOfHead + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfTail  = VolumesOfMolecules[2] * CorrectionToVolumeOfTail;

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
    Constraints[2] = Concentration;
    Constraints[3] = VolumeOfHead;
    Constraints[4] = VolumeOfTail;
}
