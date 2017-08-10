void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Variables describing scattering lengths
    double ScatteringLengthDensityOfHeads;
    double ScatteringLengthDensityOfTails;
    double ScatteringLengthOfHeads;
    double ScatteringLengthOfTails;
    double Radius;

    // Variables describing the volumes
    double VolumeOfHead;
    double VolumeOfTail;
    double CorrectionToVolumeOfHead;
    double CorrectionToVolumeOfTail;
    double ThicknessOfHeads;
    double ThicknessOfTails;

    // Variables describing water
    double VolumeOfWater;
    double CorrectionToVolumeOfWater;
    double NumberOfWaterAtHead;
    double ScatteringLengthOfWater;
    double ScatteringLengthDensityOfWater;

    /// Get parameters from Parameters
    CorrectionToVolumeOfHead  = Parameters[5];
    CorrectionToVolumeOfTail  = Parameters[5];
    CorrectionToVolumeOfWater = Parameters[11];
    NumberOfWaterAtHead       = fabs(Parameters[3]);
    Radius                    = Parameters[4];

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater = VolumesOfMolecules[0];
    VolumeOfHead  = VolumesOfMolecules[1] * CorrectionToVolumeOfHead + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfTail  = (VolumesOfMolecules[2] + VolumesOfMolecules[3]) * CorrectionToVolumeOfTail;

    /// Compute bilayer thickness
	ThicknessOfTails = Parameters[19] * VolumeOfTail / (VolumeOfTail + VolumeOfHead);
	ThicknessOfHeads = Parameters[19] * VolumeOfHead / (VolumeOfTail + VolumeOfHead);

    /// Obtain scattering lengths
    ScatteringLengthOfWater = ScatteringLengths[0];
    ScatteringLengthOfHeads = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    ScatteringLengthOfTails = ScatteringLengths[2] + ScatteringLengths[3];

    /// Derive scattering length densities
    ScatteringLengthDensityOfWater = ScatteringLengthOfWater / VolumeOfWater;
    ScatteringLengthDensityOfHeads = ScatteringLengthOfHeads / VolumeOfHead - ScatteringLengthDensityOfWater;
    ScatteringLengthDensityOfTails = ScatteringLengthOfTails / VolumeOfTail - ScatteringLengthDensityOfWater;

    /// Assign values of the constraints
    Constraints[0] = ScatteringLengthDensityOfHeads;
    Constraints[1] = ScatteringLengthDensityOfTails;
    Constraints[2] = ThicknessOfHeads;
    Constraints[3] = ThicknessOfTails;
    Constraints[4] = Concentration;
    Constraints[5] = Radius;
}
