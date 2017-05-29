void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Variables describing scattering lengths
    double ScatteringLengthDensity;
    double ScatteringLength;
    double Radius;

    // Variables describing the volumes
    double VolumeOfLipid;
    double CorrectionToVolume;
    double Thickness;

    // Variables describing water
    double VolumeOfWater;
    double CorrectionToVolumeOfWater;
    double NumberOfWaterAtHead;
    double ScatteringLengthOfWater;
    double ScatteringLengthDensityOfWater;

    /// Get parameters from Parameters
    CorrectionToVolume = Parameters[5];
    CorrectionToVolumeOfWater = Parameters[11];
    NumberOfWaterAtHead = fabs(Parameters[3]);
    Radius = Parameters[4];

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater = VolumesOfMolecules[0];
    VolumeOfLipid = (VolumesOfMolecules[1] + VolumesOfMolecules[2] + VolumesOfMolecules[3]) * CorrectionToVolume +
                    NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;

	Thickness = Parameters[15];

    /// Obtain scattering length densities
    ScatteringLengthOfWater = ScatteringLengths[0];
    ScatteringLength = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0] +
                       ScatteringLengths[2] + ScatteringLengths[3];

    /// Derive scattering lengths
    ScatteringLengthDensityOfWater = ScatteringLengthOfWater / VolumeOfWater;
    ScatteringLengthDensity = ScatteringLength / VolumeOfLipid - ScatteringLengthDensityOfWater;

    /// Assign values of the constraints
    Constraints[0] = ScatteringLengthDensity;
    Constraints[1] = 0.0;
    Constraints[2] = Thickness;
    Constraints[3] = 0.0;
    Constraints[4] = Concentration;
    Constraints[5] = Radius;
}
