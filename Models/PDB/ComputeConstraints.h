void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    // Declare dummy variables needed in function
    double ScatteringLengthDensityOfWater;
    double ScatteringLengthOfWater;
    double VolumeOfWater;

    // Volume of water
    VolumeOfWater = VolumesOfMolecules[0];

    // Obtain scattering length densities
    ScatteringLengthOfWater = ScatteringLengths[0];

    // Derive scattering lengths
    ScatteringLengthDensityOfWater = ScatteringLengthOfWater / VolumeOfWater;
        
    // Assign values of the constraints
    Constraints[0] = ScatteringLengthDensityOfWater;
    Constraints[1] = Concentration;
}
