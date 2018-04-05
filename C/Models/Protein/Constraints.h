// create int synonymes of parameter names
enum CONSTRAINTS {
    SLDWater, // 0
    Conc, // 1
};

enum PARAMETERS {
    BACKN100, // 0
    SCALEN100, // 1
    BACKN0, // 2
    SCALEN0, // 3
    BACKX, // 4
    SCALEX, // 5
    SCALECONC, // 6
    HYDR, // 7
    GLYCV // 8
};


void ComputeConstraints(double * Parameters, double * Volumes, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
  // Declare dummy variables needed in function
  double ScatteringLengthDensityOfWater;
  double ScatteringLengthOfWater;
  double VolumeOfWater;

  // Volume of water
  VolumeOfWater = Volumes[0]; // [Å^3]

  // Obtain scattering length densities
  ScatteringLengthOfWater  = ScatteringLengths[0]; // [cm]

  // Derive scattering lengths
  ScatteringLengthDensityOfWater = ScatteringLengthOfWater / VolumeOfWater; // [cm/Å^3]

  // Assign values of the constraints
  Constraints[SLDWater] = ScatteringLengthDensityOfWater; // [cm/Å^3]
  Constraints[Conc] = Concentration; // [1/cm^3]
}
